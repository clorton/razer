// ════════════════════════════════════════════════════════════════════════════
// column.rs — Rust-owned, dtype-tagged 1-D property arrays exported to R.
//
// A `Column` is the backing store for one agent property (e.g. `state`). The
// data lives in a Rust `Vec<T>` that R never sees directly: R holds only an
// opaque external-pointer handle (the same mechanism as `Distribution`). This buys us:
//
//   * Every integer/float width — i8/u8/i16/u16/i32/u32/f32/f64 — even though R's
//     own atomic vectors only cover i32 (integer), f64 (double), and u8 (raw).
//   * Exact memory: `sizeof(T)` per element, no per-element R header or boxing.
//   * Zero-copy, in-place mutation: the step kernels borrow `&mut [T]` straight
//     into the Vec — no R copy-on-modify, and the slices are `Send`/`Sync` so
//     Rayon can split them across worker threads.
//
// Copying back into an R vector happens ONLY on explicit inspection (`$values()`),
// widening to the nearest native R type. That copy is the accepted, secondary
// cost; it never touches the compute path.
//
// Orientation for readers coming from C / C++ / C# / Python:
//   * `enum Storage` is a tagged union (like a C `union` + discriminant, or a
//     numpy dtype). Each variant owns a `Vec<T>` (a growable array / std::vector).
//   * `match &mut self.data { Storage::U8(v) => ... }` destructures the active
//     variant and borrows its Vec mutably — the compiler proves only one variant
//     is live, so this is the safe, checked equivalent of a union switch.
//   * `x as i8` is an explicit, possibly-truncating numeric cast (never implicit).
//   * `#[extendr]` on the struct generates the external-pointer wrapper; on the
//     `impl` block it exposes each method to R as `column$method(...)`.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;

// The closed set of supported element types. Parsed from the R-side dtype string.
#[derive(Clone, Copy, PartialEq, Eq)]
enum DType {
    I8, U8, I16, U16, I32, U32, F32, F64,
}

impl DType {
    // Map an R-side dtype string (with common aliases) to a DType, or abort.
    fn parse(s: &str) -> DType {
        match s {
            "i8"  | "int8"            => DType::I8,
            "u8"  | "uint8"  | "raw"  => DType::U8,
            "i16" | "int16"           => DType::I16,
            "u16" | "uint16"          => DType::U16,
            "i32" | "int32"  | "integer" => DType::I32,
            "u32" | "uint32"          => DType::U32,
            "f32" | "float32" | "single" => DType::F32,
            "f64" | "float64" | "double" | "real" => DType::F64,
            other => panic!(
                "unknown dtype '{other}'; use one of i8, u8, i16, u16, i32, u32, f32, f64"
            ),
        }
    }

    fn name(self) -> &'static str {
        match self {
            DType::I8 => "i8",   DType::U8 => "u8",
            DType::I16 => "i16", DType::U16 => "u16",
            DType::I32 => "i32", DType::U32 => "u32",
            DType::F32 => "f32", DType::F64 => "f64",
        }
    }
}

// The typed, owned backing buffer. One variant is live per Column. `pub(crate)`
// so sibling modules (e.g. `bincount`) can match on the element type to dispatch
// to a generic kernel; it stays invisible to R (R only sees the `Column` handle).
pub(crate) enum Storage {
    I8(Vec<i8>),   U8(Vec<u8>),
    I16(Vec<i16>), U16(Vec<u16>),
    I32(Vec<i32>), U32(Vec<u32>),
    F32(Vec<f32>), F64(Vec<f64>),
}

// Allocate a zero-filled `Storage` of the given dtype and element count.
fn zeroed_storage(dtype: DType, n: usize) -> Storage {
    match dtype {
        DType::I8  => Storage::I8(vec![0; n]),
        DType::U8  => Storage::U8(vec![0; n]),
        DType::I16 => Storage::I16(vec![0; n]),
        DType::U16 => Storage::U16(vec![0; n]),
        DType::I32 => Storage::I32(vec![0; n]),
        DType::U32 => Storage::U32(vec![0; n]),
        DType::F32 => Storage::F32(vec![0.0; n]),
        DType::F64 => Storage::F64(vec![0.0; n]),
    }
}

/// A Rust-owned, dtype-tagged property array (1-D scalar, or a 2-D vector report).
///
/// Allocate one with [allocate_scalar()] (`nrows` elements, `ncols == 1`) or
/// [allocate_vector()] (`nrows`-per-slot × `ncols` slots, stored COLUMN-MAJOR so a
/// whole slot/time-slice is contiguous). The data is held in Rust and exposed to
/// R only as an opaque handle; use `$values()` to copy a snapshot back into an R
/// vector (or matrix, when `ncols > 1`) for inspection, `$fill()` / `$set()` to
/// write, and `$length()` / `$dtype()` to query. The simulation step kernels
/// operate on the buffer in place with no copies.
///
/// @export
#[extendr]
pub struct Column {
    data: Storage,
    // Logical shape, SLICE-MAJOR (a.k.a. row-major): `n_slices` slices, each of
    // `slice_len` contiguous elements. The backing Vec holds `n_slices * slice_len`
    // elements; slice `s` is the contiguous range `s*slice_len .. (s+1)*slice_len`,
    // so indexing the FIRST (slice) dimension returns a contiguous block. A scalar
    // column has `n_slices == 1`. For a report buffer the first/slice dimension is
    // TIME and the inner dimension is NODE — shape (n_ticks, n_nodes) — so a whole
    // tick's per-node values are contiguous.
    n_slices: usize,
    slice_len: usize,
}

impl Column {
    // Total element count (`n_slices * slice_len`). `pub(crate)` for sibling modules.
    pub(crate) fn len(&self) -> usize {
        match &self.data {
            Storage::I8(v) => v.len(),   Storage::U8(v) => v.len(),
            Storage::I16(v) => v.len(),  Storage::U16(v) => v.len(),
            Storage::I32(v) => v.len(),  Storage::U32(v) => v.len(),
            Storage::F32(v) => v.len(),  Storage::F64(v) => v.len(),
        }
    }

    // Shape accessors for crate-internal kernels: number of slices (e.g. ticks),
    // and the contiguous length of each slice (e.g. node count).
    pub(crate) fn n_slices(&self) -> usize { self.n_slices }
    pub(crate) fn slice_len(&self) -> usize { self.slice_len }

    // Copy the data into an R vector, mapping each output position k through `idx`
    // (a buffer index). Used by `values()`: identity for a scalar, transposed for a
    // 2-D column. Integer widths widen to R `integer`, u32/f32/f64 to R `double`.
    fn gather_to_robj(&self, n: usize, idx: impl Fn(usize) -> usize) -> Robj {
        match &self.data {
            Storage::I8(v)  => (0..n).map(|k| v[idx(k)] as i32).collect::<Vec<i32>>().iter().collect_robj(),
            Storage::U8(v)  => (0..n).map(|k| v[idx(k)] as i32).collect::<Vec<i32>>().iter().collect_robj(),
            Storage::I16(v) => (0..n).map(|k| v[idx(k)] as i32).collect::<Vec<i32>>().iter().collect_robj(),
            Storage::U16(v) => (0..n).map(|k| v[idx(k)] as i32).collect::<Vec<i32>>().iter().collect_robj(),
            Storage::I32(v) => (0..n).map(|k| v[idx(k)]).collect::<Vec<i32>>().iter().collect_robj(),
            Storage::U32(v) => (0..n).map(|k| v[idx(k)] as f64).collect::<Vec<f64>>().iter().collect_robj(),
            Storage::F32(v) => (0..n).map(|k| v[idx(k)] as f64).collect::<Vec<f64>>().iter().collect_robj(),
            Storage::F64(v) => (0..n).map(|k| v[idx(k)]).collect::<Vec<f64>>().iter().collect_robj(),
        }
    }

    // Borrow the typed backing store (for crate-internal kernels to dispatch on).
    pub(crate) fn storage(&self) -> &Storage {
        &self.data
    }

    pub(crate) fn storage_mut(&mut self) -> &mut Storage {
        &mut self.data
    }

    // Typed slice accessors used by the step kernels; panic on a dtype mismatch.
    pub(crate) fn as_u8(&self) -> &[u8] {
        match &self.data {
            Storage::U8(v) => v,
            _ => panic!("expected a u8 Column"),
        }
    }

    pub(crate) fn as_u8_mut(&mut self) -> &mut [u8] {
        match &mut self.data {
            Storage::U8(v) => v,
            _ => panic!("expected a u8 Column"),
        }
    }

    pub(crate) fn as_u16(&self) -> &[u16] {
        match &self.data {
            Storage::U16(v) => v,
            _ => panic!("expected a u16 Column"),
        }
    }

    pub(crate) fn as_u16_mut(&mut self) -> &mut [u16] {
        match &mut self.data {
            Storage::U16(v) => v,
            _ => panic!("expected a u16 Column"),
        }
    }

    pub(crate) fn as_u32(&self) -> &[u32] {
        match &self.data {
            Storage::U32(v) => v,
            _ => panic!("expected a u32 Column"),
        }
    }

    pub(crate) fn as_u32_mut(&mut self) -> &mut [u32] {
        match &mut self.data {
            Storage::U32(v) => v,
            _ => panic!("expected a u32 Column"),
        }
    }

    pub(crate) fn as_i32_mut(&mut self) -> &mut [i32] {
        match &mut self.data {
            Storage::I32(v) => v,
            _ => panic!("expected an i32 Column"),
        }
    }

    pub(crate) fn as_f64(&self) -> &[f64] {
        match &self.data {
            Storage::F64(v) => v,
            _ => panic!("expected an f64 Column"),
        }
    }

    pub(crate) fn as_f64_mut(&mut self) -> &mut [f64] {
        match &mut self.data {
            Storage::F64(v) => v,
            _ => panic!("expected an f64 Column"),
        }
    }

    // Copy the whole buffer into an owned Vec<f64> (any numeric dtype widens to f64).
    pub(crate) fn to_f64(&self) -> Vec<f64> {
        match &self.data {
            Storage::I8(v)  => v.iter().map(|&x| x as f64).collect(),
            Storage::U8(v)  => v.iter().map(|&x| x as f64).collect(),
            Storage::I16(v) => v.iter().map(|&x| x as f64).collect(),
            Storage::U16(v) => v.iter().map(|&x| x as f64).collect(),
            Storage::I32(v) => v.iter().map(|&x| x as f64).collect(),
            Storage::U32(v) => v.iter().map(|&x| x as f64).collect(),
            Storage::F32(v) => v.iter().map(|&x| x as f64).collect(),
            Storage::F64(v) => v.clone(),
        }
    }

    fn dtype_enum(&self) -> DType {
        match &self.data {
            Storage::I8(_) => DType::I8,   Storage::U8(_) => DType::U8,
            Storage::I16(_) => DType::I16, Storage::U16(_) => DType::U16,
            Storage::I32(_) => DType::I32, Storage::U32(_) => DType::U32,
            Storage::F32(_) => DType::F32, Storage::F64(_) => DType::F64,
        }
    }
}

#[extendr]
impl Column {
    /// Number of elements in the array.
    /// @return An integer length.
    fn length(&self) -> i32 {
        self.len() as i32
    }

    /// The element data type as a string (e.g. `"u8"`, `"f32"`).
    /// @return A length-1 character vector.
    fn dtype(&self) -> String {
        self.dtype_enum().name().to_string()
    }

    /// Copy the array into an R vector for inspection (NOT a view — a snapshot).
    ///
    /// Integer-width types (i8, u8, i16, u16, i32) widen to R `integer`; `u32`,
    /// `f32`, and `f64` widen to R `double` (since `u32` overflows R's signed
    /// 32-bit integer). This O(n) copy is the only place data leaves Rust.
    ///
    /// For a 2-D column (`n_slices > 1`, from [allocate_vector()]) the result carries
    /// a `dim` attribute and reads back as an `n_slices × slice_len` R matrix (e.g.
    /// `n_ticks × n_nodes`), so row `t` is tick `t`'s per-node vector; otherwise a
    /// plain vector. The snapshot is transposed during the copy (our buffer is
    /// slice-major, R matrices are column-major) — inexpensive, inspection-only.
    ///
    /// @return A numeric vector (or matrix) — integer or double — of `length()` elements.
    fn values(&self) -> Robj {
        let (ns, sl) = (self.n_slices, self.slice_len);
        // Scalar / single-slice: copy in natural order, return a plain vector.
        if ns <= 1 {
            return self.gather_to_robj(ns * sl, |k| k);
        }
        // 2-D: emit column-major over (n_slices, slice_len). R reads output position
        // k as element [s, i] with s = k % ns, i = k / ns; that element lives at
        // buffer index s*sl + i (slice-major), so gather through that mapping.
        let mut robj = self.gather_to_robj(ns * sl, |k| (k % ns) * sl + k / ns);
        robj.set_attrib("dim", [ns as i32, sl as i32]).expect("set dim");
        robj
    }

    /// Set every element to `value`, cast to the array's data type.
    ///
    /// For integer-typed arrays the value is truncated toward zero (e.g. `2.9`
    /// becomes `2`); out-of-range values wrap per Rust's `as` cast.
    ///
    /// @param value A single numeric value to broadcast across the array.
    fn fill(&mut self, value: f64) {
        match &mut self.data {
            Storage::I8(v)  => v.iter_mut().for_each(|e| *e = value as i8),
            Storage::U8(v)  => v.iter_mut().for_each(|e| *e = value as u8),
            Storage::I16(v) => v.iter_mut().for_each(|e| *e = value as i16),
            Storage::U16(v) => v.iter_mut().for_each(|e| *e = value as u16),
            Storage::I32(v) => v.iter_mut().for_each(|e| *e = value as i32),
            Storage::U32(v) => v.iter_mut().for_each(|e| *e = value as u32),
            Storage::F32(v) => v.iter_mut().for_each(|e| *e = value as f32),
            Storage::F64(v) => v.iter_mut().for_each(|e| *e = value),
        }
    }

    /// Overwrite the array from an R numeric vector (integer or double).
    ///
    /// The input length must equal `length()`. Each element is cast to the array's
    /// data type (integer-typed arrays truncate toward zero). Useful for setup —
    /// e.g. writing per-agent node assignments or seeding initial states from R.
    ///
    /// @param values A numeric vector of length `length()`.
    fn set(&mut self, values: Robj) {
        // Accept either an R integer or double vector; read both as f64 then cast.
        let vals: Vec<f64> = if let Some(s) = values.as_real_slice() {
            s.to_vec()
        } else if let Some(s) = values.as_integer_slice() {
            s.iter().map(|&i| i as f64).collect()
        } else {
            panic!("`values` must be a numeric (integer or double) vector");
        };
        assert_eq!(
            vals.len(), self.len(),
            "`values` length ({}) must equal the column length ({})",
            vals.len(), self.len()
        );
        match &mut self.data {
            Storage::I8(v)  => for (e, &x) in v.iter_mut().zip(&vals) { *e = x as i8 },
            Storage::U8(v)  => for (e, &x) in v.iter_mut().zip(&vals) { *e = x as u8 },
            Storage::I16(v) => for (e, &x) in v.iter_mut().zip(&vals) { *e = x as i16 },
            Storage::U16(v) => for (e, &x) in v.iter_mut().zip(&vals) { *e = x as u16 },
            Storage::I32(v) => for (e, &x) in v.iter_mut().zip(&vals) { *e = x as i32 },
            Storage::U32(v) => for (e, &x) in v.iter_mut().zip(&vals) { *e = x as u32 },
            Storage::F32(v) => for (e, &x) in v.iter_mut().zip(&vals) { *e = x as f32 },
            Storage::F64(v) => for (e, &x) in v.iter_mut().zip(&vals) { *e = x },
        }
    }

    /// Read one column (`slot`) of a 2-D Column as an R vector snapshot.
    ///
    /// Returns the `slice_len` values in column `slot` (e.g. all nodes for one tick),
    /// widened to R `integer`/`double` like `values()`. For a scalar column the only
    /// valid `slot` is 0 (the whole vector).
    ///
    /// @param slot 0-based column index, less than the number of columns.
    /// @return A numeric vector of length `slice_len`.
    fn col(&self, slot: i32) -> Robj {
        assert!(slot >= 0, "`slot` must be non-negative, got {slot}");
        let slot = slot as usize;
        assert!(slot < self.n_slices, "`slot` ({slot}) out of range for {} columns", self.n_slices);
        let n = self.slice_len;
        self.gather_to_robj(n, move |k| slot * n + k)
    }

    /// Write `values` into one column (`slot`) of a 2-D Column, in place.
    ///
    /// `values` must have length `slice_len`; each element is cast to the column's
    /// data type. Lets the caller update a single column (e.g. one tick's per-node
    /// slice — a derived population total) without rewriting the whole buffer.
    ///
    /// @param slot   0-based column index, less than the number of columns.
    /// @param values A numeric vector of length `slice_len`.
    fn set_col(&mut self, slot: i32, values: Robj) {
        assert!(slot >= 0, "`slot` must be non-negative, got {slot}");
        let slot = slot as usize;
        assert!(slot < self.n_slices, "`slot` ({slot}) out of range for {} columns", self.n_slices);
        let n = self.slice_len;
        let vals: Vec<f64> = if let Some(s) = values.as_real_slice() {
            s.to_vec()
        } else if let Some(s) = values.as_integer_slice() {
            s.iter().map(|&i| i as f64).collect()
        } else {
            panic!("`values` must be a numeric (integer or double) vector");
        };
        assert_eq!(vals.len(), n, "`values` length ({}) must equal slice_len ({n})", vals.len());
        let start = slot * n;
        macro_rules! write_slice { ($v:expr, $t:ty) => {
            for (e, &x) in $v[start..start + n].iter_mut().zip(&vals) { *e = x as $t; }
        } }
        match &mut self.data {
            Storage::I8(v)  => write_slice!(v, i8),
            Storage::U8(v)  => write_slice!(v, u8),
            Storage::I16(v) => write_slice!(v, i16),
            Storage::U16(v) => write_slice!(v, u16),
            Storage::I32(v) => write_slice!(v, i32),
            Storage::U32(v) => write_slice!(v, u32),
            Storage::F32(v) => write_slice!(v, f32),
            Storage::F64(v) => write_slice!(v, f64),
        }
    }

    /// Compact the first `length(keep)` elements in place, keeping those flagged `TRUE`.
    ///
    /// Drops elements where `keep` is `FALSE`/`NA`, shifting the survivors to the front
    /// (order preserved), and returns the number kept. Use it to reclaim the slots of
    /// deceased agents: apply the SAME `keep` mask to every per-agent Column (so they stay
    /// aligned) and set the active count to the returned value. The R [squash()] helper
    /// does exactly this across a people environment. Only valid for a 1-D (scalar) Column.
    ///
    /// @param keep A logical vector whose length is at most the column length (typically
    ///   the active agent count).
    /// @return The number of kept elements (an integer); elements past it are left as-is.
    fn squash(&mut self, keep: Robj) -> i32 {
        assert_eq!(self.n_slices, 1, "squash is only valid for a 1-D (scalar) Column");
        let k = keep.as_logical_slice().expect("`keep` must be a logical vector");
        let l = k.len();
        assert!(l <= self.len(), "`keep` length ({l}) exceeds column length ({})", self.len());
        // In-place stable compaction: copy each kept element down to the write cursor `w`.
        macro_rules! sq { ($v:expr) => {{
            let mut w = 0usize;
            for i in 0..l { if k[i].is_true() { $v[w] = $v[i]; w += 1; } }
            w
        }} }
        let w = match &mut self.data {
            Storage::I8(v)  => sq!(v), Storage::U8(v)  => sq!(v),
            Storage::I16(v) => sq!(v), Storage::U16(v) => sq!(v),
            Storage::I32(v) => sq!(v), Storage::U32(v) => sq!(v),
            Storage::F32(v) => sq!(v), Storage::F64(v) => sq!(v),
        };
        w as i32
    }
}

/// Allocate a fresh, zero-filled property array of a given type and length.
///
/// Returns an opaque [Column] handle backed by a Rust-owned buffer. Choose the
/// narrowest type that fits the data to minimize memory (e.g. `"u8"` for a small
/// set of disease-state codes).
///
/// @param dtype Element type, one of `"i8"`, `"u8"`, `"i16"`, `"u16"`, `"i32"`,
///   `"u32"`, `"f32"`, `"f64"` (aliases: `"int8"`, `"uint8"`/`"raw"`, …,
///   `"integer"` = i32, `"double"`/`"real"` = f64, `"single"` = f32).
/// @param count Array length (number of elements); a non-negative integer.
/// @return A [Column] object whose elements are all zero.
/// @examples
/// state <- allocate_scalar("u8", 5L)
/// state$dtype()    # "u8"
/// state$length()   # 5
/// state$values()   # 0 0 0 0 0
/// @export
#[extendr]
fn allocate_scalar(dtype: &str, count: i32) -> Column {
    assert!(count >= 0, "`count` must be non-negative, got {count}");
    let n = count as usize;
    Column { data: zeroed_storage(DType::parse(dtype), n), n_slices: 1, slice_len: n }
}

/// Allocate a fresh, zero-filled 2-D property array (a per-slot report buffer).
///
/// Returns an opaque [Column] of `n_slices * slice_len` elements laid out
/// SLICE-MAJOR (row-major): `n_slices` slices, each a contiguous run of `slice_len`
/// elements. Slice `s` is the block `s*slice_len .. (s+1)*slice_len`, so indexing
/// the FIRST dimension yields a contiguous array. The conventional use is a
/// time-series report with the **first dimension time and the second node** —
/// `allocate_vector(dtype, n_ticks, n_nodes)` — so each tick's per-node values are
/// contiguous (cache-friendly for the step kernels that fill one tick at a time).
/// `$values()` reads it back as an `n_slices × slice_len` (e.g. `n_ticks × n_nodes`)
/// R matrix, so row `t` is tick `t`'s vector.
///
/// @param dtype     Element type (see [allocate_scalar()] for the accepted names).
/// @param n_slices  Number of slices — the first/outer dimension (e.g. the tick
///   count). A non-negative integer.
/// @param slice_len Contiguous length of each slice — the inner dimension (e.g. the
///   node count). A non-negative integer.
/// @return A [Column] of shape `n_slices × slice_len`, all elements zero.
/// @examples
/// recoveries <- allocate_vector("u32", 4L, 3L)   # 4 ticks x 3 nodes
/// dim(recoveries$values())                        # 4 3
/// @export
#[extendr]
fn allocate_vector(dtype: &str, n_slices: i32, slice_len: i32) -> Column {
    assert!(n_slices >= 0, "`n_slices` must be non-negative, got {n_slices}");
    assert!(slice_len >= 0, "`slice_len` must be non-negative, got {slice_len}");
    let (n_slices, slice_len) = (n_slices as usize, slice_len as usize);
    Column {
        data: zeroed_storage(DType::parse(dtype), n_slices * slice_len),
        n_slices,
        slice_len,
    }
}

/// Carry a per-node counter forward one tick: copy column `tick` onto `tick + 1`.
///
/// For a 2-D report [Column] (e.g. `n_ticks+1 × n_nodes`), copies the contiguous
/// slice for `tick` onto the slice for `tick + 1`. This seeds the next tick with the
/// current counts so a dynamics kernel can then update it in place — keeping the
/// census invariant `count[t+1] = count[t] ± delta`. Call it once per state that
/// must persist across ticks (e.g. S, I, R for an SIR model; add E for SEIR, or a
/// user-defined "V" vaccinated count). Works for any element type.
///
/// @param counter A 2-D Column to carry forward.
/// @param tick    0-based source tick; column `tick` is copied onto column `tick+1`.
/// @return `NULL` (invisibly); `counter` is modified in place.
/// @examples
/// counts <- allocate_vector("i32", 3L, 2L)   # 3 ticks x 2 nodes
/// counts$set(c(5L, 7L, 0L, 0L, 0L, 0L))      # tick 0 = (5, 7)
/// carry_forward(counts, 0L)
/// counts$values()[2L, ]                       # 5 7  (carried to tick 1)
/// @export
#[extendr]
fn carry_forward(counter: &mut Column, tick: i32) {
    assert!(tick >= 0, "`tick` must be non-negative, got {tick}");
    let tick = tick as usize;
    let n = counter.slice_len();
    assert!(
        tick + 1 < counter.n_slices(),
        "`tick`+1 ({}) out of range for a Column with {} slices",
        tick + 1, counter.n_slices()
    );
    let (src, dst) = (tick * n, (tick + 1) * n);
    // `copy_within(src_range, dest)` is an in-place memmove within the one buffer.
    match counter.storage_mut() {
        Storage::I8(v)  => v.copy_within(src..src + n, dst),
        Storage::U8(v)  => v.copy_within(src..src + n, dst),
        Storage::I16(v) => v.copy_within(src..src + n, dst),
        Storage::U16(v) => v.copy_within(src..src + n, dst),
        Storage::I32(v) => v.copy_within(src..src + n, dst),
        Storage::U32(v) => v.copy_within(src..src + n, dst),
        Storage::F32(v) => v.copy_within(src..src + n, dst),
        Storage::F64(v) => v.copy_within(src..src + n, dst),
    }
}

extendr_module! {
    mod column;
    impl Column;
    fn allocate_scalar;
    fn allocate_vector;
    fn carry_forward;
}
