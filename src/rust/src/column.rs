// ════════════════════════════════════════════════════════════════════════════
// column.rs — Rust-owned, dtype-tagged 1-D property arrays exported to R.
//
// A `Column` is the backing store for one agent property (e.g. `state`). The
// data lives in a Rust `Vec<T>` that R never sees directly: R holds only an
// opaque external-pointer handle (the same mechanism as `LaserFrame` /
// `Distribution`). This buys us:
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

/// A Rust-owned, dtype-tagged 1-D property array.
///
/// Allocate one with [allocate_scalar()]. The data is held in Rust and exposed to
/// R only as an opaque handle; use `$values()` to copy a snapshot back into an R
/// vector for inspection, `$fill()` / `$set()` to write, and `$length()` /
/// `$dtype()` to query. The simulation step kernels operate on the buffer in
/// place with no copies.
///
/// @export
#[extendr]
pub struct Column {
    data: Storage,
}

impl Column {
    // Element count of the live variant. `pub(crate)` for sibling modules.
    pub(crate) fn len(&self) -> usize {
        match &self.data {
            Storage::I8(v) => v.len(),   Storage::U8(v) => v.len(),
            Storage::I16(v) => v.len(),  Storage::U16(v) => v.len(),
            Storage::I32(v) => v.len(),  Storage::U32(v) => v.len(),
            Storage::F32(v) => v.len(),  Storage::F64(v) => v.len(),
        }
    }

    // Borrow the typed backing store (for crate-internal kernels to dispatch on).
    pub(crate) fn storage(&self) -> &Storage {
        &self.data
    }

    pub(crate) fn storage_mut(&mut self) -> &mut Storage {
        &mut self.data
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
    /// @export
    fn length(&self) -> i32 {
        self.len() as i32
    }

    /// The element data type as a string (e.g. `"u8"`, `"f32"`).
    /// @return A length-1 character vector.
    /// @export
    fn dtype(&self) -> String {
        self.dtype_enum().name().to_string()
    }

    /// Copy the array into an R vector for inspection (NOT a view — a snapshot).
    ///
    /// Integer-width types (i8, u8, i16, u16, i32) widen to R `integer`; `u32`,
    /// `f32`, and `f64` widen to R `double` (since `u32` overflows R's signed
    /// 32-bit integer). This O(n) copy is the only place data leaves Rust.
    ///
    /// @return A numeric vector (integer or double) of length `length()`.
    /// @export
    fn values(&self) -> Robj {
        match &self.data {
            Storage::I8(v)  => v.iter().map(|&x| x as i32).collect::<Vec<i32>>().iter().collect_robj(),
            Storage::U8(v)  => v.iter().map(|&x| x as i32).collect::<Vec<i32>>().iter().collect_robj(),
            Storage::I16(v) => v.iter().map(|&x| x as i32).collect::<Vec<i32>>().iter().collect_robj(),
            Storage::U16(v) => v.iter().map(|&x| x as i32).collect::<Vec<i32>>().iter().collect_robj(),
            Storage::I32(v) => v.iter().collect_robj(),
            Storage::U32(v) => v.iter().map(|&x| x as f64).collect::<Vec<f64>>().iter().collect_robj(),
            Storage::F32(v) => v.iter().map(|&x| x as f64).collect::<Vec<f64>>().iter().collect_robj(),
            Storage::F64(v) => v.iter().collect_robj(),
        }
    }

    /// Set every element to `value`, cast to the array's data type.
    ///
    /// For integer-typed arrays the value is truncated toward zero (e.g. `2.9`
    /// becomes `2`); out-of-range values wrap per Rust's `as` cast.
    ///
    /// @param value A single numeric value to broadcast across the array.
    /// @export
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
    /// @export
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
    let data = match DType::parse(dtype) {
        DType::I8  => Storage::I8(vec![0; n]),
        DType::U8  => Storage::U8(vec![0; n]),
        DType::I16 => Storage::I16(vec![0; n]),
        DType::U16 => Storage::U16(vec![0; n]),
        DType::I32 => Storage::I32(vec![0; n]),
        DType::U32 => Storage::U32(vec![0; n]),
        DType::F32 => Storage::F32(vec![0.0; n]),
        DType::F64 => Storage::F64(vec![0.0; n]),
    };
    Column { data }
}

extendr_module! {
    mod column;
    impl Column;
    fn allocate_scalar;
}
