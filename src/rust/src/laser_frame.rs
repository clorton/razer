use extendr_api::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

// ── helpers ───────────────────────────────────────────────────────────────────

/// In-place compaction: slide kept rows down toward index 0. Returns new count.
fn squash_slice<T: Copy>(v: &mut [T], mask: &[bool], count: usize) -> usize {
    let mut dst = 0usize;
    for src in 0..count {
        if mask[src] {
            v[dst] = v[src];
            dst += 1;
        }
    }
    dst
}

fn robj_to_integers(src: &Robj, expected_len: usize) -> Vec<i32> {
    let vals: Vec<i32> = if let Some(v) = src.as_integer_vector() {
        v
    } else if let Some(v) = src.as_real_vector() {
        v.iter().map(|&r| r as i32).collect()
    } else {
        panic!("expected integer or real vector")
    };
    assert_eq!(
        vals.len(),
        expected_len,
        "length mismatch: expected {expected_len}, got {}",
        vals.len()
    );
    vals
}

fn robj_to_reals(src: &Robj, expected_len: usize) -> Vec<f64> {
    let vals: Vec<f64> = if let Some(v) = src.as_real_vector() {
        v
    } else if let Some(v) = src.as_integer_vector() {
        v.iter().map(|&i| i as f64).collect()
    } else {
        panic!("expected real or integer vector")
    };
    assert_eq!(
        vals.len(),
        expected_len,
        "length mismatch: expected {expected_len}, got {}",
        vals.len()
    );
    vals
}

fn robj_to_logicals(src: &Robj, expected_len: usize) -> Vec<bool> {
    let vals: Vec<bool> = src
        .as_logical_vector()
        .unwrap_or_else(|| panic!("expected logical vector"))
        .iter()
        .map(|b| b.is_true())
        .collect();
    assert_eq!(
        vals.len(),
        expected_len,
        "length mismatch: expected {expected_len}, got {}",
        vals.len()
    );
    vals
}

// ── PropData ──────────────────────────────────────────────────────────────────

/// Typed backing store for one property, mirroring R's three native vector types.
pub(crate) enum PropData {
    Integer(Vec<i32>),
    Real(Vec<f64>),
    Logical(Vec<bool>),
}

impl PropData {
    pub fn fill(n: usize, dtype: &str, default: &Robj) -> Self {
        match dtype {
            "integer" => {
                let d = default.as_integer().unwrap_or(0);
                PropData::Integer(vec![d; n])
            }
            "real" | "double" => {
                let d = default.as_real().unwrap_or(0.0);
                PropData::Real(vec![d; n])
            }
            "logical" => {
                let d = default.as_logical().map_or(false, |b| b.is_true());
                PropData::Logical(vec![d; n])
            }
            other => panic!(
                "unknown dtype '{other}'; use \"integer\", \"real\", or \"logical\""
            ),
        }
    }

    pub fn dtype_name(&self) -> &'static str {
        match self {
            PropData::Integer(_) => "integer",
            PropData::Real(_) => "real",
            PropData::Logical(_) => "logical",
        }
    }

    /// Reorder `data[0..count]` by the 0-based permutation `perm`.
    pub fn permute_rows(&mut self, perm: &[usize], count: usize) {
        match self {
            PropData::Integer(v) => {
                let tmp: Vec<i32> = perm.iter().map(|&i| v[i]).collect();
                v[..count].copy_from_slice(&tmp);
            }
            PropData::Real(v) => {
                let tmp: Vec<f64> = perm.iter().map(|&i| v[i]).collect();
                v[..count].copy_from_slice(&tmp);
            }
            PropData::Logical(v) => {
                let tmp: Vec<bool> = perm.iter().map(|&i| v[i]).collect();
                v[..count].copy_from_slice(&tmp);
            }
        }
    }

    /// Return `data[start..end]` as an R vector.
    pub fn slice_to_robj(&self, start: usize, end: usize) -> Robj {
        match self {
            PropData::Integer(v) => v[start..end].iter().collect_robj(),
            PropData::Real(v) => v[start..end].iter().collect_robj(),
            PropData::Logical(v) => v[start..end]
                .iter()
                .map(|&b| Rbool::from(b))
                .collect_robj(),
        }
    }

    /// Write `src` into `data[start..end]`, coercing integer↔real as needed.
    pub fn set_from_robj(&mut self, start: usize, end: usize, src: &Robj) {
        let n = end - start;
        match self {
            PropData::Integer(v) => {
                let vals = robj_to_integers(src, n);
                v[start..end].copy_from_slice(&vals);
            }
            PropData::Real(v) => {
                let vals = robj_to_reals(src, n);
                v[start..end].copy_from_slice(&vals);
            }
            PropData::Logical(v) => {
                let vals = robj_to_logicals(src, n);
                for (i, b) in vals.iter().enumerate() {
                    v[start + i] = *b;
                }
            }
        }
    }

    /// Compact `data[0..count]`, keeping elements where `mask[i]` is true. Returns new count.
    pub fn squash_rows(&mut self, mask: &[bool], count: usize) -> usize {
        match self {
            PropData::Integer(v) => squash_slice(v, mask, count),
            PropData::Real(v) => squash_slice(v, mask, count),
            PropData::Logical(v) => squash_slice(v, mask, count),
        }
    }
}

// ── ScalarProp ────────────────────────────────────────────────────────────────

/// A 1-D per-entry property; backing array has length `capacity`.
pub(crate) struct ScalarProp {
    pub data: PropData,
}

// ── VectorProp ────────────────────────────────────────────────────────────────

/// A 2-D per-entry property, stored **column-major** (R / Fortran order):
///
/// - `nrows` = `capacity` of the parent `LaserFrame` (agents or nodes)
/// - `ncols` = `length` argument (e.g., `nticks + 1` for a time-series)
///
/// Element `[row, col]` lives at `data[row + col * nrows]`.
/// A full column `[:, col]` occupies `data[col*nrows .. (col+1)*nrows]` — contiguous.
pub(crate) struct VectorProp {
    pub data: PropData,
    pub nrows: usize,
    pub ncols: usize,
}

impl VectorProp {
    /// Return column `col` (0-based), active rows only, as an R vector.
    pub fn get_col(&self, col: usize, active_rows: usize) -> Robj {
        let start = col * self.nrows;
        self.data.slice_to_robj(start, start + active_rows)
    }

    /// Overwrite column `col` (0-based), active rows only, from an R vector.
    pub fn set_col(&mut self, col: usize, active_rows: usize, src: &Robj) {
        let start = col * self.nrows;
        self.data.set_from_robj(start, start + active_rows, src);
    }

    /// Return the active portion as an R matrix of shape `(active_rows, ncols)`.
    ///
    /// Because R matrices are column-major and our backing store is also column-major,
    /// each output column is a single contiguous copy from the backing store.
    pub fn to_rmatrix(&self, active_rows: usize) -> Robj {
        match &self.data {
            PropData::Integer(v) => {
                let mut out = Vec::with_capacity(active_rows * self.ncols);
                for col in 0..self.ncols {
                    let start = col * self.nrows;
                    out.extend_from_slice(&v[start..start + active_rows]);
                }
                let mut robj = out.iter().collect_robj();
                robj.set_attrib("dim", vec![active_rows as i32, self.ncols as i32])
                    .expect("set dim");
                robj
            }
            PropData::Real(v) => {
                let mut out = Vec::with_capacity(active_rows * self.ncols);
                for col in 0..self.ncols {
                    let start = col * self.nrows;
                    out.extend_from_slice(&v[start..start + active_rows]);
                }
                let mut robj = out.iter().collect_robj();
                robj.set_attrib("dim", vec![active_rows as i32, self.ncols as i32])
                    .expect("set dim");
                robj
            }
            PropData::Logical(v) => {
                let mut out: Vec<Rbool> = Vec::with_capacity(active_rows * self.ncols);
                for col in 0..self.ncols {
                    let start = col * self.nrows;
                    out.extend(
                        v[start..start + active_rows]
                            .iter()
                            .map(|&b| Rbool::from(b)),
                    );
                }
                let mut robj = out.iter().collect_robj();
                robj.set_attrib("dim", vec![active_rows as i32, self.ncols as i32])
                    .expect("set dim");
                robj
            }
        }
    }

    /// Compact each column's active rows, keeping rows where `mask[i]` is true.
    pub fn squash_rows(&mut self, mask: &[bool], count: usize) {
        let nrows = self.nrows;
        match &mut self.data {
            PropData::Integer(v) => {
                for col in 0..self.ncols {
                    squash_slice(&mut v[col * nrows..], mask, count);
                }
            }
            PropData::Real(v) => {
                for col in 0..self.ncols {
                    squash_slice(&mut v[col * nrows..], mask, count);
                }
            }
            PropData::Logical(v) => {
                for col in 0..self.ncols {
                    squash_slice(&mut v[col * nrows..], mask, count);
                }
            }
        }
    }
}

// ── LaserFrame ────────────────────────────────────────────────────────────────

/// Fixed-capacity struct-of-arrays population or patch data store.
///
/// Mirrors `laser.core.LaserFrame` from Python. Each property occupies a
/// pre-allocated backing array; `count` tracks how many entries are active.
///
/// **Scalar properties** store one value per entry (length = `capacity`).
///
/// **Vector properties** store `length` values per entry, laid out column-major
/// so that `get_col(name, col)` returns a contiguous memory slice. This is
/// R's native matrix storage order and makes per-tick node-state updates fast.
///
/// @export
#[extendr]
pub struct LaserFrame {
    pub(crate) capacity: usize,
    pub(crate) count: usize,
    pub(crate) scalars: HashMap<String, ScalarProp>,
    pub(crate) vectors: HashMap<String, VectorProp>,
}

impl LaserFrame {
    fn assert_name_free(&self, name: &str) {
        assert!(
            !self.scalars.contains_key(name) && !self.vectors.contains_key(name),
            "property '{name}' already exists on this LaserFrame"
        );
    }
}

#[extendr]
impl LaserFrame {
    // ── construction ──────────────────────────────────────────────────────────

    /// Create a new `LaserFrame`.
    ///
    /// @param capacity Maximum number of entries (agents or nodes). Must be positive.
    /// @param initial_count Number of entries active at construction.
    ///   Pass `-1` (the default) to set active count equal to `capacity`.
    /// @return A new `LaserFrame` object.
    /// @export
    fn new(capacity: i64, initial_count: i64) -> Self {
        assert!(capacity > 0, "capacity must be positive, got {capacity}");
        let capacity = capacity as usize;
        let count = if initial_count < 0 {
            capacity
        } else {
            assert!(
                initial_count as usize <= capacity,
                "initial_count ({initial_count}) must not exceed capacity ({capacity})"
            );
            initial_count as usize
        };
        LaserFrame {
            capacity,
            count,
            scalars: HashMap::new(),
            vectors: HashMap::new(),
        }
    }

    // ── metadata ──────────────────────────────────────────────────────────────

    /// Number of currently active entries.
    /// @export
    fn count(&self) -> i64 {
        self.count as i64
    }

    /// Total capacity (fixed at construction).
    /// @export
    fn capacity(&self) -> i64 {
        self.capacity as i64
    }

    /// Names of all scalar properties, sorted alphabetically.
    /// @export
    fn scalar_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self.scalars.keys().cloned().collect();
        names.sort();
        names
    }

    /// Names of all vector properties, sorted alphabetically.
    /// @export
    fn vector_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self.vectors.keys().cloned().collect();
        names.sort();
        names
    }

    /// Number of columns in a named vector property.
    ///
    /// @param name Vector property name.
    /// @export
    fn vector_ncols(&self, name: &str) -> i64 {
        self.vectors
            .get(name)
            .unwrap_or_else(|| panic!("no vector property '{name}'"))
            .ncols as i64
    }

    /// Human-readable summary of the frame's properties and memory layout.
    /// @export
    fn describe(&self) -> String {
        let mut s = format!(
            "LaserFrame  capacity={}, count={}\n",
            self.capacity, self.count
        );

        let mut snames: Vec<&str> = self.scalars.keys().map(String::as_str).collect();
        snames.sort();
        if !snames.is_empty() {
            s.push_str("  Scalars:\n");
            for n in snames {
                s.push_str(&format!(
                    "    {n}  [{}]\n",
                    self.scalars[n].data.dtype_name()
                ));
            }
        }

        let mut vnames: Vec<&str> = self.vectors.keys().map(String::as_str).collect();
        vnames.sort();
        if !vnames.is_empty() {
            s.push_str("  Vectors:\n");
            for n in vnames {
                let p = &self.vectors[n];
                s.push_str(&format!(
                    "    {n}  [{}, {} cols]\n",
                    p.data.dtype_name(),
                    p.ncols
                ));
            }
        }
        s
    }

    // ── property registration ─────────────────────────────────────────────────

    /// Add a scalar property (one value per entry).
    ///
    /// @param name Property name. Must not already exist.
    /// @param dtype `"integer"`, `"real"`, or `"logical"`.
    /// @param default Fill value for the backing array.
    /// @export
    fn add_scalar_property(&mut self, name: &str, dtype: &str, default: Robj) {
        self.assert_name_free(name);
        let data = PropData::fill(self.capacity, dtype, &default);
        self.scalars.insert(name.to_string(), ScalarProp { data });
    }

    /// Add a vector property (`capacity × length`, column-major).
    ///
    /// The backing array has `capacity * length` elements.
    /// Element `[entry, col]` is stored at offset `entry + col * capacity`.
    /// Column `col` (all active entries for one time-step) is contiguous.
    ///
    /// @param name Property name. Must not already exist.
    /// @param length Number of columns (e.g., `nticks + 1` for a time-series).
    /// @param dtype `"integer"`, `"real"`, or `"logical"`.
    /// @param default Fill value.
    /// @export
    fn add_vector_property(&mut self, name: &str, length: i64, dtype: &str, default: Robj) {
        self.assert_name_free(name);
        let ncols = length as usize;
        let data = PropData::fill(self.capacity * ncols, dtype, &default);
        self.vectors.insert(
            name.to_string(),
            VectorProp {
                data,
                nrows: self.capacity,
                ncols,
            },
        );
    }

    // ── lifecycle ─────────────────────────────────────────────────────────────

    /// Activate `n` additional entries and return their 1-based index range.
    ///
    /// Returns `c(start, end)` (both inclusive, 1-based) so that
    /// `frame$get("prop")[start:end]` addresses the newly activated entries.
    ///
    /// @param n Number of entries to activate.
    /// @export
    fn add(&mut self, n: i64) -> Vec<i32> {
        let n = n as usize;
        assert!(
            self.count + n <= self.capacity,
            "add({n}) would exceed capacity: count={}, capacity={}",
            self.count,
            self.capacity
        );
        let start = self.count + 1; // 1-based
        self.count += n;
        vec![start as i32, self.count as i32]
    }

    /// Compact active entries, keeping only those where `mask` is `TRUE`.
    ///
    /// All scalar and vector properties are squashed in place.
    /// `count` is updated to the number of kept entries.
    ///
    /// @param mask Logical vector of length `count`.
    /// @export
    fn squash(&mut self, mask: Robj) {
        let mask_vec = robj_to_logicals(&mask, self.count);
        let count = self.count;
        let new_count = mask_vec.iter().filter(|&&b| b).count();

        for prop in self.scalars.values_mut() {
            prop.data.squash_rows(&mask_vec, count);
        }
        for prop in self.vectors.values_mut() {
            prop.squash_rows(&mask_vec, count);
        }
        self.count = new_count;
    }

    /// Reorder active scalar properties by `perm` (1-based permutation of length `count`).
    ///
    /// Each property is permuted in parallel across available CPU threads via Rayon.
    /// Only scalar properties are reordered; vector properties (time-series) are left unchanged.
    ///
    /// @param perm Integer vector of length `count`. Must be a valid permutation of `1:count`.
    /// @export
    fn sort_by(&mut self, perm: Robj) {
        let perm_1based = perm
            .as_integer_vector()
            .expect("perm must be an integer vector");
        assert_eq!(
            perm_1based.len(),
            self.count,
            "perm must have length count={}, got {}",
            self.count,
            perm_1based.len()
        );

        // Convert 1-based R indices to 0-based Rust indices
        let perm_0based: Vec<usize> = perm_1based
            .iter()
            .map(|&i| {
                assert!(
                    i >= 1 && i as usize <= self.count,
                    "perm value {i} is out of range (must be in 1..={})",
                    self.count
                );
                (i - 1) as usize
            })
            .collect();

        let count = self.count;
        // Permute each scalar property independently, in parallel across properties.
        self.scalars
            .par_iter_mut()
            .for_each(|(_, prop)| prop.data.permute_rows(&perm_0based, count));
    }

    // ── scalar access ─────────────────────────────────────────────────────────

    /// Return the active slice of a scalar property as an R vector.
    ///
    /// @param name Scalar property name.
    /// @export
    fn get(&self, name: &str) -> Robj {
        let prop = self
            .scalars
            .get(name)
            .unwrap_or_else(|| panic!("no scalar property '{name}'"));
        prop.data.slice_to_robj(0, self.count)
    }

    /// Overwrite the active slice of a scalar property from an R vector.
    ///
    /// Integer and real vectors are accepted for integer and real properties
    /// respectively; integer values are accepted for real properties (coerced).
    ///
    /// @param name Scalar property name.
    /// @param values R vector of length `count`.
    /// @export
    fn set(&mut self, name: &str, values: Robj) {
        let count = self.count;
        let prop = self
            .scalars
            .get_mut(name)
            .unwrap_or_else(|| panic!("no scalar property '{name}'"));
        prop.data.set_from_robj(0, count, &values);
    }

    // ── vector property access ────────────────────────────────────────────────

    /// Return one column of a vector property as an R vector.
    ///
    /// The column is returned for active entries only (length = `count`).
    /// Columns are **1-based**: column 1 is the first column.
    ///
    /// @param name Vector property name.
    /// @param col Column index (1-based).
    /// @export
    fn get_col(&self, name: &str, col: i64) -> Robj {
        let prop = self
            .vectors
            .get(name)
            .unwrap_or_else(|| panic!("no vector property '{name}'"));
        assert!(
            col >= 1 && col as usize <= prop.ncols,
            "col {col} out of range (1..={})",
            prop.ncols
        );
        prop.get_col((col - 1) as usize, self.count)
    }

    /// Overwrite one column of a vector property from an R vector.
    ///
    /// @param name Vector property name.
    /// @param col Column index (1-based).
    /// @param values R vector of length `count`.
    /// @export
    fn set_col(&mut self, name: &str, col: i64, values: Robj) {
        let count = self.count;
        let prop = self
            .vectors
            .get_mut(name)
            .unwrap_or_else(|| panic!("no vector property '{name}'"));
        assert!(
            col >= 1 && col as usize <= prop.ncols,
            "col {col} out of range (1..={})",
            prop.ncols
        );
        prop.set_col((col - 1) as usize, count, &values);
    }

    /// Return a vector property as an R matrix of shape `(count, ncols)`.
    ///
    /// Because both R and the backing store use column-major layout, each
    /// column of the returned matrix is a direct contiguous copy.
    ///
    /// @param name Vector property name.
    /// @export
    fn get_matrix(&self, name: &str) -> Robj {
        let prop = self
            .vectors
            .get(name)
            .unwrap_or_else(|| panic!("no vector property '{name}'"));
        prop.to_rmatrix(self.count)
    }
}

extendr_module! {
    mod laser_frame;
    impl LaserFrame;
}
