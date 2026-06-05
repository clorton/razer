// ════════════════════════════════════════════════════════════════════════════
// migration.rs — spatial helpers exported to R via extendr.
//
// Currently: `distances`, a port of laser-core's `distance` (laser/core/
// migration.py) for the common all-pairs case. It builds the N×N great-circle
// distance matrix (kilometres) between geographic points using the haversine
// formula.
//
// Orientation for readers coming from C / C++ / C# / Python (see epidemic.rs for
// the fuller tour of the Rust idioms):
//
//   * `&[f64]` is a *slice*: a (pointer, length) view into a contiguous f64 array
//     that the caller owns — like a `std::span<double>` / `ReadOnlySpan<double>`.
//     We never copy the input; we only read through the borrow.
//   * `x.to_radians()`, `.sin()`, `.cos()`, `.sqrt()`, `.asin()`, `.powi(2)` are
//     inherent methods on `f64` (Rust has no free-standing `sin(x)`; it is
//     `x.sin()`).
//   * `par_chunks_mut(n)` (from Rayon) splits a mutable slice into consecutive
//     non-overlapping `n`-element chunks and hands each to a worker thread — the
//     parallel analogue of an OpenMP `parallel for` over columns. Disjoint chunks
//     mean no data races, so the borrow checker allows the concurrent writes.
//   * `panic!(...)` aborts like a C++ `throw`; extendr catches it at the FFI
//     boundary and re-raises it as an R `stop()` error.
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;   // `#[extendr]`, `Robj`, R<->Rust conversions
use rayon::prelude::*;         // data-parallel iterators (par_chunks_mut)

// Mean Earth radius in kilometres — the same constant laser-core uses, so the
// distances match its reference implementation to floating-point precision.
const EARTH_RADIUS_KM: f64 = 6371.0;

// Read an R numeric/integer vector argument into an owned `Vec<f64>`.
//
// `as_real_slice()` borrows the f64 buffer of an R "double" vector; `as_integer_slice()`
// borrows an R "integer" vector (which we widen to f64). `if let Some(s) = ...` is
// Rust's match-on-`Option` shorthand (Some = present, None = absent). We copy into
// an owned Vec so no borrow on `obj` outlives this call.
fn as_f64_vec(obj: &Robj, name: &str) -> Vec<f64> {
    if let Some(s) = obj.as_real_slice() {
        s.to_vec()
    } else if let Some(s) = obj.as_integer_slice() {
        s.iter().map(|&x| x as f64).collect()
    } else {
        panic!("`{name}` must be a numeric vector");
    }
}

/// Great-circle distance matrix between geographic points.
///
/// Ports laser-core's `distance` (all-pairs case): given the latitudes and
/// longitudes of `N` points in decimal degrees, returns the symmetric `N × N`
/// matrix whose `[i, j]` entry is the haversine great-circle distance, in
/// **kilometres**, between point `i` and point `j`. The diagonal is zero.
///
/// The haversine formula (see <https://en.wikipedia.org/wiki/Haversine_formula>)
/// is evaluated with a mean Earth radius of 6371 km, matching laser-core.
///
/// @param latitude  Numeric vector of latitudes in decimal degrees, in `[-90, 90]`.
/// @param longitude Numeric vector of longitudes in decimal degrees, in `[-180, 180]`.
///   Must be the same length as `latitude`.
/// @return An `N × N` numeric matrix of pairwise distances in kilometres
///   (symmetric, zero diagonal), where `N` is the number of points.
/// @examples
/// # London and Paris, ~344 km apart:
/// d <- distances(c(51.5074, 48.8566), c(-0.1278, 2.3522))
/// d[1, 2]
/// @export
#[extendr]
fn distances(latitude: Robj, longitude: Robj) -> Robj {
    let lat = as_f64_vec(&latitude, "latitude");
    let lon = as_f64_vec(&longitude, "longitude");

    let n = lat.len();
    // Contract checks. `assert!`/`assert_eq!` panic with the message on failure;
    // extendr surfaces the panic as an R error before any work is done.
    assert_eq!(
        lon.len(), n,
        "`latitude` and `longitude` must have the same length ({} vs {})",
        n, lon.len()
    );
    assert!(
        lat.iter().all(|&v| (-90.0..=90.0).contains(&v)),
        "`latitude` values must be in [-90, 90]"
    );
    assert!(
        lon.iter().all(|&v| (-180.0..=180.0).contains(&v)),
        "`longitude` values must be in [-180, 180]"
    );

    // Pre-convert to radians once (haversine works in radians). `.collect()` runs
    // the lazy `.map(...)` into an owned Vec, like a Python list comprehension.
    let lat_rad: Vec<f64> = lat.iter().map(|d| d.to_radians()).collect();
    let lon_rad: Vec<f64> = lon.iter().map(|d| d.to_radians()).collect();

    // Allocate the flat backing buffer for the matrix. R matrices are COLUMN-MAJOR:
    // element `[i, j]` lives at flat index `i + j*n`, so column `j` is the
    // contiguous slice `out[j*n .. (j+1)*n]`. We fill one column per Rayon task.
    let mut out = vec![0.0_f64; n * n];
    out.par_chunks_mut(n).enumerate().for_each(|(j, col)| {
        let (lat_j, lon_j) = (lat_rad[j], lon_rad[j]);
        // `col[i]` is the distance from point i to this column's point j.
        for i in 0..n {
            let dlat = lat_rad[i] - lat_j;
            let dlon = lon_rad[i] - lon_j;
            // haversine: a = sin²(Δφ/2) + cos φ_i · cos φ_j · sin²(Δλ/2)
            let a = (dlat * 0.5).sin().powi(2)
                + lat_rad[i].cos() * lat_j.cos() * (dlon * 0.5).sin().powi(2);
            // central angle c = 2·asin(√a); arc length = R·c.
            col[i] = 2.0 * a.sqrt().asin() * EARTH_RADIUS_KM;
        }
    });

    // Turn the flat Vec into an R matrix by attaching a `dim` attribute (the only
    // thing that distinguishes an R matrix from a plain vector).
    let mut robj: Robj = out.into();
    robj.set_attrib("dim", [n as i32, n as i32]).expect("set dim");
    robj
}

// ════════════════════════════════════════════════════════════════════════════
// Migration-network models — ports of laser-core's `gravity`, `radiation`,
// `competing_destinations`, `stouffer`, and `row_normalizer` (laser/core/
// migration.py).
//
// Every model takes a 1-D `pops` (one population per node) and a symmetric N×N
// `distances` matrix, and returns an N×N `network` matrix whose `[i, j]` entry is
// the (un-normalized) migration weight from node i to node j. The diagonal is 0
// (no self-migration). The result is generally NOT symmetric.
//
// R matrices are COLUMN-MAJOR: element `[row, col]` lives at flat index
// `row + col*N`. So `network[i, j]` (flow i→j) is stored at `out[i + j*N]`, and
// `distances[i, j]` is read at `dist[i + j*N]` (== `dist[j + i*N]` since distances
// are symmetric).
// ════════════════════════════════════════════════════════════════════════════

// Read the `distances` argument into a column-major `Vec<f64>` and validate it is
// an N×N symmetric, non-negative matrix matching `pops`. `n` comes from `pops`.
fn read_distances(distances: &Robj, n: usize) -> Vec<f64> {
    let dist = as_f64_vec(distances, "distances");
    assert_eq!(
        dist.len(), n * n,
        "`distances` must be an {n} x {n} matrix matching `pops` ({} elements), got {}",
        n * n, dist.len()
    );
    // Non-negativity and symmetry (laser-core requires `distances == distances.T`).
    // A small relative tolerance absorbs any floating-point asymmetry.
    for j in 0..n {
        for i in 0..n {
            let d = dist[i + j * n];
            assert!(d >= 0.0, "`distances` must contain only non-negative values");
            let dt = dist[j + i * n];
            assert!(
                (d - dt).abs() <= 1e-9 * (1.0 + d.abs()),
                "`distances` must be a symmetric matrix"
            );
        }
    }
    dist
}

// Validate `pops` is non-negative; reused by every model.
fn check_pops(pops: &[f64]) {
    assert!(
        pops.iter().all(|&p| p >= 0.0),
        "`pops` must contain only non-negative values"
    );
}

// Wrap a column-major flat buffer as an R N×N matrix (attach the `dim` attribute).
fn as_square_matrix(data: Vec<f64>, n: usize) -> Robj {
    let mut robj: Robj = data.into();
    robj.set_attrib("dim", [n as i32, n as i32]).expect("set dim");
    robj
}

// Stable argsort of source row `i`'s distances (ascending). Returns the node
// indices in nearest-to-farthest order. `sort_by` is a stable sort (ties keep
// input order), matching numpy's `argsort(..., kind="stable")`.
fn argsort_row(i: usize, dist: &[f64], n: usize) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| dist[i + a * n].partial_cmp(&dist[i + b * n]).unwrap());
    idx
}

// Cumulative population of all nodes "as close or closer" to the source, given a
// distance-sorted population row and its sorted distances. Port of laser-core's
// `sum_populations_as_close_or_closer`.
//
// The base case is a plain prefix sum (cumulative sum) of the sorted populations.
// The subtlety: when several destinations are *equidistant* from the source, the
// model sums over "all k with d_ik <= d_ij", so every node in a tie group must
// see the SAME cumulative total — the value at the end of the group. We therefore
// flatten each run of equal distances to the run's final prefix-sum value.
fn sum_close_or_closer(sorted_pops: &[f64], sorted_dist: &[f64]) -> Vec<f64> {
    let n = sorted_pops.len();
    let mut cum = vec![0.0_f64; n];
    let mut acc = 0.0;
    for m in 0..n {
        acc += sorted_pops[m];
        cum[m] = acc;
    }
    // Flatten ties: scan consecutive runs of equal distance and set every entry in
    // a run to the run's last cumulative value.
    let mut start = 0;
    while start < n {
        let mut end = start;
        while end + 1 < n && sorted_dist[end + 1] == sorted_dist[start] {
            end += 1;
        }
        let v = cum[end];
        for c in cum.iter_mut().take(end + 1).skip(start) {
            *c = v;
        }
        start = end + 1;
    }
    cum
}

// Core gravity computation, returned as a column-major flat buffer (so it can be
// reused by `competing_destinations` without round-tripping through R).
fn gravity_impl(pops: &[f64], dist: &[f64], n: usize, k: f64, a: f64, b: f64, c: f64) -> Vec<f64> {
    let mut out = vec![0.0_f64; n * n];
    // `network[i, j] = k * p_i^a * p_j^b * d_ij^(-c)`; the diagonal stays 0
    // (laser-core sets the diagonal distance to 1 to avoid div-by-zero, then zeros
    // the diagonal network entry — skipping i==j is equivalent).
    for j in 0..n {
        for i in 0..n {
            if i == j {
                continue;
            }
            let d = dist[i + j * n];
            out[i + j * n] = k * pops[i].powf(a) * pops[j].powf(b) * d.powf(-c);
        }
    }
    out
}

/// Gravity migration-network model.
///
/// `network[i, j] = k * pops[i]^a * pops[j]^b / distance[i, j]^c`, with a zero
/// diagonal (no self-migration). Port of laser-core's `gravity`.
///
/// @param pops      Numeric vector of node populations (length N, non-negative).
/// @param distances Symmetric `N × N` numeric distance matrix (e.g. from
///   [distances()]).
/// @param k Scaling constant for the overall flow magnitude (non-negative).
/// @param a Exponent on the origin population.
/// @param b Exponent on the destination population.
/// @param c Exponent on the distance (larger `c` = stronger distance decay).
/// @return An `N × N` numeric matrix; `[i, j]` is the migration weight from node
///   i to node j.
/// @export
#[extendr]
fn gravity(pops: Robj, distances: Robj, k: f64, a: f64, b: f64, c: f64) -> Robj {
    let pops = as_f64_vec(&pops, "pops");
    check_pops(&pops);
    let n = pops.len();
    let dist = read_distances(&distances, n);
    assert!(k >= 0.0 && a >= 0.0 && b >= 0.0 && c >= 0.0, "`k`, `a`, `b`, `c` must be non-negative");
    as_square_matrix(gravity_impl(&pops, &dist, n, k, a, b, c), n)
}

/// Radiation migration-network model (Simini et al., Nature 2012).
///
/// For each source `i`, destinations are ranked by distance and the weight to
/// destination `j` is
/// `k * p_i * p_j / (p_i + s) / (p_i + p_j + s)`, where `s` is the total
/// population of all nodes as close or closer to `i` than `j` (excluding the home
/// population `p_i` when `include_home = FALSE`). Port of laser-core's
/// `radiation`. The diagonal is 0.
///
/// @param pops      Numeric vector of node populations (length N, non-negative).
/// @param distances Symmetric `N × N` numeric distance matrix.
/// @param k Scaling constant for the flow magnitude (non-negative).
/// @param include_home Logical; whether the home (source) population is included
///   in the "as close or closer" cumulative sum.
/// @return An `N × N` numeric migration-weight matrix (generally not symmetric).
/// @export
#[extendr]
fn radiation(pops: Robj, distances: Robj, k: f64, include_home: bool) -> Robj {
    let pops = as_f64_vec(&pops, "pops");
    check_pops(&pops);
    let n = pops.len();
    let dist = read_distances(&distances, n);
    assert!(k >= 0.0, "`k` must be non-negative");

    // Compute each source row independently in parallel, then assemble. Each row
    // is returned in ORIGINAL column order (the sort is undone by writing each
    // weight straight to `row[order[m]]`).
    let rows: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let order = argsort_row(i, &dist, n);
            let sorted_pops: Vec<f64> = order.iter().map(|&j| pops[j]).collect();
            let sorted_dist: Vec<f64> = order.iter().map(|&j| dist[i + j * n]).collect();
            let mut cum = sum_close_or_closer(&sorted_pops, &sorted_dist);
            if !include_home {
                let home = sorted_pops[0]; // closest node is the source itself (d = 0)
                for c in cum.iter_mut() {
                    *c -= home;
                }
            }
            let pi = pops[i];
            let mut row = vec![0.0_f64; n];
            for m in 0..n {
                let pj = sorted_pops[m];
                let s = cum[m];
                row[order[m]] = k * pi * pj / (pi + s) / (pi + pj + s);
            }
            row
        })
        .collect();

    let mut out = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            out[i + j * n] = if i == j { 0.0 } else { rows[i][j] };
        }
    }
    as_square_matrix(out, n)
}

/// Stouffer's intervening-opportunities migration model (Stouffer, 1940).
///
/// For each source `i` and destination `j`,
/// `network[i, j] = k * p_i^a * (p_j / s)^b`, where `s` is the cumulative
/// population as close or closer to `i` than `j` (excluding the home population
/// when `include_home = FALSE`); the nearest node (the source) gets weight 0.
/// Port of laser-core's `stouffer`. The diagonal is 0.
///
/// @param pops      Numeric vector of node populations (length N, non-negative).
/// @param distances Symmetric `N × N` numeric distance matrix.
/// @param k Scaling constant for the flow magnitude (non-negative).
/// @param a Exponent on the origin population.
/// @param b Exponent on the destination/cumulative-population ratio.
/// @param include_home Logical; whether the home population is included in the
///   cumulative sum.
/// @return An `N × N` numeric migration-weight matrix (generally not symmetric).
/// @export
#[extendr]
fn stouffer(pops: Robj, distances: Robj, k: f64, a: f64, b: f64, include_home: bool) -> Robj {
    let pops = as_f64_vec(&pops, "pops");
    check_pops(&pops);
    let n = pops.len();
    let dist = read_distances(&distances, n);
    assert!(k >= 0.0 && a >= 0.0 && b >= 0.0, "`k`, `a`, `b` must be non-negative");

    let rows: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let order = argsort_row(i, &dist, n);
            let sorted_pops: Vec<f64> = order.iter().map(|&j| pops[j]).collect();
            let sorted_dist: Vec<f64> = order.iter().map(|&j| dist[i + j * n]).collect();
            let mut cum = sum_close_or_closer(&sorted_pops, &sorted_dist);
            if !include_home {
                let home = sorted_pops[0];
                for c in cum.iter_mut() {
                    *c -= home;
                }
            }
            let pi = pops[i];
            let mut row = vec![0.0_f64; n];
            // Skip m == 0 (the source/nearest node), matching numpy's `network[i, 1:]`.
            for m in 1..n {
                let ratio = sorted_pops[m] / cum[m];
                row[order[m]] = k * pi.powf(a) * ratio.powf(b);
            }
            row
        })
        .collect();

    let mut out = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            out[i + j * n] = if i == j { 0.0 } else { rows[i][j] };
        }
    }
    as_square_matrix(out, n)
}

/// Competing-destinations migration model (Fotheringham, 1984).
///
/// Starts from the [gravity()] weights and multiplies each `[i, j]` by an
/// accessibility term `(Σ_k p_k^b / d_jk^c)^delta`, summed over competing
/// destinations `k ∉ {i, j}`. Port of laser-core's `competing_destinations`. The
/// diagonal is 0.
///
/// @param pops      Numeric vector of node populations (length N, non-negative).
/// @param distances Symmetric `N × N` numeric distance matrix.
/// @param k Scaling constant for the flow magnitude (non-negative).
/// @param a Exponent on the origin population (gravity term).
/// @param b Exponent on the destination population (gravity and competition terms).
/// @param c Exponent on distance (gravity and competition terms).
/// @param delta Exponent on the competing-destinations accessibility term.
/// @return An `N × N` numeric migration-weight matrix (generally not symmetric).
/// @export
#[extendr]
fn competing_destinations(pops: Robj, distances: Robj, k: f64, a: f64, b: f64, c: f64, delta: f64) -> Robj {
    let pops = as_f64_vec(&pops, "pops");
    check_pops(&pops);
    let n = pops.len();
    let dist = read_distances(&distances, n);
    assert!(k >= 0.0 && a >= 0.0 && b >= 0.0 && c >= 0.0, "`k`, `a`, `b`, `c` must be non-negative");

    let mut network = gravity_impl(&pops, &dist, n, k, a, b, c);

    // competition[j, kk] = p_kk^b * d_jkk^(-c), with a zero diagonal (the j==kk
    // term is dropped to avoid the zero-distance singularity).
    let mut comp = vec![0.0_f64; n * n];
    for kk in 0..n {
        for j in 0..n {
            if j == kk {
                continue;
            }
            comp[j + kk * n] = pops[kk].powf(b) * dist[j + kk * n].powf(-c);
        }
    }
    // Row sums over all competing destinations kk for each j.
    let mut row_sums = vec![0.0_f64; n];
    for j in 0..n {
        let mut s = 0.0;
        for kk in 0..n {
            s += comp[j + kk * n];
        }
        row_sums[j] = s;
    }
    // Multiply each gravity weight by (Σ_kk≠i competition[j, kk])^delta. Subtracting
    // comp[j, i] removes the k==i term; the k==j term is already 0 in `comp`.
    for j in 0..n {
        for i in 0..n {
            if i == j {
                continue;
            }
            let access = row_sums[j] - comp[j + i * n];
            network[i + j * n] *= access.powf(delta);
        }
    }
    as_square_matrix(network, n)
}

/// Cap each row sum of a network matrix at `max_rowsum`.
///
/// Rows whose total exceeds `max_rowsum` are scaled down proportionally so they
/// sum to exactly `max_rowsum`; smaller rows are left unchanged. Port of
/// laser-core's `row_normalizer`. Useful to bound the fraction of a node's force
/// of infection that is exported before passing the matrix as the transmission
/// `network`.
///
/// @param network    A square non-negative numeric matrix.
/// @param max_rowsum Maximum allowed row sum, in `[0, 1]`.
/// @return The row-capped matrix, same shape as `network`.
/// @export
#[extendr]
fn row_normalizer(network: Robj, max_rowsum: f64) -> Robj {
    let data = as_f64_vec(&network, "network");
    // Recover N from the flat length (the matrix is square).
    let n = (data.len() as f64).sqrt() as usize;
    assert_eq!(n * n, data.len(), "`network` must be a square matrix");
    assert!(data.iter().all(|&v| v >= 0.0), "`network` must contain only non-negative values");
    assert!((0.0..=1.0).contains(&max_rowsum), "`max_rowsum` must be in [0, 1]");

    let mut out = data;
    for i in 0..n {
        let mut rowsum = 0.0;
        for j in 0..n {
            rowsum += out[i + j * n];
        }
        if rowsum > max_rowsum {
            let scale = max_rowsum / rowsum;
            for j in 0..n {
                out[i + j * n] *= scale;
            }
        }
    }
    as_square_matrix(out, n)
}

extendr_module! {
    mod migration;
    fn distances;
    fn gravity;
    fn radiation;
    fn competing_destinations;
    fn stouffer;
    fn row_normalizer;
}
