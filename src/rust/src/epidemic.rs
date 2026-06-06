// ════════════════════════════════════════════════════════════════════════════
// epidemic.rs — disease-compartment state codes.
//
// The per-tick dynamics live in the Column-based kernels (sir.rs, measles.rs,
// vitals.rs, mortality.rs, births.rs); those kernels and the R model scripts share
// the compartment codes defined here. `laser_states()` exposes the same codes to R.
//
// Orientation for readers coming from C / C++ / C# / Python:
//   * `pub const` is a compile-time constant (like a C `#define` / C# `const`).
//   * `#[extendr]` is a procedural macro (an attribute, like a C# attribute or a
//     Python decorator) that generates the C-ABI shim and the R wrapper so the
//     annotated Rust fn is callable from R. `Robj` is an owning handle to an R
//     object (SEXP).
// ════════════════════════════════════════════════════════════════════════════

use extendr_api::prelude::*;

// ── State constants ───────────────────────────────────────────────────────────
// -1 for D keeps the alive compartments in the non-negative range, so a
// `state != D` test is an unambiguous "is alive?" check. The u8 `state` Column
// stores D (= -1) as 255.

pub const STATE_S: i32 =  0;
pub const STATE_E: i32 =  1;
pub const STATE_I: i32 =  2;
pub const STATE_R: i32 =  3;
pub const STATE_M: i32 =  4;   // maternal immunity (newborns); wanes to S
pub const STATE_D: i32 = -1;

/// Named integer vector of epidemic compartment state codes.
///
/// Returns `c(S=0L, E=1L, I=2L, R=3L, M=4L, D=-1L)`. Use these constants to set
/// and test the `state` property on a people frame. `M` is maternal immunity
/// (newborns protected by maternal antibodies, waning to `S`); `D` is deceased.
///
/// @return Named integer vector with elements S, E, I, R, M, D.
/// @export
#[extendr]
fn laser_states() -> Robj {
    // `vec![...]` is the Vec (growable array, ~ std::vector / List<T>) literal.
    let vals: Vec<i32> = vec![STATE_S, STATE_E, STATE_I, STATE_R, STATE_M, STATE_D];
    // `.iter()` borrows each element; `collect_robj()` materializes the iterator
    // into an R integer vector. `mut` marks the binding as reassignable/mutable
    // (Rust bindings are immutable by default, the opposite of C/C#).
    let mut robj = vals.iter().collect_robj();
    // Build the names vector. `["S",..]` is a fixed-size array; `.map(...)` is the
    // lazy transform (LINQ `Select` / Python `map`); `.collect()` runs it into a
    // `Vec<String>`. `s.to_string()` copies the `&str` literal into an owned String.
    let names: Vec<String> = ["S", "E", "I", "R", "M", "D"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let names_robj = names.iter().map(|s| s.as_str()).collect_robj();
    // `set_attrib` returns a `Result` (Ok/Err, Rust's checked-error type — there
    // are no exceptions). `.expect(msg)` unwraps Ok or panics with `msg` on Err.
    robj.set_attrib("names", names_robj).expect("set names");
    robj            // trailing expression with no `;` is the return value
}

extendr_module! {
    mod epidemic;
    fn laser_states;
}
