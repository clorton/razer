# razer — project notes for Claude

## Modeling convention: downstream-first transition ordering

When composing a model from the per-tick step kernels, **call the transitions in
downstream-first order** — i.e. move agents *out* of each timed compartment before
moving agents *into* it, walking the disease-progression chain backwards:

| Model  | Per-tick call order                                                        |
|--------|----------------------------------------------------------------------------|
| SIR    | `step_infectious_ir` → `step_transmission_si`                              |
| SEIR   | `step_infectious_ir` → `step_exposed_ei` → `step_transmission_se`         |
| SEIRS  | `step_recovered_rs` → `step_infectious_ir` → `step_exposed_ei` → `step_transmission_se` |
| SIS    | `step_infectious_is` → `step_transmission_si`                             |

**Why.** Each timed transition sets a fresh duration timer when an agent *enters*
a state, and the step for the *next* transition decrements that timer. If the
inflow step runs before the outflow step within the same tick, a just-arrived
agent is decremented immediately and loses one tick of residency — its effective
exposed/infectious period is `duration − 1` instead of `duration`. Running
downstream-first (outflow before inflow) prevents this: a newly exposed or newly
infectious agent is not touched again until the next tick.

This is validated by the attack-fraction examples: with downstream-first ordering
the **SEIR** model matches the Kermack–McKendrick final size `A = 1 − exp(−R0·A)`
with `R0 = beta · D` (full infectious period `D`).

**No `beta · (D − 1)` models — always realize the full `beta · D`.** With direct
`S→I` transmission a newly infectious agent enters `I` *after* the tick's
force-of-infection tally, so it is never counted on its entry tick. If recoveries are
*also* excluded from the tally, the agent contributes on only `D − 1` ticks and the
effective reproduction number collapses to `R0 = beta · (D − 1)`. **This is a bug, not
an accepted convention; no model in this package should exhibit it.**

The fix depends on which kernels you use:

- **Column kernels (`calc_foi` / `sir_step` / `transmission`) — the standard path.**
  Run `calc_foi` *before* the recovery step (`sir_step`), so an agent is still counted
  on its recovery tick. The directly-infected agent loses its entry tick but gains its
  recovery tick → the full `D`. Per-tick order: `carry_forward` → **`calc_foi`** →
  `sir_step` → `transmission` → (vitals). (SEIR-style models — where agents enter `I`
  via a step run *before* the tally, e.g. `measles_step`'s E→I — instead run `calc_foi`
  *after* that step: the new infectious are counted on entry and recoveries excluded,
  also `beta · D`.) `calc_foi`'s docstring spells out both orderings.
- **Monolithic LaserFrame `step_transmission_si`.** Because it computes its infectious
  tally and applies the S→I infection in the *same* call, neither reordering can
  recover the full `D` — it structurally yields `beta · (D − 1)`. **Do not use it for
  new SIR models; prefer the Column kernels** (or SEIR-style entry). `run_sir()` and
  `examples/sir_attack_fraction.R` still ride this legacy kernel — see their notes.

See `examples/simple_sir.R` / `endemic_sir.R` (Column SIR, `beta · D`) and
`examples/seir_attack_fraction.R` (SEIR, `beta · D`).

Steps that transition into untimed states (`step_infectious_is`,
`step_recovered_rs` into S; `step_mortality_cdr` into D; `step_births_cbr`) do not
themselves suffer the artifact, but `step_recovered_rs` is still an *outflow* of R
and so leads in SEIRS.

## Commenting conventions

The audience for this codebase is a programmer fluent in **C / C++ / C# /
Python** but **not** in Rust or R. Comment for that reader.

- **Always comment Rust code (`.rs`) for this audience.** Explain Rust-specific
  idioms a C/C#/Python programmer would not immediately recognize — ownership and
  borrowing (`&`, `&mut`), slices vs. raw pointers, `match`/`enum` destructuring,
  closures, lazy and parallel iterators (`iter`/`map`/`fold`/`par_iter`), `Option`
  / `Result` / `?`, the turbofish (`::<T>`), `unsafe` blocks and their SAFETY
  contracts, and what attribute macros like `#[extendr]` generate. Don't explain
  general programming; explain what differs from the languages above.
  `src/rust/src/epidemic.rs` is the reference exemplar for the level and style.
- **Always comment R code (`.R`) for this audience.** Explain R idioms that trip
  up non-R programmers — `<-` assignment, S3 dispatch (`generic.class`,
  `` `$<-` `` replacement methods), `.Call` into compiled code, environments and
  closure rebinding, vectors / `c()` / `%in%`, `[[` vs `$`, `invisible()`, `NULL`,
  column-major matrices, and vectorized operations. `R/laser_frame.R` is the
  reference exemplar. Note: `R/extendr-wrappers.R` and `tools/*.R` are
  auto-generated — do not comment them (regeneration overwrites edits).

## Build / dev notes

- `cargo` is at `~/.cargo/bin` (not on the default PATH); prepend
  `export PATH="$HOME/.cargo/bin:$PATH"` for cargo / rextendr / R CMD INSTALL.
- After changing Rust signatures, regenerate the R wrappers, NAMESPACE, and man
  pages with `Rscript -e 'devtools::document()'` from the package root. This
  recompiles the Rust, regenerates `R/extendr-wrappers.R` (via the `document`
  cargo bin), and writes NAMESPACE + `man/*.Rd`. (Do **not** use
  `rextendr::document()` — it was deprecated in rextendr 0.4.0 and merely wraps
  `devtools::document()` after a config check.)
- Distribution constructors use a `dist_` prefix (`dist_normal`, `dist_gamma`, …)
  to avoid masking base/stats functions such as `base::gamma` / `stats::poisson`.
- Run tests with `Rscript -e 'devtools::test()'`.

## Git workflow

- **Integrate with fast-forward only — never create a merge commit.** Merge a
  feature or worktree branch into `main` with `git merge --ff-only <branch>`;
  if the branch has diverged, rebase it onto `main` first, then fast-forward.
  Enforced repo-locally via `merge.ff = only` and `pull.ff = only` in git config.
