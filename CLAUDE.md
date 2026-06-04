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

**SIR caveat.** For direct `S→I` transmission the offset persists regardless of
order: an agent enters `I` *inside* `step_transmission_si`, after that step's
infectious tally is computed, so it first contributes to the force of infection on
the next tick. The SIR effective reproduction number is therefore
`R0 = beta · (D − 1)`, not `beta · D`. (SEIR avoids this because agents enter `I`
via `step_exposed_ei`, a separate step that runs before the transmission tally.)
See `examples/sir_attack_fraction.R` and `examples/seir_attack_fraction.R`.

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
  pages with `Rscript -e 'rextendr::document()'` from the package root.
- Distribution constructors use a `dist_` prefix (`dist_normal`, `dist_gamma`, …)
  to avoid masking base/stats functions such as `base::gamma` / `stats::poisson`.
- Run tests with `Rscript -e 'devtools::test()'`.

## Git workflow

- **Integrate with fast-forward only — never create a merge commit.** Merge a
  feature or worktree branch into `main` with `git merge --ff-only <branch>`;
  if the branch has diverged, rebase it onto `main` first, then fast-forward.
  Enforced repo-locally via `merge.ff = only` and `pull.ff = only` in git config.
