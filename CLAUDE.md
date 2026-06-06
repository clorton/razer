# razer — project notes for Claude

## Modeling convention: the Column kernels and per-tick ordering

Models are composed from the **Column-based** per-tick kernels operating on
`Column` buffers (allocated with `allocate_scalar` / `allocate_vector`):

- `carry_forward(col, t)` / `carry_forward_states(list, t, total)` — copy each census
  counter from column `t` to `t+1` (and optionally total them into `N`); every kernel
  then applies only its *delta* to column `t+1`, so `count[t+1] = count[t] ± delta`.
- `calc_foi(I, N, beta, seasonality, network, foi, t)` — per-node force of infection.
- `sir_step(...)` — I→R recovery (u8 timer); `measles_step(...)` — M→S, E→I, I→R in one
  u16-timer pass (use it for SEIR/SEIRS too, M left empty).
- `transmission(...)` — S→I writing a u8 timer; `transmission_u16(...)` — S→E (or S→I)
  writing a u16 timer.
- `constant_pop_vitals_sir(...)`, `mortality(...)`, `births(...)` — vital dynamics.

**No `beta · (D − 1)` models — always realize the full `beta · D`.** An agent must
contribute to the force of infection on exactly the `D` ticks it is infectious. The
trap is the placement of `calc_foi` relative to the entry into / exit from `I`:

- **Direct S→I (SIR), Column kernels.** A directly-infected agent enters `I` *after*
  the tick's FOI tally, so it is never counted on its entry tick. Run **`calc_foi`
  before `sir_step`** so it is still counted on its *recovery* tick — losing the entry
  tick but gaining the recovery tick nets the full `D`. Per-tick order:
  `carry_forward` → **`calc_foi`** → `sir_step` → `transmission` → (vitals).
- **SEIR-style entry (`measles_step`'s E→I).** Agents enter `I` via `measles_step`,
  which runs *before* the tally, so run **`calc_foi` after `measles_step`**: new
  infectious are counted on entry, recoveries excluded — also `beta · D`. Per-tick
  order: `carry_forward` → `measles_step` → **`calc_foi`** → `transmission_u16` →
  (births, mortality). `calc_foi`'s docstring spells out both orderings.

This is validated by `examples/sir_attack_fraction.R` and `seir_attack_fraction.R`:
both match the Kermack–McKendrick final size `A = 1 − exp(−R0·A)` with `R0 = beta · D`.

**Combined timed transitions avoid the residency artifact structurally.**
`measles_step` advances M→S, E→I, and I→R in a *single* pass, branching on each agent's
state at the start of the pass, so an agent that does E→I this tick is not also recovered
this tick. (A newly entered timed state must not be decremented again in the same tick,
or its residency would be `duration − 1`.) See `examples/engwal_measles.R` (full M-S-E-I-R
model), `simple_sir.R` / `endemic_sir.R` (Column SIR, `beta · D`).

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
  `src/rust/src/sir.rs` is the reference exemplar for the level and style.
- **Always comment R code (`.R`) for this audience.** Explain R idioms that trip
  up non-R programmers — `<-` assignment, S3 dispatch (`generic.class`,
  `` `$<-` `` replacement methods), `.Call` into compiled code, environments and
  closure rebinding, vectors / `c()` / `%in%`, `[[` vs `$`, `invisible()`, `NULL`,
  column-major matrices, and vectorized operations. `R/pyramid.R` is the
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
