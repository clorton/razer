# razer — project notes for Claude

## Modeling convention: the Column kernels and per-tick ordering

Models are composed from the **Column-based** per-tick kernels over `Column` buffers
(`allocate_scalar` / `allocate_vector`). Agent `state` is a `u8`; the per-agent `timer`
is a **`u16`** everywhere (maternal / immunity periods exceed a `u8`'s 255).

**Kernels mutate the per-agent arrays and RETURN per-node counts; the model applies the
deltas to the census it maintains.** So no kernel takes census/flow buffers, and a model
allocates only the compartments it has. The flow is: `carry_forward(_states)` copies each
census column `t → t+1`, then the model applies each kernel's returned counts to column
`t+1` with `move_count(from, to, counts, t)` (`from`/`to` may be `NULL` for one-sided
moves — a death decrements only, a birth increments only).

- `calc_foi(I, N, beta, seasonality, network, foi, t)` — per-node force of infection
  (writes the `foi` report Column; the one kernel that still reads/writes Columns).
- `transmission(state, timer, nodeid, count, foi, t, to_state, duration) → counts` —
  S→`to_state` (E or I), sets the u16 timer; returns new infections per node.
- `transmission_si(state, nodeid, count, foi, t) → counts` — S→I **absorbing** (no
  timer); the SI model.
- `step_si(…, inf_dur) → list(waned, onset)` — M→S, E→I.
- `step_sir(…, inf_dur, absorbing_state) → list(waned, onset, cleared)` — adds
  I→`absorbing_state` (S or R).
- `step_sirs(…, inf_dur, imm_dur) → list(waned, onset, recovered, waned_r)` — adds I→R
  (sets an immunity timer) and R→S.
- `births(…) → list(count, born)`, `mortality(…) → list(m, s, e, i, r)`,
  `constant_pop_vitals_sir(…)` (constant-pop convenience; still writes its census),
  `import_infections(…)`.

**The eight-model menagerie** — every model is a transmission + a step kernel:

| Model | transmission | step kernel |
|---|---|---|
| SI | `transmission_si` (S→I absorbing) | `step_si` |
| SEI | `transmission` (S→E) | `step_si` |
| SIS / SIR | `transmission` (S→I) | `step_sir`, absorbing = S / R |
| SEIS / SEIR | `transmission` (S→E) | `step_sir`, absorbing = S / R |
| SIRS | `transmission` (S→I) | `step_sirs` |
| SEIRS | `transmission` (S→E) | `step_sirs` |

Each step kernel is a **single pass branching on the agent's entry state** (and leads
with M→S, so any model can add a maternal compartment), so a just-entered timed state is
never decremented again the same tick — the residency artifact (`duration − 1`) is
avoided structurally.

**One ordering for every model — always realize the full `beta · D`, never `beta · (D −
1)`.** Each tick runs, for every model:

> `carry_forward(_states)` → **`step_*`** (apply counts) → **`calc_foi`** → **`transmission`** (apply counts)

with `calc_foi` placed **immediately before `transmission`** and **no step kernel between
them**. This works because `calc_foi` reads the **settled start-of-interval** infectious
census `I[t]` (not the working column `t+1` that this tick's step and transmission build),
so its value does NOT depend on where the step kernel runs. An agent enters `I` one census
column after it is infected (via either the step kernel's E→I or `transmission`'s S→I) and
recovers `D` columns later, so it contributes to the FOI on exactly the `D` columns it
occupies — `R0 = beta · D` for both direct S→I and E-entry, with no special-casing and no
residency artifact. (Historical note: `calc_foi` used to read `I[t+1]`, which forced a
different ordering per family; reading `I[t]` removed that footgun.)

Validated by `examples/sir_attack_fraction.R` / `seir_attack_fraction.R` (both match the
Kermack–McKendrick final size with `R0 = beta · D`). See `examples/engwal_measles.R` for
a full M-S-E-I-R model and `simple_sir.R` / `endemic_sir.R` for spatial Column SIR.

**Higher-level helpers.** `run_model(scenario, model, …)` wires the whole menagerie in the
correct order and returns a **`model` environment** bundling `$people`, `$nodes`,
`$network`, and `$carry`; it seeds only the states the model has (any `E`/`I`/`R` column
the scenario supplies), records the per-node flows `incidence` / `onset` / `recovery` /
`waning` as applicable, and takes optional `init(model)` / `step_enter(model)` /
`step_exit(model)` callbacks for user extension. `squash(people)` compacts dead agents
out; `calc_capacity(...)` bounds cumulative-births capacity (no reclaim) while
`calc_capacity_cdr(...)` bounds peak-living capacity for squash-reclaimed long runs. The
binning family is `bincount()` / `bincount_wt()` (weighted) / `bincount_where()` (predicate-
filtered, count-aware) / `bincount_where_wt()` (weighted + filtered).

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
