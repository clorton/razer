# AGENTS.md — operating manual for coding agents in the Razer repo

You are working in **Razer**, an R package whose performance kernels are **Rust**, exposed
to R via [extendr](https://extendr.github.io/). This file is the quick, actionable contract
for making changes here without breaking things. It is deliberately short; the deeper
references are:

- **`CLAUDE.md`** — the full modeling convention (the per-tick ordering, the eight-model
  menagerie, `run_model`) and the commenting/build/git conventions. Treat it as the source
  of truth for anything not spelled out here.
- **The "Using and extending Razer" article** (`vignettes/articles/extending_razer.Rmd`,
  rendered at <https://clorton.github.io/razer/articles/extending_razer.html>) — the kernel
  ABI in full, the R-vs-Rust decision, and an LLM prompt template for writing a kernel.
- **`writeup.md`** — the design retrospective vs. Python LASER (`laser-core`/`laser-generic`).

## What Razer is (one paragraph)

Agent-based spatial disease models. Each agent property is a Rust-owned, dtype-tagged array
(a `Column`, e.g. `u8` `state`, `u16` `timer`, `u16` `nodeid`) held in R as an opaque handle;
the per-tick kernels (Rust + Rayon) mutate them in place and **return per-node counts** that
the model applies to a per-node census. `run_model()` wires the SI/SEI/SIS/SEIS/SIR/SEIR/
SIRS/SEIRS menagerie in the correct order; users extend via R callbacks or new Rust kernels.

## Repo map

```
R/                      R wrappers + helpers (run_model.R, bincount.R, calc_capacity.R, …)
R/extendr-wrappers.R    AUTO-GENERATED — never edit (devtools::document() overwrites it)
src/rust/src/           the Rust kernels:
  lib.rs                module registration (extendr_module!)
  column.rs             the Column typed-array store
  steps.rs              EXEMPLAR step kernels + the generic step_timer_expire(_set)
  transmission.rs       calc_foi + transmission/transmission_si
  rng.rs                seedable, thread-count-independent RNG
  vitals.rs mortality.rs births.rs migration.rs pyramid.rs kmestimator.rs bincount.rs …
man/                    AUTO-GENERATED Rd help — never edit by hand
examples/               runnable scripts (+ examples/data, examples/output)
vignettes/articles/     pkgdown teaching articles (website-only, .Rbuildignore'd)
tests/testthat/         tests
NEWS.md                 the changelog (NOT CHANGELOG.md — see below)
```

## Environment & build

- **`cargo` is at `~/.cargo/bin`, not on the default PATH.** Prepend it for any cargo /
  document / install step:
  ```bash
  export PATH="$HOME/.cargo/bin:$PATH"
  ```
- **After ANY change to Rust signatures, run `devtools::document()`** — it recompiles the
  Rust, regenerates `R/extendr-wrappers.R`, and writes `NAMESPACE` + `man/*.Rd`:
  ```bash
  Rscript -e 'devtools::document()'
  ```
  Do **not** use `rextendr::document()` (deprecated). Do **not** hand-edit
  `R/extendr-wrappers.R` or `man/*.Rd`.
- Fast compile-only check while iterating: `cd src/rust && cargo check`.
- Run tests: `Rscript -e 'devtools::test()'`.
- **`load_all()`/`test()` do not rebuild a stale *installed* package.** When verifying an
  example end to end, install first (`R CMD INSTALL .` or `devtools::install()`) and then run
  the `examples/*.R` script — that is the only way to catch install-time/wrapper issues.
- Rust **panics** (`assert!`, `panic!`, OOB index) surface in R as `stop()` errors — so
  `assert!` is the idiomatic input check.

## Writing a Rust kernel — the hard rules

Read `src/rust/src/steps.rs` first; mirror it. Full detail + an LLM prompt template are in
the extending article (§5–§7). The contract:

1. `#[extendr] fn`; args = the agent `Column`s + `count: i32` + `n_nodes: i32` +
   any `&Distribution`s. Take `Column`s by `&` (read) / `&mut` (write); slice to `[..count]`.
2. Storage types: `state` `u8`, `timer` `u16`, `nodeid` `u16`. Use the matching accessor —
   `as_u8_mut` / `as_u16_mut` / `as_u16` (also `as_u32` / `as_i32_mut` / `as_f64`). They
   panic on a dtype mismatch.
3. **Return per-node counts** — a `Vec<i32>` or a named `List` of `integer[n_nodes]`. Never
   write the node census inside a kernel; the R caller applies counts with `move_count`.
4. Parallelize with `par_chunks_mut(rng::RNG_CHUNK)` + a private per-node `Vec<i64>` reduced
   by summing (no shared writes).
5. **RNG:** `let base = rng::next_call_base();` once before the parallel fan-out, then
   `rng::chunk_rng(base, chunk_index)` inside each chunk. This is what keeps a seeded run
   reproducible **independent of thread count** — never use `thread_rng()` or seed per-thread.
   Draw a uniform with `rng.gen::<f64>()` (`use rand::Rng;`).
6. Guard `u16` timer decrements: `if t[j] > 0 { t[j] -= 1; }` before testing `== 0`.
7. Validate inputs with `assert!`.
8. **Register** the function in the file's `extendr_module! { … }`, and ensure the module is
   listed in `src/rust/src/lib.rs` (both `mod x;` and `use x;`). New file ⇒ add both.
9. Ship a test (see Testing) and run `devtools::document()` then `devtools::test()`.

**Do you even need Rust?** Only if the change requires a *novel per-agent computation every
tick* that is not a composition of existing kernels (`step_timer_expire(_set)`, `move_count`,
`births`/`mortality`/`import_infections`, the `bincount` family, `squash`). Vaccination,
waning, quarantine, vital dynamics, imports, age-targeted reports, and time/space-varying
drivers are all **pure-R callbacks** — no Rust. If you want to write a per-agent `for` loop in
R, or round-trip a whole agent `Column` through R every tick, that is the signal to write a
kernel instead. See the extending article §3.

## House conventions (do not violate)

- **Comment for a C/C++/C#/Python reader who does NOT know Rust or R.** In `.rs`, explain
  Rust-specific idioms (`&`/`&mut`, slices, `match`, closures, iterators, `Option`/`Result`/
  `?`, turbofish, `unsafe`, what `#[extendr]` generates). In `.R`, explain R idioms (`<-`, S3
  dispatch, `.Call`, environments, `[[` vs `$`, vectorization). Exemplars: `steps.rs` /
  `transmission.rs` (`.rs`), `R/pyramid.R` (`.R`). Do NOT comment auto-generated files
  (`R/extendr-wrappers.R`, `tools/*.R`).
- **Changelog is `NEWS.md`, not `CHANGELOG.md`** (a repo-specific exception). Record accepted
  edits there. **Keep the `# razer (development version)` H1 lowercase** — pkgdown's news
  parser matches it against the lowercase `Package:` name; changing it re-breaks the
  Changelog page.
- **The project is branded "Razer" in prose, but the package/code name is lowercase `razer`.**
  Never change `library(razer)`, `razer::`, the `clorton.github.io/razer` / `github.com/
  clorton/razer` URLs, file paths, or the `Package:` field.
- **Node ids are 0-based on the Rust side**; convert to 1-based only at the R boundary.
- Distribution constructors use a **`dist_` prefix** (`dist_normal`, `dist_gamma`, …) to avoid
  masking base/stats functions.
- **Articles** live in `vignettes/articles/` (website-only, `.Rbuildignore`d, not run by
  `R CMD check`). In docs, link example scripts with **absolute GitHub URLs**
  (`https://github.com/clorton/razer/blob/main/examples/…`) and cross-link articles with
  **relative** `name.html`.
- Examples are **device-aware** (write PNGs to `examples/output/` under `Rscript`, draw to the
  Plots pane when `source()`d in RStudio) and should be **liberal with informative plots**.
- Internal root docs (`CLAUDE.md`, `AGENTS.md`, `writeup.md`) are `.Rbuildignore`d and dropped
  from the published pkgdown site in `.github/workflows/pkgdown.yaml`. If you add another,
  do the same.

## Testing

- Write tests in **given–when–then** style with a docstring stating the purpose. New behavior
  ⇒ a test in the same change (tests are not optional follow-up).
- For a new kernel, the standard check is **parallel-vs-serial census**: run it on ~1e6 agents
  and assert the returned per-node counts equal a direct R tally of the resulting agent
  states; add a reproducibility check (run twice under `set_seed()`, assert identical).
- **Run the tests before declaring done.** Report failures honestly with output.

## Git & finishing

- **Fast-forward-only integration** — merge into `main` with `git merge --ff-only`; never a
  merge commit (rebase a diverged branch first). Enforced via repo git config.
- Commit messages use **emoji prefixes**: ✨ feature, 🧪 tests, 📝 docs, ⚡️ perf, 🐛 fix,
  🚧 wip, 🔧 infra, 🚚 rename, 🦺 validation (concatenate if a commit does several). If you
  are an AI agent, add a `Co-Authored-By:` trailer identifying the assistant.
- **Stage explicit paths, not `git add -A`** (avoids sweeping in build artifacts or stray
  files). Commit only when asked; **do not push** — the human pushes.

### Before you finish, confirm
- [ ] If Rust signatures changed: `devtools::document()` ran (wrappers + `man/` regenerated).
- [ ] `devtools::test()` passes (and examples installed+run if behavior changed).
- [ ] New/changed behavior has tests; `NEWS.md` updated.
- [ ] If docs changed: `pkgdown::check_pkgdown()` is clean (articles/reference indexed).
- [ ] Comments suit the C/C#/Python reader; no auto-generated files hand-edited.
