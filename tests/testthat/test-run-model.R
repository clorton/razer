# Tests for run_model(): the high-level runner for the closed-population menagerie. It
# builds people/nodes, seeds, runs the correct per-tick loop, and applies the census
# deltas. Closed population (no births/deaths), so the living states sum to N at
# every tick. Written given-when-then.

# Run each model with the durations it needs; returns the result. We pass every period for
# convenience, so models lacking E/waning warn that those are ignored (see the dedicated
# warning test) — suppress that expected noise here.
run_one <- function(model, nticks = 60L, n = 20000L, seed = 1L) {
  scenario <- data.frame(population = n, I = 50L)
  suppressWarnings(
    run_model(scenario, model, nticks = nticks, r0 = 2.5,
              infectious_period = 8, incubation_period = 5, immunity_period = 60, seed = seed))
}

# Sum the present census states at each recorded tick (an n_ticks-length vector).
living_by_tick <- function(nodes, model) {
  cols <- list(nodes$S, nodes$I)
  if (grepl("E", model)) cols <- c(cols, list(nodes$E))
  if (grepl("R", model)) cols <- c(cols, list(nodes$R))
  Reduce(`+`, lapply(cols, function(c) rowSums(c$values())))
}

test_that("run_model runs every model and conserves the closed population", {
  # Given each of the eight models
  # When run for 60 ticks at R0 = 2.5
  # Then the present states sum to N at every recorded tick (no births/deaths), and
  #      some transmission has occurred (final S < initial S). Failure would mean a
  #      mis-applied census delta (the desync the runner exists to prevent) or a model
  #      mis-wired.
  for (model in c("SI", "SEI", "SIS", "SEIS", "SIR", "SEIR", "SIRS", "SEIRS")) {
    res <- run_one(model)
    living <- living_by_tick(res$nodes, model)
    expect_true(all(living == 20000L), info = paste(model, "conserves population"))
    s <- rowSums(res$nodes$S$values())
    expect_lt(s[length(s)], s[1L])             # susceptibles fell -> transmission happened
  }
})

test_that("run_model is reproducible under a seed", {
  # Given the same seed
  # When an SEIR model is run twice
  # Then the infectious trajectories are identical. Failure would mean the runner does not
  #      thread the seed through deterministically.
  a <- run_one("SEIR", seed = 7L)
  b <- run_one("SEIR", seed = 7L)
  expect_identical(a$nodes$I$values(), b$nodes$I$values())
})

test_that("SIS returns to susceptible (no R state); SIR accumulates R", {
  # SIS: I clears back to S, so there is no R census and S recovers after the peak.
  sis <- run_one("SIS")
  expect_null(sis$nodes$R)
  # SIR: recovereds accumulate monotonically (lifelong immunity, closed population).
  sir <- run_one("SIR")
  r <- rowSums(sir$nodes$R$values())
  expect_true(all(diff(r) >= 0))
  expect_gt(r[length(r)], 0)
})

test_that("run_model validates its inputs", {
  scen <- data.frame(population = 1000L, I = 10L)
  expect_error(run_model(scen, "SXYZ", nticks = 10L, r0 = 2, infectious_period = 5), "must be one of")
  expect_error(run_model(scen, "SEIR", nticks = 10L, r0 = 2, infectious_period = 5), "incubation_period")
  expect_error(run_model(scen, "SIRS", nticks = 10L, r0 = 2, infectious_period = 5), "immunity_period")
  expect_error(run_model(scen, "SIR", nticks = 1L, r0 = 2, infectious_period = 5), "nticks")
  expect_error(run_model(scen, "SIR", nticks = 10L, r0 = 2, infectious_period = 5, init = 7),
               "init.*function")
  expect_error(run_model(scen, "SIR", nticks = 10L, r0 = -2, infectious_period = 5), "r0")
  expect_error(run_model(scen, "SIR", nticks = 10L, r0 = NA, infectious_period = 5), "r0")
  expect_error(run_model(data.frame(population = c(100, NA)), "SIR", nticks = 10L, r0 = 2,
                         infectious_period = 5), "population.*NA|NA.*population|finite")
  expect_error(run_model(data.frame(population = 100.5), "SIR", nticks = 10L, r0 = 2,
                         infectious_period = 5), "whole")
  expect_error(run_model(data.frame(population = 1000L, I = NA), "SIR", nticks = 10L, r0 = 2,
                         infectious_period = 5), "I.*finite|finite")
})

test_that("run_model returns a model environment bundling people/nodes/network/carry", {
  # Given a run
  # When we inspect the result
  # Then it is an environment exposing $people, $nodes, $network, and $carry — the bundle
  #      the callbacks receive. Failure would mean the model wasn't packaged as intended.
  m <- run_one("SIR")
  expect_true(is.environment(m))
  expect_true(all(c("people", "nodes", "network", "carry") %in% ls(m)))
  expect_true(is.environment(m$people) && is.environment(m$nodes))
  expect_true(inherits(m$network, "Column") && is.list(m$carry))   # network as a 2-D Column
  expect_equal(m$network$length(), m$nodes$count^2)
})

test_that("run_model records the per-model node flows", {
  # Given each model
  # When run
  # Then exactly the right flows exist: incidence (all), onset (E models), recovery (I-exit
  #      models), waning (waning models); and incidence is non-trivial. Failure would mean a
  #      flow is missing, spurious, or not recorded.
  expect_setequal(intersect(c("incidence","onset","recovery","waning"), ls(run_one("SI")$nodes)),
                  "incidence")
  expect_setequal(intersect(c("incidence","onset","recovery","waning"), ls(run_one("SEI")$nodes)),
                  c("incidence","onset"))
  expect_setequal(intersect(c("incidence","onset","recovery","waning"), ls(run_one("SIR")$nodes)),
                  c("incidence","recovery"))
  expect_setequal(intersect(c("incidence","onset","recovery","waning"), ls(run_one("SEIR")$nodes)),
                  c("incidence","onset","recovery"))
  expect_setequal(intersect(c("incidence","onset","recovery","waning"), ls(run_one("SIRS")$nodes)),
                  c("incidence","recovery","waning"))
  expect_setequal(intersect(c("incidence","onset","recovery","waning"), ls(run_one("SEIRS")$nodes)),
                  c("incidence","onset","recovery","waning"))
  expect_gt(sum(run_one("SEIR")$nodes$incidence$values()), 0)
})

test_that("flow accounting matches the census changes (SEIR)", {
  # Given an SEIR run
  # When we compare flows to census deltas
  # Then incidence equals the rise in (E entries) and onset equals E->I transitions: each
  #      tick, S[t+1]-S[t] = -incidence (the only way to leave S), and the running R rise
  #      equals cumulative recovery. Failure would mean a flow and its delta disagree.
  m <- run_one("SEIR")
  S <- rowSums(m$nodes$S$values()); R <- rowSums(m$nodes$R$values())
  inc <- rowSums(m$nodes$incidence$values()); rec <- rowSums(m$nodes$recovery$values())
  expect_equal(-diff(S), inc)                       # S only ever leaves via new infection
  expect_equal(diff(R), rec)                         # R only ever grows via recovery
})

test_that("run_model fires the init and per-tick callbacks in order", {
  # Given init and all three per-tick callbacks (step_enter, step_update, step_exit)
  # When an SIR model is run for nticks
  # Then init runs once, each per-tick callback fires once per interval (nticks - 1), and
  #      within a tick they fire in the order enter -> update -> exit (so update sits
  #      between the step kernel and calc_foi). Failure would mean a callback is not invoked
  #      at the documented point or in the documented order.
  scen <- data.frame(population = 5000L, I = 20L)
  counters <- new.env()
  counters$enter <- 0L; counters$update <- 0L; counters$exit <- 0L
  counters$order <- character(0)
  m <- run_model(scen, "SIR", nticks = 12L, r0 = 2, infectious_period = 6,
                 init = function(model) model$nodes$inited <- TRUE,
                 step_enter  = function(model) { counters$enter  <- counters$enter  + 1L
                                                 counters$order <- c(counters$order, "enter") },
                 step_update = function(model) { counters$update <- counters$update + 1L
                                                 counters$order <- c(counters$order, "update") },
                 step_exit   = function(model) { counters$exit   <- counters$exit   + 1L
                                                 counters$order <- c(counters$order, "exit") })
  expect_true(isTRUE(m$nodes$inited))
  expect_equal(counters$enter, 11L)
  expect_equal(counters$update, 11L)
  expect_equal(counters$exit, 11L)
  expect_equal(counters$order[1:3], c("enter", "update", "exit"))   # within-tick order
})

test_that("run_model seeds only states the model has (item 3)", {
  # Given a scenario carrying an E column
  # When run as SEIR (has E) vs SIR (no E)
  # Then SEIR seeds the E state from the column, while SIR ignores it (no E census,
  #      and those agents stay susceptible). Failure would mean states absent from the
  #      model get initialized, or present ones are skipped.
  scen <- data.frame(population = 1000L, I = 10L, E = 30L, R = 5L)
  seir <- run_model(scen, "SEIR", nticks = 5L, r0 = 2, infectious_period = 6,
                    incubation_period = 4, seed = 1L)
  expect_equal(seir$nodes$E$values()[1L, 1L], 30L)   # E seeded
  expect_equal(seir$nodes$I$values()[1L, 1L], 10L)
  expect_equal(seir$nodes$R$values()[1L, 1L], 5L)
  expect_equal(seir$nodes$S$values()[1L, 1L], 1000L - 10L - 30L - 5L)

  sir <- suppressWarnings(                            # SIR ignores the E column (warns; see below)
    run_model(scen, "SIR", nticks = 5L, r0 = 2, infectious_period = 6, seed = 1L))
  expect_null(sir$nodes$E)                            # SIR has no E census
  expect_equal(sir$nodes$S$values()[1L, 1L], 1000L - 10L - 5L)  # the 30 "E" stay S
})

test_that("run_model warns about inputs the chosen model does not use", {
  # Given a model lacking a state/parameter that the caller still supplies
  # When run_model runs
  # Then it warns (rather than silently ignoring the input) — catching typos and wrong
  #      expectations. Failure would let an ignored E column or unused period pass silently.
  scen <- data.frame(population = 1000L, I = 10L, E = 30L, R = 5L)
  expect_warning(run_model(scen, "SIR", nticks = 5L, r0 = 2, infectious_period = 6), "no E state")
  expect_warning(run_model(data.frame(population = 1000L, I = 10L), "SIS",
                           nticks = 5L, r0 = 2, infectious_period = 6, immunity_period = 30),
                 "immunity_period")
  expect_warning(run_model(data.frame(population = 1000L, I = 10L), "SIR",
                           nticks = 5L, r0 = 2, infectious_period = 6, incubation_period = 4),
                 "incubation_period")
})

test_that("run_model reserves agent-array capacity for population growth", {
  # Given a capacity larger than the initial population
  # When run_model builds the model
  # Then the people arrays are allocated to capacity with count == n (reserved slots a
  #      callback can later activate); capacity < n is rejected. Failure would mean the
  #      closed-population allocation can't grow.
  m <- run_model(data.frame(population = 1000L, I = 10L), "SIR", nticks = 5L, r0 = 2,
                 infectious_period = 6, capacity = 4000L)
  expect_equal(m$people$capacity, 4000L)
  expect_equal(m$people$count, 1000L)
  expect_equal(m$people$state$length(), 4000L)          # arrays sized to capacity
  expect_error(run_model(data.frame(population = 1000L), "SIR", nticks = 5L, r0 = 2,
                         infectious_period = 6, capacity = 500L), "capacity")
})

test_that("run_model tracks an extra M state and applies its M->S waning", {
  # Given an SIR model with M registered and 100 agents seeded into M (maternal timer 5)
  # When run for 20 ticks (closed population)
  # Then nodes$M and the waning_m flow exist, M is seeded at tick 0, the maternally immune
  #      wane to S (waning_m > 0), and the living total S+I+R+M is conserved AND equals N
  #      (no vital dynamics, so M->S is purely internal). Failure would mean the extra
  #      state isn't carried/applied — the desync the menagerie step kernels' `waned`
  #      leg would otherwise cause.
  st <- laser_states()
  scen <- data.frame(population = 1000L, I = 10L)
  m <- run_model(scen, "SIR", nticks = 20L, r0 = 2, infectious_period = 6, seed = 1L,
                 extra_states = "M",
                 init = function(model) {
                   s <- model$people$state$values(); s[901:1000] <- st[["M"]]; model$people$state$set(s)
                   tm <- model$people$timer$values(); tm[901:1000] <- 5L; model$people$timer$set(tm)
                 })
  expect_false(is.null(m$nodes$M)); expect_false(is.null(m$nodes$waning_m))
  expect_equal(m$nodes$M$values()[1L, 1L], 100L)         # 100 seeded into M at tick 0
  expect_equal(m$nodes$S$values()[1L, 1L], 1000L - 10L - 100L)
  expect_gt(sum(m$nodes$waning_m$values()), 0)           # maternal waning M->S occurred
  living <- rowSums(m$nodes$S$values()) + rowSums(m$nodes$I$values()) +
            rowSums(m$nodes$R$values()) + rowSums(m$nodes$M$values())
  expect_true(all(living == 1000L))                      # conserved (closed pop)
  expect_true(all(living == rowSums(m$nodes$N$values())))# N agrees (no within-tick vitals)

  expect_error(run_model(scen, "SIR", nticks = 5L, r0 = 2, infectious_period = 6,
                         extra_states = "I"), "states")   # can't repeat a model state
})

test_that("run_model supports a user-defined vaccinated state V (no waning)", {
  # Given an SIR model with a NEW state "V" registered via extra_states, and a step_exit
  #       callback that moves 1,000 susceptibles to V at tick 5 (no timer)
  # When run
  # Then V is assigned a fresh state code (exposed in model$states), V is tracked/carried,
  #      the 1,000 vaccinated are permanent (no waning) and never infected (the disease
  #      kernels leave state V untouched), and S+I+R+V is conserved. Failure would mean a
  #      new state isn't really inert/protected, or the census desyncs.
  scen <- data.frame(population = 10000L, I = 50L)
  m <- run_model(scen, "SIR", nticks = 30L, r0 = 2.5, infectious_period = 6, seed = 1L,
                 extra_states = "V",
                 step_exit = function(model) {
                   if (model$tick != 5L) return(invisible())
                   V <- model$states[["V"]]; S <- model$states[["S"]]
                   s <- model$people$state$values()
                   take <- head(which(s == S), 1000L)        # single node: all active
                   s[take] <- V; model$people$state$set(s)
                   move_count(model$nodes$S, model$nodes$V, length(take), model$tick)  # S->V at slice t+1
                 })
  expect_true(!is.null(m$states[["V"]]) && m$states[["V"]] == 5L)   # new code assigned (5)
  expect_false(is.null(m$nodes$V))                                  # V census tracked
  Vtraj <- rowSums(m$nodes$V$values())
  expect_equal(Vtraj[7L], 1000L)                                   # 1,000 vaccinated (slice 6 -> row 7)
  expect_true(all(Vtraj[7:30] == 1000L))                           # permanent (no waning), never infected
  living <- rowSums(m$nodes$S$values()) + rowSums(m$nodes$I$values()) +
            rowSums(m$nodes$R$values()) + rowSums(m$nodes$V$values())
  expect_true(all(living == 10000L))                               # conserved
})

test_that("a vaccinated state V can wane back to S via step_timer_expire in step_update", {
  # Given a population vaccinated into V with a finite immunity timer, and a step_update
  #       callback running the generic step_timer_expire(V -> S) kernel each tick
  # When run with no disease (r0 = 0) to isolate the V dynamics
  # Then V fills at vaccination then drains back to S as timers expire — waning needs no new
  #      kernel and no step-kernel change, just the existing step_timer_expire. Failure would
  #      mean a user-defined state can't be given timed transitions.
  scen <- data.frame(population = 10000L, I = 0L)
  m <- run_model(scen, "SIR", nticks = 40L, r0 = 0, infectious_period = 6, seed = 1L,
                 extra_states = "V",
                 step_exit = function(model) {           # vaccinate 4,000 at tick 2, timer 8
                   if (model$tick != 2L) return(invisible())
                   V <- model$states[["V"]]; S <- model$states[["S"]]
                   s  <- model$people$state$values(); tm <- model$people$timer$values()
                   take <- head(which(s == S), 4000L)
                   s[take] <- V; tm[take] <- 8L
                   model$people$state$set(s); model$people$timer$set(tm)
                   move_count(model$nodes$S, model$nodes$V, length(take), model$tick)
                 },
                 step_update = function(model) {         # V -> S waning (generic kernel)
                   V <- model$states[["V"]]; S <- model$states[["S"]]
                   waned <- step_timer_expire(model$people$state, model$people$timer,
                              model$people$nodeid, model$people$count, model$nodes$count, V, S)
                   move_count(model$nodes$V, model$nodes$S, waned, model$tick)
                 })
  Vtraj <- rowSums(m$nodes$V$values())
  expect_equal(max(Vtraj), 4000L)                        # V filled to the vaccinated count
  expect_equal(Vtraj[nrow(m$nodes$V$values())], 0L)      # all waned back to S by the end
  expect_true(any(Vtraj == 0L) && which.max(Vtraj) < 40L)# rose then fell
  living <- rowSums(m$nodes$S$values()) + rowSums(m$nodes$I$values()) +
            rowSums(m$nodes$R$values()) + rowSums(m$nodes$V$values())
  expect_true(all(living == 10000L))                     # conserved through vaccinate + wane
})
