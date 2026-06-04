# Generalized timer-expiry transition into an *absorbing* (untimed) state.

For each agent currently in `from_state`, decrements `timer` by 1; when
the timer reaches 0 the agent moves to `to_state` and its timer is left
at 0. The destination carries no duration of its own — this is the
"transition to an absorbing state" generalization (laser-generic's
`nb_timer_update`), e.g. I→S in SIS, I→R in SIR, or R→S waning.

## Usage

``` r
step_timer_expire(people, from_state, to_state)
```

## Arguments

- people:

  LaserFrame of agents.

- from_state:

  Integer state code an agent must currently occupy to be eligible.

- to_state:

  Integer state code an agent moves to when its timer expires.

## Details

This is the engine behind
[`step_infectious_is()`](https://clorton.github.io/razer/reference/step_infectious_is.md)
(I→S) and
[`step_recovered_rs()`](https://clorton.github.io/razer/reference/step_recovered_rs.md)
(R→S):
`step_timer_expire(people, laser_states()[["I"]], laser_states()[["S"]])`
is exactly `step_infectious_is(people)`.

**Required people properties:** `state`, `timer` (both integer).
