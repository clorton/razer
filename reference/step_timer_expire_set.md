# Generic timed transition `from_state -> to_state` into a TIMED destination.

Like
[`step_timer_expire()`](https://clorton.github.io/razer/reference/step_timer_expire.md)
but on expiry the agent's `timer` is reset to a fresh draw from
`duration` (the destination state's own clock — e.g. E→I sets the
infectious period, I→R sets a waning-immunity period). Returns per-node
transition counts.

## Usage

``` r
step_timer_expire_set(
  state,
  timer,
  nodeid,
  count,
  n_nodes,
  from_state,
  to_state,
  duration
)
```

## Arguments

- state, timer, nodeid, count, n_nodes:

  As in
  [`step_si()`](https://clorton.github.io/razer/reference/step_si.md).

- from_state:

  Integer state code an agent must occupy to be eligible.

- to_state:

  Integer state code an agent moves to on expiry.

- duration:

  A Distribution for the destination state's timer.

## Value

An integer vector of per-node transition counts (length `n_nodes`).
