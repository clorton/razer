# Generic timed transition `from_state -> to_state` into an UNTIMED destination.

For each agent in `from_state`, decrements its u16 `timer`; on expiry
the agent moves to `to_state` (timer left at 0). Returns the per-node
count of transitions. Compose these (downstream-first) to build models
beyond the named menagerie; apply the counts with `move_count`.
Generalizes the M→S / R→S / I→S or I→R legs.

## Usage

``` r
step_timer_expire(state, timer, nodeid, count, n_nodes, from_state, to_state)
```

## Arguments

- state, timer, nodeid, count, n_nodes:

  As in
  [`step_si()`](https://clorton.github.io/razer/reference/step_si.md).

- from_state:

  Integer state code an agent must occupy to be eligible.

- to_state:

  Integer (untimed) state code an agent moves to on expiry.

## Value

An integer vector of per-node transition counts (length `n_nodes`).
