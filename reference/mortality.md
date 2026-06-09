# Apply natural mortality for one tick, returning deaths per node by state.

For each of the first `count` living agents (state != D) whose `dod` (an
absolute tick) is `<= tick`, sets the agent's `state` to D and tallies
the death against the state it occupied. Returns `list(m, s, e, i, r)`
of per-node death counts; the caller decrements those census states (and
records the total deaths flow).

## Usage

``` r
mortality(state, dod, nodeid, count, n_nodes, tick)
```

## Arguments

- state:

  Per-agent `u8` state Column (mutated; the deceased become D = 255).

- dod:

  Per-agent `u32` date-of-death Column (an absolute tick index).

- nodeid:

  Per-agent `u16` 0-based node-id Column.

- count:

  Number of active agents to process.

- n_nodes:

  Number of nodes (the length of each returned vector).

- tick:

  0-based tick index; agents with `dod <= tick` die.

## Value

`list(m, s, e, i, r)` of `integer[n_nodes]` death counts by source
state.
