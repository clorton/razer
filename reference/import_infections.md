# Import new infectious cases from a schedule, activating reserved agent slots.

For the given `tick`, scans the schedule (parallel `sched_tick` /
`sched_node` / `sched_count` vectors) and, for every entry whose
`sched_tick == tick`, activates `sched_count` new agents in node
`sched_node`: each takes the next free slot after the active `count`, is
set Infectious with a `timer` drawn from `duration`, and has its
`nodeid` set. The new agents must fit in the reserved capacity (the
`state`/`timer`/`nodeid` Columns are allocated larger than the initial
`count`). Per-node import counts are added to the I census at column
`tick + 1` and written to the `importations` flow at column `tick`.

## Usage

``` r
import_infections(
  state,
  timer,
  nodeid,
  count,
  i_count,
  importations,
  sched_tick,
  sched_node,
  sched_count,
  duration,
  tick
)
```

## Arguments

- state:

  Per-agent `u8` state Column (capacity-sized; imported slots set to I).

- timer:

  Per-agent `u16` timer Column (imported slots set from `duration`).

- nodeid:

  Per-agent `u16` node-id Column (imported slots set to their node).

- count:

  Current active agent count (the first free slot).

- i_count:

  `n_ticks x n_nodes` i32 I census Column (mutated at `tick + 1`).

- importations:

  `(n_ticks-1) x n_nodes` i32 flow Column (set at `tick`).

- sched_tick, sched_node, sched_count:

  Equal-length integer schedule vectors.

- duration:

  A Distribution for the imported cases' infectious timer.

- tick:

  0-based tick index.

## Value

The new active agent count (`count` plus the number imported this tick).

## Details

Returns the new active agent count (the caller stores it back into
`people$count`). Sequential — it touches only the handful of imported
slots, not all agents.
