# Apply crude-birth-rate births for one tick; newborns enter maternal immunity (M).

For each of the first `count` living agents, draws a birth with per-node
probability `1 - exp(-birth_rate[node, tick])`. Each birth activates the
next reserved slot as a new agent: state M, `timer` from
`maternal_duration`, `dob = tick`, and `dod = tick` plus a Kaplan–Meier
age at death (from `km`). Returns `list(count, born)`: the new active
count (store it back into `people$count`) and the per-node birth count
(add it to the `M` census and a birth report). Capped at the allocated
capacity.

## Usage

``` r
births(
  state,
  timer,
  nodeid,
  dob,
  dod,
  count,
  n_nodes,
  birth_rate,
  maternal_duration,
  km,
  tick
)
```

## Arguments

- state:

  Per-agent `u8` state Column (capacity-sized; newborn slots set to M).

- timer:

  Per-agent `u16` timer Column (newborn slots set from
  `maternal_duration`).

- nodeid:

  Per-agent `u16` node-id Column (newborn slots set to the parent's
  node).

- dob:

  Per-agent `i32` date-of-birth Column (newborn slots set to `tick`).

- dod:

  Per-agent `u32` date-of-death Column (newborn slots set via `km`).

- count:

  Current active agent count (the first free slot).

- n_nodes:

  Number of nodes (the length of the returned `born` vector).

- birth_rate:

  `n_ticks x n_nodes` f64 daily-birth-rate grid; column `tick` is read.

- maternal_duration:

  A Distribution for the newborns' maternal-immunity timer.

- km:

  A KaplanMeierEstimator giving each newborn its age at death.

- tick:

  0-based tick index.

## Value

`list(count = <new active count>, born = integer[n_nodes])`.
