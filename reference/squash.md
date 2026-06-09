# Compact a people environment, dropping excluded agents and reclaiming their slots.

Applies a logical `keep` mask (length `people$count`) to every per-agent
[Column](https://clorton.github.io/razer/reference/Column.md) in the
`people` environment — shifting the survivors to the front of each
array, in order — and sets `people$count` to the number kept. Reuse
frees the slots of agents that have left the simulation (e.g. the
deceased) so the per-tick kernels stop iterating over them and `births`
can refill the slots. All Columns are compacted by the SAME mask, so
they remain row-aligned.

## Usage

``` r
squash(people, keep = NULL)
```

## Arguments

- people:

  An environment whose per-agent properties are scalar
  [Column](https://clorton.github.io/razer/reference/Column.md)s (e.g.
  `state`, `timer`, `nodeid`, `dob`, `dod`) plus an integer `count`.

- keep:

  Optional logical vector of length `people$count`; `TRUE` keeps the
  agent. Defaults to "still alive" — every agent whose `state` is not
  `D` (stored as 255 in the u8 state Column).

## Value

The new active count (invisibly); `people$count` is updated in place.

## Examples

``` r
if (FALSE) { # \dontrun{
# Periodically reclaim dead agents during a long run with mortality:
people$count <- squash(people)
} # }
```
