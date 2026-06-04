# Timer-based I→S transition for SIS models (no immunity).

For each infectious agent, decrements `timer` by 1. When `timer` reaches
0 the agent transitions directly back to state S.

## Usage

``` r
step_infectious_is(people)
```

## Arguments

- people:

  LaserFrame of agents.
