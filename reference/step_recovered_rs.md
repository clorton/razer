# Timer-based R→S transition for waning-immunity models.

For each recovered agent, decrements `timer` by 1. When `timer` reaches
0 the agent becomes susceptible again (state S).

## Usage

``` r
step_recovered_rs(people)
```

## Arguments

- people:

  LaserFrame of agents.

## Details

A fixed-state shorthand for
[`step_timer_expire()`](https://clorton.github.io/razer/reference/step_timer_expire.md)`(people, R, S)`.
