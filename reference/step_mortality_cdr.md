# Stochastic mortality step using a crude death rate.

For each living (non-D) agent, draws Bernoulli(`cdr`) for death. Agents
that die have their state set to D (-1) and `timer` set to 0. Dead
agents remain in the frame — call `$squash(people$state >= 0L)`
periodically to compact.

## Usage

``` r
step_mortality_cdr(people, cdr)
```

## Arguments

- people:

  LaserFrame of agents.

- cdr:

  Crude death rate per agent per tick (probability in \[0, 1\]).

## Details

Each agent's draw is independent; Rayon assigns a fixed slice to each
worker thread. RNG is thread-local (no locking).
