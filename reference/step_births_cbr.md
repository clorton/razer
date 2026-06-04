# Stochastic birth step using a crude birth rate.

**Parallelism:** the Bernoulli draw (does this agent give birth?) is
parallelised across all active agents. The subsequent bookkeeping —
incrementing `count` and writing new-agent properties — is serial
because it modifies the frame's active count.

## Usage

``` r
step_births_cbr(people, cbr)
```

## Arguments

- people:

  LaserFrame of agents.

- cbr:

  Crude birth rate per agent per tick (probability in \[0, 1\]).

## Details

For each active (non-D) agent, draws Bernoulli(`cbr`) for a birth event.
Each birth creates one new agent that inherits the parent's `node` and
starts in state S with `timer = 0`. Other scalar properties of new
agents take the default value set at `add_scalar_property` time.

Excess births beyond `capacity - count` are silently dropped.

**Required people properties:** `state`, `node`, `timer` (all integer).
