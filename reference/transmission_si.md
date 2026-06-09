# Stochastic transmission S→I into an ABSORBING `I` (the SI model), returning new infections per node.

Like
[`transmission()`](https://clorton.github.io/razer/reference/transmission.md)
but the agent enters `I` permanently — no `timer` is set (`I` is
terminal in SI). Returns the per-node count of new infections; the
caller applies the `S` ↓ / `I` ↑ delta.

## Usage

``` r
transmission_si(state, nodeid, count, foi, tick)
```

## Arguments

- state:

  Per-agent `u8` state Column (mutated; S→I).

- nodeid:

  Per-agent `u16` 0-based node-id Column.

- count:

  Number of active agents to process.

- foi:

  `n_ticks x n_nodes` f64 FOI Column; column `tick` read.

- tick:

  0-based tick index.

## Value

An integer vector of new infections per node (length `n_nodes`).
