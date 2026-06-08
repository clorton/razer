"""Run the worldwide agent-based SEIR model.

Loads the global admin-2 node table, builds the sparse intra-country coupling (and optionally
the cross-border air coupling for a chosen top-N), runs the model, and writes a time-series
CSV plus plots into ``output/seir/``.

Examples:
    # 1/200-scale world, seed in China, intra-country spread only, 180 days:
    python scripts/10_run_seir.py --subsample 200 --seed-iso3 CHN

    # add the cross-border air network (top-250 airports) and compare arrival timing:
    python scripts/10_run_seir.py --subsample 200 --seed-iso3 CHN --air --air-top-n 250

    # full resolution (needs >32 GB RAM):
    python scripts/10_run_seir.py --subsample 1 --seed-iso3 CHN --nticks 365
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from wwsim import covid_reference as cov  # noqa: E402
from wwsim.abm import SEIRParams, WorldSEIR  # noqa: E402
from wwsim.abm import networks as abm_networks  # noqa: E402
from wwsim.abm import plots as abm_plots  # noqa: E402
from wwsim.config import load_config  # noqa: E402
from wwsim.logging import logger  # noqa: E402
from wwsim.nodes import load_global_nodes  # noqa: E402


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--subsample", type=int, default=200, help="divide each node's pop by this")
    p.add_argument("--nticks", type=int, default=180)
    p.add_argument("--beta", type=float, default=0.35)
    p.add_argument("--incubation", type=int, default=4)
    p.add_argument("--infectious", type=int, default=6)
    p.add_argument("--waning", type=int, default=0, help="0 = SEIR; >0 days = SEIRS")
    p.add_argument("--intra-fraction", type=float, default=0.1, help="intra-country FOI export")
    p.add_argument("--air", action="store_true", help="enable the cross-border air coupling")
    p.add_argument("--air-top-n", type=int, default=None, help="airport cut-off for the air layer")
    p.add_argument("--air-fraction", type=float, default=0.02, help="international FOI export")
    p.add_argument("--seed-count", type=int, default=100)
    p.add_argument("--seed-iso3", default=None, help="seed the most populous node of this country")
    p.add_argument("--seed-nodeid", type=int, default=None)
    p.add_argument("--prng-seed", type=int, default=20260101)
    p.add_argument("--rt-forcing", action="store_true",
                   help="apply a per-country, per-tick beta multiplier from real COVID Rt "
                        "(OWID reproduction_rate); each admin-2 node uses its country's factor")
    p.add_argument("--rt-anchor", default="2020-01-22", help="calendar date for sim day 0")
    p.add_argument("--rt-baseline", type=float, default=None,
                   help="Rt divisor for the multiplier (default = model R0 = beta*infectious)")
    p.add_argument("--rt-clip-max", type=float, default=5.0)
    p.add_argument("--tag", default=None, help="output filename tag")
    p.add_argument("--save-history", action="store_true",
                   help="write per-node, per-tick S/E/I/R to output/seir/history_<tag>.npz "
                        "(input for the animated choropleth, scripts/12_animate_choropleth.py)")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    cfg = load_config()
    out_dir = cfg.output_dir / "seir"
    out_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_global_nodes(cfg).sort_values("global_nodeid").reset_index(drop=True)

    # Resolve the seed node (a specific country's most populous admin-2, if requested).
    seed_nodeid = args.seed_nodeid
    if seed_nodeid is None and args.seed_iso3:
        sub = nodes[nodes["iso3"] == args.seed_iso3.upper()]
        if len(sub) == 0:
            raise SystemExit(f"no nodes for ISO3 {args.seed_iso3!r}")
        seed_nodeid = int(sub.loc[sub["population"].idxmax(), "global_nodeid"])

    params = SEIRParams(
        nticks=args.nticks, beta=args.beta, incubation_period=args.incubation,
        infectious_period=args.infectious, waning_period=args.waning, subsample=args.subsample,
        intra_out_fraction=args.intra_fraction, use_air=args.air, air_top_n=args.air_top_n,
        air_out_fraction=args.air_fraction, seed_count=args.seed_count, seed_nodeid=seed_nodeid,
        prng_seed=args.prng_seed,
    )

    logger.info("building intra-country coupling...")
    intra = abm_networks.build_intra_coupling(nodes, cfg, params.intra_out_fraction)
    air = None
    if params.use_air:
        logger.info("building air coupling (top-%s)...", params.air_top_n)
        air = abm_networks.build_air_coupling(cfg, params.air_top_n, len(nodes), params.air_out_fraction)

    forcing = None
    if args.rt_forcing:
        iso3_order = sorted(nodes["iso3"].unique())
        baseline = args.rt_baseline if args.rt_baseline is not None else args.beta * args.infectious
        logger.info("building Rt forcing (anchor=%s, baseline=%.3g)...", args.rt_anchor, baseline)
        m_ct = cov.rt_forcing(cfg, iso3_order, pd.Timestamp(args.rt_anchor), params.nticks,
                              baseline, clip_max=args.rt_clip_max)
        forcing = (m_ct, iso3_order)

    model = WorldSEIR(nodes, params, intra, air, forcing=forcing)
    model.run()

    totals = model.totals()
    tag = args.tag or f"sub{args.subsample}_beta{args.beta}{'_air' if args.air else ''}{'_rt' if args.rt_forcing else ''}"
    totals.to_csv(out_dir / f"totals_{tag}.csv")
    if args.save_history:
        model.save_history(out_dir / f"history_{tag}.npz")

    # Final attack rate per node (use the SUBSAMPLED denominator the model actually ran on,
    # not the full census population -- R is in subsampled counts).
    attack = np.clip(model.nodes.R[-1] / model.node_pop, 0, 1)
    final = totals[["S", "E", "I", "R"]].iloc[-1]
    logger.info("SEIR done [%s]: peak I=%s (day %d), final attack=%.3f",
                tag, f"{int(totals['I'].max()):,}", int(totals["I"].idxmax()),
                final["R"] / final.sum())

    abm_plots.plot_epi_curve(totals, out_dir / f"epicurve_{tag}.png", title=f"Global SEIR ({tag})")
    abm_plots.plot_attack_choropleth(nodes, attack, out_dir / f"attack_{tag}.png")
    sample_isos = ["CHN", "USA", "GBR", "ITA", "BRA", "ZAF", "IND", "AUS"]
    abm_plots.plot_country_curves(model, sample_isos, out_dir / f"countries_{tag}.png")
    logger.info("SEIR outputs -> %s", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
