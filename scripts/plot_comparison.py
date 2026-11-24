#!/usr/bin/env python3
"""Plot a hydro_run dump against an exact solution.

Draws density, x-velocity, pressure and specific internal energy against x:
every cell of each dump as a point (the cells of a 1D problem on a slab
overlap), over the meshes' current x extent, the oracle as a line, drawn as
a vertical step at each of the oracle's "jumps". Several dumps (e.g. order 1
and order 2 of one deck) overlay on the same axes, labeled by --labels. The exact internal energy is the oracle's
specific_internal_energy, else P / ((gamma - 1) rho) with gamma from
--gamma, the oracle, or 1.4, in that order. Run on demand, not from ctest.

Usage:
  python3 scripts/plot_comparison.py DUMP [DUMP ...]
                                     [--oracle tests/oracle_data/sod.json]
                                     [--labels L ...] [--gamma G]
                                     [-o comparison.png]
"""
import argparse
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = pathlib.Path(__file__).parent.parent
# Categorical series colors in fixed order (dataviz default palette, light).
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def with_jumps(x, f, jumps):
    """Samples (x, f) with a vertical step at each jump between two samples,
    instead of the slanted segment a line plot would draw."""
    x, f = list(x), list(f)
    for jump in sorted(jumps, reverse=True):
        i = next((k for k in range(1, len(x)) if x[k - 1] < jump <= x[k]), None)
        if i is not None:
            x[i:i] = [jump, jump]
            f[i:i] = [f[i - 1], f[i]]
    return x, f


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("dumps", type=pathlib.Path, nargs="+",
                        help="hydro_run dumps (JSON), overlaid")
    parser.add_argument("--labels", nargs="+", default=None,
                        help="one legend label per dump (default: the file names)")
    parser.add_argument("--oracle", type=pathlib.Path,
                        default=REPO / "tests" / "oracle_data" / "sod.json",
                        help="exact solution (default: Sod)")
    parser.add_argument("--gamma", type=float, default=None,
                        help="ideal-gas gamma, for the exact internal energy when "
                             "the oracle has none (default: the oracle's, else 1.4)")
    parser.add_argument("-o", "--output", type=pathlib.Path,
                        default=pathlib.Path("comparison.png"))
    args = parser.parse_args()

    if len(args.dumps) > len(SERIES):
        parser.error(f"at most {len(SERIES)} dumps")
    labels = args.labels or [d.stem for d in args.dumps]
    if len(labels) != len(args.dumps):
        parser.error("--labels needs one label per dump")
    dumps = [json.loads(d.read_text()) for d in args.dumps]
    oracle = json.loads(args.oracle.read_text())
    for path, dump in zip(args.dumps, dumps):
        if abs(dump["t"] - oracle["t"]) > 1e-12 * max(1.0, abs(oracle["t"])):
            print(f"warning: {path} at t = {dump['t']}, oracle at t = {oracle['t']}")

    x_ex = oracle["x"]
    jumps = oracle.get("jumps", [])
    rho_ex = np.asarray(oracle["density"])
    p_ex = np.asarray(oracle["pressure"])
    if "specific_internal_energy" in oracle:
        eps_ex = oracle["specific_internal_energy"]
    else:
        gamma = args.gamma or oracle.get("gamma", 1.4)
        eps_ex = p_ex / ((gamma - 1.0) * rho_ex)
    panels = [("density", "density", rho_ex),
              ("velocity", "velocity_x", oracle["velocity"]),
              ("pressure", "pressure", p_ex),
              ("internal energy", "internal_energy", eps_ex)]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for ax, (name, key, exact) in zip(axes.flat, panels):
        ax.plot(*with_jumps(x_ex, exact, jumps), color="0.2", lw=1.2, label="exact",
                zorder=3)
        for k, (dump, label) in enumerate(zip(dumps, labels)):
            ax.plot(dump["cells"]["centroid_x"], dump["cells"][key], "o", ms=3,
                    mfc=SERIES[k], mec="white", mew=0.4, label=label, zorder=2)
        ax.set_ylabel(name)
        ax.grid(alpha=0.3)
    # The mesh moves: show where the gas is now, not the oracle's whole range.
    nodes_x = [x for dump in dumps for x in dump["nodes"]["x"]]
    axes[0, 0].set_xlim(min(nodes_x), max(nodes_x))
    for ax in axes[1]:
        ax.set_xlabel("x")
    axes[0, 0].legend()
    first = dumps[0]
    fig.suptitle(f"{oracle['case']} at t = {first['t']:g}, "
                 f"{len(first['cells']['density'])} cells")
    fig.tight_layout()
    fig.savefig(args.output, dpi=120)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
