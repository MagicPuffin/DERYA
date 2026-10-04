#!/usr/bin/env python3
"""Plot a hydro_run dump against an exact solution.

Draws density, x-velocity, pressure and specific internal energy against x:
every cell of the dump as a point (the cells of a 1D problem on a slab
overlap), the oracle as a line. The oracle's internal energy is
P / ((gamma - 1) rho). Run on demand, not from ctest.

Usage:
  python3 scripts/plot_comparison.py DUMP [--oracle tests/oracle_data/sod.json]
                                     [--gamma 1.4] [-o comparison.png]
"""
import argparse
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = pathlib.Path(__file__).parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("dump", type=pathlib.Path, help="hydro_run dump (JSON)")
    parser.add_argument("--oracle", type=pathlib.Path,
                        default=REPO / "tests" / "oracle_data" / "sod.json",
                        help="exact solution (default: Sod)")
    parser.add_argument("--gamma", type=float, default=1.4,
                        help="ideal-gas gamma, for the exact internal energy")
    parser.add_argument("-o", "--output", type=pathlib.Path,
                        default=pathlib.Path("comparison.png"))
    args = parser.parse_args()

    dump = json.loads(args.dump.read_text())
    oracle = json.loads(args.oracle.read_text())
    if abs(dump["t"] - oracle["t"]) > 1e-12 * max(1.0, abs(oracle["t"])):
        print(f"warning: dump at t = {dump['t']}, oracle at t = {oracle['t']}")

    cells = dump["cells"]
    x = np.asarray(cells["centroid_x"])
    x_ex = np.asarray(oracle["x"])
    rho_ex = np.asarray(oracle["density"])
    p_ex = np.asarray(oracle["pressure"])
    panels = [
        ("density", cells["density"], rho_ex),
        ("velocity", cells["velocity_x"], oracle["velocity"]),
        ("pressure", cells["pressure"], p_ex),
        ("internal energy", cells["internal_energy"], p_ex / ((args.gamma - 1.0) * rho_ex)),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for ax, (name, numerical, exact) in zip(axes.flat, panels):
        ax.plot(x_ex, exact, color="0.3", lw=1.2, label="exact")
        ax.plot(x, numerical, ".", ms=3, color="C0", label="hydro_run")
        ax.set_ylabel(name)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("x")
    axes[0, 0].legend()
    fig.suptitle(f"{oracle['case']} at t = {dump['t']:g}, "
                 f"{len(x)} cells, {dump['steps']} steps")
    fig.tight_layout()
    fig.savefig(args.output, dpi=120)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
