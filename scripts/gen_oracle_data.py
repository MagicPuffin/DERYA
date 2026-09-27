#!/usr/bin/env python3
"""Generate exact-solution oracle data for the V&V regression suite.

Writes tests/oracle_data/<case>.json. Each entry is checked into the repo
(small, decouples C++ test runs from needing Python+ExactPack installed) but
this script is kept so the data is reproducible/regenerable.

Usage: python3 scripts/gen_oracle_data.py
"""
import json
import pathlib

import numpy as np
from exactpack.solvers.riemann.ep_riemann import IGEOS_Solver

OUT_DIR = pathlib.Path(__file__).parent.parent / "tests" / "oracle_data"


def gen_sod():
    """Standard Sod shocktube. IGEOS_Solver's defaults ARE the Sod problem
    (xd0=0.5, rho_l=1/p_l=1 vs rho_r=0.125/p_r=0.1, gamma=1.4 both sides) —
    note the class is *not* named `Sod`; it's the general ideal-gas Riemann
    solver, Sod-configured by default."""
    solver = IGEOS_Solver()
    x = np.linspace(0.0, 1.0, 201)
    result = solver(x, t=0.2)
    return {
        "case": "sod",
        "description": "Sod shocktube, ExactPack IGEOS_Solver defaults",
        "t": 0.2,
        "x": x.tolist(),
        "density": result["density"].tolist(),
        "pressure": result["pressure"].tolist(),
        "velocity": result["velocity"].tolist(),
    }


CASES = {
    "sod": gen_sod,
    # "noh": ..., "saltzman": ..., "kidder": ..., "guderley": ... — added as
    # each test lands in the v0.1/v0.2/v0.3 backlog, one CASES entry at a time.
}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, gen in CASES.items():
        data = gen()
        out_path = OUT_DIR / f"{name}.json"
        out_path.write_text(json.dumps(data, indent=2))
        print(f"wrote {out_path} ({len(data['x'])} points)")


if __name__ == "__main__":
    main()
