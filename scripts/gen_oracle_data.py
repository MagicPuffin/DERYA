#!/usr/bin/env python3
"""Generate exact-solution oracle data for the V&V regression suite.

Writes tests/oracle_data/<case>.json. Each entry is checked into the repo
(small, decouples C++ test runs from needing Python+ExactPack installed) but
this script is kept so the data is reproducible/regenerable.

Usage: python3 scripts/gen_oracle_data.py [case ...]   (default: all cases)

Regenerating a case on another machine can change values by 1 ulp (see
docs/SETUP.md), so name only the cases you mean to regenerate.
"""
import json
import pathlib
import sys
import warnings

import numpy as np
from exactpack.solvers.noh import PlanarNoh
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


def gen_noh():
    """Planar Noh: cold gas, rho0 = 1, u0 = -1 toward a wall at x = 0,
    gamma = 5/3. The shock is at x = t/3 = 0.2 at t = 0.6, with rho = 4,
    eps = 1/2, P = 4/3 and u = 0 behind it. ExactPack's pressure formula
    hard-codes gamma = 5/3 (it is not (gamma-1) rho eps for other gamma), so
    keep gamma = 5/3 here."""
    gamma = 5.0 / 3.0
    # PlanarNoh takes only gamma; rho0 = 1 and u0 = -1 are fixed defaults.
    solver = PlanarNoh(gamma=gamma)
    assert solver.rho0 == 1.0 and solver.u0 == -1.0
    x = np.linspace(0.0, 1.0, 201)
    with warnings.catch_warnings():
        # Its unused spherical-convergence branch divides by x = 0.
        warnings.simplefilter("ignore", RuntimeWarning)
        result = solver(x, t=0.6)
    return {
        "case": "noh",
        "description": "Planar Noh, ExactPack PlanarNoh, gamma = 5/3, rho0 = 1, u0 = -1",
        "t": 0.6,
        "gamma": gamma,
        "x": x.tolist(),
        "density": result["density"].tolist(),
        "pressure": result["pressure"].tolist(),
        "velocity": result["velocity"].tolist(),
        "specific_internal_energy": result["specific_internal_energy"].tolist(),
        # Discontinuity positions: comparisons must not interpolate across them.
        "jumps": [float(j) for j in result.jumps],
    }


CASES = {
    "sod": gen_sod,
    "noh": gen_noh,
    # "saltzman": ..., "kidder": ..., "guderley": ... — added as
    # each test lands in the v0.1/v0.2/v0.3 backlog, one CASES entry at a time.
}


def main():
    names = sys.argv[1:] or list(CASES)
    unknown = [n for n in names if n not in CASES]
    if unknown:
        sys.exit(f"unknown case(s) {unknown}; known: {list(CASES)}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name in names:
        data = CASES[name]()
        out_path = OUT_DIR / f"{name}.json"
        out_path.write_text(json.dumps(data, indent=2))
        print(f"wrote {out_path} ({len(data['x'])} points)")


if __name__ == "__main__":
    main()
