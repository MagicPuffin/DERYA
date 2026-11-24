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


def gen_leblanc():
    """LeBlanc shock tube (decisions.md, "v0.3 test problems"): gamma = 5/3,
    (rho, eps) = (1, 0.1) for x < 3 and (1e-3, 1e-7) beyond, at rest, x in
    [0, 9], t = 6. Sampled every 0.005 (the 900-cell mesh's cell is 0.01).
    IGEOS_Solver gives no jump positions, so they come from its own star
    state: the contact moves at u*, the shock at rho2 u* / (rho2 - rho_r)
    (mass conservation across it, the gas ahead at rest). The rarefaction's
    ends are continuous, not jumps."""
    gamma, t, x0 = 5.0 / 3.0, 6.0, 3.0
    rho_l, eps_l, rho_r, eps_r = 1.0, 0.1, 1e-3, 1e-7
    solver = IGEOS_Solver(xmin=0.0, xd0=x0, xmax=9.0, gl=gamma, gr=gamma,
                          rl=rho_l, ul=0.0, pl=(gamma - 1.0) * rho_l * eps_l,
                          rr=rho_r, ur=0.0, pr=(gamma - 1.0) * rho_r * eps_r)
    x = np.linspace(0.0, 9.0, 1801)
    with warnings.catch_warnings():
        # Its JWL branch divides by zero for an ideal gas.
        warnings.simplefilter("ignore", RuntimeWarning)
        result = solver(x, t=t)
        # Between the contact and the shock at t = 6 (6.731, 7.975).
        star = solver(np.array([7.5]), t=t)
    u_star, rho_2 = float(star["velocity"][0]), float(star["density"][0])
    jumps = [x0 + u_star * t, x0 + rho_2 * u_star / (rho_2 - rho_r) * t]
    return {
        "case": "leblanc",
        "description": "LeBlanc shock tube, ExactPack IGEOS_Solver, gamma = 5/3, "
                       "(rho, eps) = (1, 0.1) | (1e-3, 1e-7) at x = 3",
        "t": t,
        "gamma": gamma,
        "x": x.tolist(),
        "density": result["density"].tolist(),
        "pressure": result["pressure"].tolist(),
        "velocity": result["velocity"].tolist(),
        "specific_internal_energy": result["specific_internal_energy"].tolist(),
        "jumps": jumps,
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


def gen_saltzman():
    """Saltzman piston: cold gas, rho0 = 1, at rest, gamma = 5/3, pushed by a
    piston at x = t (velocity 1). It is planar Noh in the piston's frame
    (x' = x - t, u' = u - 1), so it comes from PlanarNoh shifted back: shock at
    x = 4t/3 = 0.8 at t = 0.6, rho = 4, u = 1, eps = 1/2, P = 4/3 behind it.
    Sampled only where there is gas, x in [t, 1]."""
    gamma = 5.0 / 3.0
    t = 0.6
    solver = PlanarNoh(gamma=gamma)
    assert solver.rho0 == 1.0 and solver.u0 == -1.0
    x = np.linspace(t, 1.0, 81)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = solver(x - t, t=t)
    return {
        "case": "saltzman",
        "description": "Saltzman piston (u = 1, gamma = 5/3, rho0 = 1), from ExactPack "
                       "PlanarNoh in the piston frame",
        "t": t,
        "gamma": gamma,
        "x": x.tolist(),
        "density": result["density"].tolist(),
        "pressure": result["pressure"].tolist(),
        "velocity": (result["velocity"] + 1.0).tolist(),
        "specific_internal_energy": result["specific_internal_energy"].tolist(),
        "jumps": [float(j) + t for j in result.jumps],
    }


CASES = {
    "sod": gen_sod,
    "noh": gen_noh,
    "saltzman": gen_saltzman,
    "leblanc": gen_leblanc,
    # "sedov": ..., "guderley": ... — added as
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
