#!/usr/bin/env python3
"""LeBlanc on an independent 1D model of hydro_run's first-order step (Task 19b).

The model: cells between nodes x_j, each with mass m, velocity u and total
energy E; node velocity u*_j and pressure p*_j from the two-shock (Dukowicz)
interface solve hydro_run uses, Z = rho (a + (gamma + 1)/2 |u* - u|), or
with --solver exact the exact Riemann solver (Lagrangian Godunov); forward
Euler, m du = -(p*_R - p*_L) dt, m dE = -(p*_R u*_R - p*_L u*_L) dt,
x += u* dt; walls at both ends; dt from hydro_run's controller (CFL 0.25 on
min(dx, h) / a with h the slab width, volume fraction 0.1, growth 1.01). On
the uniform 225- and 450-cell meshes it puts the shock where hydro_run does,
to 2e-4 (--compare-dump checks a dump).

Prints the experiments of vv_catalog.md, "LeBlanc": the shock position
(last crossing of rho = 2.5e-3, exact 7.975) on uniform meshes, with only
the heavy or only the light gas refined, with a refined layer next to the
contact, and the entropy P / rho^gamma of the heavy gas near the contact.
About 2 minutes. Run on demand, not from ctest.

Usage:
  python3 scripts/leblanc_1d.py [--solver dukowicz|exact] [--max-cells 3600]
      [--compare-dump DUMP]
"""
import argparse
import json
import pathlib
import time

import numpy as np

GAMMA = 5.0 / 3.0
EXACT_SHOCK = 7.974710


def leblanc_state(nodes):
    """(rho, eps) per cell: (1, 0.1) left of x = 3, (1e-3, 1e-7) right."""
    xc = 0.5 * (nodes[1:] + nodes[:-1])
    return np.where(xc < 3.0, 1.0, 1e-3), np.where(xc < 3.0, 0.1, 1e-7)


def segments(*parts):
    """Nodes from (x0, x1, n) segments, n uniform cells each."""
    nodes = [np.linspace(x0, x1, n + 1) for x0, x1, n in parts]
    return np.concatenate([nodes[0]] + [s[1:] for s in nodes[1:]])


def layer_mesh(n, x_fine, k):
    """n-cell uniform mesh of [0, 9] with [x_fine, 3] refined k times."""
    d = 9.0 / n
    return segments((0.0, x_fine, round(x_fine / d)), (x_fine, 3.0, round((3.0 - x_fine) / d) * k),
                    (3.0, 9.0, round(6.0 / d)))


def star_two_shock(rl, ul, pl, al, rr, ur, pr, ar):
    G = 0.5 * (GAMMA + 1.0)
    u = (rl * al * ul + rr * ar * ur + pl - pr) / (rl * al + rr * ar)
    for _ in range(60):  # Newton on the monotone momentum balance
        dl, dr = u - ul, u - ur
        f = pl - rl * (al + G * np.abs(dl)) * dl - pr - rr * (ar + G * np.abs(dr)) * dr
        df = -rl * (al + 2 * G * np.abs(dl)) - rr * (ar + 2 * G * np.abs(dr))
        step = f / df
        u = u - step
        if np.all(np.abs(step) <= 1e-15 * (1.0 + np.abs(u))):
            break
    dl = u - ul
    return u, pl - rl * (al + G * np.abs(dl)) * dl


def star_exact(rl, ul, pl, al, rr, ur, pr, ar):
    g = GAMMA

    def f(p, r, pk, a):
        shock = (p - pk) * np.sqrt(2.0 / ((g + 1) * r) / (p + (g - 1) / (g + 1) * pk))
        rare = 2 * a / (g - 1) * ((p / pk) ** ((g - 1) / (2 * g)) - 1)
        return np.where(p > pk, shock, rare)

    lo = np.zeros_like(pl)
    hi = 10.0 * np.maximum(pl, pr) + 1e3 * (rl + rr) * (np.abs(ul - ur) + al + ar) ** 2
    for _ in range(120):  # bisection: robust at a 1e9 pressure ratio
        p = 0.5 * (lo + hi)
        above = f(p, rl, pl, al) + f(p, rr, pr, ar) + ur - ul > 0
        hi, lo = np.where(above, p, hi), np.where(above, lo, p)
    p = 0.5 * (lo + hi)
    return 0.5 * (ul + ur) + 0.5 * (f(p, rr, pr, ar) - f(p, rl, pl, al)), p


def run(nodes, h, solver="dukowicz", t_end=6.0, cfl=0.25, volume_fraction=0.1, growth=1.01):
    rho, eps = leblanc_state(nodes)
    x = nodes.astype(float)
    m = rho * np.diff(x)
    u = np.zeros_like(m)
    E = eps.copy()
    t, dt_prev, steps = 0.0, np.inf, 0
    while t < t_end:
        dx = np.diff(x)
        r = m / dx
        p = (GAMMA - 1.0) * r * (E - 0.5 * u * u)
        a = np.sqrt(GAMMA * p / r)
        left, right = (r[:-1], u[:-1], p[:-1], a[:-1]), (r[1:], u[1:], p[1:], a[1:])
        if solver == "exact":
            us, ps = star_exact(*left, *right)
            p0 = star_exact(r[0], -u[0], p[0], a[0], r[0], u[0], p[0], a[0])[1]
            p1 = star_exact(r[-1], u[-1], p[-1], a[-1], r[-1], -u[-1], p[-1], a[-1])[1]
        else:
            us, ps = star_two_shock(*left, *right)
            G = 0.5 * (GAMMA + 1.0)  # walls: u* = 0
            p0 = p[0] - r[0] * (a[0] + G * abs(u[0])) * u[0]
            p1 = p[-1] + r[-1] * (a[-1] + G * abs(u[-1])) * u[-1]
        U = np.concatenate([[0.0], us, [0.0]])
        P = np.concatenate([[p0], ps, [p1]])
        dt = min(growth * dt_prev, np.min(cfl * np.minimum(dx, h) / a))
        rate = np.abs(np.diff(U))
        moving = rate > 0.0
        if moving.any():
            dt = min(dt, np.min(volume_fraction * dx[moving] / rate[moving]))
        dt = min(dt, t_end - t)
        u = u - dt / m * (P[1:] - P[:-1])
        E = E - dt / m * (P[1:] * U[1:] - P[:-1] * U[:-1])
        x = x + dt * U
        t, dt_prev, steps = t + dt, dt, steps + 1
    r = m / np.diff(x)
    return {"xc": 0.5 * (x[1:] + x[:-1]), "rho": r, "u": u, "m": m, "steps": steps,
            "p": (GAMMA - 1.0) * r * (E - 0.5 * u * u), "heavy": nodes[:-1] < 3.0 - 1e-12}


def shock_position(xc, rho, level=2.5e-3):
    last = np.nonzero(rho > level)[0][-1]
    if last + 1 == len(rho):
        return xc[last]
    return xc[last] + (rho[last] - level) / (rho[last] - rho[last + 1]) * (xc[last + 1] - xc[last])


def report(label, nodes, h, solver):
    start = time.time()
    R = run(nodes, h, solver)
    s = shock_position(R["xc"], R["rho"])
    behind = (R["xc"] > s - 1.0) & (R["xc"] < s - 0.2)
    print(f"  {label:<34} {len(nodes) - 1:>6} cells {R['steps']:>6} steps  shock {s:.4f} "
          f"({s - EXACT_SHOCK:+.3f})  u {R['u'][behind].min():.4f}-{R['u'][behind].max():.4f}"
          f"  {time.time() - start:.0f} s", flush=True)
    return R


def entropy_near_contact(R, masses=(0.001, 0.005, 0.01, 0.02, 0.04)):
    """P / rho^gamma over its initial value, at heavy-gas mass m from the
    contact (the heavy plateau holds 0.041; exact: 1 everywhere)."""
    heavy = R["heavy"]
    mass = np.cumsum(R["m"][heavy][::-1])
    S = (R["p"] / R["rho"] ** GAMMA)[heavy][::-1] / ((GAMMA - 1.0) * 0.1)
    return "  ".join(f"m={m:g}: {S[min(np.searchsorted(mass, m), len(S) - 1)]:.2f}"
                     for m in masses)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--solver", choices=["dukowicz", "exact"], default="dukowicz")
    parser.add_argument("--max-cells", type=int, default=3600,
                        help="finest uniform mesh (14400 takes ~5 min)")
    parser.add_argument("--compare-dump", type=pathlib.Path, default=None,
                        help="a hydro_run order-1 dump of a uniform LeBlanc deck: compare")
    args = parser.parse_args()

    if args.compare_dump:
        dump = json.loads(args.compare_dump.read_text())
        xs, rs = np.array(dump["cells"]["centroid_x"]), np.array(dump["cells"]["density"])
        n = len(xs)
        R = run(np.linspace(0.0, 9.0, n + 1), 9.0 / n, args.solver)
        print(f"{n} cells: shock model {shock_position(R['xc'], R['rho']):.4f}, "
              f"hydro_run {shock_position(xs, rs):.4f}; steps {R['steps']} vs {dump['steps']}")
        return

    print(f"Uniform meshes ({args.solver}):")
    n = 225
    while n <= args.max_cells:
        R = report(f"uniform {n}", np.linspace(0.0, 9.0, n + 1), 9.0 / n, args.solver)
        n *= 2
    print("Heavy or light gas refined alone (0.01 cells elsewhere):")
    for nl, nr in ((1200, 600), (4800, 600), (300, 2400)):
        report(f"heavy {nl}, light {nr}", segments((0.0, 3.0, nl), (3.0, 9.0, nr)), 0.01,
               args.solver)
    print("A refined layer [x, 3] next to the contact:")
    for n, x_fine, k in ((900, 2.99, 64), (900, 2.95, 16), (900, 2.95, 64), (450, 2.96, 64)):
        report(f"{n} cells, [{x_fine}, 3] x{k}", layer_mesh(n, x_fine, k), 9.0 / n, args.solver)
    print("Entropy over its initial value, heavy gas at mass m from the contact:")
    for label, nodes in (("uniform 900", np.linspace(0.0, 9.0, 901)),
                         ("900, [2.95, 3] x256", layer_mesh(900, 2.95, 256))):
        print(f"  {label:<20} {entropy_near_contact(run(nodes, 0.01, args.solver))}")


if __name__ == "__main__":
    main()
