#!/usr/bin/env python3
"""Convergence of Sod, Noh, Saltzman and LeBlanc: L1 density errors and rates.

Runs hydro_run on each deck family decks/convergence/<problem>_<n>.json
(n cells along the flow) at each scheme order, and prints the L1 density
error against the exact solution, sum_c |rho_c - rho(x_c)| V_c / A (A the
fixed y-z cross-section, so V_c / A is the cell's length; rho at the cell's
centroid), with the observed rate between successive resolutions. Checks
the v0.2 gate (roadmap.md): the second-order error below the first-order
one at every resolution, and for LeBlanc the shock-position bound of the
v0.3 gate. Optionally plots error against resolution. Runs are independent
and go in parallel (--jobs). Run on demand, not from ctest; the ctest covers
the cheaper families.

The exact solutions are evaluated here, not interpolated from
tests/oracle_data/ (sampled every 0.005, too coarse at 200 cells, and the
Sod file marks no jumps): Sod and LeBlanc by an exact Riemann solver, planar
Noh and the Saltzman piston in closed form. They agree with the oracle files
at their sample points (--check-oracles); for LeBlanc that compares this
solver with ExactPack's, an independent evaluation.

Usage:
  python3 scripts/convergence.py [--problems sod noh saltzman leblanc]
      [--sizes N ...] [--orders 1 2] [--jobs 8] [--plot convergence.png]
      [--markdown] [--hydro-run build/src/app/hydro_run] [--check-oracles]
"""
import argparse
import concurrent.futures
import json
import math
import pathlib
import subprocess
import sys
import tempfile
import time

REPO = pathlib.Path(__file__).parent.parent
DECKS = REPO / "decks" / "convergence"


def sod_exact(x, t=0.2, gamma=1.4, left=(1.0, 0.0, 1.0), right=(0.125, 0.0, 0.1),
              x0=0.5):
    """Density of the exact Riemann solution (rho, u, P states) at x."""
    rl, ul, pl = left
    rr, ur, pr = right
    al, ar = math.sqrt(gamma * pl / rl), math.sqrt(gamma * pr / rr)

    def f(p, r, pk, a):  # velocity change across a wave to pressure p
        if p > pk:  # shock
            A, B = 2.0 / ((gamma + 1.0) * r), (gamma - 1.0) / (gamma + 1.0) * pk
            return (p - pk) * math.sqrt(A / (p + B))
        return 2.0 * a / (gamma - 1.0) * ((p / pk) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)

    lo, hi = 1e-12, 10.0 * max(pl, pr)  # bisection on f_l + f_r + du = 0
    for _ in range(200):
        p = 0.5 * (lo + hi)
        if f(p, rl, pl, al) + f(p, rr, pr, ar) + ur - ul > 0.0:
            hi = p
        else:
            lo = p
    p = 0.5 * (lo + hi)
    u = 0.5 * (ul + ur) + 0.5 * (f(p, rr, pr, ar) - f(p, rl, pl, al))
    s = (x - x0) / t
    g = (gamma - 1.0) / (gamma + 1.0)
    if s < u:  # left of the contact
        if p > pl:
            shock = ul - al * math.sqrt((gamma + 1.0) / (2 * gamma) * p / pl
                                        + (gamma - 1.0) / (2 * gamma))
            return rl if s < shock else rl * (p / pl + g) / (g * p / pl + 1.0)
        a_star = al * (p / pl) ** ((gamma - 1.0) / (2.0 * gamma))
        if s < ul - al:
            return rl
        if s > u - a_star:
            return rl * (p / pl) ** (1.0 / gamma)
        return rl * (2.0 / (gamma + 1.0) + g / al * (ul - s)) ** (2.0 / (gamma - 1.0))
    if p > pr:  # right of the contact
        shock = ur + ar * math.sqrt((gamma + 1.0) / (2 * gamma) * p / pr
                                    + (gamma - 1.0) / (2 * gamma))
        return rr if s > shock else rr * (p / pr + g) / (g * p / pr + 1.0)
    a_star = ar * (p / pr) ** ((gamma - 1.0) / (2.0 * gamma))
    if s > ur + ar:
        return rr
    if s < u + a_star:
        return rr * (p / pr) ** (1.0 / gamma)
    return rr * (2.0 / (gamma + 1.0) - g / ar * (ur - s)) ** (2.0 / (gamma - 1.0))


def noh_exact(x, t=0.6, gamma=5.0 / 3.0):
    """Planar Noh, rho0 = 1, u0 = -1 into a wall at x = 0: rho = (g+1)/(g-1)
    behind the shock at x = t (g - 1) / 2, 1 ahead."""
    return (gamma + 1.0) / (gamma - 1.0) if x < 0.5 * (gamma - 1.0) * t else 1.0


def saltzman_exact(x, t=0.6, gamma=5.0 / 3.0):
    """Saltzman piston at u = 1 from x = 0 into cold gas at rest: Noh in the
    piston frame, shock at x = t (g + 1) / 2."""
    return (gamma + 1.0) / (gamma - 1.0) if x < 0.5 * (gamma + 1.0) * t else 1.0


def leblanc_exact(x, t=6.0, gamma=5.0 / 3.0):
    """LeBlanc shock tube (decisions.md): (rho, eps) = (1, 0.1) | (1e-3, 1e-7)
    at x = 3, at rest; P = (gamma - 1) rho eps."""
    return sod_exact(x, t, gamma, left=(1.0, 0.0, (gamma - 1.0) * 0.1),
                     right=(1e-3, 0.0, (gamma - 1.0) * 1e-10), x0=3.0)


EXACT = {"sod": sod_exact, "noh": noh_exact, "saltzman": saltzman_exact,
         "leblanc": leblanc_exact}
SIZES = {"leblanc": [225, 450, 900]}  # the rest: --sizes

# v0.3 gate (roadmap.md): LeBlanc's shock position within 0.1 of 7.975 at
# the finest resolution.
LEBLANC_SHOCK, LEBLANC_SHOCK_TOL = 7.974710, 0.1


def check_oracles():
    """Max |exact - oracle| density over each oracle file's samples, skipping
    samples within 0.006 of a discontinuity (where they are on either side)."""
    jumps = {"sod": [0.5 + 0.92745 * 0.2, 0.5 + 1.75216 * 0.2], "noh": [0.2],
             "saltzman": [0.8], "leblanc": [6.731032, 7.974710]}
    for name, exact in EXACT.items():
        oracle = json.loads((REPO / "tests" / "oracle_data" / f"{name}.json").read_text())
        worst = max(abs(exact(x) - r) for x, r in zip(oracle["x"], oracle["density"])
                    if all(abs(x - j) > 0.006 for j in jumps[name]))
        print(f"{name}: max |exact - oracle| density {worst:.2e}")


def l1_density_error(dump, exact):
    cells = dump["cells"]
    nodes = dump["nodes"]
    area = ((max(nodes["y"]) - min(nodes["y"])) * (max(nodes["z"]) - min(nodes["z"])))
    return sum(abs(r - exact(x)) * v / area for x, v, r in
               zip(cells["centroid_x"], cells["volume"], cells["density"]))


def shock_position(dump, level=2.5e-3):
    """The last crossing of rho = level along x, interpolated linearly between
    the centroids of the last cell above it and the next (LeBlanc: halfway
    between the shocked 4e-3 and the 1e-3 ahead, by the v0.3 gate)."""
    cells = sorted(zip(dump["cells"]["centroid_x"], dump["cells"]["density"]))
    last = max(i for i, (_, r) in enumerate(cells) if r > level)
    if last + 1 == len(cells):
        return cells[last][0]
    (x0, r0), (x1, r1) = cells[last], cells[last + 1]
    return x0 + (r0 - level) / (r0 - r1) * (x1 - x0)


def run_one(hydro_run, problem, n, order, work):
    deck = json.loads((DECKS / f"{problem}_{n}.json").read_text())
    deck.setdefault("scheme", {})["order"] = order
    deck_path = work / f"{problem}_{n}_o{order}.json"
    dump_path = work / f"{problem}_{n}_o{order}_dump.json"
    deck_path.write_text(json.dumps(deck))
    start = time.time()
    proc = subprocess.run([str(hydro_run), str(deck_path), str(dump_path)],
                          capture_output=True, text=True)
    elapsed = time.time() - start
    if proc.returncode != 0:
        return {"error": proc.stderr.strip() or proc.stdout.strip(), "time": elapsed}
    dump = json.loads(dump_path.read_text())
    result = {"l1": l1_density_error(dump, EXACT[problem]), "steps": dump["steps"],
              "time": elapsed}
    if problem == "leblanc":
        result["shock"] = shock_position(dump)
    return result


def plot(results, problems, sizes, orders, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, LogLocator, NullLocator

    # Reference palette (dataviz skill), light mode: categorical slots 1-2,
    # recessive axes and grid, text in ink tokens.
    colors = {1: "#2a78d6", 2: "#eb6834"}
    markers = {1: "o", 2: "s"}
    ink, muted, grid, surface = "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"
    fig, axes = plt.subplots(1, len(problems), figsize=(4.2 * len(problems), 3.6),
                             facecolor=surface, squeeze=False)
    for ax, problem in zip(axes[0], problems):
        ax.set_facecolor(surface)
        for order in orders:
            pts = [(n, results[(problem, n, order)]["l1"]) for n in sizes[problem]
                   if "l1" in results[(problem, n, order)]]
            if not pts:
                continue
            ax.plot(*zip(*pts), color=colors[order], marker=markers[order],
                    markersize=8, linewidth=2, label=f"order {order}")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(sizes[problem])
        ax.set_xticklabels([str(n) for n in sizes[problem]])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_title(problem.capitalize(), color=ink, fontsize=11, loc="left")
        ax.set_xlabel("cells along the flow", color=muted)
        ax.grid(True, color=grid, linewidth=0.8)
        ax.tick_params(colors=muted)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(muted)
    axes[0][0].set_ylabel("L1 density error", color=muted)
    axes[0][-1].legend(frameon=False, labelcolor=ink)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=surface)
    print(f"wrote {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--problems", nargs="+", choices=list(EXACT),
                        default=list(EXACT))
    parser.add_argument("--sizes", type=int, nargs="+", default=None,
                        help="cells along the flow (default 50 100 200, LeBlanc 225 450 900)")
    parser.add_argument("--orders", type=int, nargs="+", choices=[1, 2], default=[1, 2])
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--hydro-run", type=pathlib.Path,
                        default=REPO / "build" / "src" / "app" / "hydro_run")
    parser.add_argument("--plot", type=pathlib.Path, default=None,
                        help="write a PNG of error against resolution")
    parser.add_argument("--markdown", action="store_true",
                        help="print the table as Markdown (for vv_catalog.md)")
    parser.add_argument("--check-oracles", action="store_true",
                        help="compare the exact solutions with tests/oracle_data and exit")
    args = parser.parse_args()
    if args.check_oracles:
        check_oracles()
        return

    sizes = {p: args.sizes or SIZES.get(p, [50, 100, 200]) for p in args.problems}
    keys = [(p, n, o) for p in args.problems for n in sizes[p] for o in args.orders]
    results = {}
    with tempfile.TemporaryDirectory() as tmp, \
            concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        # Largest runs first, so they don't finish last.
        order = sorted(keys, key=lambda k: -(k[1] ** 2 * k[2]))
        futures = {pool.submit(run_one, args.hydro_run, *k, pathlib.Path(tmp)): k
                   for k in order}
        for future in concurrent.futures.as_completed(futures):
            results[futures[future]] = future.result()

    gate_ok = True
    if args.markdown:
        print("| Problem | Cells | " + " | ".join(
            f"Order {o} L1 | rate" for o in args.orders) + " |")
        print("|---|---|" + "---|---|" * len(args.orders))
    for problem in args.problems:
        if not args.markdown:
            print(f"{problem}: L1 density error (rate), steps, time")
        for i, n in enumerate(sizes[problem]):
            cols = []
            for o in args.orders:
                r = results[(problem, n, o)]
                if "l1" not in r:
                    cols.append(("FAILED", "", r["error"].splitlines()[-1][:80]))
                    gate_ok = False
                    continue
                rate = ""
                prev = results[(problem, sizes[problem][i - 1], o)] if i else None
                if prev and "l1" in prev:
                    rate = f"{math.log(prev['l1'] / r['l1']) / math.log(n / sizes[problem][i - 1]):.2f}"
                extra = f", shock {r['shock']:.4f}" if "shock" in r else ""
                cols.append((f"{r['l1']:.3e}", rate,
                             f"{r['steps']} steps {r['time']:.1f} s{extra}"))
            if args.markdown:
                print(f"| {problem} | {n} | " + " | ".join(f"{c[0]} | {c[1]}" for c in cols) + " |")
            else:
                print(f"  {n:>4}  " + "   ".join(
                    f"order {o}: {c[0]} ({c[1] or '-'}), {c[2]}"
                    for o, c in zip(args.orders, cols)))
            if 1 in args.orders and 2 in args.orders:
                e1, e2 = (results[(problem, n, o)].get("l1") for o in (1, 2))
                if e1 is None or e2 is None or not e2 < e1:
                    gate_ok = False
                    print(f"  gate: {problem} at {n}: order 2 not below order 1")
        if problem == "leblanc":
            finest = sizes[problem][-1]
            for o in args.orders:
                shock = results[(problem, finest, o)].get("shock")
                if shock is None or abs(shock - LEBLANC_SHOCK) > LEBLANC_SHOCK_TOL:
                    gate_ok = False
                    print(f"  gate: leblanc at {finest}, order {o}: shock at {shock}, "
                          f"not within {LEBLANC_SHOCK_TOL} of {LEBLANC_SHOCK}")
    if 1 in args.orders and 2 in args.orders:
        print("v0.2 gate (order 2 below order 1 everywhere):", "PASS" if gate_ok else "FAIL")
    if args.plot:
        plot(results, args.problems, sizes, args.orders, args.plot)
    sys.exit(0 if gate_ok else 1)


if __name__ == "__main__":
    main()
