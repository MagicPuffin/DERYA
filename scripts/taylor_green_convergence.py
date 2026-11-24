#!/usr/bin/env python3
"""Taylor-Green convergence table: pressure errors and orders at t_end.

Runs hydro_run on decks/taylor_green.json at each n x n x 1 mesh, reads the
dumps and prints the L1, L2 and Linf pressure errors of Georges' thesis
Eq. 3.7 against the exact (stationary) pressure at each cell's centroid,
Eq. 3.5, with the observed orders between successive meshes (Eq. 3.8),
next to the thesis's Table 3.1(b) (predictor-corrector, discrete
gradients, unlimited). The deck's constants (density, c1, c2) and scheme
are used as they are, except for --order and --limiter. Run on demand, not
from ctest (ctest covers 10x10 to 40x40): at order 2, optimized, 80x80
takes ~3 min and 160x160 ~18 min (4x the cells, 2x the steps; Task 18);
320x320 (~3 h) has not been run.

Usage:
  python3 scripts/taylor_green_convergence.py [--sizes 10 20 40 80]
      [--order 2] [--limiter none] [--hydro-run build/src/app/hydro_run]
      [--workdir DIR]
"""
import argparse
import json
import math
import pathlib
import subprocess
import tempfile
import time

REPO = pathlib.Path(__file__).parent.parent

# Thesis Table 3.1(b): mesh -> (E_L1, E_L2, E_Linf). The table prints
# 2.11e-3 for E_L1 at 20x20, a typo: the order it lists, 0.82, needs 2.11e-2.
THESIS = {
    10: (3.72e-2, 4.96e-2, 1.17e-1),
    20: (2.11e-2, 2.37e-2, 4.64e-2),
    40: (8.29e-3, 9.61e-3, 1.99e-2),
    80: (2.39e-3, 2.83e-3, 7.46e-3),
    160: (6.21e-4, 7.45e-4, 2.19e-3),
    320: (1.57e-4, 1.89e-4, 5.69e-4),
}


def exact_pressure(tg, x, y):
    return (0.25 * tg["density"] * tg["c1"] ** 2
            * (math.cos(2 * math.pi * x) + math.cos(2 * math.pi * y)) + tg["c2"])


def errors(dump, tg):
    cells = dump["cells"]
    volume = l1 = l2 = linf = 0.0
    for x, y, v, p in zip(cells["centroid_x"], cells["centroid_y"],
                          cells["volume"], cells["pressure"]):
        dp = abs(p - exact_pressure(tg, x, y))
        volume += v
        l1 += v * dp
        l2 += v * dp * dp
        linf = max(linf, dp)
    return l1 / volume, math.sqrt(l2 / volume), linf


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sizes", type=int, nargs="+", default=[10, 20, 40, 80])
    parser.add_argument("--order", type=int, choices=[1, 2], default=2)
    parser.add_argument("--limiter", choices=["none", "sp-v"], default="none")
    parser.add_argument("--deck", type=pathlib.Path,
                        default=REPO / "decks" / "taylor_green.json")
    parser.add_argument("--hydro-run", type=pathlib.Path,
                        default=REPO / "build" / "src" / "app" / "hydro_run")
    parser.add_argument("--workdir", type=pathlib.Path, default=None,
                        help="where decks and dumps go (default: a temporary "
                             "directory, removed at the end)")
    args = parser.parse_args()

    base = json.loads(args.deck.read_text())
    tg = {"density": 1.0, "c1": 1.0, "c2": 1.0}
    tg.update(base["initial"]["taylor_green"])
    with tempfile.TemporaryDirectory() as tmp:
        work = args.workdir or pathlib.Path(tmp)
        work.mkdir(parents=True, exist_ok=True)
        print(f"order {args.order}, limiter {args.limiter}, t_end "
              f"{base['time']['t_end']}")
        print(f"{'mesh':>9} {'steps':>6} {'time':>7}  "
              f"{'E_L1':>9} {'O_L1':>5}  {'E_L2':>9} {'O_L2':>5}  "
              f"{'E_Linf':>9} {'O_Linf':>6}  | thesis E_L1, E_L2, E_Linf")
        previous, previous_n = None, None
        for n in args.sizes:
            deck = json.loads(json.dumps(base))
            deck["mesh"]["cells"] = [n, n, 1]
            deck.setdefault("scheme", {})
            deck["scheme"]["order"] = args.order
            deck["scheme"]["limiter"] = args.limiter
            deck.setdefault("output", {})["log_every"] = 0
            deck_path = work / f"taylor_green_{n}.json"
            dump_path = work / f"taylor_green_{n}_dump.json"
            deck_path.write_text(json.dumps(deck, indent=2))
            start = time.time()
            subprocess.run([str(args.hydro_run), str(deck_path), str(dump_path)],
                           check=True, stdout=subprocess.DEVNULL)
            elapsed = time.time() - start
            dump = json.loads(dump_path.read_text())
            e = errors(dump, tg)
            rates = ([math.log(a / b) / math.log(n / previous_n)
                      for a, b in zip(previous, e)] if previous else [None] * 3)
            cols = "  ".join(f"{err:9.3e} {'' if r is None else f'{r:.2f}':>5}"
                             for err, r in zip(e, rates))
            ref = THESIS.get(n)
            ref_text = ", ".join(f"{x:.2e}" for x in ref) if ref else "-"
            print(f"{n:>4}x{n:<4} {dump['steps']:>6} {elapsed:6.1f}s  {cols}  | {ref_text}",
                  flush=True)
            previous, previous_n = e, n


if __name__ == "__main__":
    main()
