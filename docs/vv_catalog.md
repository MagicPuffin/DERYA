# V&V Test Catalog (Hydro Operator)

Oracle sources: **ExactPack** (`pip install git+https://github.com/lanl/ExactPack.git`)
and the **Kamm-Timmes Sedov solver** (their method, LA-UR-07-2849). Generate via
`scripts/gen_oracle_data.py`, one `CASES` entry per test, output to
`tests/oracle_data/<case>.json`.

## Already in the reference paper's suite (parity checkpoint)
| Test | Stresses | Status |
|---|---|---|
| Sod | shock/contact/rarefaction | passes at 1 resolution via `hydro_run` (L1 density 0.0195, 100×2×2) |
| Noh | wall-heating, stagnation shock | planar passes at 1 resolution via `hydro_run` (L1 density 0.0194, 100 cells) |
| Woodward-Colella blast wave | shock-shock interaction | Tier A/AMR-adjacent, defer |
| Sedov | spherical/cylindrical symmetry | not started |
| Triple point | multi-material vortex | Tier B, defer |

## Additional standard tests (not in the paper)
| Test | Stresses | Priority |
|---|---|---|
| Saltzman piston | mesh-misalignment robustness | High — v0.1/v0.2 |
| **Kidder ball** | no spurious entropy under isentropic compression | **High** — most relevant to FLARE's actual drive regime |
| **Guderley converging shock** | converging-shock focusing, cylindrical symmetry | **High** — most relevant to cylindrical implosion geometry |
| LeBlanc | severe density/pressure jump, Riemann robustness | Medium |
| Coggeshall (pure-hydro subset) | adiabatic self-similar compression | Low — Kidder covers similar ground |
| Gresho vortex | low-Mach vorticity preservation | Low for v1.0 |

## Methodology (run alongside every test, not a test itself)
- **Grid Convergence Index / Richardson extrapolation** (Roache) — run every
  test at ≥3 resolutions, assert observed order matches formal order, not
  just "passes at one resolution."
- **Method of Manufactured Solutions** — the only rigorous way to verify the
  AMR remap/refine/coarsen operator in isolation (relevant from v0.4 on);
  none of the exact-solution tests above isolate AMR from the base hydro
  scheme.

## v0.1 priority order
Sod → Noh → Saltzman → Kidder → Guderley → (Sedov, LeBlanc as time allows)
