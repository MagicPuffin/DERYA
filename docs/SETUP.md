# Setup

Verified working sequence (Ubuntu 24.04, GCC 13.3.0, CMake 3.28.3, Python 3.12).

## C++ toolchain

```
apt-get install cmake   # 3.20+ required; project uses FetchContent
cmake -S . -B build
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

First configure fetches Catch2 (v3.7.1) and nlohmann::json (v3.11.3) from
GitHub via `FetchContent` — needs network access to `github.com`. No p4est
dependency until v0.4 (see `docs/NOTES.md`).

## Python tooling (oracle data generation)

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 scripts/gen_oracle_data.py
```

## Gotchas hit during bootstrap (Task 0), recorded so they aren't rediscovered

- **`pip install exactpack` fails — it's not on PyPI under that name.** Install
  from GitHub instead: `pip install git+https://github.com/lanl/ExactPack.git`.
  `requirements.txt` already does this.
- **The Sod solver class is not called `Sod`.** It's
  `exactpack.solvers.riemann.ep_riemann.IGEOS_Solver` — a general ideal-gas
  Riemann solver whose *default parameters* happen to be the Sod shocktube
  (`xd0=0.5`, ρ_l=1/p_l=1 vs ρ_r=0.125/p_r=0.1, γ=1.4 both sides). Verified
  against the known analytic Sod profile at t=0.2.
- **`target_link_libraries` must list `nlohmann_json::nlohmann_json`
  explicitly** — linking only `Catch2::Catch2WithMain` builds but fails at the
  `#include <nlohmann/json.hpp>` compile step in any file that needs it.
- **`requirements.txt` needs `matplotlib` explicitly, even though nothing here
  plots anything.** `exactpack/solvers/riemann/ep_riemann.py` does
  `import matplotlib.pyplot as plt` at module load time, so importing
  `IGEOS_Solver` fails with `ModuleNotFoundError: No module named
  'matplotlib'` on a clean venv that doesn't happen to have it already
  installed. This didn't surface during the original bootstrap because the
  sandbox environment used to verify Task 0 had matplotlib preinstalled for
  unrelated reasons — a real gap in that verification, not just a missing
  requirements line. `requirements.txt` now lists it explicitly.
- **Regenerated oracle data is not bit-for-bit identical across machines.**
  Running `scripts/gen_oracle_data.py` on Luca's machine (numpy 2.5.3,
  scipy-openblas 0.3.34, ExactPack 1.7.11) changes 6 of the 804 values in
  `sod.json` by 1 ulp (max relative difference 2.1e-16) compared with the
  committed file from the bootstrap environment. Repeat runs on the same
  machine are identical. The cause is floating-point rounding differences
  between library builds; `requirements.txt` doesn't pin versions. Don't
  commit a regenerated oracle file whose only changes are at this level: the
  committed file is the reference, and regression tests compare against it
  with a tolerance, never exact equality.
