# HYDRO

A 3D-native, cell-centred Lagrangian hydrodynamics code in C++17, following
the EUCCLHYD-type scheme of Colaïtis, Guisset & Breil, *"A cell-centered
AMR-ALE framework for 3D multi-material hydrodynamics, Part I"* (SSRN
preprint 5167115).

The aim is a verified hydro operator, checked against exact solutions at
several resolutions, that can run CINDER's canonical 150 µg DT point design
as a 2D cylindrical hydro-only problem. That run tests whether the resolved
trajectory departs from CINDER's 0D assumption in the way CINDER's own
qualification report predicts. This is a personal project, not a First Light
Fusion deliverable.

## Status

**v0.1 in progress** (first-order Lagrangian hydro, no AMR). The EOS, mesh,
boundary conditions, the first-order Lagrangian step and the timestep
controller are done and unit-tested. The `hydro_run` driver runs JSON decks;
Sod through it matches the exact solution. Noh and Saltzman decks are next.

- Milestones and Definition of Done: [docs/roadmap.md](docs/roadmap.md)
- Task-level progress and what's next: [docs/backlog.md](docs/backlog.md)

## What it does today

- Ideal-gas EOS
- Structured hexahedral mesh (`nz = 1` is the "2D" thin-slab mode), with
  nodal face-area vectors (paper Eq. 5)
- Boundary conditions on the nodal solver: symmetry/reflecting walls and
  outflow
- First-order Lagrangian step (paper Eqs. 1–4): nodal solver, corner forces,
  cell and node updates. It conserves mass, energy and volume to round-off
  in a closed box, and Sod on a 100×2×2 slab matches the exact solution
  (L1 density error 0.020) with a fixed timestep
- Two-shock (Dukowicz) impedance in the nodal solver instead of the paper's
  acoustic one, so cold gases work: planar Noh on 100 cells gives the
  post-shock plateau within 0.2% and converges at first order (see
  [docs/decisions.md](docs/decisions.md))
- Timestep controller: acoustic CFL, relative volume change (which bounds
  dt in a cold gas), growth limit and `dt_max`
- `hydro_run <deck.json> [dump.json]`: runs a JSON deck (mesh, ideal gas,
  boundaries, initial regions, time control) and writes the final state as
  JSON. `decks/sod.json` gives L1 density error 0.0195 against the exact
  solution in 211 steps. `scripts/plot_comparison.py` plots a dump against
  the exact solution. Formats: [docs/deck_format.md](docs/deck_format.md)

## Out of scope for v1.0

Multi-material/ALE remap, real MPI, radiation, drive coupling and burn
physics. See [docs/roadmap.md](docs/roadmap.md) for the reasons.

## Build and test

Requires CMake ≥ 3.20 and a C++17 compiler. The first configure fetches
Catch2 and nlohmann::json from GitHub.

```
cmake -S . -B build
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

To run Sod and plot it (the plot needs the Python environment below):

```
./build/src/app/hydro_run decks/sod.json sod_final.json
python3 scripts/plot_comparison.py sod_final.json -o sod.png
```

Oracle data (exact solutions from ExactPack) is committed under
`tests/oracle_data/`. To regenerate it:

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 scripts/gen_oracle_data.py
```

[docs/SETUP.md](docs/SETUP.md) lists the verified toolchain versions and the
known setup problems.

## Layout

| Path | Contents |
|---|---|
| `src/<module>/` | One CMake target per module: `eos`, `mesh`, `bc`, `hydro`, `timestep`, `sync` (stub), `io`, `app` (`hydro_run`); `fields` (placeholder) |
| `decks/` | Example input decks |
| `tests/unit/` | Catch2 unit and integration tests |
| `tests/oracle_data/` | Exact-solution reference data (JSON) |
| `scripts/` | Oracle data generation, dump plotting |
| `docs/` | Project documentation, see below |

## Documentation

| File | Read it when |
|---|---|
| [roadmap.md](docs/roadmap.md) | You want the milestones, Definition of Done and known design gaps |
| [backlog.md](docs/backlog.md) | You're picking up the next task |
| [task_workflow.md](docs/task_workflow.md) | Before starting any task: sizing, stop condition, review checklist |
| [decisions.md](docs/decisions.md) | Before changing architecture: locked technical decisions |
| [deck_format.md](docs/deck_format.md) | You're writing a deck or reading a dump |
| [vv_catalog.md](docs/vv_catalog.md) | You're adding or checking a verification test |
| [SETUP.md](docs/SETUP.md) | Your build or oracle setup fails |
| [NOTES.md](docs/NOTES.md) | You want the dev log |
