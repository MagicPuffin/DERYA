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

**v0.1 done** (first-order Lagrangian hydro, no AMR): Sod, planar Noh and
the Saltzman piston, run through the `hydro_run` driver, match their exact
solutions at one resolution. **v0.2 done** (second-order extension, paper
Sec. 2.3, selectable per deck): Sod, Noh and Saltzman at 3 resolutions have
second-order errors below first-order ones at every resolution, and the
Taylor-Green vortex converges at order ~2
([docs/vv_catalog.md](docs/vv_catalog.md)).

- Milestones and Definition of Done: [docs/roadmap.md](docs/roadmap.md)
- Task-level progress and what's next: [docs/backlog.md](docs/backlog.md)

## What it does today

- Ideal-gas EOS
- Structured hexahedral mesh (`nz = 1` is the "2D" thin-slab mode), with
  nodal face-area vectors (paper Eq. 5)
- Saltzman's skewed mesh as an optional perturbation of the structured mesh
- Polar mapping of the structured mesh (`"mesh": {"mapping": "polar"}`):
  x is the radius, y the angle in degrees, for quarter shells and quarter
  discs. Down to r = 0 the axis is one node per level and the cells there
  are triangular prisms; walls on the two radial planes pin it
- Boundary conditions on the nodal solver: symmetry/reflecting walls,
  outflow, pistons (walls moving with a prescribed velocity), and a
  prescribed pressure, constant or Kidder's time law, evaluated at each
  stage's time
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
  solution in 211 steps, `decks/noh.json` 0.0194 in 856 steps, and
  `decks/saltzman.json` 0.0454 in 1052 steps. `scripts/plot_comparison.py` plots a dump against
  the exact solution. `--probe cell_id` and `--probe-node node_id` trace
  cells and nodes step by step as JSON lines. Formats:
  [docs/deck_format.md](docs/deck_format.md)
- Second-order scheme, selected with `"scheme": {"order": 2}` in a deck
  (experimental, first order stays the default): Green-formula gradients,
  the SP-V limiter (or none) and a predictor-corrector step, ~4.5x the cost
  of first order. Sod on the same slab: L1 density error 0.0072 (first
  order: 0.0195), flow still exactly 1D. Noh and Saltzman run at second
  order too, with smaller errors than first order
- `scripts/convergence.py`: runs Sod, Noh and Saltzman at 50/100/200 cells
  and LeBlanc at 225/450/900 (`decks/convergence/`) at both orders,
  tabulates L1 density errors and rates, checks the v0.2 gate (and
  LeBlanc's shock position) and plots error against resolution
- LeBlanc shock tube (`decks/leblanc.json`, v0.3): runs at both orders
  through the 1e9 pressure ratio. On a uniform mesh the shock lands
  0.3-0.4 ahead of the exact one at 900 cells: start-up entropy in the few
  heavy cells at the contact, which a uniform mesh resolves only at ~10^4
  cells. With that layer refined (`mesh.x_segments`,
  `decks/convergence/leblanc_layer_*.json`) the shock is within 0.013 at
  order 2 and the v0.3 gate passes. `scripts/leblanc_1d.py` (an independent
  1D model of the first-order step) shows it (see
  [docs/vv_catalog.md](docs/vv_catalog.md))
- Taylor-Green vortex (`decks/taylor_green.json`, with the prescribed
  energy source that keeps it stationary): the second-order scheme
  converges at order ~2 on this smooth flow (L1 pressure order 1.95 at
  160×160), with errors within 25% of Georges' thesis's Table 3.1(b).
  `scripts/taylor_green_convergence.py` tabulates errors and orders

## Out of scope for v1.0

Multi-material/ALE remap, real MPI, radiation, drive coupling and burn
physics. See [docs/roadmap.md](docs/roadmap.md) for the reasons.

## Build and test

Requires CMake ≥ 3.20 and a C++17 compiler. The first configure fetches
Catch2 and nlohmann::json from GitHub.

```
cmake --preset default        # build/: optimized, with debug symbols
cmake --build --preset default
ctest --preset default
```

For stepping through code in a debugger, use the separate unoptimized tree
`build-debug/`: the same three commands with `debug` in place of `default`.
It is about 12 times slower, so the suite takes several minutes there
(almost all of it the convergence tests) instead of about a minute.

To run Sod and plot it (the plot needs the Python environment below):

```
./build/src/app/hydro_run decks/sod.json sod_final.json
python3 scripts/plot_comparison.py sod_final.json -o sod.png
```

For Noh, Saltzman or LeBlanc, use `decks/<case>.json` and add
`--oracle tests/oracle_data/<case>.json` (for LeBlanc also
`--log density pressure`).

Oracle data (exact solutions from ExactPack) is committed under
`tests/oracle_data/`. To regenerate it:

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 scripts/gen_oracle_data.py [case ...]    # e.g. noh; default: all
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
| `scripts/` | Oracle data generation, dump plotting, convergence tables |
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
