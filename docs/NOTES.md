# Dev log

## Task 0 — repo & toolchain bootstrap

- Repo skeleton created per the agreed module layout (`src/eos`, `src/mesh`,
  `src/bc`, `src/timestep`, `src/hydro`, `src/sync`, `src/fields`, `src/io`,
  `src/app`) — directories exist, modules themselves are empty until Task 1+.
- CMake project configures and builds cleanly with Catch2 (v3.7.1) and
  nlohmann::json (v3.11.3) via `FetchContent`. Toolchain smoke test
  (`tests/unit/test_toolchain_smoke.cpp`) passes both cases via CTest.
- Deliberately no p4est/libsc dependency yet — v0.1-v0.3 use a minimal
  non-p4est structured mesh per the v0.1 backlog; p4est enters at v0.4.
- Python oracle pipeline verified: ExactPack installs from GitHub (not on
  PyPI), `IGEOS_Solver` (Sod-by-default) generates a correct Sod profile at
  t=0.2, `scripts/gen_oracle_data.py` writes `tests/oracle_data/sod.json`.
- See `docs/SETUP.md` for the exact verified command sequence and the three
  gotchas hit along the way (ExactPack install source, solver class name,
  missing json link target).
- Next: Task 1 (EOS module) per the v0.1 backlog.

## Task 18 — cost scaling check (2026-11-24)

Question: Task 14 saw Taylor-Green at order 2 cost 8.3x from 40x40 to 80x80
but 24x from 80x80 to 160x160 (78 min). Is any part of the step superlinear?
**No.** The 24x was the shared machine. Measured again on the i7-7500U
(2 cores, 4 MB L3), optimized build (RelWithDebInfo), with a desktop
(Chrome) running alongside, so expect ~10% noise.

Full runs, `decks/taylor_green.json` (order 2, unlimited, t = 0.5):

| Mesh | Steps | Wall time | Per step | Per cell-step |
|---|---|---|---|---|
| 40x40 | 280 | 19.0 s | 68 ms | 42 us |
| 80x80 | 562 | 169 s | 300 ms | 47 us |
| 160x160 | 1124 | 1094 s (18 min) | 973 ms | 38 us |

So 80 -> 160 cost 6.5x, under the 8x of 4x the cells and 2x the steps. The
160x160 errors are Task 14's (E_L1 5.98e-4, E_L2 7.23e-4), so it is the same
run. 320x320 would take ~2250 steps x ~4.5 s, about 3 h; not run.

The steps double per refinement, as they should. CFL binds throughout, at
the corner cells: the domain corners are stagnation points of strain rate
pi, so the corner cells are squeezed by exp(-pi t), and the shortest cell
length is 0.208 h at t = 0.5 on every mesh (exp(-pi/2) = 0.208). That is why
dt falls ~5x over the run, at every resolution.

One step's parts, at t = 0, best of 3, in ms (scratchpad benchmark: a copy
of `src/hydro/*.cpp` with a Newton iteration counter, linked against the
build's mesh/bc/timestep):

| Part | 40x40 | 80x80 | 160x160 | 320x320 |
|---|---|---|---|---|
| geometry (`cell_face_vectors`) | 0.59 | 2.6 | 12.4 | 84 |
| thermo (`cell_thermo`) | 0.17 | 0.90 | 3.8 | 11 |
| `wall_constraints` | 0.78 | 4.2 | 11.8 | 64 |
| first-order nodal solve | 12.4 | 61 | 250 | 1036 |
| Z + nodal pressures | 0.7 | 3.0 | 10.6 | 72 |
| gradients | 0.54 | 2.5 | 9.9 | 49 |
| SP-V limiters (P + V) | 2.1 | 8.7 | 54 | 229 |
| reconstruction | 0.12 | 0.59 | 2.9 | 24 |
| second-order nodal solve | 12.3 | 37 | 180 | 793 |
| Z + corner forces | 0.53 | 3.8 | 11.3 | 48 |
| `stable_dt` | 0.92 | 4.9 | 18.5 | 79 |
| **step, order 1** | 15.2 | 70 | 301 | 1332 |
| **step, order 2 unlimited** | 58.6 | 243 | 1010 | 4425 |
| **step, order 2 SP-V** | 69.6 | 252 | 1099 | 5105 |

- Per 4x cells the order-2 step costs 4.15x, 4.15x, 4.4x; order 1 4.6x,
  4.3x, 4.4x. The 320x320 excess is in the memory-bound parts (geometry
  6.8x, reconstruction 8x): its per-cell arrays (CellFaceVectors alone are
  576 B/cell, 59 MB) are far beyond L3. Not a superlinear algorithm.
- Newton: 4 residual evaluations per solve (3 for the second-order solve at
  t = 0 from 80x80 up) at every resolution, and 3.96-4.00 per solve over the
  whole 80x80 run. Iteration counts do not grow with refinement or in time.
- The nodal solves are ~85% of an order-2 step (4 solves per step, two per
  stage). Each Newton evaluation makes four passes over the corners
  (impedances, assembly, flux scale, Jacobian impedances and assembly). That
  is where an optimization would go. It is not needed for scaling and was
  not done (Task 18 is not an optimization pass).

Cost model for the v0.3 resolutions (Task 16, provisional): ~40 us per
cell-step at order 2, ~13 us at order 1, on this machine, alone. A run's
cost is cells x steps x that; the step counts depend on each problem's
dt history.

No follow-up fix task: nothing is superlinear.
