# v0.1 Backlog

v0.1 acceptance bar (from `roadmap.md`): Sod, Noh, Saltzman pass at 1
resolution, no AMR, no p4est. Update the checkboxes below as tasks land —
this file is the source of truth for "what's next," not the chat history.

## Task 0 — Repo & toolchain bootstrap — ✅ DONE
- [x] CMake skeleton, Catch2 v3.7.1 + nlohmann::json v3.11.3 via FetchContent
- [x] Toolchain smoke test green
- [x] ExactPack installed from GitHub (not on PyPI), `matplotlib` added to
      `requirements.txt` (ep_riemann imports pyplot at module load)
- [x] `scripts/gen_oracle_data.py` → `tests/oracle_data/sod.json`,
      reproducible to within 1 ulp across machines, **not** bit-for-bit:
      regenerating on Luca's machine changes 6 of 804 values by 1 ulp (max
      relative difference 2.1e-16); repeat runs on one machine are identical.
      The committed file stays the reference; regression tests must compare
      with a tolerance (see `docs/SETUP.md`)
- [x] Pushed to GitHub (`DERYA`)

## Task 1 — EOS module (`src/eos/`) — ✅ DONE
- [x] `EOS` interface: `pressure(rho, eps)`, `sound_speed(rho, eps)`
- [x] `IdealGasEOS`, header-only, in `src/eos/eos.hpp`
- [x] Test passes: p=(γ-1)ρε=1.0, a=√1.4 for ρ=1.0/ε=2.5/γ=1.4 (note: the
      original task prompt had a wrong expected value, 1.5 — correctly
      caught and flagged rather than "fixed" by editing the test to match)

## Task 2 — Mesh module (`src/mesh/`) — ✅ DONE
- [x] `Mesh` type + test for a structured mesh: node/cell/face counts and
      connectivity (interface + test only). First built as a 2D quad mesh,
      then generalized to 3D-native hexahedra per `decisions.md` ("Mesh is
      3D-native"): `(i,j,k)` indexing, `nz` defaults to 1. Test asserts the
      2×2×1 case (18 nodes, 4 cells, 20 faces, cell 0 corners
      {0,1,4,3,9,10,13,12}) plus a 3×2×2 case and face-orientation checks
- [x] Implement `generate_structured_mesh()` to pass the above
- [x] Test: face-area vector (paper Eq. 5) on one known cell of the generated
      mesh. Eq. 5 is per node of each face (S_pf n_pf, face split into
      triangles around its barycenter p*_f), so a rectangular face alone
      can't tell it from an even A/4 split; the test adds a trapezoidal face
      (hand-derived 5/12, 1/3 weights) plus closure and GCL volume checks
- [x] Implement it: `face_area_vectors(mesh, f)`, oriented out of the owner

## Task 3 — Boundary conditions (`src/bc/`) — ✅ DONE
- [x] Interface. The placeholder `BoundaryCondition::apply(mesh, fields)`
      was replaced before starting: EUCCLHYD applies BCs at boundary nodes,
      on the nodal system `M_p V_p = B` (Eq. 4), not on cell fields (see
      `decisions.md`, "Boundary conditions act on the nodal solver"). One
      `BoundaryType` per domain side (`BoundarySet`), keyed by the new
      `Mesh::face_boundary` side tag (flagged mesh edit, tested)
- [x] Reflecting/symmetry BC: `wall_normals()` (per-node orthonormal wall
      normals; coplanar faces count once) + `solve_nodal_velocity()`
      (Lagrange-multiplier solve with `V_p . n = 0`, 0-3 walls). Needed for Noh
- [x] Outflow BC: `apply_pressure_bcs()`, prescribed pressure P* = owner
      cell pressure, `B_p -= P* S_pf n_pf`. Test includes a uniform state at
      rest staying at rest under outflow + symmetry. Needed for Sod
- [ ] (Kidder's prescribed-motion BC waits for its own task batch; a fixed
      pressure P* (free surface) is a one-line extension of the outflow BC
      when a test needs it)

## Task 4 — Core Lagrangian solver (`src/hydro/`) — ✅ DONE
First-order EUCCLHYD step, paper Sec. 2.2, in `src/hydro/hydro.{hpp,cpp}`
(namespace `hydro::lagrangian`; `hydro::hydro` was ambiguous under `using
namespace hydro`). Links `mesh`, `bc`, `eos`. Tests in `tests/unit/test_hydro.cpp`.
- [x] Vec3: `mesh::Vec3` shared by all modules (`bc::Vec3` aliases it). Index
      stays a plain `int32`; the ghost-flag addressing is deferred to v0.4,
      when there are ghosts to address
- [x] Face-area vector (Eq. 5): already `mesh::face_area_vectors` (Task 2)
- [x] Corner geometry: `cell_face_vectors()` (outward S_pf n_pf per cell,
      local face, face node), `corner_vectors()`, `cell_volume()`. Tests:
      owner/neighbor signs, closure, GCL (n_cp = dV_c/dx_p by finite
      differences on a distorted mesh), box and non-planar volumes
- [x] State + impedance (Eq. 3): `HydroState` (moving mesh, SoA m_c, V_c,
      E_c), `make_state()`, `cell_thermo()` (rho, P, Z = rho a; throws on
      non-positive volume)
- [x] Nodal solver (Eq. 4): `assemble_nodal_system()` over all faces, then
      `nodal_velocities()` = `bc::apply_pressure_bcs` + `bc::solve_nodal_velocity`.
      Tests: hand-computed corner M/B; uniform flow gives V_p = V_c at every
      node on a distorted mesh (symmetry + outflow, and all outflow)
- [x] Corner forces with P_cfp (Eq. 3), momentum + energy update (Eq. 1),
      node update (Eq. 2): `corner_forces()`, `update_cells()`,
      `move_nodes()`, all forward Euler, composed in `lagrangian_step(s, bcs,
      eos, dt)` (dt comes from Task 5)
- [x] Integration tests: uniform rest state stays at rest; closed
      symmetry box conserves mass, total energy and volume to round-off
      (1e-13) over 50 steps of a non-uniform state; Sod on a 100x2x2 slab to
      t = 0.2 with fixed dt stays exactly 1D and has L1 density error 0.0202
      vs `sod.json` (converges at order ~0.6 over nx = 50-400; bound 0.025)
- [x] Before Noh: a cold gas (eps = 0) gives a = 0, Z = 0 and a singular M_p.
      Decided: two-shock impedance, solved by per-node Newton (Task 4b,
      `decisions.md`)

## Task 4b — Two-shock impedance (`src/eos/`, `src/hydro/`) — ✅ DONE
Files: `src/eos/eos.hpp`, `src/hydro/hydro.{hpp,cpp}`,
`tests/unit/test_{eos,hydro}.cpp`.
- [x] `EOS::shock_coefficient(rho, eps)`; ideal gas `(gamma+1)/2`
- [x] `CellThermo` carries `sound_speed` and `shock_coefficient` instead of
      a per-cell impedance; `corner_impedances(s, faces, thermo, V_p)` gives
      `Z_cfp`, which `assemble_nodal_system` and `corner_forces` take
- [x] `nodal_velocities(s, faces, thermo, bcs)`: per-node Newton solve of
      Eq. 4. Tests: two cold slabs colliding (density ratio 4) give the
      two-shock interface velocity 2/3; uniform cold flow gives V_p = V_c
- [x] Existing tests still hold (Sod L1 0.0202 at nx = 100, unchanged to 3
      digits); planar Noh smoke test (100 cells, gamma = 5/3, t = 0.6):
      plateau rho within 1% of 4, energy conserved to 1e-12, L1 density
      error 0.0188, first order over nx = 50-400
- [x] Flagged addition: `cell_thermo` clamps eps = E - |V|^2/2 to 0 when it
      is negative by round-off (1e-13 of the kinetic energy) and throws
      beyond that. A cold, moving gas gave eps = -5.6e-17, so a = NaN; the
      acoustic solver would have hit it too
- Cost: Newton takes ~3-4 iterations per step, so Sod runs ~3.8x slower
  (15 s vs 4 s in the test suite). Lagging Z from the previous step is the
  fallback if this matters later

## Task 5 — Timestep/CFL controller (`src/timestep/`) — ✅ DONE
Files: new `src/timestep/`, new `src/sync/` (header-only `global_reduce`),
`src/CMakeLists.txt`, `tests/unit/CMakeLists.txt`, new
`tests/unit/test_timestep.cpp`. Flagged: `src/hydro/hydro.{hpp,cpp}`
(`lagrangian_step` returns the V_p it used, for the volume criterion).
Criteria in `decisions.md` ("Lagrangian timestep").
- [x] CFL-limited dt test on a trivial single-cell case, known analytic value
      (box cell 2x1x0.5, a = 1: dt = C_cfl * 0.5); volume criterion on an
      expanding cold cell (V_p = x_p: dV/dt = 3V, dt = C_V / 3); growth
      limit and dt_max; cold uniform flow has no CFL/volume limit
- [x] Implement, wired through a (no-op for now) `global_reduce`:
      `timestep::stable_dt(s, eos, V_p, dt_prev, params)`,
      `timestep::characteristic_length()`, `sync::global_reduce()`. A rigid
      translation's dV/dt is only zero to round-off, so rates below 1e-12 of
      their terms count as zero
- Checked outside the suite (a scratch driver, as Task 6 will do): Noh from
  a cold start with dt^0 = 1e-4 takes 856 steps, dt growing to 8.2e-4 (the
  CFL limit of the compressed post-shock cells), L1 0.0194 vs 0.0188 with
  the fixed dt = 5e-4

## Task 6 — Driver & first real run (`src/app/`, minimal `src/io/`) — ✅ DONE
Files: new `src/io/` (`deck`, `dump`), new `src/app/` (`driver`, `main.cpp`
→ `hydro_run`), `src/CMakeLists.txt`, `tests/unit/CMakeLists.txt`, new
`tests/unit/test_driver.cpp`, new `decks/sod.json`, new
`scripts/plot_comparison.py`, new `docs/deck_format.md`.
- [x] Minimal JSON deck schema for Sod (`docs/deck_format.md`): structured
      mesh, ideal gas, one BC per side, background state plus boxes
      (density, velocity, pressure), time control, output. Unknown keys are
      errors, naming the key. Tests: valid deck and defaults, 12 rejections
- [x] `hydro_run <deck.json> [dump.json]`: `app::initial_state()` +
      `app::run()` (stable_dt → lagrangian_step, last step cut to land on
      t_end, `max_steps` cap, failures rethrown with step and time). Tests:
      region assignment, exact t_end, max_steps, a failing first step
- [x] Minimal dump (`io::write_dump`, JSON "hydro-dump-0"): per-cell
      centroid, volume, mass, rho, V, P, eps, plus node coordinates. Not a
      restart file (no connectivity or BCs)
- [x] **First integration milestone:** the test runs the `hydro_run`
      executable on `decks/sod.json` (100x2x2) and compares the dump with
      `sod.json`: 211 steps, L1 density error 0.0195 (bound 0.021; 0.0202
      with the fixed dt of the `lagrangian_step` test), flow exactly 1D.
      ~11 s in the Debug build
- [x] Thin-slab guard: `io::thin_slab_warnings()`, run by `parse_deck`,
      printed by `hydro_run` to stderr. Test: a 0.005-thick z extrusion
      with dx = 0.01 warns; cubic cells, the 100x2x2 slab and an all-1 mesh
      do not
- [x] `scripts/plot_comparison.py DUMP [--oracle] [--gamma] [-o]`: rho, u,
      P, eps vs x, cells as points, oracle as a line. On Sod it shows the
      expected first-order smearing and the eps overshoot at the contact
- Failure policy (roadmap gap) is unchanged: a failed step stops the run
  with exit status 1 and no dump

## Task 7 — Noh through the driver — ✅ DONE
Files: `scripts/gen_oracle_data.py`, new `tests/oracle_data/noh.json`, new
`decks/noh.json`, `tests/unit/test_driver.cpp`, `scripts/plot_comparison.py`.
- [x] Oracle: planar Noh from ExactPack's `PlanarNoh` (gamma = 5/3, rho0 = 1,
      u0 = -1, t = 0.6) in `gen_oracle_data.py`, written without
      regenerating `sod.json` (`gen_oracle_data.py [case ...]`). `PlanarNoh`
      takes only gamma (rho0, u0 are fixed defaults, asserted), and its
      pressure formula hard-codes gamma = 5/3. `noh.json` also stores gamma,
      eps and the shock position (`jumps`)
- [x] `decks/noh.json`: 100x1x1 cubic cells of [0,1]x[0,0.01]x[0,0.01], cold
      inflow, wall at x = 0, outflow at x = 1, `dt_initial` = 1e-4
- [x] Acceptance test running `hydro_run` on it: plateau rho and eps within
      1% of 4 and 1/2, untouched inflow ahead of the shock, total energy
      conserved to 1e-12, L1 density error vs `noh.json` 0.0194 in 856 steps
      (bound 0.021). Comparisons no longer interpolate the oracle across a
      listed jump: plain interpolation over the 0.005 sample spacing gave
      0.0219, because cells 0.0025 wide next to the shock got invented
      intermediate "exact" values. The run-and-compare code of the Sod test
      moved into shared helpers (Sod unchanged at 0.0195)
- [x] `plot_comparison.py` handles Noh: exact eps from the oracle (or its
      gamma) rather than the Sod default gamma = 1.4, jumps drawn as
      vertical steps, x axis limited to the mesh's current extent. The Noh
      plot shows the classic wall heating (eps 0.80, rho 2.51 at the wall)
      and a shock ~3 cells wide

## Task 8 — Saltzman piston
Needs two new features before the test itself; split further when started.
- [ ] Moving-piston BC: a wall with prescribed normal velocity,
      `V_p . n = u_piston` (extends `bc::solve_nodal_velocity`). Kidder will
      need prescribed motion too
- [ ] Skewed initial mesh: Saltzman's node perturbation, as a deck option
      and in the mesh generator (or a separate perturbation function)
- [ ] Oracle: the planar piston shock (exact, simple)
- [ ] Deck + acceptance test through `hydro_run`. First-order schemes are
      known to distort on this mesh; if it tangles, that forces the
      numerical failure policy (`roadmap.md`) to be decided

## After Task 8
v0.1 is done when Sod, Noh and Saltzman pass at one resolution. The v0.2
backlog (second-order MUSCL/SP-V reconstruction, paper Sec. 2.3) gets
written then, not before.
