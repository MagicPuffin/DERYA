# Backlog

This file is the source of truth for "what's next," not the chat history.
Update the checkboxes as tasks land.

- **v0.1** (Tasks 0-8, done): first-order Lagrangian hydro. Acceptance bar
  (from `roadmap.md`): Sod, Noh, Saltzman pass at 1 resolution, no AMR, no
  p4est.
- **v0.2** (Tasks 9-15, below): second-order extension, paper Sec. 2.3.

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

## Task 8 — Saltzman piston — ✅ DONE (v0.1 gate met)
Files: `src/bc/bc.{hpp,cpp}`, `src/mesh/mesh.{hpp,cpp}`, `src/io/deck.{hpp,cpp}`,
`src/app/driver.{hpp,cpp}`, `scripts/gen_oracle_data.py`, new
`tests/oracle_data/saltzman.json`, new `decks/saltzman.json`,
`tests/unit/test_{bc,mesh,hydro,driver}.cpp`. Flagged: `src/hydro/hydro.{hpp,cpp}`
(piston walls in the Newton initial guess; underflow guard, below).
- [x] Moving-piston BC: `BoundaryType::Piston`, a wall moving with velocity
      u_w, `V_p . n = u_w . n`, tangential motion free. `BoundarySet` is now
      an array of `bc::Boundary {type, velocity}`, implicitly constructible
      from a `BoundaryType`, so existing `{kOut, kSym, ...}` initializers are
      unchanged. `bc::wall_constraints()` gives each node's orthonormal wall
      normals plus the prescribed normal part of V_p (Gram-Schmidt carries
      the values); `nodal_velocities` adds it to the initial guess and the
      Newton steps stay tangent, so `solve_nodal_velocity` is unchanged.
      Tests: prescribed part on a slab; piston meeting a wall at 63 degrees;
      one step of a piston into cold gas gives the exact post-shock pressure
      4/3 (P = rho Gamma s^2) and the matching work
- [x] Skewed initial mesh: `mesh::apply_saltzman_skew()`, x += (y_max - y)
      sin(pi (x - x_min) / L); deck key `mesh.perturbation: "saltzman"`.
      Test: node positions of the standard 100x10 mesh, positive cell areas
      summing to the domain
- [x] Oracle: planar Noh in the piston frame (ExactPack has no gas piston;
      its `ep_piston` is elastic-plastic), shifted back: shock at x = 0.8,
      piston at 0.6 at t = 0.6. Sampled on [0.6, 1] only
- [x] Deck + acceptance test through `hydro_run` (100x10x1, standard
      Saltzman setup to t = 0.6): no tangling, 1052 steps. Shock in every
      row within x = 0.789-0.811 (tilted along the initial skew, as for
      first-order schemes), plateau rho 3.87-4.22 per cell with mass-weighted
      means rho 4.009, u 0.998, eps 0.504, an overshoot (4.75) where the
      shock meets the top wall, energy -1.2% from the exact piston work, gas
      ahead untouched, L1 density error 0.0454 (bound 0.05)
- [x] Found on the way, in the nodal solver: on the skewed mesh, cold gas
      ahead of the shock carries velocity precursors that decay to
      subnormal values (1e-315 and below). Their Jacobian's 1e-8 tr(J)
      regularization underflowed to zero and the elimination gave 0/0, NaN
      velocities at step 12. Nodes whose trace is below DBL_MIN / 1e-8 now
      count as singular and keep their V_p (an underflow guard, not a
      physical scale). Test: 3 cold cells moving at 1, 1e-320, 0, red without
      the guard. The numerical failure policy (`roadmap.md`) was not forced:
      nothing tangled
- Cost: the test takes ~2 min 20 s in the default Debug build (11 s
  optimized, bit-identical), so it dominates the suite

# v0.2 — second-order extension (paper Sec. 2.3)

The paper gives the gradients (Eq. 6), the SP-V limiter (Eq. 7) and
Barth-Jespersen (Eq. 8), and defers the rest to ref. [30] (Georges' PhD
thesis, Bordeaux 2016, in French, not in `docs/`). Task 9 settles what the
paper leaves open before any code depends on it, as Task 4b did for the
impedance. Every task keeps the first-order path working and selectable:
`order = 1` must reproduce today's results exactly, which is the regression
check for every v0.2 task.

## Task 9 — Decide the unspecified parts of the second-order scheme
Files: `docs/decisions.md`, `docs/backlog.md` (Tasks 10-15 adjusted to the
outcome), `docs/roadmap.md` (gate wording, if agreed). No code.
- [ ] Try to get ref. [30] (theses.fr / HAL may have it open access) and
      Georges, Breil & Maire [29]; read their second-order sections. What
      can't be found gets decided here and marked "our choice"
- [ ] Time integration. The paper does not say. Candidates: the EUCCLHYD
      predictor-corrector (predict geometry and state to t^{n+1/2}, solve
      there, then a full step with those fluxes) or two-stage RK2 (Heun).
      Either needs two nodal solves per step, so roughly twice the cost
- [ ] Where Eq. 6's nodal values come from. (grad P)_c uses P_cfp and
      (grad V)_c uses V_p, which are outputs of a nodal solve: a first-order
      solve at the start of each step (or stage), or lagged values from the
      previous step
- [ ] How the reconstruction enters the solver: P_c and V_c in Eqs. 3-4
      and in the corner forces replaced, per corner, by
      P_c + (grad P)_c . (x_p - x_c) and V_c + (grad V)_c (x_p - x_c) ("the
      reconstruction is evaluated at the nodes of the cell"). Is rho_c a_c
      in the two-shock Z still the cell value? (Expected: yes)
- [ ] x_c: volume centroid or node average (the dump uses the node
      average; a linear field is reproduced exactly only with the right one)
- [ ] Limiter edge cases: Eq. 7's velocity limiter divides by |V_c| (a gas
      at rest, or the Noh wall); stencil C(P(c)) at domain boundaries
- [ ] Timestep: is the current controller (`decisions.md`) still right, or
      does the second-order step need a smaller C_cfl?
- [ ] **v0.2 gate wording** (needs Luca's agreement). The roadmap says "same
      tests, 3 resolutions, correct convergence order", but on solutions
      with shocks every scheme converges at order <= 1; second order shows
      as a smaller error constant there, not a higher rate. Proposed gate:
      Sod/Noh/Saltzman at 3 resolutions with second-order errors below
      first-order ones and the observed rates documented, plus order ~2 on
      one smooth problem (Task 14)
- [ ] Smooth test problem for Task 14. Candidates: a 1D isentropic simple
      wave before it steepens (exact solution by characteristics, oracle
      from a short script), or Kidder's isentropic compression (smooth and
      exact, but brought forward from v0.3)

## Task 10 — Unlimited cell gradients (Eq. 6)
Files: `src/hydro/` (new `gradients.{hpp,cpp}` or added to `hydro`),
`tests/unit/test_gradients.cpp`, the test CMakeLists.
- [ ] Cell centre x_c as decided in Task 9
- [ ] `cell_gradients()`: (grad P)_c = 1/V_c sum_p sum_f S_pf P_cfp n_pf
      and (grad V)_c = 1/V_c sum_p sum_f S_pf V_p (x) n_pf (a 3x3 tensor),
      from given nodal values
- [ ] Tests: a constant field gives zero gradient (closure); a linear field
      sampled at the nodes gives its exact gradient on the distorted test
      mesh (by the GCL); hand-computed values on a single box cell

## Task 11 — Barth-Jespersen and the SP-V limiter (Eqs. 7-8)
Files: `src/hydro/` (limiter), `tests/unit/test_limiter.cpp`.
- [ ] `barth_jespersen(controls, extrapolated, reference)` (Eq. 8)
- [ ] chi_c(P): controls = P of the cells sharing a node with c (C(P(c))),
      extrapolated values at c's own nodes (the SP-V difference from [30])
- [ ] chi_c(V): the same on velocities projected onto V_c / |V_c|, with the
      |V_c| = 0 rule from Task 9
- [ ] Tests: Eq. 8 by hand (above, below, inside the bounds); a linear
      field away from boundaries gives chi = 1; a cell at a local extremum
      gives chi = 0; on a slab with a 1D field every cell of a column gets
      the same chi (1D symmetry, the property the paper wants)

## Task 12 — Second-order nodal solve and corner forces
Files: `src/hydro/hydro.{hpp,cpp}`, `tests/unit/test_hydro.cpp`.
- [ ] Reconstructed corner values (Task 9) used in the assembly of Eq. 4,
      in the two-shock Newton solve and in the corner forces (Eq. 3)
- [ ] Tests: zero gradients reproduce the first-order V_p and forces
      exactly; a uniform flow still gives V_p = V_c; momentum and total
      energy are still conserved in a closed box (conservation does not
      depend on the reconstruction)

## Task 13 — Second-order time step and the scheme option
Files: `src/hydro/`, `src/io/deck.{hpp,cpp}`, `src/app/driver.cpp`,
`docs/deck_format.md`, tests.
- [ ] The time integration from Task 9, as `lagrangian_step_2()` (or an
      order argument), with gradients, limiting and the reconstructed solve
      in each stage
- [ ] Deck keys `scheme.order` (1 or 2, default 1 until v0.2 is accepted)
      and `scheme.limiter` (`"sp-v"` or `"none"`, the latter for smooth
      convergence tests)
- [ ] Tests: `order = 1` decks give results identical to today's;
      `order = 2` conserves mass, momentum and energy to round-off in a
      closed box; Sod at `order = 2` stays exactly 1D on the slab

## Task 14 — Smooth convergence test
Files: `scripts/gen_oracle_data.py`, new oracle and decks, tests.
- [ ] Oracle for the smooth problem chosen in Task 9
- [ ] Order ~1 at `order = 1` and ~2 at `order = 2` over 3+ resolutions
      (L1, unlimited and with SP-V; the limiter may cost some order at
      smooth extrema, which is expected and gets documented)
- [ ] Runs as a ctest; the resolutions are kept small enough for the
      optimized build

## Task 15 — v0.2 acceptance: Sod, Noh, Saltzman at 3 resolutions
Files: decks, tests, `scripts/convergence.py` (new), `docs/vv_catalog.md`.
- [ ] Decks at 3 resolutions per problem (e.g. 50/100/200 cells along the
      flow), both orders
- [ ] `scripts/convergence.py`: runs a deck family, tabulates L1 errors and
      observed rates, plots error vs resolution. The first piece of
      `make_report.py` (roadmap DoD item 6)
- [ ] Acceptance per the gate agreed in Task 9: second-order errors below
      first-order ones at every resolution, rates documented in
      `vv_catalog.md`. Saltzman at second order is where a limiter problem
      would show (the paper's own triple-point bug); if it tangles, the
      numerical failure policy (`roadmap.md`) gets decided then
- [ ] Optional, cheap, recommended before this task: the `--probe cell_id`
      debug flag from `roadmap.md`'s design gaps

## After Task 15
v0.2 is done when the gate agreed in Task 9 is met. The v0.3 backlog
(Sedov, Kidder, Guderley, LeBlanc) gets written then.
