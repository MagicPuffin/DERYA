# Backlog

This file is the source of truth for "what's next," not the chat history.
Update the checkboxes as tasks land.

- **v0.1** (Tasks 0-8, done): first-order Lagrangian hydro. Acceptance bar
  (from `roadmap.md`): Sod, Noh, Saltzman pass at 1 resolution, no AMR, no
  p4est.
- **v0.2** (Tasks 9-15, done): second-order extension, paper Sec. 2.3.
- **v0.3** (Tasks 16-25, below): Sedov, Kidder, Guderley, LeBlanc, in 2D.

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
- [x] (Kidder's boundary waited for its own task batch: a prescribed,
      time-dependent pressure, an extension of the outflow BC; done in
      Task 21)

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

## Task 9 — Decide the unspecified parts of the second-order scheme — ✅ DONE
Files: `docs/decisions.md` ("Second-order extension (v0.2)"), this backlog,
`docs/roadmap.md` (gate). No code.
- [x] Ref. [30] found open access (HAL tel-01376682), kept locally as
      `docs/GEORGES_GABRIEL_2016.pdf` (excluded from git like the paper).
      Its Chapter 2 (English) settles most of the open points
- [x] Time integration: the thesis's predictor-corrector (Sec. 2.6.1); GRP
      deferred
- [x] Eq. 6's nodal values: from a first-order solve at each stage (our
      choice; the thesis says only that the gradients "depend on the first
      order solver")
- [x] Reconstruction in the solver: thesis Eqs. 2.43-2.44; cell rho, a; our
      choice: reconstructed V in the two-shock jump too
- [x] x_c = node average (thesis Eq. 2.3)
- [x] Limiter edge cases: `V_c = 0` falls back to the pressure-gradient
      direction, then e_x; no ghost cells in the stencil
- [x] Timestep: unchanged
- [x] **v0.2 gate wording**, agreed with Luca 2026-10-25 and now in
      `roadmap.md`. The old "same tests, 3 resolutions, correct convergence
      order" can't be met: on solutions with shocks every scheme converges
      at order <= 1, and second order shows as a smaller error constant
      there, not a higher rate. New gate: Sod/Noh/Saltzman at 3 resolutions
      with second-order errors below first-order ones and the observed
      rates documented, plus order ~2 on one smooth problem (Task 14)
- [x] Smooth test problem: the Taylor-Green vortex, with the thesis's Table
      3.1(b) as reference numbers

## Task 10 — Unlimited cell gradients (Eq. 6) — ✅ DONE
Files: `src/hydro/` (new `gradients.{hpp,cpp}` or added to `hydro`),
`tests/unit/test_gradients.cpp`, the test CMakeLists.
- [x] `cell_centers()`: x_c = mean of the cell's nodes
- [x] `cell_gradients()`: (grad P)_c = 1/V_c sum_p sum_f S_pf P_cfp n_pf
      and (grad V)_c = 1/V_c sum_p sum_f S_pf V_p (x) n_pf (a 3x3 tensor),
      from given nodal values; plus a helper giving the P_cfp of a
      first-order solve (Eq. 3 at the solved V_p)
- [x] Tests: a constant field gives zero gradient (closure); a linear field
      sampled at the nodes gives its exact gradient on the distorted test
      mesh (thesis Appendix A.2); hand-computed values on a single box cell

## Task 11 — Barth-Jespersen and the SP-V limiter (Eqs. 7-8) — ✅ DONE
Files: `src/hydro/` (limiter), `tests/unit/test_limiter.cpp`.
- [x] `barth_jespersen(controls, extrapolated, reference)` (Eq. 8; the
      preprint's `max(., 1)` in the `b_j < d` branch is a typo for `min`)
- [x] chi_c(P): controls = P of the cells sharing a node with c (C(P(c)),
      c included), extrapolated to c's own nodes
- [x] chi_c(V): the same on velocities projected onto V_c / |V_c|; for
      V_c = 0 onto (grad P)_c / |(grad P)_c|, else e_x
- [x] Tests: Eq. 8 by hand (above, below, inside the bounds); a linear
      field away from boundaries gives chi = 1; a cell at a local extremum
      gives chi = 0; on a slab with a 1D field every cell of a column gets
      the same chi (1D symmetry, the property the paper wants); the V_c = 0
      fallback

## Task 12 — Second-order nodal solve and corner forces — ✅ DONE
Files: `src/hydro/hydro.{hpp,cpp}`, `tests/unit/test_hydro.cpp`. Flagged:
`src/hydro/gradients.{hpp,cpp}` (`reconstructed_corners` needs
`CellGradients` and `cell_centers`, and `gradients.hpp` includes `hydro.hpp`).
- [x] `CornerState`: P~_c(x_p), V~_c(x_p) per cell and local corner;
      `first_order_corners()` (P_c, V_c) and `reconstructed_corners(s,
      thermo, gradients, chi_P, chi_V)` (thesis Eq. 2.1, limited).
      Overloads of `corner_impedances`, `assemble_nodal_system`,
      `nodal_velocities`, `nodal_pressures` and `corner_forces` take one and
      use it in B (thesis Eq. 2.43), in the two-shock jump (against
      V~_c(x_p)), in the Newton initial guess and flux scale, and in Eq. 3 /
      the corner forces (thesis Eq. 2.44). The old signatures are the
      first-order case and delegate. Outflow boundaries keep P* = P_c
      (`decisions.md`)
- [x] Tests: reconstruction by hand on a box cell, zero factors/gradients
      give `first_order_corners` exactly; P_cfp = P~_c(x_p) and Z = rho a
      when V_p = V~_c(x_p); first-order corners reproduce V_p, Z, P_cfp and
      the forces bit for bit (and nonzero ones change V_p); a uniform warm
      flow through the full limited stage keeps V_p = V_c; a closed box
      conserves mass, volume and energy to 1e-13 over 50 forward-Euler steps,
      with unlimited arbitrary and with limited corner values, and the
      corner forces balance at every node
- [x] Sod, Noh and Saltzman dumps through `hydro_run` byte-identical to
      before the task
- Known limitation, left as is (agreed with Luca 2026-10-25, to keep
  `order = 1` exact): a gas that is cold *and* uniform everywhere fails the
  second-order Newton. SP-V limits velocity only along V_c, so round-off
  (1e-16) differences in V~_c(x_p) pass through; the two-shock residual and
  the global flux scale it is judged against are then both ~1e-33, a
  degenerate root where Newton only halves its step. Adding rho_c |V~|^2
  to the flux scale fixes it but changes first-order Sod at round-off
  (<= 2e-13). Noh and Saltzman have a shock from the first step that sets
  the scale; revisit if a real deck hits it

## Task 13 — Predictor-corrector step and the scheme option — ✅ DONE
Files: new `src/hydro/scheme.{hpp,cpp}`, `src/hydro/CMakeLists.txt`,
`src/io/deck.{hpp,cpp}`, `src/app/driver.{hpp,cpp}`, `docs/deck_format.md`,
`tests/unit/test_{hydro,driver,limiter}.cpp`. Flagged:
`src/hydro/limiter.{hpp,cpp}` (round-off V_c in SP-V, below), `README.md`.
- [x] `second_order_fluxes(s, bcs, eos, limiter)`: one stage (first-order
      solve, gradients, SP-V or none, second-order solve, corner forces at
      the converged V_p). `predictor_corrector_step(s, bcs, eos, dt,
      limiter)`: thesis Eqs. 2.45-2.49, the corrector applying each stage's
      fluxes and V_p over dt / 2 from s^n; returns the mean V_p for the
      volume criterion. `scheme_step(..., Scheme{order, limiter})`
      dispatches; order 1 is `lagrangian_step`. dt once per step, before
      the predictor
- [x] Deck block `scheme`: `order` (1 or 2, default 1) and `limiter`
      (`"sp-v"`, default, or `"none"`), in `docs/deck_format.md`
- [x] Tests: `scheme_step` at order 1 is `lagrangian_step` bit for bit;
      the step equals the two stages composed by hand (red when the
      corrector uses only one stage); a closed box conserves mass, energy
      and volume to 1e-13 over 50 steps with both limiters, and its
      momentum changes by exactly the wall nodes' impulse; warm uniform flow
      stays uniform and the mesh translates; deck parsing and 5 rejections;
      Sod at order 2 through `hydro_run` stays exactly 1D (transverse
      velocity < 1e-12, all 4 cells of a column equal to 1e-12), 212 steps,
      L1 density error 0.0072 (first order 0.0195, bound 0.008)
- [x] Sod, Noh and Saltzman dumps at order 1 byte-identical to before
- Found on the way, fixed (flagged, `decisions.md`): Sod at order 2 was
  *not* 1D at first, columns 1.4e-4 apart in density. SP-V's "V_c = 0
  exactly" fallback is decided by round-off: at rest ahead of a wave one
  mirror cell has V_c = 0 and projects on grad P (chi_V = 0, correct), the
  other has a 1e-17 V_c and projects on its arbitrary direction (chi_V
  ~1). V_c now counts as zero below 1e-8 of the largest |V_c'| over
  C(P(c)). Test: a 1e-17 V_c falls back like V_c = 0 (red before). Each
  limiter alone, or none, had stayed 1D
- Measured, not tested (Task 15 material): Noh at order 2, 999 steps, L1
  0.0156 (first order 0.0194), plateau peak rho 4.81 vs 4.19 at first order.
  Cost ~4.5x first order (Sod 3.5 s vs 0.74 s optimized)

## Task 14 — Smooth convergence: the Taylor-Green vortex — ✅ DONE
Files: `src/hydro/scheme.{hpp,cpp}` (energy source), `src/io/deck.{hpp,cpp}`,
`src/app/driver.{hpp,cpp}`, new `decks/taylor_green.json`, new
`scripts/taylor_green_convergence.py`, `tests/unit/test_{hydro,driver}.cpp`,
`docs/{decisions,deck_format,vv_catalog}.md`, `README.md`.
- [x] Prescribed energy source: `EnergySource` (power per unit volume
      S(x)), `source_power` (V_c S(x_c)) and `apply_source`; an optional
      argument of `scheme_step`, `predictor_corrector_step` and
      `second_order_fluxes`, entering each stage like the fluxes (thesis
      Appendix B only adapts GRP; our choice, `decisions.md`). Order 1:
      forward Euler at t^n. Deck: `"source": "taylor_green"` (thesis Eq. 3.6)
- [x] Taylor-Green setup: `initial.taylor_green` {density, c1, c2}
      (default 1, 1, 1) in place of background/regions;
      `taylor_green_{velocity,pressure,source}` in the driver;
      `decks/taylor_green.json` (gamma 1.4, symmetry walls, [0,1]^3 at
      20x20x1, t = 0.5, order 2 unlimited)
- [x] Pressure errors (thesis Eq. 3.7) at order 2, unlimited, against Table
      3.1(b) (whose 20x20 E_L1 "2.11e-3" is a typo for 2.11e-2):

      | Mesh | E_L1 | O_L1 | E_L2 | O_L2 | thesis E_L1 / E_L2 |
      |---|---|---|---|---|---|
      | 10x10 | 2.81e-2 | - | 3.98e-2 | - | 3.72e-2 / 4.96e-2 |
      | 20x20 | 2.03e-2 | 0.47 | 2.28e-2 | 0.80 | 2.11e-2 / 2.37e-2 |
      | 40x40 | 8.08e-3 | 1.33 | 9.44e-3 | 1.27 | 8.29e-3 / 9.61e-3 |
      | 80x80 | 2.31e-3 | 1.81 | 2.76e-3 | 1.77 | 2.39e-3 / 2.83e-3 |
      | 160x160 | 5.98e-4 | 1.95 | 7.23e-4 | 1.93 | 6.21e-4 / 7.45e-4 |

      Thesis orders at 160x160: 1.94 / 1.92. 320x320 not run (160x160
      took 78 min optimized on a shared machine; 18 min alone, Task 18)

      Order 1: E_L1 1.37e-1, 9.09e-2, 5.29e-2, 2.86e-2 (orders 0.59, 0.78,
      0.89)
- [x] ctest `hydro_run: Taylor-Green vortex converges at order ~2`: 10x10,
      20x20, 40x40 at both orders, in-process; errors decrease, order 2
      below order 1 everywhere, the order-2 rate rises (> 1.2 at 40x40,
      0.3 above order 1's), order 1's in [0.6, 1.2], order 2 within a factor
      1.5 of the thesis. ~26 s optimized. Red without the source: the 40x40
      error grows (rate -0.9). The full table is
      `scripts/taylor_green_convergence.py`, on demand (80x80 ~3 min,
      160x160 18 min, Task 18)
- [x] Unit tests: a uniform source heats a gas at rest uniformly (both
      orders); in a closed box the total energy changes by dt sum V_c S(x_c)
      at order 1 and by the mean over the two stages' geometries at order 2
      (the stages differ by 1.5e-5 relative, the check is 1e-10); an empty
      source is bit for bit the step without one; the Taylor-Green fields
      solve steady Euler with the source (central differences, Appendix B);
      the initial state; deck parsing and 5 rejections
- [x] Sod, Noh and Saltzman dumps byte-identical to before the task

## Task 15 — v0.2 acceptance: Sod, Noh, Saltzman at 3 resolutions — ✅ DONE
Files: new `decks/convergence/{sod,noh,saltzman}_{50,100,200}.json`, new
`scripts/convergence.py`, new `docs/figures/convergence_v0.2.png`,
`tests/unit/test_{driver,hydro}.cpp`, `docs/{vv_catalog,decisions,roadmap}.md`,
`README.md`. Flagged: `src/hydro/hydro.{hpp,cpp}` (the negative-eps
tolerance, below; the task named the fix but not the file).
- [x] Decks at 50/100/200 cells along the flow: Sod n×2×2 and Noh n×1×1
      with cubic cells, Saltzman n×(n/10)×1 with z = 1/n; the 100 members
      equal `decks/{sod,noh,saltzman}.json` but for the output block. The
      order is set by the script or test, not the deck
- [x] `scripts/convergence.py`: runs the families in parallel at both
      orders, L1 density error against exact solutions evaluated in the
      script (exact Riemann solver for Sod, closed form for Noh and
      Saltzman; they match `tests/oracle_data/` to 2e-13,
      `--check-oracles`), rates, the gate check (exit 1 on failure),
      `--plot`, `--markdown`. Not yet `make_report.py`: it has no report
      output beyond the table and the PNG
- [x] Gate passes: order 2 below order 1 at every resolution of all three.
      Rates and the figure in `vv_catalog.md`. L1 at 200 cells, order 1 /
      order 2: Sod 1.24e-2 / 3.53e-3, Noh 9.62e-3 / 7.91e-3, Saltzman
      2.19e-2 / 1.56e-2. Rates ~1 on Noh and Saltzman at both orders; Sod
      0.63 / 0.96 (the contact)
- [x] ctest: `v0.2 gate: Sod` and `Noh` at 50/100/200, `Saltzman` at 50×5,
      both orders, in-process, with the oracle interpolation of the older
      tests (~50 s, ~60 s and ~30 s in parallel; the suite now takes about
      a minute). Saltzman 100×10 and 200×20 at order 2 (1 and 8 min) are
      script-only
- [x] Saltzman at order 2 (failed at step 537 since Task 13): the
      negative-eps tolerance in `cell_thermo` is now 1e-13 of the larger
      of the cell's kinetic energy and the mesh's largest |E_c|, so the
      round-off a cold cell picks up from the shocked flow's energy fluxes
      (-2.1e-14) is clamped. It only widens the tolerance, so a run that
      did not throw before is unchanged: order-1 Sod, Noh and Saltzman
      dumps byte-identical. Unit test, red before: eps = -2e-14 in a cold
      cell at rest is clamped next to a cell with E = 1, throws when no
      cell has energy, and -1e-10 always throws. Saltzman 100×10 at order
      2: 1121 steps, L1 3.10e-2 (first order 4.54e-2)
- [x] Noh order-2 overshoot (peak rho 4.81) checked: it is in the 2nd-3rd
      cells from the wall, the start-up/wall-heating error (cell 1 has rho
      2.50, eps 0.80, as at first order); beyond x = 0.02 the plateau is
      within 0.6% of 4 (first order: 1.6%), the shock is two cells wide
      with no oscillation train. SP-V does not show the thesis's SP
      failure on Noh, so the MM limiter is not needed. Saltzman at order 2
      did not tangle at any resolution; the numerical failure policy stays
      undecided
- [x] Not done here (optional): the `--probe cell_id` flag; nothing in
      this task needed it. Done in Task 17

# v0.3 — Sedov, Kidder, Guderley, LeBlanc (paper Sec. 6.4, thesis Ch. 3)

Gate (`roadmap.md`): "all pass, convergence documented". As in v0.2, each
problem runs at 3 resolutions at both orders through
`scripts/convergence.py`, and Task 16 turns "pass" into a bar per problem.

Scope, proposed here and settled in Task 16: the **cylindrical** versions
in 2D (thin slab, one cell in z), the DoD's "2D at 3+ resolutions". The
thesis runs Sedov on a 3D cube and Kidder on an eighth of a spherical
shell; spherical and genuinely 3D runs belong to v0.5 (DoD item 2), and at
today's cost (Taylor-Green 160x160 at order 2: 18 min, Task 18) a 3D cube at 3
resolutions would take hours. Every task keeps the v0.2 decks' results
exactly (`order = 1` and `order = 2` dumps of Sod, Noh, Saltzman and
Taylor-Green byte-identical), unless a task says otherwise and why.

Every problem task also delivers a **profile overlay** (asked by Luca,
2026-11-24): the computed cells over the exact solution, both orders at
the finest resolution, committed to `docs/figures/` and shown in
`vv_catalog.md`. 1D problems against x, as `scripts/plot_comparison.py`
does today; radial problems as the thesis's scatter plots, every cell's
density against its radius over the exact radial profile, which shows
the symmetry spread too. Where the thesis or the paper has a comparable
figure, the published curve goes on the same axes. Those exist only as
images, so they are digitized (e.g. WebPlotDigitizer), the points kept in
`tests/oracle_data/digitized/` with the source figure and a note that
they are approximate, and drawn as a separate, labeled series, never used
in a gate. Candidates: thesis Figs. 3.13-3.15 (Sedov, but in 3D), 3.17c
and 3.19c (Kidder, 3D), the paper's Fig. 29 (spherical Sedov). None is the
same 2D problem, so they are qualitative comparisons only.

## Task 16 — Decide the v0.3 setups — ✅ DONE (gate wording pending Luca)
Files: `docs/decisions.md` (new "v0.3 test problems" entry), this backlog,
`docs/roadmap.md` (proposed gate). No code. Read: thesis Secs. 3.5-3.6,
paper Sec. 6.4 (spherical Sedov on an eighth-sphere AMR mesh, so v0.4/v0.5
parity, not v0.3), ExactPack's `sedov`, `guderley`, `riemann`. Not read:
Maire et al. 2007 (the thesis's [93]), not available here; the 2D Kidder
constants were derived instead (below).
- [x] Per problem: geometry, domain, constants, final time, boundaries,
      mesh and reference, all in `decisions.md`. Corrections to the
      proposals here:
  - *LeBlanc*: constants as proposed, checked: ExactPack and
    `sod_exact` agree to 1e-11; contact at 6.731, shock at 7.975, rho =
    4e-3 and eps = 0.1933 behind it, at t = 6
  - *Sedov*: ExactPack's E0 = 0.311357 puts the shock at r = **0.75**,
    not 1, at t = 1. For r = 1 the quarter plane needs **E = 0.246019 per
    unit length** (alpha / 4, alpha = 0.984074). Cartesian, not polar
  - *Kidder*: **gamma = 2** in 2D (derived: h'' h^(d(gamma-1)+1) must be
    constant, so gamma = 1 + 2/d), t_foc = 0.217945, tau_1 = 0.188746
    (h = 1/2). The same formula at gamma = 5/3 gives the thesis's 3D
    t_foc, 0.254345
  - *Guderley*: gamma = 1.4, ExactPack covers cylindrical (lambda =
    1.197141); Caramana-Whalen time, shock at r = 1 at t = 0, collapse at
    t = 0.750024. Start at t = 0 from the exact profile, polar mesh r in
    [0, 2] with degenerate axis cells, outflow at r = 2, compared on r <=
    1 at t = 0.6 and t = 0.9
- [x] Mesh: Cartesian for LeBlanc and Sedov, polar for Kidder and
      Guderley; theta in degrees in the deck; curved boundaries only as
      pressure boundaries, walls only on planes
- [x] Error metric: radial L1 of density per unit length of slab, at the
      node-average radius; Kidder also the radii and the entropy; symmetry
      as the density spread per radial bin
- [x] Oracles: Sedov and Guderley from `gen_oracle_data.py` (Guderley takes
      ~1 min per time level in ExactPack, so a table is required); LeBlanc
      and Kidder in closed form, in place
- [ ] Gate wording: **proposed** in `roadmap.md` ("v0.3 gate"), not yet
      agreed with Luca. Tasks 19-24 can proceed; the bounds are final once
      agreed, and are not to be changed after seeing results
- [ ] Resolutions, provisional until Task 18 measures the cost: LeBlanc
      225/450/900 x 1 x 1; Sedov 30/60/120 squared; Kidder n_r x n_theta =
      10x40, 20x80, 40x160; Guderley 50x20, 100x40, 200x80 (r in [0, 2]).
      Budget: the finest order-2 run of each family within ~15 min alone,
      optimized; ctest runs the cheap resolutions only (each test <~60 s)
- [x] Known risks recorded in `decisions.md`: SP limiter on Sedov (node
      crossing, bowtie hourglass) and Kidder (inner-face oscillations at
      second order), Kidder's mesh-defect growth, LeBlanc's eps overshoot
      and cold right state, and the zero-area faces of the polar axis cells
      (`src/bc/bc.cpp` normalizes face normals unguarded)

## Task 17 — `--probe` debug flag — ✅ DONE
Files: `src/app/{main,driver}.{hpp,cpp}`, `docs/deck_format.md`,
`tests/unit/test_driver.cpp`. The `roadmap.md` design gap. Also changed,
outside the list: `src/app/CMakeLists.txt` (`app` links nlohmann_json for
the lines), `README.md` and the `roadmap.md` gap entry.
- [x] `hydro_run deck.json [dump.json] --probe cell_id` (repeatable) and
      `--probe-node node_id`, to `--probe-out PATH` or stderr, as JSON
      lines: per probed cell (step, t, dt, centroid, V, m, rho, velocity,
      P, eps, a), then per probed node (position and the step's V_p), at
      step 0 and after every step. `app::Probes` in `driver.hpp`, keys in
      `deck_format.md` ("Probes")
- [x] Tests: the probed values equal the dump's at the last step, exactly;
      the dump is byte-identical with and without probes; an id outside the
      mesh fails the run before the first step, a non-id is a usage error.
      The v0.2 decks' dumps (Sod, Noh, Saltzman, Taylor-Green, both
      orders) are byte-identical to the previous build's, with and without
      probes

## Task 18 — Cost scaling check — ✅ DONE
Files: none required (a report in `NOTES.md` and this backlog); any fix
goes to its own task. Measured in Task 14: Taylor-Green at order 2, 40x40
-> 80x80 cost 8.3x (4x cells, 2x steps: expected 8x), 80x80 -> 160x160
cost 24x. Part of that run shared the machine, so measure again first.
Also changed: the "78 min" in Task 14, the v0.3 scope note and
`scripts/taylor_green_convergence.py`'s docstring.
- [x] Time one step (both orders) at 40x40, 80x80, 160x160 and 320x320
      Taylor-Green, and the parts of it (geometry, thermo, nodal solves and
      Newton iteration counts, gradients, limiter, corner forces): is any
      part superlinear in the number of cells, or does Newton take more
      iterations on finer meshes? **No** (`NOTES.md`, "Task 18"). Per 4x
      cells an order-2 step costs 4.15x, 4.15x, 4.4x (59 ms, 243 ms,
      1.01 s, 4.4 s); the 320x320 excess is cache misses in the memory-bound
      parts. Newton takes 4 evaluations per solve at every mesh and over
      the whole run. Full runs alone: 19 s, 169 s, 18 min (80 -> 160 is
      6.5x); the 78 min was the shared machine. The steps double per
      refinement: CFL binds at the corner stagnation cells, squeezed by
      exp(-pi t) on every mesh. The nodal solves are ~85% of a step.
      Cost model: ~40 us per cell-step at order 2, ~13 us at order 1
- [x] If something is superlinear, a follow-up task to fix it, with results
      bit-identical. Not a general optimization pass. Nothing is: no
      follow-up

## Task 19 — LeBlanc shock tube — ✅ DONE (shock position fails the proposed gate)
Files: `decks/leblanc.json` and `decks/convergence/leblanc_*.json`,
`scripts/gen_oracle_data.py` (+ `tests/oracle_data/leblanc.json`) or an
exact solver in `scripts/convergence.py`, `tests/unit/test_driver.cpp`,
`docs/vv_catalog.md`. No new physics: the first v0.3 problem because it
only needs a deck. Both of the exact-solution options were taken (the
C++ test needs the oracle file). Also changed: `scripts/plot_comparison.py`
(`--log`, `--xlim`), `README.md`, `docs/deck_format.md`, `decisions.md`.
- [x] Decks at the Task 16 resolutions (cubic slab cells), both orders run
      by the script. `decks/leblanc.json` is the 225-cell one
- [x] Exact solution (ExactPack `riemann`), checked against an independent
      evaluation; `convergence.py` family `leblanc` (sizes 225/450/900 by
      default). ExactPack's `IGEOS_Solver` in `tests/oracle_data/leblanc.json`
      (every 0.005, jumps from its star state) and `sod_exact` with the
      LeBlanc states agree to 1.1e-11 (`--check-oracles`)
- [x] Tests: runs to t = 6 at both orders; shock position and the plateau
      density within Task 16's bounds; the gate on the cheap resolutions.
      ctest (`hydro_run: LeBlanc shock tube at both orders`, ~30 s) runs
      225 cells at both orders: mass and energy conserved, flow 1D, no eps
      below the initial 1e-7, the gas ahead untouched, shock and plateau at
      their measured values, order 2's L1 below order 1's. Task 16's
      shock bound applies at 900 cells, in `convergence.py`, and **fails**:
      shock at 8.383 (order 1) and 8.311 (order 2), not within 0.1 of
      7.975 (8.398 at 1800 cells, order 1). Gate items 1-3 pass (L1 order
      1 8.86e-2, 6.73e-2, 4.53e-2; order 2 4.13e-2, 3.18e-2, 1.87e-2)
- [x] Risks: the 1e9 pressure ratio and the cold right state (the
      negative-eps tolerance, the Newton flux scale, SP-V near zero
      velocity); the contact's eps overshoot. Report, don't tune. None of
      the robustness risks was hit; the eps overshoot is there (0.242 and
      0.225 at 900 cells, exact 0.1933). The shocked gas is overdriven (u
      0.67-0.69 at 900 cells, order 1, exact 0.622), converging slowly;
      the few heavy cells behind the contact are the likely cause
      (`vv_catalog.md`, "LeBlanc")
- [x] Profile overlay: rho, u, P, eps against x over the exact solution,
      both orders at 900 cells (log scale for rho and P):
      `docs/figures/profiles_leblanc_900.png`

## Task 19b — LeBlanc shock position — ✅ DONE (understood; gate decision open)
Files: `scripts/leblanc_1d.py` (new), `docs/vv_catalog.md` ("Why the shock
is ahead"), `docs/NOTES.md`, this backlog, `roadmap.md`. No change to the
scheme.
- [x] Is the 0.34-0.41 offset at 900 cells a defect? **No.** An
      independent 1D model of the order-1 step (`leblanc_1d.py`) puts the
      shock where hydro_run does to 2e-4, and on uniform meshes converges:
      +0.42 at 1800 cells, +0.32, +0.21, +0.13 at 3600, 7200, 14400
- [x] Cause: start-up entropy in the few heavy cells at the contact (the
      first one expands against a near-vacuum and ends with 5.3x the
      isentropic entropy at 900 cells). They are the whole heavy plateau,
      so the contact pressure is 6.2e-4 instead of 5.16e-4. Refining the
      light gas changes nothing; refining only [2.95, 3] fixes it (+0.035
      at x64 on 1215 cells, +0.007 at x256). The exact Riemann solver is
      worse (+0.84 at 450), so not the two-shock approximation
- [x] **Decision** (Luca, 2026-11-26): **(c)**. Task 19c builds it. The
      options were:
  - (a) Keep the uniform 225/450/900 meshes and replace the shock bound
        with what they can meet (e.g. the shock error falls from 450 to
        900, or a bound of 0.5), with the 1D model's convergence as the
        evidence
  - (b) Refine the heavy layer at the contact: a graded-mesh deck key
        (e.g. `mesh.x_segments`, cells per segment), decks n = 225/450/900
        with [2.95, 3] refined x64. The 1D model gives +0.035 at 900
        (1215 cells, ~7600 steps: ~2.5 min at order 1, maybe ~10 min at
        order 2 in hydro_run) and the 0.1 bound would hold as written.
        Its own task (mesh + deck + tests), before LeBlanc's gate run
  - (c) Both: (b) for the gate, the uniform family kept as the reported
        literature-style setup

## Task 19c — LeBlanc with a refined layer at the contact — ✅ DONE (gate met)
Decided in Task 19b, option (c): the gate runs on meshes with the heavy
gas next to the contact refined; the uniform family stays, reported, not
gated. Files: `src/mesh/mesh.{hpp,cpp}` (if the structured mesh needs
non-uniform x), `src/io/deck.{hpp,cpp}`, `docs/deck_format.md`,
`decks/convergence/leblanc_layer_*.json`, `scripts/convergence.py`,
`tests/unit/test_{mesh,driver}.cpp`, `docs/{vv_catalog,roadmap}.md`.
Also changed: `README.md`, `decisions.md`, the figure.
- [x] A graded-mesh deck key, `mesh.x_segments: [[x0, x1, n], ...]`
      (uniform cells per segment, contiguous, covering [min, max] in x),
      filling a new `StructuredMeshSpec::x_nodes`. Tests: node positions
      and connectivity, cell masses by region, the rejections, a state at
      rest stays at rest (cells 0.01 to 0.24 wide, walls). Decks without
      the key: the generator's path is unchanged and the Sod, Noh and
      LeBlanc dumps are byte-identical to the previous build's. A deck
      with one segment over the box gives the uniform nodes to 1e-15, not
      bit for bit (segment ends are exact as written). The thin-slab guard
      uses the smallest x cell; the timestep's length is min(dx, slab
      width) on a slab, as it should be
- [x] Decks: base n = 225/450/900 (slab width 9/n), **[2.96, 3]** refined
      x64 (288/576/1152 cells): [2.95, 3] does not fit the 225 mesh. The
      1D model gives +0.051/+0.035/+0.031; hydro_run's order-1 runs match
      it exactly (same steps, shocks to 1e-4)
- [x] `convergence.py`: family `leblanc_layer`, gated; `leblanc`
      (uniform) reported, not gated; `--keep-dumps DIR`
- [x] Tests: `v0.3 gate: LeBlanc with a refined layer, 225 base` at both
      orders (~45 s alone; shared helper with the uniform ctest, which
      stays)
- [x] Order 2 on the layered mesh: shock +0.021, +0.021, +0.013; L1
      1.27e-2, 5.73e-3, 2.12e-3 (rates 1.15, 1.44), 5-13 times below
      order 1. **The gate passes** (`vv_catalog.md`, "Gated: a refined
      layer at the contact")
- [x] Profile overlay: `docs/figures/profiles_leblanc_layer_900.png`

## Task 20 — Polar mesh mapping — ✅ DONE
Files: `src/mesh/mesh.{hpp,cpp}`, `src/io/deck.{hpp,cpp}`,
`src/app/driver.cpp` (if the mapping needs it), `docs/deck_format.md`,
`tests/unit/test_{mesh,bc,driver}.cpp`. Also changed: `src/bc/bc.{hpp,cpp}`
(named in the third item), `tests/unit/test_hydro.cpp` (the axis cells'
face vectors), `README.md`, `decisions.md`. The driver needed nothing.
- [x] `mesh.mapping: "polar"`: the structured mesh's x range is r, its y
      range is theta (degrees, Task 16), z unchanged; nodes at
      (r cos theta, r sin theta, z). Sides keep their tags (x_min = inner
      surface, x_max = outer, y_min/y_max the two radial planes).
      `StructuredMeshSpec::mapping`; theta spans at most 180 degrees, each
      cell less than 180, r >= 0. cos/sin exact at multiples of 90 and
      mirrored about 45 (theta = 30 and 60 swap x and y bit for bit). The
      thin-slab guard measures theta cells as arcs at the innermost ring.
      Cartesian dumps (Sod, Noh, Saltzman, LeBlanc, LeBlanc layer 225,
      Taylor-Green) byte-identical to the previous build's
- [x] Tests: total volume = the annular sector's polygonal approximation
      (to 1e-15), face vectors close per cell (GCL), wall normals on the
      theta = 0 and 90 planes are exactly (0, -1, 0) and (-1, 0, 0), a
      uniform state at rest stays at rest with walls and outflow (both
      orders, velocities below 1e-14)
- [x] r_min = 0 (Guderley): the axis nodes are **one node per level**
      (separate copies would neither be pinned nor have a regular M_p), the
      axis cells triangular prisms, their x_min faces of zero area.
      `wall_constraints` skips zero-area faces (the nodal solver already
      did). The axis node is pinned by the two walls and stays at r = 0.
      Tests: node numbering and merged corners, zero area exactly on the
      axis faces only, no NaN in the wall constraints, corner vectors of
      the prisms (closure, dV/dx_p), volumes positive, a uniform state at
      rest stays at rest at both orders. No Cartesian fallback needed
      (`decisions.md`)

## Task 21 — Prescribed time-dependent pressure boundary — ✅ DONE
Files: `src/bc/bc.{hpp,cpp}`, `src/hydro/{hydro,scheme}.{hpp,cpp}`,
`src/app/driver.cpp`, `src/io/deck.{hpp,cpp}`, `docs/deck_format.md`,
tests. Also changed: `README.md`, `decisions.md`. `hydro.cpp` and
`deck.hpp` needed nothing (`hydro.hpp` only comments).
- [x] `BoundaryType::Pressure`: P* a constant (`bc::pressure_boundary(p)`)
      or a law (`pressure_boundary(PressureLaw)`, `bc::kidder_pressure(p0,
      t_focus, gamma)`, Eqs. 3.13-3.14, throwing at t >= t_focus), where
      outflow has P* = P_c; not a wall. Deck: `{"type": "pressure",
      "value": P}` or `{"type": "pressure", "law": "kidder", "p0": ..,
      "t_focus": ..}`, gamma from `eos`. Time: `bc::at_time(bcs, t)`
      evaluates the laws once per stage, the nodal solver stays time-free
      and `apply_pressure_bcs` throws `std::logic_error` on an unevaluated
      law. `scheme_step` and `predictor_corrector_step` take t^n (before
      dt); the driver passes `r.t`. Stage times from thesis Eq. 2.42
      (checked in Sec. 2.6.1): order 1 at t^n, order 2's predictor at t^n
      and corrector stage at t^n + dt (`decisions.md`). Sod, Noh,
      Saltzman, Taylor-Green and LeBlanc layer 225 dumps byte-identical to
      the previous build's
- [x] Tests: P* equal to the gas pressure (a constant and a law) keeps a
      state at rest on the distorted mesh at both orders; P* = P + 0.1
      pushes the x_max nodes in at the two-shock speed solving
      rho (a + Gamma |V|) |V| = dP (the acoustic rate dP / (rho a) for small
      dP), with the last cell's dV and dE by hand (first step);
      a recording law is called at {t^n} at order 1 and {t^n, t^n + dt} at
      order 2, and the step equals the hand-built one with each stage's
      constant P* (red when the corrector uses t^n: checked by mutation);
      `apply_pressure_bcs` with a law, `at_time`, `kidder_pressure` values
      (h = 1/2 at tau_1: x16 at gamma = 2, x32 at 5/3); deck parsing and
      rejections; through the driver, a Kidder law with t_focus inside the
      run fails one step earlier at order 2 (its corrector reaches t_focus
      first)

## Task 22 — Sedov (cylindrical) — ✅ DONE (order 2 in Task 22b)
Files: decks, `src/io/deck.{hpp,cpp}` and `src/app/driver.{hpp,cpp}` (the
energy-spike initialization), `scripts/{gen_oracle_data,convergence}.py`,
oracle data, tests, `docs/vv_catalog.md`. Also changed:
`scripts/plot_comparison.py` (`--radial`), `docs/{deck_format,decisions}.md`,
`README.md`, the figure.
- [x] Initialization: `initial.energy_spike: {"at": [x, y, z],
      "energy_per_length": E}`; E L_z (slab thickness) goes into the cells
      with a node at `at`, split by mass, replacing their background eps.
      Decks `decks/convergence/sedov_{30,60,120}.json`, E = 0.246019.
      Tests: parsing and rejections; the corner cell holds E L_z, an
      interior node's four cells share it, a point with no node throws
- [x] Oracle: ExactPack `Sedov(geometry=2, eblast=4 x 0.246019)` at t = 1,
      every 0.001 on r in [0, 1.7] (the box corner is at 1.697), shock at
      1.0000005 as the jump (`tests/oracle_data/sedov.json`). Family
      `sedov` in `convergence.py`: radial L1 per unit length, shock = the
      densest cell's radius (gated within 0.05 of 1 at the finest), peak
      rho and the symmetry spread reported
- [x] Tests (order 1): `hydro_run: cylindrical Sedov at order 1, 30x30`
      (~7 s): energy and mass to 1e-12, v_z ~ 0, diagonal mirror symmetry,
      the gas ahead untouched, shock within 0.05 of 1 and at the measured
      0.9677, peak rho 5.007, radial L1 0.1664. Order 1 at 30/60/120: L1
      0.166, 0.114 (listed as 0.116 at first: a slip, Task 22b), 0.052; shock 0.968, 0.981, 0.992 (`vv_catalog.md`,
      "Sedov")
- [x] **Order 2 with SP-V fails** at every resolution (fixed in Task 22b) (negative eps in a
      cold cell just ahead of the shock, step 85 at n = 30). Cause: SP-V
      limits V only along u_c = V_c/|V_c|, and a ~1e-7 tangential
      precursor velocity passes the 1e-8 "negligible speed" threshold
      (Task 13), so the radial gradient is left unlimited. With no limiter
      order 2 runs; with the threshold at 1e-2 (throwaway build) it runs
      and beats order 1 (L1 0.110 and 0.048 at 30 and 60, shock 0.981 and
      0.991). **Decision** (Luca, 2026-11-26): **(a)**, built in Task
      22b. The options were:
  - (a) Fix the SP-V direction rule (e.g. a relative threshold like 1e-2,
        or u_c from the pressure gradient when |V_c| is small against its
        neighbours), re-checking the v0.2 dumps and gate. Smallest change;
        the rule for V_c ~ 0 is already our choice, not the paper's
  - (b) The MM limiter (thesis Sec. 2.4), as Task 16 planned for an SP-V
        failure on Sedov
- [x] Cost (taken up in Task 22b): order 1 at 120 took 5185 steps, 32 min with the machine
      shared (steps x3.2 per doubling: dt is set by the hot core). Order 2
      at 120 would be ~1 h, above the ~15 min budget; the 120 size (or the
      budget) needs revisiting with the decision above
- [x] Risk (thesis Sec. 3.5): at order 1 no node crossing and no
      tangling up to 120 (it runs; mirror symmetric to 7e-14)
- [x] Profile overlay, order 1 at 120: `docs/figures/profiles_sedov_120.png`
      (density, radial u, P, eps against r). Order 2 waits for the fix. The
      thesis's Figs. 3.13-3.15 are the spherical problem, so not overlaid

## Task 22b — SP-V velocity direction for near-rest cells; Sedov at order 2 — ✅ DONE
Decided in Task 22, option (a): the SP-V rule for the direction u_c is
fixed rather than replaced by the MM limiter. Files:
`src/hydro/limiter.{hpp,cpp}`, `tests/unit/test_limiter.cpp`,
`tests/unit/test_driver.cpp`, `scripts/convergence.py` (if needed),
`docs/{decisions,vv_catalog,roadmap}.md`, `README.md`, the Sedov figure.
- [x] The rule: u_c = V_c / |V_c| only when |V_c| > 1e-2 of the largest
      |V_c'| over C(P(c)); otherwise the pressure-gradient direction, else
      e_x (`kNegligibleSpeed` in `limiter.cpp`, was 1e-8). Justified by a
      sweep with the threshold read at run time (throwaway build): Sedov
      30 at order 2 fails at 1e-8 and 1e-6, and every value from 1e-5 to
      3e-1 runs with the same L1 (1.0429e-1, 1.0432e-1 at 3e-1), so 1e-2
      is inside a flat band, not tuned to the gate. Recorded in
      `decisions.md`, "Limiter: SP-V", with the rejected alternative
      (|V_c| against |grad V| h)
- [x] Unit test first: `velocity_limiters: projection direction and the
      V_c = 0 fallback`, section "V_c near rest, tangential to the
      gradient": V_c = 2e-8 e_x next to neighbours at 0.1 with dV_y/dy,
      chi_V = 2/3 (red at 1.0 before). The "V_c != 0" section's V_c =
      1e-3 e_y sat at the new threshold: now 2e-3 (limited along V_c) and
      5e-4 (counts as at rest)
- [x] Regressions, order 2 against the previous build: Noh and
      Taylor-Green (unlimited and SP-V) byte-identical; Sod 1.6e-15 and
      LeBlanc layer 225 3e-13 in density (round-off); Saltzman changed in
      its shock cells only (x in [0.60, 0.83], up to 1e-2 mid-jump, 1122
      steps instead of 1121; the shock enters P = 0 gas on the skewed
      mesh). v0.2 gate passes; Saltzman's order-2 L1 7.041e-2, 3.104e-2,
      1.557e-2 (was 7.030e-2, 3.102e-2, 1.559e-2), Sod, Noh, LeBlanc layer
      unchanged to 4 digits. Sod's mirror columns 1.2e-14 apart. Full
      ctest green (111)
- [x] Sedov at order 2 on 30/60/120: L1 1.043e-1, 4.842e-2, 1.896e-02
      (order 1: 1.664e-1, 1.139e-1, 5.227e-2); shock 0.981, 0.991,
      0.999. The ctest is now `hydro_run: cylindrical Sedov at both
      orders, 30x30` (order 2: shock 0.9807, peak 5.650, L1 0.1043,
      energy to 2e-14, mirror 8e-14). Overlay with both orders:
      `docs/figures/profiles_sedov_120.png`
- [ ] Cost, **open for Luca**: order 2 at 120 took 90.6 min of CPU time (5233 steps; 9.0 h wall clock on the shared desktop, so the CPU time is the measure). Order
      1 at 120 took 32 min in Task 22. Both are above the ~15 min budget,
      so: keep 120 as the finest Sedov mesh (and the budget for it), or
      drop to 30/60/(90)? Not decided here

## Task 23 — Kidder isentropic compression (cylindrical shell)
Files: decks, deck/driver support for the initial profile (as
`initial.taylor_green` was done), the Task 21 pressure law,
`scripts/convergence.py`, tests, `docs/vv_catalog.md`.
- [ ] Initial profile (thesis Eq. 3.10, gamma = 2 in 2D) and boundary
      pressures (Eqs. 3.13-3.14, h^-4), on the polar mesh (Task 20), to
      tau_1 = 0.188746
- [ ] Exact solution in closed form: R(r, t) = h(t) r, rho and P from the
      isentrope; errors on density and on the inner/outer radii
- [ ] Tests: entropy P / rho^gamma stays 1 to a stated bound (isentropic);
      radii within bounds at tau_1; the gate on the cheap resolutions
- [ ] Risk (thesis Sec. 3.6): exponential growth of mesh defects near the
      outer surface. A quarter shell on a polar mesh has no triple point,
      unlike the thesis's 3D mesh, so tau_1 should hold
- [ ] Profile overlay: density scatter against r at tau_1 over the exact
      profile, and the inner/outer radii against t over h(t) r (the
      thesis's Figs. 3.17c-d and 3.19c-d), both orders

## Task 24 — Guderley converging shock (cylindrical)
Files: decks, deck/driver support for an initial radial profile from a
table (e.g. `initial.radial_profile: "<oracle file>"`), oracle data,
`scripts/convergence.py`, tests, `docs/vv_catalog.md`.
- [ ] Setup per Task 16: start at t = 0 (Caramana-Whalen time) from the
      exact profile sampled at cell radii, polar r in [0, 2], outflow at
      r = 2, dumps at t = 0.6 and 0.9, compared on r <= 1. The oracle holds
      the t = 0 profile (initial data) and the two comparison profiles
- [ ] Tests: shock radius against the exact trajectory over the window,
      density behind it within bounds, symmetry spread; the gate on the
      cheap resolutions
- [ ] Risk: the degenerate axis cells (Task 20), and the cold gas ahead
      of the shock (needs `dt_initial`, as Noh)
- [ ] Profile overlay: density, u, P scatter against r at t = 0.6 and
      t = 0.9 over the exact profiles, and the shock radius against t over
      r_s(t), both orders at the finest mesh

## Task 25 — v0.3 acceptance
Files: `scripts/convergence.py` (or its successor `make_report.py`),
`docs/{vv_catalog,roadmap}.md`, `README.md`, figures.
- [ ] All four families in the convergence table and figure, the gate as
      worded in Task 16, rates documented in `vv_catalog.md`
- [ ] `make_report.py` (DoD item 6): one command for every V&V table and
      plot so far (v0.2 shock problems, Taylor-Green, v0.3), the
      expensive runs behind a flag
- [ ] Profile overlays for the v0.2 problems too, regenerated by
      `make_report.py` with the v0.3 ones. Sod, Noh, Saltzman at 200 cells
      already in `vv_catalog.md` (made by hand with `plot_comparison.py`,
      which now overlays several dumps); Taylor-Green still to do
- [ ] v0.4's backlog (p4est, AMR) gets written after this task

## After Task 25
v0.3 is done when the gate agreed in Task 16 is met.

## Deferred: performance (not scheduled)
Luca, 2026-11-24: correctness first; speed waits until the physics is
established. Not before v0.3 is accepted, and only as its own task. From
Task 18 (`NOTES.md`): ~40 us per cell-step at order 2, ~3.4 KB peak memory
per cell (a 100^3 3D mesh would need ~3.4 GB), the nodal solves ~85% of a
step. Candidates, each with dumps bit-identical to before:
- The Newton loop: face areas and per-step quantities computed once, the
  four corner passes per evaluation fused, no allocation per iteration,
  `wall_constraints` once per step instead of per solve
- Fewer temporaries per cell (CornerValues arrays, the predictor's state
  copy), for the memory per cell
- Threads (OpenMP): node assembly as a gather in a fixed order, so the sums
  keep their order
