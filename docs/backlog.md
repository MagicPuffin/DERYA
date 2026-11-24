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
- [ ] (Kidder's boundary waits for its own task batch, now Task 21: a
      prescribed, time-dependent pressure, an extension of the outflow BC)

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
      took 78 min optimized, ~24x 80x80)

      Order 1: E_L1 1.37e-1, 9.09e-2, 5.29e-2, 2.86e-2 (orders 0.59, 0.78,
      0.89)
- [x] ctest `hydro_run: Taylor-Green vortex converges at order ~2`: 10x10,
      20x20, 40x40 at both orders, in-process; errors decrease, order 2
      below order 1 everywhere, the order-2 rate rises (> 1.2 at 40x40,
      0.3 above order 1's), order 1's in [0.6, 1.2], order 2 within a factor
      1.5 of the thesis. ~26 s optimized. Red without the source: the 40x40
      error grows (rate -0.9). The full table is
      `scripts/taylor_green_convergence.py`, on demand (80x80 ~3 min,
      160x160 78 min)
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
- [ ] Not done (optional): the `--probe cell_id` flag; nothing in this task
      needed it. Still a `roadmap.md` design gap

# v0.3 — Sedov, Kidder, Guderley, LeBlanc (paper Sec. 6.4, thesis Ch. 3)

Gate (`roadmap.md`): "all pass, convergence documented". As in v0.2, each
problem runs at 3 resolutions at both orders through
`scripts/convergence.py`, and Task 16 turns "pass" into a bar per problem.

Scope, proposed here and settled in Task 16: the **cylindrical** versions
in 2D (thin slab, one cell in z), the DoD's "2D at 3+ resolutions". The
thesis runs Sedov on a 3D cube and Kidder on an eighth of a spherical
shell; spherical and genuinely 3D runs belong to v0.5 (DoD item 2), and at
today's cost (Taylor-Green 160x160 at order 2: 78 min) a 3D cube at 3
resolutions would take hours. Every task keeps the v0.2 decks' results
exactly (`order = 1` and `order = 2` dumps of Sod, Noh, Saltzman and
Taylor-Green byte-identical), unless a task says otherwise and why.

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

## Task 17 — `--probe` debug flag
Files: `src/app/{main,driver}.{hpp,cpp}`, `docs/deck_format.md`,
`tests/unit/test_driver.cpp`. The `roadmap.md` design gap.
- [ ] `hydro_run deck.json [dump.json] --probe cell_id` (repeatable; maybe
      `--probe-node node_id`): one line per step per probed cell (step, t,
      dt, x_c, rho, V, eps, P, a) to a file or stderr, JSON lines or CSV
- [ ] Tests: the probed values equal the dump's at the last step; no probe
      changes nothing (dump byte-identical); a bad id is rejected

## Task 18 — Cost scaling check
Files: none required (a report in `NOTES.md` and this backlog); any fix
goes to its own task. Measured in Task 14: Taylor-Green at order 2, 40x40
-> 80x80 cost 8.3x (4x cells, 2x steps: expected 8x), 80x80 -> 160x160
cost 24x. Part of that run shared the machine, so measure again first.
- [ ] Time one step (both orders) at 40x40, 80x80, 160x160 and 320x320
      Taylor-Green, and the parts of it (geometry, thermo, nodal solves and
      Newton iteration counts, gradients, limiter, corner forces): is any
      part superlinear in the number of cells, or does Newton take more
      iterations on finer meshes?
- [ ] If something is superlinear, a follow-up task to fix it, with results
      bit-identical. Not a general optimization pass

## Task 19 — LeBlanc shock tube
Files: `decks/leblanc.json` and `decks/convergence/leblanc_*.json`,
`scripts/gen_oracle_data.py` (+ `tests/oracle_data/leblanc.json`) or an
exact solver in `scripts/convergence.py`, `tests/unit/test_driver.cpp`,
`docs/vv_catalog.md`. No new physics: the first v0.3 problem because it
only needs a deck.
- [ ] Decks at the Task 16 resolutions (cubic slab cells), both orders run
      by the script
- [ ] Exact solution (ExactPack `riemann`), checked against an independent
      evaluation; `convergence.py` family `leblanc`
- [ ] Tests: runs to t = 6 at both orders; shock position and the plateau
      density within Task 16's bounds; the gate on the cheap resolutions
- [ ] Risks: the 1e9 pressure ratio and the cold right state (the
      negative-eps tolerance, the Newton flux scale, SP-V near zero
      velocity); the contact's eps overshoot. Report, don't tune

## Task 20 — Polar mesh mapping
Files: `src/mesh/mesh.{hpp,cpp}`, `src/io/deck.{hpp,cpp}`,
`src/app/driver.cpp` (if the mapping needs it), `docs/deck_format.md`,
`tests/unit/test_{mesh,bc,driver}.cpp`.
- [ ] `mesh.mapping: "polar"`: the structured mesh's x range is r, its y
      range is theta (unit to decide in Task 16), z unchanged; nodes at
      (r cos theta, r sin theta, z). Sides keep their tags (x_min = inner
      surface, x_max = outer, y_min/y_max the two radial planes)
- [ ] Tests: total volume = the annular sector's (to the polygonal
      approximation), face vectors close per cell (GCL), wall normals on
      the theta = 0 and theta = pi/2 planes are the plane normals, a
      uniform state at rest stays at rest with walls and outflow
- [ ] r_min = 0 (Guderley): the axis cells are degenerate, their x_min
      faces have zero area. Zero-area faces skipped in the wall
      constraints (`src/bc/bc.cpp` normalizes without a check today) and
      wherever a face normal is normalized; the axis node is pinned by the
      two walls. Tests: no NaN on a mesh with r_min = 0, volumes positive,
      a uniform state at rest stays at rest. If this cannot be made clean,
      Guderley falls back to a Cartesian quarter plane (`decisions.md`)

## Task 21 — Prescribed time-dependent pressure boundary
Files: `src/bc/bc.{hpp,cpp}`, `src/hydro/{hydro,scheme}.{hpp,cpp}`,
`src/app/driver.cpp`, `src/io/deck.{hpp,cpp}`, `docs/deck_format.md`,
tests.
- [ ] A `Pressure` boundary type with P*(t) given (a constant, or a named
      law such as Kidder's Eqs. 3.13-3.14), where today's outflow has P* =
      P_c. The steps then need the time: each stage at its own time
      (predictor at t^n, corrector stage at t^n + dt; check against thesis
      Sec. 2.6.1), decided in the task and recorded
- [ ] Tests: a constant P* equal to the gas pressure keeps a state at rest;
      a constant P* above it pushes the boundary inward at the acoustic
      rate (first step, by hand); P*(t) evaluated at the right stage times
      (red when the corrector uses t^n)

## Task 22 — Sedov (cylindrical)
Files: decks, `src/io/deck.{hpp,cpp}` and `src/app/driver.cpp` (the
energy-spike initialization), `scripts/{gen_oracle_data,convergence}.py`,
oracle data, tests, `docs/vv_catalog.md`.
- [ ] Initialization: total energy E deposited in the cell(s) touching the
      origin (a deck key, e.g. `initial.energy_spike: {"at": [0,0,0],
      "energy_per_length": E}`; the driver multiplies by the slab
      thickness), the rest from the background. E = 0.246019 for the
      quarter plane (shock at r = 1 at t = 1; `decisions.md`)
- [ ] Oracle: ExactPack `Sedov(geometry=2, eblast=4 x 0.246019)` (not
      `CylindricalSedov`'s default energy, whose shock is at 0.75) at
      t = 1, fine-sampled, with
      the shock as a jump; radial L1 metric in `convergence.py`
- [ ] Tests: total energy conserved to round-off; shock radius and peak
      density within Task 16's bounds; symmetry spread; the gate on the
      cheap resolutions
- [ ] Risk (thesis Sec. 3.5): node crossing along the axes and non-convex
      cells at 45 degrees. If SP-V does not hold, the MM limiter (thesis
      Sec. 2.4) becomes its own task

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

## Task 25 — v0.3 acceptance
Files: `scripts/convergence.py` (or its successor `make_report.py`),
`docs/{vv_catalog,roadmap}.md`, `README.md`, figures.
- [ ] All four families in the convergence table and figure, the gate as
      worded in Task 16, rates documented in `vv_catalog.md`
- [ ] `make_report.py` (DoD item 6): one command for every V&V table and
      plot so far (v0.2 shock problems, Taylor-Green, v0.3), the
      expensive runs behind a flag
- [ ] v0.4's backlog (p4est, AMR) gets written after this task

## After Task 25
v0.3 is done when the gate agreed in Task 16 is met.
