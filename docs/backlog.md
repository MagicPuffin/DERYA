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

## Task 4 — Core Lagrangian solver (`src/hydro/`)
19-task breakdown, paper Sec. 2 — Vec3/index types → face-area vector (Eq. 5;
already done in Task 2 as `mesh::face_area_vectors`) → corner-area assembly → acoustic impedance Z_c (Eq. 3, needs EOS) → nodal
solver M_p/B assembly (Eq. 4, over all faces incl. boundary faces, then
`bc::apply_pressure_bcs`) → velocity solve (`bc::solve_nodal_velocity`,
already done in Task 3) → momentum update (Eq. 1) →
energy update (Eq. 1) → node position update (Eq. 2) → single-patch
integration test. Each implementation task preceded by its own test task.
(Full table in the original chat planning doc if the detail is needed again
— reconstruct the same granularity here as each sub-task starts.)

## Task 5 — Timestep/CFL controller (`src/timestep/`)
- [ ] CFL-limited dt test on a trivial single-cell case, known analytic value
- [ ] Implement, wired through a (no-op for now) `global_reduce`

## Task 6 — Driver & first real run (`src/app/`, minimal `src/io/`)
- [ ] Minimal JSON deck schema for Sod
- [ ] `hydro_run <deck.json>` wiring Tasks 1-5 into a time-stepping loop
- [ ] Minimal checkpoint output (CSV/JSON dump, full schema deferred)
- [ ] **First integration milestone:** Sod via `hydro_run`, checked against
      `tests/oracle_data/sod.json`

## After Task 6
Extend the same pattern to Noh and Saltzman (v0.1's actual acceptance bar).
v0.2 backlog (second-order MUSCL/SP-V reconstruction, paper Sec. 2.3) gets
written once v0.1 is genuinely done, not before.
