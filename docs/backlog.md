# v0.1 Backlog

v0.1 acceptance bar (from `roadmap.md`): Sod, Noh, Saltzman pass at 1
resolution, no AMR, no p4est. Update the checkboxes below as tasks land —
this file is the source of truth for "what's next," not the chat history.

## Task 0 — Repo & toolchain bootstrap — ✅ DONE
- [x] CMake skeleton, Catch2 v3.7.1 + nlohmann::json v3.11.3 via FetchContent
- [x] Toolchain smoke test green
- [x] ExactPack installed from GitHub (not on PyPI), `matplotlib` added to
      `requirements.txt` (ep_riemann imports pyplot at module load)
- [x] `scripts/gen_oracle_data.py` → `tests/oracle_data/sod.json`, verified
      reproducible bit-for-bit on Luca's machine via SSH
- [x] Pushed to GitHub (`DERYA`)

## Task 1 — EOS module (`src/eos/`) — ✅ DONE
- [x] `EOS` interface: `pressure(rho, eps)`, `sound_speed(rho, eps)`
- [x] `IdealGasEOS`, header-only, in `src/eos/eos.hpp`
- [x] Test passes: p=(γ-1)ρε=1.0, a=√1.4 for ρ=1.0/ε=2.5/γ=1.4 (note: the
      original task prompt had a wrong expected value, 1.5 — correctly
      caught and flagged rather than "fixed" by editing the test to match)

## Task 2 — Mesh module (`src/mesh/`) — IN PROGRESS
- [x] `Mesh` type + test for a structured mesh: node/cell/face counts and
      connectivity (interface + test only). First built as a 2D quad mesh,
      then generalized to 3D-native hexahedra per `decisions.md` ("Mesh is
      3D-native"): `(i,j,k)` indexing, `nz` defaults to 1. Test asserts the
      2×2×1 case (18 nodes, 4 cells, 20 faces, cell 0 corners
      {0,1,4,3,9,10,13,12}) plus a 3×2×2 case and face-orientation checks
- [x] Implement `generate_structured_mesh()` to pass the above
- [ ] Test: face-area vector (paper Eq. 5) on one known cell of the generated
      mesh
- [ ] Implement it

## Task 3 — Boundary conditions (`src/bc/`)
- [ ] `BoundaryCondition::apply(mesh, fields)` interface
- [ ] Reflecting/symmetry BC (test + implement) — needed for Noh
- [ ] Outflow BC (test + implement) — needed for Sod
- [ ] (Kidder's prescribed-motion BC waits for its own task batch)

## Task 4 — Core Lagrangian solver (`src/hydro/`)
19-task breakdown, paper Sec. 2 — Vec3/index types → face-area vector (Eq. 5)
→ corner-area assembly → acoustic impedance Z_c (Eq. 3, needs EOS) → nodal
solver M_p/B assembly (Eq. 4) → velocity solve → momentum update (Eq. 1) →
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
