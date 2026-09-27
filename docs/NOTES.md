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
