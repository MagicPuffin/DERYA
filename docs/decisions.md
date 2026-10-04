# Locked Technical Decisions

Read this before starting any task — these are settled, not open for
re-litigation mid-task. If a task seems to require violating one of these,
stop and report rather than deciding to override it.

## Architecture
- **Language/build:** C++17, CMake with `FetchContent` (no vendored deps, no
  system-package assumptions beyond the compiler/CMake itself).
- **AMR library:** p4est — but NOT a dependency until v0.4. v0.1-v0.3 use a
  minimal hand-rolled structured mesh (`src/mesh/`).
- **Parallelism:** built against p4est/libsc's real parallel API from day one,
  compiled with `--enable-mpi=no` (a genuine libsc build mode, single-rank
  stub, not a hack) — so real parallelism is a build-flag flip later, not a
  rewrite. `sync_ghosts`/`global_reduce` will be the sole MPI-aware layer
  above the p4est wrapper, once p4est lands at v0.4.
- **Field layout:** Structure-of-Arrays, local-index-plus-ghost-flag
  addressing throughout (the ghost half is a no-op until v0.4, but the
  addressing scheme is designed in from the start).
- **Module boundaries are enforced by the build system**, not convention —
  each `src/<module>/` is its own CMake target; a module linking against
  something it shouldn't is a build failure.
- **Mesh is 3D-native from Task 2 onward** — hexahedral cells, `(i,j,k)`
  indexing, `nz` defaults to 1. There is no separate 2D mesh implementation.
  "2D" test setups (Sod, etc.) run as a thin extrusion (`nz=1`, later
  possibly `nz=2` for symmetry checks) with symmetry/periodic BCs on the two
  thin faces — this matches the reference paper's own approach exactly
  (Sec. 6.2: their "Sod test" runs on a `4×4×200` 3D slab, not a true 1D
  mesh; their AMR implementation is "purely 3D" throughout). Caught and
  corrected after the first mesh task was built genuinely 2D (quad cells,
  `(i,j)` indexing) — that version is being generalized before anything
  downstream (face-area vectors, the Lagrangian solver) depends on it.
- **Boundary conditions act on the nodal solver**, not via ghost cells or
  cell fields (EUCCLHYD, paper refs [29, 30]; the paper itself doesn't spell
  them out). Each domain side gets one type. Symmetry/wall = `V_p . n = 0`,
  enforced by a Lagrange-multiplier solve of Eq. 4 (handles edges and corners
  where 2-3 walls meet, e.g. thin-slab z faces + Sod end walls). Outflow =
  prescribed pressure `P*` = owner cell pressure, added to `B`. Eq. 4's
  `M_p`/`B` are assembled over all faces including boundary faces. `bc`
  depends on `mesh` only. Decided at the start of Task 3, replacing the
  backlog's placeholder `apply(mesh, fields)` interface. Piston = a wall
  moving with velocity `u_w`, `V_p . n = u_w . n` (Task 8): the nodal solver
  starts from a guess with that normal part and takes Newton steps tangent
  to the walls, so the wall solve itself is unchanged.
- **Two-shock (Dukowicz) impedance in the nodal solver**, not the paper's
  acoustic `Z_c = rho_c a_c`: `Z_cfp = rho_c (a_c + Gamma_c |(V_p - V_c) . n_pf|)`
  per cell, face and node, with `Gamma = (gamma+1)/2` for an ideal gas (from
  `EOS::shock_coefficient`). The acoustic Z vanishes in a cold gas (Noh,
  Saltzman) and makes `M_p` singular. A floor `Z >= rho a_min` was rejected:
  it has units, and since Z is the scheme's only dissipation, a small floor
  under-heats strong shocks while a large one turns into a tuned artificial
  viscosity. The two-shock Z tends to the Rankine-Hugoniot shock impedance
  `rho D` for strong shocks and needs no parameters. Conservation and
  entropy production hold for any `Z >= 0`.
  Since `Z_cfp` depends only on its own node's `V_p`, Eq. 4 becomes a
  per-node nonlinear system `sum S [rho (a + Gamma|s|) s - P] n = 0`
  (`s = (V_p - V_c) . n`), the gradient of a convex function. It is solved by
  Newton's method, Jacobian `sum S rho (a + 2 Gamma|s|) n x n`, starting from
  the mean of the surrounding cell velocities. Plain fixed-point iteration
  on Z was rejected: in the cold limit its slope at the solution is
  `-(1+r)^2/4` for a density ratio `r^2`, so it does not converge. Where every
  surrounding cell is cold and has the same velocity, `M_p` is singular. The
  forces there do not depend on `V_p`, and the Newton step uses a
  `1e-8 tr(J)` regularization, so V_p keeps its initial value. The corner
  forces use Z evaluated at the converged `V_p`; conservation is then exact
  up to the Newton residual. A cold, moving gas also needs
  `eps = E - |V|^2/2` clamped to 0 when it is negative by round-off (it
  throws beyond 1e-13 of the kinetic energy). Decided before Task 5 (Task 4b).
  A Jacobian whose trace is below `DBL_MIN / 1e-8` also counts as singular,
  because its regularization would underflow (subnormal precursors ahead of
  the Saltzman shock gave 0/0; Task 8).
- **Lagrangian timestep** (the paper gives none; EUCCLHYD practice, after
  Maire et al.): `dt^{n+1} = min(C_cfl min_c lambda_c / a_c,
  C_V min_c V_c / |dV_c/dt|, C_M dt^n, dt_max)`, with `lambda_c = V_c /
  (largest face area)` (the shortest edge of a box cell) and `dV_c/dt =
  sum_p n_cp . V_p` (exact by the GCL) from the step just taken. Defaults
  `C_cfl = 0.25`, `C_V = 0.1`, `C_M = 1.01`, `dt_max = inf`. The volume
  criterion is what bounds dt in a cold gas (a = 0); a cold cell in
  uniform motion has no limit. The first step has no V_p, so the volume
  criterion is skipped and `dt^0 = min(dt_initial, CFL)`. The global min goes
  through `sync::global_reduce` (a no-op until v0.4). `lambda_c` uses every
  direction, so a thin-slab cell thinner than its in-plane size limits dt
  more than needed; slab tests keep the extrusion at least one cell wide.
  Decided at the start of Task 5.

## Testing
- **Tier 1 — unit tests (Catch2):** one function/equation at a time, test
  written before or alongside implementation, never after.
- **Tier 2 — oracle regression tests:** exact solutions generated by
  ExactPack (`pip install git+https://github.com/lanl/ExactPack.git` — NOT on
  PyPI under `exactpack`) and the Kamm-Timmes Sedov solver, checked into
  `tests/oracle_data/`, regenerated via `scripts/gen_oracle_data.py`.
- **Tier 3 — parallel structural smoke tests:** deferred until v0.4 (no
  parallelism to smoke-test before then).
- Deck format: `nlohmann::json`.

## Scope
- **Hydro operator only.** No multi-material/ALE remap (Tier B), no MPI, no
  radiation, no drive coupling, no burn physics for v1.0. See `roadmap.md`
  for the explicit out-of-scope list and why.
- **Personal project**, not a First Light Fusion deliverable or publication
  target — no formal review process, but the engineering discipline here
  (atomic tasks, test-first, locked decisions) is deliberately as strict as
  if it were.
- Priority V&V tests beyond the reference paper's own suite: **Kidder ball**
  and **Guderley converging shock** — higher validation-relevance to
  cylindrical hybrid-ignition physics than Sod/triple-point.

## Reference material
- Colaïtis, Guisset & Breil, "A cell-centered AMR-ALE framework for 3D
  multi-material hydrodynamics, Part I" (SSRN preprint 5167115) — motivating
  paper; unpublished/non-peer-reviewed, DOE FES-funded. Cite properly in code
  comments where an equation is implemented directly from it.
- CINDER / `simple_model` (separate repo) — cross-validation anchor. Canonical
  target: 150 µg DT, 100 g/cm³, 40×40 µm cylinder, 10 µm gold sidewall.