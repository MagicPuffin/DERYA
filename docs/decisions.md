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
  Widened in Task 15 to 1e-13 of the larger of the cell's kinetic energy
  and the mesh's largest |E_c|: a cold cell at rest ahead of a shock gets
  its E_c as a sum of energy fluxes of the shocked flow's size, so it holds
  round-off of that size (Saltzman at second order: eps = -2.1e-14 with no
  kinetic energy, which used to throw).
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

- **Second-order extension (v0.2)**, decided in Task 9. The paper (Sec. 2.3)
  gives the Green-formula gradients (Eq. 6), SP-V (Eq. 7) and
  Barth-Jespersen (Eq. 8) and defers the rest to Georges' thesis [30]
  (HAL tel-01376682, kept locally as `docs/GEORGES_GABRIEL_2016.pdf`, not
  committed; Chapter 2 is in English). From the thesis unless marked "our
  choice":
  - *Reconstruction* (thesis Eq. 2.1): pressure and velocity only,
    `P~_c(x) = P_c + (grad P)_c . (x - x_c)`, `V~_c(x) = V_c + (grad V)_c
    (x - x_c)`, evaluated at the cell's nodes. `rho_c`, `a_c` (so the
    impedance's cell part) stay cell values.
  - *x_c = node average* (thesis Eq. 2.3), not the volume centroid. The
    reconstruction is then not conservative, which does not matter: the
    nodal solver guarantees momentum and energy conservation for any
    corner values.
  - *Gradients* (paper Eq. 6, thesis Eqs. 2.14-2.15) use the nodal values
    `P_cfp`, `V_p` of a *first-order* nodal solve at the same state ("these
    gradients depend on the first order solver"). Our choice where the
    thesis is silent: that solve is done afresh at each stage, not lagged
    from the previous step, so each stage costs two nodal solves. Lagging
    is a possible later optimization, to be measured against this.
  - *Second-order nodal solve* (thesis Eqs. 2.43-2.44): `P_c`, `V_c`
    replaced per corner by `P~_c(x_p)`, `V~_c(x_p)` in `B`, in Eq. 3 and in
    the corner forces. Our choice: the two-shock impedance's velocity jump
    is `s = (V_p - V~_c(x_p)) . n_pf` too (the thesis has the acoustic Z).
    Our choice (Task 12): outflow boundaries keep the owner cell's `P_c` as
    their prescribed pressure `P*`, not `P~_c(x_p)`; the boundary face then
    sees the reconstruction's jump `P~_c(x_p) - P_c`. Irrelevant where the
    state near the boundary is uniform (Sod); revisit if an outflow
    boundary misbehaves at second order.
  - *Limiter: SP-V* as in the paper (Eqs. 7-8): Barth-Jespersen with
    controls on the cells sharing a node with c (`C(P(c))`, c included, no
    ghost cells at boundaries: our choice) and extrapolated values at c's
    own nodes; velocity limited as a whole by one factor, from its
    projection on `V_c/|V_c|`. Our choice for `V_c = 0`: project on the
    pressure-gradient direction, else `e_x` (the thesis SP rule). Changed
    in Task 13 from "exactly 0" to "`|V_c|` <= 1e-8 of the largest `|V_c'|`
    over `C(P(c))`": at rest ahead of a wave, mirror cells of a 1D flow
    hold V_c = 0 in one and a 1e-17 round-off vector in the other, which
    projects on an arbitrary direction. On Sod at second order the factors
    then differed by up to 1 between mirror cells and the columns drifted
    1e-4 apart (v_y 2e-4); 1e-12 left 2e-9, 1e-8 leaves round-off. The
    pressure-gradient fallback has the same exact-zero test; no run has
    needed a tolerance there yet.
  - *Time integration: predictor-corrector* (thesis Sec. 2.6.1, Eqs.
    2.42-2.49, alpha = 1): an Euler predictor to state (1), then the
    corrector averages the predictor's and its own fluxes and nodal
    velocities. Each stage = first-order solve, gradients, SP-V,
    second-order solve, so ~4 nodal solves per step (~4x the first-order
    cost; measured ~4.5x on Sod and Noh). dt is computed once per step,
    before the predictor (the thesis computes it inside the predictor,
    after the nodal solve); the step hands the timestep controller the
    mean nodal velocity the mesh moved with. The thesis's
    GRP alternative (Sec. 2.6.2) is cheaper and less dissipative but needs
    time derivatives of the nodal velocities and face vectors; deferred
    until cost matters.
  - *Timestep:* the existing controller unchanged. The thesis (Sec. 1.5)
    also has an internal-energy positivity criterion; added only if a
    second-order run needs it.
  - *Smooth convergence problem:* the Taylor-Green vortex (thesis Sec. 3.2,
    after Dobrev et al.): stationary, smooth, exact, and the thesis's Table
    3.1(b) gives reference errors for this exact scheme (PC, Green
    gradients, unlimited: L1 order 1.98 at 320x320). It needs a prescribed
    energy source term (thesis Eq. 3.6, Appendix B). Rejected: a 1D simple
    wave (needs a characteristics solver for the oracle and has no
    reference numbers); Kidder (stays in v0.3).
    Settled in Task 14, our choices where the thesis is silent: the source
    adds `V_c S(x_c)` (x_c the node average, the cell's current volume and
    position) to `m_c dE_c/dt`, at each predictor-corrector stage like the
    fluxes (the predictor's at s^n, the corrector the mean of s^n's and
    s^(1)'s; Appendix B only adapts the GRP step); at first order it is
    forward Euler at t^n. The initial state and the exact pressure are
    point values at each cell's node average, the exact one at the cell's
    *current* position (the profile is stationary in the Eulerian sense,
    the cells move). Domain [0,1]^3, one cell in z, as the thesis. Eq. 3.5
    is `rho0 C1^2 / 4 (cos 2 pi x + cos 2 pi y) + C2`, so (1, 1, 1) gives
    P in [0.5, 1.5]; a deck must keep C2 > rho0 C1^2 / 2. Table 3.1(b)'s
    20x20 E_L1 "2.11e-3" is a typo for 2.11e-2 (its order, 0.82, needs it).
    Our errors are within 25% of the thesis's at every mesh to 160x160,
    and the orders match (L1 1.95 vs 1.94 at 160x160).
  - *Known risk:* the thesis reports its SP limiter failing on Noh and
    Sedov (non-monotone, oscillating), where only its MM limiter (Sec. 2.4)
    works. The paper's SP-V reportedly fixes this. If Noh at second order
    oscillates, the MM limiter is the fallback, decided then. Checked in
    Task 15: planar Noh at second order with SP-V does not oscillate (the
    peak rho 4.81 is the wall-heating start-up error, 2nd-3rd cell from the
    wall), so no MM limiter for now; Sedov (v0.3) is the next check.

- **v0.3 test problems**, decided in Task 16. All four run as the
  *cylindrical* (or planar) problem on a thin slab, one cell in z, with
  symmetry walls on z; spherical and 3D runs are v0.5. The paper's own
  Sedov (Sec. 6.4) is spherical on an eighth-sphere AMR mesh, so it is a
  v0.4/v0.5 parity check, not v0.3. Constants below were checked against
  ExactPack (`.venv`) or derived here where marked.
  - *LeBlanc* (planar shock tube): gamma = 5/3, (rho, eps) = (1, 0.1) for
    x < 3 and (1e-3, 1e-7) beyond, so P = 2/3 x 1e-1 and 2/3 x 1e-10,
    x in [0, 9], t = 6, at rest, walls at both ends (no wave reaches them
    by t = 6). The Riemann solution, from ExactPack's `IGEOS_Solver` and
    from `convergence.py`'s `sod_exact` with these states (they agree to
    1e-11): p* = 5.156e-4, u* = 0.6218, rarefaction head at x = 3 - 6 a_L
    = 1 (a_L = 1/3), contact at x = 6.731, shock at x = 7.975 with
    rho = 4e-3 (the strong-shock limit (gamma+1)/(gamma-1) x 1e-3, to 5e-7)
    and eps = 0.1933 behind it. Solved in place like Sod in
    `convergence.py`; Task 19 also wrote an oracle file
    (`tests/oracle_data/leblanc.json`, ExactPack, every 0.005 with the two
    jumps) for the ctest and the profile plots, and the two agree to
    1.1e-11. Mesh n x 1 x 1, cubic cells.
    *Gated setup changed* (Task 19b, agreed with Luca 2026-11-26): the
    uniform mesh leaves the shock 0.3-0.4 ahead at 900 cells from
    start-up entropy in the few heavy cells at the contact (`vv_catalog.md`,
    "Why the shock is ahead"), so the gate runs on meshes with [2.96, 3]
    refined x64 (Task 19c; [2.95, 3] does not fit the 225 mesh); the
    uniform family stays, reported.
  - *Sedov* (cylindrical): Cartesian quarter plane [0, 1.2]^2, one cell in
    z, symmetry walls on all sides (Maire et al.'s 2D setup and the
    thesis's cube cut to 2D), gamma = 1.4, rho = 1, P = 1e-6, at rest; the
    blast energy goes into the one cell touching the origin as internal
    energy. ExactPack's `CylindricalSedov` default E0 = 0.311357 is per
    unit length *for the full plane* and puts the shock at r = 0.75 at
    t = 1, not 1 (its docstring; checked). For r = 1 at t = 1 the full
    plane needs E = alpha rho = 0.984074 (ExactPack's alpha for gamma =
    1.4), so the quarter plane gets **E = 0.246019 per unit length**. The
    literature's 0.244816 (Maire et al.) puts the shock at r = 0.9988.
    The deck gives the energy per unit length in z; the driver multiplies
    by the slab thickness, so the decks differ only in cell counts. Final
    time t = 1, peak density 6 at the shock (the strong-shock limit).
    Polar mesh rejected for Sedov: it would need degenerate cells at the
    origin, where the energy goes, and the Cartesian mesh is the standard
    2D test of grid imprinting (the thesis's node crossing on the axes).
  - *Kidder* (cylindrical shell): quarter shell r in [0.9, 1], theta in
    [0, 90] degrees, polar mesh (Task 20), symmetry walls on the planes
    theta = 0 and 90, the exact pressure (Task 21) on the inner and outer
    surfaces: P_b(t) = P0(r_b) h(t)^(-2 gamma / (gamma - 1)) (thesis Eqs.
    3.13-3.14), with P0(0.9) = 1, P0(1) = 2^gamma. Initial state thesis
    Eq. 3.10 with rho_inn = 1, rho_out = 2, P = rho^gamma, at rest.
    *Derived here:* Kidder's homogeneous isentropic compression R = h(t) r
    in d dimensions needs h'' h^(d (gamma - 1) + 1) = const, and the
    thesis's h = sqrt(1 - (t / t_foc)^2) has h'' h^3 = -1 / t_foc^2, so
    **gamma = 1 + 2/d = 2 in 2D** (5/3 in 3D, the thesis's). Eq. 3.10 and
    Eq. 3.12 for t_foc hold for any d; density scales as h^-d and pressure
    as h^(-d gamma), which equals the boundary law's h^(-2 gamma/(gamma-1))
    only at that gamma. With gamma = 2: a_inn^2 = 2, a_out^2 = 4, t_foc =
    sqrt(0.19 / 4) = 0.217945 (Eq. 3.12 at gamma = 5/3 gives the thesis's
    0.254345, which checks the formula). Final time **tau_1 = (sqrt 3 / 2)
    t_foc = 0.188746**, h = 1/2, radii 0.45 and 0.5. tau_2 = 0.99 t_foc is
    reported if it runs, not gated. Not checked against Maire et al. 2007
    (the thesis's [93], not available here); the derivation is the check.
    Closed form evaluated in place, no oracle file.
  - *Guderley* (cylindrical): ExactPack's `Guderley(geometry=2)` covers it
    (Lazarus's Table 6.4, lambda = 1.197141 at gamma = 1.4). Its time is
    Caramana-Whalen's, t = 0.750024322 (t_L + 1): the converging shock is
    at r = 1 at t = 0 (rho = 6, u = -0.696, P = 0.581 just behind it) and
    reaches the axis at **t_c = 0.750024**, at r_s(t) = (1 - t /
    t_c)^(1/lambda) before that (0.2607 at t = 0.6). Ahead of it the gas
    is cold, rho = 1, at rest; behind it the inflow extends to infinity
    (rho 10.3, u -0.53 at r = 2). Setup: gamma = 1.4, start at t = 0 from
    the exact solution sampled at each cell's radius (a radial profile
    table from the oracle, `initial.radial_profile`, Task 24), **polar
    mesh r in [0, 2]**, theta in [0, 90], symmetry walls on the planes,
    outflow at r = 2. The outer boundary is not prescribed from the exact
    solution: that needs the boundary's particle path, and ExactPack's
    solver takes ~1 min per time level. Instead the domain is big enough:
    a disturbance from r = 2 moves in at u - c ~ -0.83, so by t = 0.9 it
    is near r = 1.25, outside the compared window **r <= 1**. Compared at
    **t = 0.6** (converging shock at 0.2607, before the collapse) and
    **t = 0.9** (after the reflection, rho ~ 58 near the axis), one run
    to t = 0.9 with an intermediate dump. The polar mesh down to r = 0 has
    degenerate cells on the axis (the x_min side collapses to a point, its
    faces have zero area), as in the classic polar-grid Noh. Today
    `wall_constraints` (`src/bc/bc.cpp`) normalizes each face normal with
    no zero check, so Task 20 must skip zero-area faces there and in the
    corner geometry, and test it. The axis node lies on both walls, so it
    is pinned, which is the symmetry the problem needs. If degenerate cells
    turn out not to work cleanly in Task 20, the fallback is a Cartesian
    quarter plane [0, 2]^2 (grid imprinting at the focus), decided there.
  - *Meshes and walls:* Cartesian for LeBlanc and Sedov, polar for Kidder
    and Guderley. The polar mapping takes the structured box's x range as
    r and its y range as theta in **degrees** (0 and 90 are exact in a
    deck), nodes at (r cos theta, r sin theta, z); sides keep their tags
    (x_min inner surface, x_max outer, y_min the theta = 0 plane, y_max
    the other). A curved boundary is only ever a pressure boundary
    (outflow or prescribed P*), never a wall: `wall_normals` counts
    non-coplanar faces at a node as separate walls, so a curved wall would
    pin its nodes like an edge. Walls only on planes.
  - *Error metric:* L1 of density against the exact profile at each
    cell's radius, r of its node average (the dump's `centroid`, as
    v0.2), `sum_c |rho_c - rho(r_c, t)| V_c / L_z` (per unit length of
    slab; for LeBlanc per unit cross-section, as Sod). Kidder adds the
    inner and outer radii (mean radius of the nodes on x_min, x_max)
    against h(t) 0.9 and h(t), and the entropy max_c |P_c / rho_c^gamma
    - 1|. Symmetry, for the radial problems: per radial bin of one cell's
    width, the spread (max - min) of cell density over the cells in it;
    reported as the largest spread behind the shock, and plotted as the
    thesis's scatter plots.
  - *Oracles:* closed forms (LeBlanc's Riemann solution, Kidder) are
    evaluated where used, as Sod is today (LeBlanc also has an oracle
    file, as Sod does). Sedov and Guderley come from
    ExactPack via `scripts/gen_oracle_data.py` as fine-sampled radial
    profiles with their jumps in `tests/oracle_data/`, so neither ctest
    nor `convergence.py` needs ExactPack at run time (ExactPack's Guderley
    takes minutes per time level).
  - *Known risks*, to report, not tune: the thesis's SP limiter kept the
    node crossing on the axes and the bowtie hourglass at 45 degrees on
    Sedov where only its MM limiter cleared them (seen at first order too,
    in 3D; the thesis says 2D EUCCLHYD does not show them); Kidder's
    exponential growth of mesh defects near the outer surface, where the
    thesis's SP limiter gave inner-face oscillations at second order (a
    2D polar mesh is equi-angular, no triple point, so tau_1 should hold);
    LeBlanc's 1e9 pressure ratio and cold right state (the negative-eps
    tolerance, the Newton flux scale, SP-V near zero velocity) and its
    known eps overshoot at the contact. If SP-V fails on Sedov or Kidder,
    the MM limiter (thesis Sec. 2.4) becomes its own task; the gate is not
    loosened to pass.

- **Build trees:** `build/` is optimized with debug symbols
  (`RelWithDebInfo`, the default build type and the `default` preset) and is
  where tests run; `build-debug/` (`Debug`, the `debug` preset) is only for
  stepping through code. The integration tests are ~12x slower unoptimized
  (Saltzman: 2.5 min vs 11 s). Results are bit-identical between the two,
  since neither uses `-ffast-math` or FMA contraction (no `-march`); keep it
  that way, or expect last-bit differences the test tolerances must absorb.
  Decided after Task 8.

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