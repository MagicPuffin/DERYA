# Roadmap & Definition of Done

## v1.0 Definition of Done
1. All hydro-operator V&V tests (`vv_catalog.md`) pass in 2D at 3+
   resolutions with documented convergence rates matching formal order,
   first and second order.
2. At least one test (Guderley or Kidder) also verified in 3D on cylindrical
   geometry specifically.
3. Lagrange-AMR and Euler-AMR both work in 2D, reproducing the reference
   paper's own Sod/Sedov/triple-point-adjacent results.
4. **Science demo:** CINDER's canonical 150 µg DT point design run as 2D
   cylindrical hydro-only — does the resolved trajectory diverge from
   CINDER's 0D assumption the way CINDER's own qualification report predicts
   it should? This is the actual point of the project; schedule it, don't
   let it happen "incidentally."
5. Checkpoint/restart exercised at least once for a real reason.
6. One command regenerates every V&V plot/table (`scripts/make_report.py` —
   not yet started, build incrementally as tests land).

**v1.0 explicitly excludes:** multi-material/ALE remap, real MPI, radiation,
drive coupling, burn physics, any formal documentation/review process.

## Milestones
| Version | Contains | Gate |
|---|---|---|
| v0.1 | First-order Lagrangian hydro (3D-native mesh, run as thin-slab/"2D" test setups), no AMR, no p4est | Sod/Noh/Saltzman pass at 1 resolution |
| v0.2 | + second-order (MUSCL/SP-V, paper Sec. 2.3) | Sod/Noh/Saltzman at 3 resolutions, second-order errors below first-order ones, rates documented; order ~2 on one smooth problem |
| v0.3 | + Sedov, Kidder, Guderley, LeBlanc (cylindrical, thin slab) | all pass, convergence documented; "pass" per problem below (proposed in Task 16, not yet agreed) |
| v0.4 | + p4est, Lagrange-AMR, Euler-AMR, thin-slab | paper's Sod-AMR/triple-point-AMR parity |
| v0.5 | + genuinely 3D problems, cylindrical geometry | DoD item 2 |
| v1.0 | + science demo, restart, report generation | DoD fully met |

## v0.3 gate (proposed in Task 16, awaiting Luca's agreement)
LeBlanc, Sedov, Kidder and Guderley (setups in `decisions.md`, "v0.3 test
problems") at 3 resolutions, both orders (SP-V), through
`scripts/convergence.py`:
1. Every run reaches its final time.
2. The L1 density error falls with resolution at both orders, and the
   rates are documented. Rates <= 1 are expected on the three shocked
   problems; Kidder is smooth.
3. Order 2 has a smaller error than order 1 at every resolution.
4. At the finest resolution, both orders:
   - *LeBlanc*: shock position (the last crossing of rho = 2.5e-3) within
     0.1 of 7.975.
   - *Sedov*: shock radius (radius of the densest cell) within 0.05 of 1.
     The symmetry spread is reported, not gated.
   - *Kidder*: inner and outer radii at tau_1 within 1% of 0.45 and 0.5;
     max |P / rho^gamma - 1| <= 1e-2; the order-2 rate above the order-1
     rate.
   - *Guderley*: converging-shock radius at t = 0.6 within 0.02 of 0.2607;
     L1 on r <= 1 at t = 0.6 and t = 0.9.

If SP-V cannot meet this on Sedov or Kidder, that is a finding that leads
to the MM limiter (thesis Sec. 2.4) as its own task. The gate is not
loosened to pass.

## Known design gaps (not yet built, tracked so they aren't forgotten)
- **Second-order rollout** — done (Tasks 9-15 in `backlog.md`; v0.2 gate
  met, results in `vv_catalog.md`).
  The details the paper defers to Georges' thesis [30] were settled from
  the thesis in Task 9 (`decisions.md`, "Second-order extension"). The
  v0.2 gate was
  reworded (agreed 2026-10-25): rates on shocked problems are <= 1 for any
  scheme, so order 2 is shown on a smooth problem instead.
- **Debug/probe tooling** — done (`backlog.md` Task 17): `hydro_run
  --probe cell_id` and `--probe-node node_id` trace a cell's or node's full
  state at every step (`deck_format.md`, "Probes"). The reference paper's
  own authors needed exactly this to catch their triple-point oscillation bug.
- **Environment reproducibility** — `docs/SETUP.md` covers this incrementally;
  keep it current as new gotchas are hit (already has 4).
- **Checkpoint/restart schema** — deliberately minimal (CSV/JSON dump) until
  v1.0's DoD item 5 actually needs it; don't over-build early.
- **Units/normalization convention** — not yet decided against CINDER's;
  needs resolving before the v1.0 science demo, not before.
- **Numerical failure policy** — negative density/pressure, EOS
  out-of-domain, NaN propagation: no policy chosen yet. Decide deliberately
  when it's first hit rather than leaving to whatever the first crash does.
- **p4est Part II gap (rezoning/disentangling)** — NOT blocking v1.0 (Tier B
  is out of scope). Revisit only if/when a Tier B decision is made later.
  Author outreach (emailing Colaïtis/Guisset/Breil) still an open,
  low-urgency action item.

## Resolved gaps (kept here as history, not open items)
- **Mesh dimensionality** — the mesh module was initially built genuinely 2D
  (quad cells). Caught before the Lagrangian solver (Task 4, which needs the
  3D face-area-vector formula from paper Eq. 5) depended on it. Generalized
  to 3D-native hexahedral, `nz=1` as the default "2D" run mode — see
  `decisions.md`. No other module was affected since nothing downstream of
  mesh existed yet.