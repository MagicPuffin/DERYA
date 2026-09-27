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
| v0.1 | First-order 2D Lagrangian hydro, no AMR, no p4est | Sod/Noh/Saltzman pass at 1 resolution |
| v0.2 | + second-order (MUSCL/SP-V, paper Sec. 2.3) | same tests, 3 resolutions, correct convergence order |
| v0.3 | + Sedov, Kidder, Guderley, LeBlanc | all pass, convergence documented |
| v0.4 | + p4est, Lagrange-AMR, Euler-AMR, 2D | paper's Sod-AMR/triple-point-AMR parity |
| v0.5 | + 3D, cylindrical geometry | DoD item 2 |
| v1.0 | + science demo, restart, report generation | DoD fully met |

## Known design gaps (not yet built, tracked so they aren't forgotten)
- **Second-order rollout** — not scheduled in detail yet; needs its own
  atomic breakdown same as Sec. 2, lands at v0.2.
- **Debug/probe tooling** — a `--probe cell_id` flag to dump full state for
  one cell/node across timesteps. Cheap, not yet added; the reference paper's
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
