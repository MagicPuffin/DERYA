# Deck and dump formats

`hydro_run <deck.json> [dump.json]` reads a JSON deck, runs it to `t_end`
and writes the final state as a JSON dump (options: [Probes](#probes---probe---probe-node)). The dump path is the second
argument, or `output.path` from the deck. Relative paths are resolved
against the current directory. Warnings, such as the thin-slab guard, go to
stderr. A failed run exits with status 1 and reports the step and time
where it failed. Examples: [decks/sod.json](../decks/sod.json),
[decks/noh.json](../decks/noh.json),
[decks/saltzman.json](../decks/saltzman.json),
[decks/taylor_green.json](../decks/taylor_green.json),
[decks/leblanc.json](../decks/leblanc.json). The convergence families, the
same problems at 50/100/200 cells (LeBlanc 225/450/900, uniform and with
a refined layer at the contact), are in
[decks/convergence/](../decks/convergence/); run them with
`python3 scripts/convergence.py` (`--keep-dumps DIR` keeps the dumps for
plotting).

The reader rejects unknown keys, so a misspelled key is an error, not a
silently ignored setting. The parser is in `src/io/deck.cpp`.

## Deck

| Key | Type | Required | Meaning |
|---|---|---|---|
| `mesh.cells` | 3 integers ≥ 1 | yes | `nx, ny, nz` |
| `mesh.min`, `mesh.max` | 3 numbers | yes | Box corners; `max > min` in every direction |
| `mesh.x_segments` | array | no | Grades the mesh in x: `[[x0, x1, n], ...]`, `n` uniform cells on `[x0, x1]` each, contiguous (each `x0` equals the previous `x1`, exactly as written), from `mesh.min[0]` to `mesh.max[0]`, the `n` adding up to `mesh.cells[0]`. Not with `mesh.perturbation`. Example: `decks/convergence/leblanc_layer_225.json` |
| `mesh.mapping` | string | no | `"cartesian"` (default) or `"polar"`: the box's x is the radius r (`mesh.min[0]` ≥ 0) and its y the angle theta in **degrees** (span at most 180, each cell less than 180), z unchanged; nodes at `(r cos theta, r sin theta, z)`. The sides keep their names: `x_min` the inner surface, `x_max` the outer, `y_min`/`y_max` the radial planes theta = `min[1]`/`max[1]`. The planes theta = 0, 90, 180 are exactly the coordinate planes. Walls go on the planes only; a curved side is `"outflow"` (a curved wall pins its nodes like an edge). With `min[0] = 0` the axis is one node per z level, the cells there triangular prisms and their `x_min` faces of zero area (any type; it acts on nothing). Works with `mesh.x_segments` (graded in r), not with `mesh.perturbation`. Regions and centroids are in the mapped x, y, z |
| `mesh.perturbation` | string | no | Only `"saltzman"`: each node moves in x by `(y_max - y) sin(pi (x - x_min) / (x_max - x_min))` |
| `eos.type` | string | yes | Only `"ideal_gas"` |
| `eos.gamma` | number > 1 | yes | |
| `boundaries.{x,y,z}_{min,max}` | string or object | yes, all six | `"symmetry"` (wall, `V·n = 0`), `"outflow"` (prescribed pressure = owner cell pressure), `{"type": "piston", "velocity": [u, v, w]}` (wall moving with that velocity: `V·n = u_w·n`, tangential motion free), or `{"type": "pressure", ...}`: a prescribed pressure `P*` in place of outflow's owner-cell pressure, not a wall. Either `"value": P` (≥ 0, constant) or `"law": "kidder", "p0": P0, "t_focus": t_f` (Kidder's shell, thesis Eqs. 3.13-3.14: `P*(t) = P0 h^(-2 gamma/(gamma-1))`, `h = sqrt(1 - (t/t_f)^2)`, gamma from `eos`; `P0` ≥ 0, `t_f` > 0). A law is evaluated at each stage's time: `t^n` at order 1, `t^n` and `t^n + dt` at order 2. The run fails at the first stage at `t ≥ t_f` |
| `initial.background` | state | yes | State of every cell not in a region |
| `initial.regions` | array | no | Each item has `min`, `max` (3 numbers) and the state keys. A cell whose centroid (on the perturbed mesh) is in `[min, max)` takes the state. Later regions override earlier ones |
| `initial.energy_spike` | object | no | A point blast on top of `background` and `regions` (Sedov): `{"at": [x, y, z], "energy_per_length": E}` (E > 0). E times the slab thickness `mesh.max[2] - mesh.min[2]` goes, as internal energy, into the cells with a node at `at` (to 1e-12 of the box size), split by mass so they share one eps; it replaces their background internal energy. At least one cell must have a node there. Example: `decks/convergence/sedov_30.json` (the corner cell at the origin of the quarter plane) |
| `initial.taylor_green` | object | no | Instead of `background`, `regions` and `energy_spike`: the Taylor-Green vortex (Georges' thesis Eqs. 3.4-3.5) at each cell's centroid, `V = c1 (sin πx cos πy, -cos πx sin πy, 0)`, `P = density c1² / 4 (cos 2πx + cos 2πy) + c2`, `rho = density`. Keys `density` (> 0), `c1`, `c2`, each default 1; `c2` must exceed `density c1² / 2` |
| `source` | string | no | Energy source added to the total energy equation. Only `"taylor_green"` (thesis Eq. 3.6, with `initial.taylor_green`'s constants and `eos.gamma`; requires that block), which keeps the vortex stationary |
| `time.t_end` | number > 0 | yes | The last step is cut so the run ends exactly at `t_end` |
| `time.dt_initial` | number > 0 | no | `dt^0 = min(dt_initial, CFL)`. Default: CFL alone. A cold gas needs it |
| `time.cfl`, `time.volume_fraction`, `time.max_growth`, `time.dt_max` | numbers | no | Timestep controller (`decisions.md`, "Lagrangian timestep"). Defaults 0.25, 0.1, 1.01, ∞ |
| `time.max_steps` | integer ≥ 1 | no | Fail if `t_end` is not reached in this many steps. Default: no limit |
| `scheme.order` | integer, 1 or 2 | no | 1: first-order step (paper Sec. 2.2). 2: second order (paper Sec. 2.3): reconstructed pressure and velocity with a predictor-corrector step, about 4.5 times the cost. Default 1 |
| `scheme.limiter` | string | no | Order 2 only: `"sp-v"` (paper Eqs. 7-8) or `"none"` (unlimited gradients, for convergence tests on smooth flows; nothing limits them at shocks). Default `"sp-v"` |
| `output.path` | string | no | Dump path when the command line gives none |
| `output.log_every` | integer ≥ 0 | no | Steps between progress lines on stdout, 0 for none. Default 100 |

A **state** is `density` (> 0, required), `pressure` (≥ 0, required) and
`velocity` (3 numbers, default zero). The specific internal energy is
`P / ((gamma - 1) rho)`.

To tabulate Taylor-Green pressure errors and orders over several meshes
against the thesis's Table 3.1(b):
`python3 scripts/taylor_green_convergence.py [--sizes 10 20 40 80] [--order 2] [--limiter none]`.

The **thin-slab guard** warns when a direction with a single cell is
thinner than the smallest cell of the multi-cell directions, because
`lambda_c` would then limit dt for no physical reason. On a polar mesh the
theta cells count as arcs, the shortest at the middle of the innermost
ring.

## Dump (`"format": "hydro-dump-0"`)

| Key | Contents |
|---|---|
| `t`, `steps` | Time and number of steps taken |
| `cells.centroid_{x,y,z}` | Mean of the cell's 8 nodes |
| `cells.volume`, `cells.mass`, `cells.density` | |
| `cells.velocity_{x,y,z}` | |
| `cells.pressure`, `cells.internal_energy` | Specific internal energy `E - |V|^2/2` |
| `nodes.{x,y,z}` | Node coordinates |

All arrays are in mesh order (`src/mesh/mesh.hpp`). This snapshot is enough
for plotting and comparison but not for restart: it has no connectivity or
boundary types. The restart schema is deferred to v1.0 (`roadmap.md`).

## Probes (`--probe`, `--probe-node`)

`hydro_run deck.json [dump.json] --probe CELL [--probe CELL ...] [--probe-node NODE ...] [--probe-out PATH]`
traces cells and nodes through the run, for debugging. Both flags repeat.
The lines go to `PATH`, else to stderr, one JSON object per line: each
probed cell in the order given, then each probed node, at step 0 (the
initial state, `dt = 0`) and after every step. A run that fails stops
after the lines of its last good step, the state the failing step started
from.

| Line | Keys |
|---|---|
| Cell | `cell`, `step`, `t`, `dt`, `centroid` (3), `volume`, `mass`, `density`, `velocity` (3), `pressure`, `internal_energy`, `sound_speed` |
| Node | `node`, `step`, `t`, `dt`, `position` (3), `velocity` (3: the `V_p` the step moved the node with; `null` at step 0) |

The values are computed as in the dump, so the last step's lines equal the
dump's entries exactly. Probing does not change the run: the dump is
byte-identical with and without probes. Ids are in mesh order, cell
`i + nx (j + ny k)` and node `i + (nx + 1) (j + (ny + 1) k)`, except on a
polar mesh down to r = 0, where each level's axis is one node (that of
`j = 0`) and the later nodes close up (`src/mesh/mesh.hpp`). An id outside
the mesh fails the run (status 1) before the first step; a value that is
not a non-negative integer is a usage error (status 2). Each step with a
probed cell costs one more geometry and thermo pass over the mesh.

To compare a dump with an exact solution:
`python3 scripts/plot_comparison.py DUMP [DUMP ...] [--oracle FILE] [--labels L ...] [--gamma G] [--log PANEL ...] [--xlim X0 X1] [-o PNG]`
(defaults: the Sod oracle, the file names as labels, the oracle's gamma or
1.4, linear axes, the meshes' extent, `comparison.png`; `--log` takes
`density`, `pressure`, `internal_energy`). Up to four dumps overlay, e.g. order 1 and order 2
of one deck.
