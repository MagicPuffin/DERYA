# Deck and dump formats

`hydro_run <deck.json> [dump.json]` reads a JSON deck, runs it to `t_end`
and writes the final state as a JSON dump. The dump path is the second
argument, or `output.path` from the deck. Relative paths are resolved
against the current directory. Warnings, such as the thin-slab guard, go to
stderr. A failed run exits with status 1 and reports the step and time
where it failed. Examples: [decks/sod.json](../decks/sod.json),
[decks/noh.json](../decks/noh.json),
[decks/saltzman.json](../decks/saltzman.json).

The reader rejects unknown keys, so a misspelled key is an error, not a
silently ignored setting. The parser is in `src/io/deck.cpp`.

## Deck

| Key | Type | Required | Meaning |
|---|---|---|---|
| `mesh.cells` | 3 integers ≥ 1 | yes | `nx, ny, nz` |
| `mesh.min`, `mesh.max` | 3 numbers | yes | Box corners; `max > min` in every direction |
| `mesh.perturbation` | string | no | Only `"saltzman"`: each node moves in x by `(y_max - y) sin(pi (x - x_min) / (x_max - x_min))` |
| `eos.type` | string | yes | Only `"ideal_gas"` |
| `eos.gamma` | number > 1 | yes | |
| `boundaries.{x,y,z}_{min,max}` | string or object | yes, all six | `"symmetry"` (wall, `V·n = 0`), `"outflow"` (prescribed pressure = owner cell pressure), or `{"type": "piston", "velocity": [u, v, w]}` (wall moving with that velocity: `V·n = u_w·n`, tangential motion free) |
| `initial.background` | state | yes | State of every cell not in a region |
| `initial.regions` | array | no | Each item has `min`, `max` (3 numbers) and the state keys. A cell whose centroid (on the perturbed mesh) is in `[min, max)` takes the state. Later regions override earlier ones |
| `time.t_end` | number > 0 | yes | The last step is cut so the run ends exactly at `t_end` |
| `time.dt_initial` | number > 0 | no | `dt^0 = min(dt_initial, CFL)`. Default: CFL alone. A cold gas needs it |
| `time.cfl`, `time.volume_fraction`, `time.max_growth`, `time.dt_max` | numbers | no | Timestep controller (`decisions.md`, "Lagrangian timestep"). Defaults 0.25, 0.1, 1.01, ∞ |
| `time.max_steps` | integer ≥ 1 | no | Fail if `t_end` is not reached in this many steps. Default: no limit |
| `output.path` | string | no | Dump path when the command line gives none |
| `output.log_every` | integer ≥ 0 | no | Steps between progress lines on stdout, 0 for none. Default 100 |

A **state** is `density` (> 0, required), `pressure` (≥ 0, required) and
`velocity` (3 numbers, default zero). The specific internal energy is
`P / ((gamma - 1) rho)`.

The **thin-slab guard** warns when a direction with a single cell is
thinner than the smallest cell of the multi-cell directions, because
`lambda_c` would then limit dt for no physical reason.

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

To compare a dump with an exact solution:
`python3 scripts/plot_comparison.py DUMP [--oracle FILE] [--gamma G] [-o PNG]`
(defaults: the Sod oracle, the oracle's gamma or 1.4, `comparison.png`).
