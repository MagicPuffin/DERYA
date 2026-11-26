#pragma once

#include <iosfwd>
#include <vector>

#include "deck.hpp"
#include "hydro.hpp"

namespace hydro::app {

// Taylor-Green vortex (Georges' thesis Sec. 3.2): the exact (stationary)
// velocity and pressure at x (Eqs. 3.4-3.5), and the energy source per unit
// volume that keeps them stationary (Eq. 3.6, Appendix B):
//   S = pi density c1^3 / (4 (gamma - 1))
//       (cos(3 pi x) cos(pi y) - cos(3 pi y) cos(pi x)).
lagrangian::Vec3 taylor_green_velocity(const io::TaylorGreen& tg, const lagrangian::Vec3& x);
double taylor_green_pressure(const io::TaylorGreen& tg, const lagrangian::Vec3& x);
double taylor_green_source(const io::TaylorGreen& tg, double gamma, const lagrangian::Vec3& x);

// State on the deck's structured mesh, perturbed if the deck asks for it:
// each cell takes the background state, overridden by every region whose box
// [min, max) contains its centroid (the last such region wins); or, for a
// Taylor-Green deck, the profile at its centroid.
// eps = P / ((gamma - 1) rho), except in the cells of the deck's energy
// spike (io::EnergySpike), whose eps holds the spike's energy. Throws
// std::runtime_error if no cell has a node at the spike's point.
lagrangian::HydroState initial_state(const io::Deck& deck);

// The deck's energy source (empty for none).
lagrangian::EnergySource energy_source(const io::Deck& deck);

struct RunResult {
  lagrangian::HydroState state;
  double t = 0.0;
  long steps = 0;
};

// Cells and nodes to trace through a run (hydro_run --probe, --probe-node).
// With an output stream and at least one id, run writes one JSON line per
// probed cell, then per probed node, at step 0 (the initial state, dt = 0) and
// after every step, so the last lines before a failure are the state the
// failing step started from. Cell lines: {"cell", "step", "t", "dt",
// "centroid", "volume", "mass", "density", "velocity", "pressure",
// "internal_energy", "sound_speed"}, computed as in the dump
// (io::write_dump), so the last step's values equal the dump's. Node lines:
// {"node", "step", "t", "dt", "position", "velocity"}, the velocity being the
// V_p the step moved the node with (null at step 0). Probing cells costs one
// more geometry and thermo pass over the mesh per step.
struct Probes {
  std::vector<mesh::Index> cells;
  std::vector<mesh::Index> nodes;
  std::ostream* out = nullptr;
};

// Runs the deck from t = 0 to t_end: dt from timestep::stable_dt (dt^0 =
// min(dt_initial, CFL)), then a step of the deck's scheme
// (lagrangian::scheme_step: first order or predictor-corrector), the last
// step cut to land on t_end. Writes a progress line every deck.log_every steps and a
// summary at the end to log. Throws std::runtime_error, with the step and
// time, if a step fails (tangled mesh, negative internal energy, Newton
// failure) or max_steps is reached before t_end. Throws std::invalid_argument,
// before the first step, if a probed id is not a cell (node) of the mesh.
RunResult run(const io::Deck& deck, std::ostream& log, const Probes& probes = {});

}  // namespace hydro::app
