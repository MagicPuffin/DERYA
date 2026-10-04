#pragma once

#include <iosfwd>

#include "deck.hpp"
#include "hydro.hpp"

namespace hydro::app {

// State on the deck's structured mesh, perturbed if the deck asks for it:
// each cell takes the background state, overridden by every region whose box
// [min, max) contains its centroid (the last such region wins).
// eps = P / ((gamma - 1) rho).
lagrangian::HydroState initial_state(const io::Deck& deck);

struct RunResult {
  lagrangian::HydroState state;
  double t = 0.0;
  long steps = 0;
};

// Runs the deck from t = 0 to t_end: dt from timestep::stable_dt (dt^0 =
// min(dt_initial, CFL)), then lagrangian::lagrangian_step, the last step cut
// to land on t_end. Writes a progress line every deck.log_every steps and a
// summary at the end to log. Throws std::runtime_error, with the step and
// time, if a step fails (tangled mesh, negative internal energy, Newton
// failure) or max_steps is reached before t_end.
RunResult run(const io::Deck& deck, std::ostream& log);

}  // namespace hydro::app
