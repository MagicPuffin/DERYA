#pragma once

#include <iosfwd>

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
// eps = P / ((gamma - 1) rho).
lagrangian::HydroState initial_state(const io::Deck& deck);

// The deck's energy source (empty for none).
lagrangian::EnergySource energy_source(const io::Deck& deck);

struct RunResult {
  lagrangian::HydroState state;
  double t = 0.0;
  long steps = 0;
};

// Runs the deck from t = 0 to t_end: dt from timestep::stable_dt (dt^0 =
// min(dt_initial, CFL)), then a step of the deck's scheme
// (lagrangian::scheme_step: first order or predictor-corrector), the last
// step cut to land on t_end. Writes a progress line every deck.log_every steps and a
// summary at the end to log. Throws std::runtime_error, with the step and
// time, if a step fails (tangled mesh, negative internal energy, Newton
// failure) or max_steps is reached before t_end.
RunResult run(const io::Deck& deck, std::ostream& log);

}  // namespace hydro::app
