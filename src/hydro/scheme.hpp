#pragma once

#include <array>
#include <vector>

#include "bc.hpp"
#include "eos.hpp"
#include "hydro.hpp"

namespace hydro::lagrangian {

// Second-order step in time and space (Colaïtis, Guisset & Breil, SSRN
// 5167115, Sec. 2.3; Georges' thesis, Sec. 2.6.1, Eqs. 2.42-2.49) and the
// choice between it and the first-order step. See decisions.md,
// "Second-order extension".

// Slope limiter of the second-order reconstruction: SP-V (limiter.hpp), or
// none (unlimited gradients, for convergence tests on smooth flows).
enum class Limiter { SpV, None };

struct Scheme {
  int order = 1;  // 1: lagrangian_step; 2: predictor_corrector_step
  Limiter limiter = Limiter::SpV;  // order 2 only
};

// Nodal velocities V_p and corner forces F_cp (per cell and local corner) of
// one second-order stage at the state s: a first-order nodal solve, Green
// gradients from its P_cfp and V_p, the limiter, the second-order nodal
// solve on the reconstructed corner values, and the corner forces with Z at
// the converged V_p (so they balance at every node).
struct StageFluxes {
  std::vector<Vec3> node_velocity;
  std::vector<std::array<Vec3, mesh::kNodesPerCell>> forces;
};
StageFluxes second_order_fluxes(const HydroState& s, const bc::BoundarySet& bcs,
                                const eos::EOS& eos, Limiter limiter);

// One predictor-corrector step of size dt (thesis Eqs. 2.42-2.49, alpha = 1):
// the predictor is a forward-Euler step with the fluxes of s^n, to s^(1); the
// corrector moves the nodes of s^n by dt (V_p^n + V_p^(1)) / 2 and updates its
// cells with the mean of the two stages' fluxes. Returns that mean nodal
// velocity (for the timestep controller's volume criterion).
std::vector<Vec3> predictor_corrector_step(HydroState& s, const bc::BoundarySet& bcs,
                                           const eos::EOS& eos, double dt,
                                           Limiter limiter);

// One step of the given scheme: lagrangian_step at order 1,
// predictor_corrector_step at order 2.
std::vector<Vec3> scheme_step(HydroState& s, const bc::BoundarySet& bcs,
                              const eos::EOS& eos, double dt, const Scheme& scheme);

}  // namespace hydro::lagrangian
