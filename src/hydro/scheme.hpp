#pragma once

#include <array>
#include <functional>
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

// Prescribed energy source: power per unit volume S(x), added to the total
// energy equation, m_c dE_c/dt = -sum_p F_cp . V_p + V_c S(x_c), with x_c the
// node average (cell_centers). Taylor-Green needs one to stay stationary
// (thesis Eq. 3.6, Appendix B). Empty: no source.
using EnergySource = std::function<double(const Vec3& x)>;

// V_c S(x_c) per cell at the current geometry.
std::vector<double> source_power(const HydroState& s,
                                 const std::vector<CellFaceVectors>& faces,
                                 const EnergySource& source);

// Forward Euler on the source alone: E_c += dt power_c / m_c.
void apply_source(HydroState& s, const std::vector<double>& power, double dt);

// Nodal velocities V_p and corner forces F_cp (per cell and local corner) of
// one second-order stage at the state s: a first-order nodal solve, Green
// gradients from its P_cfp and V_p, the limiter, the second-order nodal
// solve on the reconstructed corner values, and the corner forces with Z at
// the converged V_p (so they balance at every node), and the source power at
// the stage's geometry (empty without a source). bcs are the boundaries at
// the stage's time (bc::at_time).
struct StageFluxes {
  std::vector<Vec3> node_velocity;
  std::vector<std::array<Vec3, mesh::kNodesPerCell>> forces;
  std::vector<double> source_power;
};
StageFluxes second_order_fluxes(const HydroState& s, const bc::BoundarySet& bcs,
                                const eos::EOS& eos, Limiter limiter,
                                const EnergySource& source = {});

// One predictor-corrector step of size dt from time t (thesis Eqs. 2.42-2.49,
// alpha = 1): the predictor is a forward-Euler step with the fluxes of s^n,
// F(t^n, U^n), to s^(1); the corrector moves the nodes of s^n by dt (V_p^n +
// V_p^(1)) / 2 and updates its cells with the mean of the two stages' fluxes,
// the second being F(t^n + dt, U^(1)): boundary pressure laws are evaluated at
// t in the predictor and at t + dt in the corrector stage. A source enters
// each stage like the fluxes: V_c S(x_c) at s^n in the predictor, the mean of
// s^n's and s^(1)'s in the corrector. Returns the mean nodal velocity (for
// the timestep controller's volume criterion).
std::vector<Vec3> predictor_corrector_step(HydroState& s, const bc::BoundarySet& bcs,
                                           const eos::EOS& eos, double t, double dt,
                                           Limiter limiter,
                                           const EnergySource& source = {});

// One step of the given scheme from time t: lagrangian_step at order 1 (the
// boundaries at t, the source, if any, at t^n's geometry, forward Euler),
// predictor_corrector_step at order 2.
std::vector<Vec3> scheme_step(HydroState& s, const bc::BoundarySet& bcs,
                              const eos::EOS& eos, double t, double dt,
                              const Scheme& scheme, const EnergySource& source = {});

}  // namespace hydro::lagrangian
