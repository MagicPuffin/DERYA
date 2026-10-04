#pragma once

#include <limits>
#include <vector>

#include "eos.hpp"
#include "hydro.hpp"

namespace hydro::timestep {

using lagrangian::Vec3;

// Timestep control of the Lagrangian step (the paper gives none; EUCCLHYD
// practice, see decisions.md, "Lagrangian timestep"):
//   dt^{n+1} = min(C_cfl min_c lambda_c / a_c, C_V min_c V_c / |dV_c/dt|,
//                  C_M dt^n, dt_max).
struct TimestepParams {
  double cfl = 0.25;             // C_cfl
  double volume_fraction = 0.1;  // C_V, the largest relative volume change
  double max_growth = 1.01;      // C_M
  double dt_max = std::numeric_limits<double>::infinity();
};

// lambda_c = V_c / (largest face area), the face area being the sum of its
// nodal areas S_pf: the shortest edge of a box cell.
double characteristic_length(const lagrangian::CellFaceVectors& faces, double volume);

// Largest stable dt for the state s, given the nodal velocities of the step
// just taken (dV_c/dt = sum_p n_cp . V_p, by the GCL) and its dt^n. Cells
// with a_c = 0 or dV_c/dt = 0 give no CFL or volume limit. node_velocity may
// be empty (the first step): the volume criterion is skipped, and dt_prev
// is then the initial dt (or infinity). Reduced over ranks with
// sync::global_reduce. Throws std::runtime_error if no criterion gives a
// finite dt (a cold gas at rest with no dt_prev or dt_max).
double stable_dt(const lagrangian::HydroState& s, const eos::EOS& eos,
                 const std::vector<Vec3>& node_velocity, double dt_prev,
                 const TimestepParams& params = {});

}  // namespace hydro::timestep
