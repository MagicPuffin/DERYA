#include "scheme.hpp"

#include <stdexcept>
#include <string>

#include "gradients.hpp"
#include "limiter.hpp"

namespace hydro::lagrangian {

StageFluxes second_order_fluxes(const HydroState& s, const bc::BoundarySet& bcs,
                                const eos::EOS& eos, Limiter limiter) {
  const mesh::Mesh& m = s.mesh;
  const auto faces = cell_face_vectors(m);
  const CellThermo thermo = cell_thermo(s, faces, eos);

  // Gradients from a first-order solve at the same state (decisions.md).
  const auto V1 = nodal_velocities(s, faces, thermo, bcs);
  const auto P1 =
      nodal_pressures(s, faces, thermo, corner_impedances(s, faces, thermo, V1), V1);
  const CellGradients g = cell_gradients(m, faces, P1, V1);

  std::vector<double> chi_P(m.num_cells(), 1.0);
  std::vector<double> chi_V(m.num_cells(), 1.0);
  if (limiter == Limiter::SpV) {
    chi_P = pressure_limiters(m, thermo.pressure, g.pressure);
    chi_V = velocity_limiters(m, s.vel_x, s.vel_y, s.vel_z, g.velocity, g.pressure);
  }
  const CornerState corners = reconstructed_corners(s, thermo, g, chi_P, chi_V);

  StageFluxes result;
  result.node_velocity = nodal_velocities(s, faces, thermo, corners, bcs);
  const auto Z = corner_impedances(s, faces, thermo, corners, result.node_velocity);
  result.forces = corner_forces(s, faces, corners, Z, result.node_velocity);
  return result;
}

std::vector<Vec3> predictor_corrector_step(HydroState& s, const bc::BoundarySet& bcs,
                                           const eos::EOS& eos, double dt,
                                           Limiter limiter) {
  // Predictor (Eqs. 2.45-2.46): forward Euler to s^(1).
  const StageFluxes stage0 = second_order_fluxes(s, bcs, eos, limiter);
  HydroState predicted = s;
  update_cells(predicted, stage0.forces, stage0.node_velocity, dt);
  move_nodes(predicted.mesh, stage0.node_velocity, dt);

  // Corrector (Eqs. 2.47-2.49): from s^n, each stage's fluxes over dt / 2.
  const StageFluxes stage1 = second_order_fluxes(predicted, bcs, eos, limiter);
  const double half = 0.5 * dt;
  update_cells(s, stage0.forces, stage0.node_velocity, half);
  update_cells(s, stage1.forces, stage1.node_velocity, half);
  move_nodes(s.mesh, stage0.node_velocity, half);
  move_nodes(s.mesh, stage1.node_velocity, half);

  std::vector<Vec3> mean(stage0.node_velocity.size());
  for (std::size_t p = 0; p < mean.size(); ++p) {
    for (int d = 0; d < 3; ++d) {
      mean[p][d] = 0.5 * (stage0.node_velocity[p][d] + stage1.node_velocity[p][d]);
    }
  }
  return mean;
}

std::vector<Vec3> scheme_step(HydroState& s, const bc::BoundarySet& bcs,
                              const eos::EOS& eos, double dt, const Scheme& scheme) {
  switch (scheme.order) {
    case 1:
      return lagrangian_step(s, bcs, eos, dt);
    case 2:
      return predictor_corrector_step(s, bcs, eos, dt, scheme.limiter);
    default:
      throw std::invalid_argument("scheme_step: order must be 1 or 2, not " +
                                  std::to_string(scheme.order));
  }
}

}  // namespace hydro::lagrangian
