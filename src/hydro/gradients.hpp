#pragma once

#include <vector>

#include "hydro.hpp"

namespace hydro::lagrangian {

// Cell gradients of the second-order extension (Colaïtis, Guisset & Breil,
// SSRN 5167115, Sec. 2.3, Eq. 6; Georges' thesis, Eqs. 2.3, 2.14-2.15). See
// decisions.md, "Second-order extension".

// Reconstruction point x_c of each cell: the mean of its nodes (thesis Eq.
// 2.3), not the volume centroid.
std::vector<Vec3> cell_centers(const mesh::Mesh& m);

// Green-formula gradients per cell, from given nodal values:
//   (grad P)_c = 1/V_c sum_p sum_f S_pf P_cfp n_pf,
//   (grad V)_c = 1/V_c sum_p sum_f S_pf V_p (x) n_pf,
// with velocity[c][i][j] = dV_i/dx_j. The sums run over every face of the
// cell, boundary faces included, so the nodal values carry the boundary
// conditions. Both are exact for fields linear in x sampled at the nodes,
// on any mesh (thesis Appendix A.2). The nodal values are those of a
// first-order solve: nodal_pressures and nodal_velocities.
struct CellGradients {
  std::vector<Vec3> pressure;
  std::vector<Mat3> velocity;
};
CellGradients cell_gradients(const mesh::Mesh& m,
                             const std::vector<CellFaceVectors>& faces,
                             const std::vector<CornerValues>& nodal_pressure,
                             const std::vector<Vec3>& node_velocity);

}  // namespace hydro::lagrangian
