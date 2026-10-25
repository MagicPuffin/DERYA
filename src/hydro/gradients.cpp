#include "gradients.hpp"

namespace hydro::lagrangian {

using mesh::Index;
using mesh::kFacesPerCell;
using mesh::kHexFaceNodes;
using mesh::kNodesPerCell;
using mesh::kNodesPerFace;

std::vector<Vec3> cell_centers(const mesh::Mesh& m) {
  std::vector<Vec3> x(m.num_cells(), Vec3{});
  for (Index c = 0; c < m.num_cells(); ++c) {
    for (const Index p : m.cell_nodes[c]) {
      x[c][0] += m.node_x[p];
      x[c][1] += m.node_y[p];
      x[c][2] += m.node_z[p];
    }
    for (int d = 0; d < 3; ++d) x[c][d] /= kNodesPerCell;
  }
  return x;
}

CellGradients cell_gradients(const mesh::Mesh& m,
                             const std::vector<CellFaceVectors>& faces,
                             const std::vector<CornerValues>& nodal_pressure,
                             const std::vector<Vec3>& node_velocity) {
  CellGradients g;
  g.pressure.assign(m.num_cells(), Vec3{});
  g.velocity.assign(m.num_cells(), Mat3{});
  for (Index c = 0; c < m.num_cells(); ++c) {
    Vec3& gp = g.pressure[c];
    Mat3& gv = g.velocity[c];
    for (int l = 0; l < kFacesPerCell; ++l) {
      for (int n = 0; n < kNodesPerFace; ++n) {
        const Vec3& Sn = faces[c][l][n];
        const double P = nodal_pressure[c][l][n];
        const Vec3& Vp = node_velocity[m.cell_nodes[c][kHexFaceNodes[l][n]]];
        for (int j = 0; j < 3; ++j) {
          gp[j] += P * Sn[j];
          for (int i = 0; i < 3; ++i) gv[i][j] += Vp[i] * Sn[j];
        }
      }
    }
    const double inv_volume = 1.0 / cell_volume(m, c, faces[c]);
    for (int j = 0; j < 3; ++j) {
      gp[j] *= inv_volume;
      for (int i = 0; i < 3; ++i) gv[i][j] *= inv_volume;
    }
  }
  return g;
}

}  // namespace hydro::lagrangian
