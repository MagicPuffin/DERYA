#include "mesh.hpp"

#include <cmath>

namespace hydro::mesh {

namespace {

enum LocalFace { kXMinus = 0, kXPlus, kYMinus, kYPlus, kZMinus, kZPlus };

}  // namespace

Mesh generate_structured_mesh(const StructuredMeshSpec& spec) {
  const Index nx = spec.nx, ny = spec.ny, nz = spec.nz;
  const double dx = (spec.x_max - spec.x_min) / nx;
  const double dy = (spec.y_max - spec.y_min) / ny;
  const double dz = (spec.z_max - spec.z_min) / nz;

  auto node_id = [&](Index i, Index j, Index k) {
    return k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
  };
  auto cell_id = [&](Index i, Index j, Index k) {
    return k * ny * nx + j * nx + i;
  };
  const Index fx_count = (nx + 1) * ny * nz;
  const Index fy_count = nx * (ny + 1) * nz;
  auto x_face_id = [&](Index i, Index j, Index k) {
    return k * ny * (nx + 1) + j * (nx + 1) + i;
  };
  auto y_face_id = [&](Index i, Index j, Index k) {
    return fx_count + k * (ny + 1) * nx + j * nx + i;
  };
  auto z_face_id = [&](Index i, Index j, Index k) {
    return fx_count + fy_count + k * ny * nx + j * nx + i;
  };

  Mesh m;

  const Index num_nodes = (nx + 1) * (ny + 1) * (nz + 1);
  m.node_x.reserve(num_nodes);
  m.node_y.reserve(num_nodes);
  m.node_z.reserve(num_nodes);
  for (Index k = 0; k <= nz; ++k) {
    for (Index j = 0; j <= ny; ++j) {
      for (Index i = 0; i <= nx; ++i) {
        m.node_x.push_back(spec.x_min + i * dx);
        m.node_y.push_back(spec.y_min + j * dy);
        m.node_z.push_back(spec.z_min + k * dz);
      }
    }
  }

  const Index num_cells = nx * ny * nz;
  m.cell_nodes.reserve(num_cells);
  m.cell_faces.reserve(num_cells);
  for (Index k = 0; k < nz; ++k) {
    for (Index j = 0; j < ny; ++j) {
      for (Index i = 0; i < nx; ++i) {
        m.cell_nodes.push_back({node_id(i, j, k), node_id(i + 1, j, k),
                                node_id(i + 1, j + 1, k), node_id(i, j + 1, k),
                                node_id(i, j, k + 1), node_id(i + 1, j, k + 1),
                                node_id(i + 1, j + 1, k + 1),
                                node_id(i, j + 1, k + 1)});
        m.cell_faces.push_back({x_face_id(i, j, k), x_face_id(i + 1, j, k),
                                y_face_id(i, j, k), y_face_id(i, j + 1, k),
                                z_face_id(i, j, k), z_face_id(i, j, k + 1)});
      }
    }
  }

  const Index num_faces = fx_count + fy_count + nx * ny * (nz + 1);
  m.face_nodes.reserve(num_faces);
  m.face_owner.reserve(num_faces);
  m.face_neighbor.reserve(num_faces);
  m.face_boundary.reserve(num_faces);

  // A face on the low side of the domain is owned by the cell above it (seen
  // through that cell's minus face); every other face is owned by the cell
  // below it (seen through its plus face), with the cell above as neighbor
  // unless the face is on the high boundary. A boundary face's side is the
  // owner's local face it was seen through.
  auto add_face = [&](Index owner, int local_face, Index neighbor) {
    std::array<Index, kNodesPerFace> nodes;
    for (int n = 0; n < kNodesPerFace; ++n) {
      nodes[n] = m.cell_nodes[owner][kHexFaceNodes[local_face][n]];
    }
    m.face_nodes.push_back(nodes);
    m.face_owner.push_back(owner);
    m.face_neighbor.push_back(neighbor);
    m.face_boundary.push_back(neighbor == kNoCell ? local_face : kInteriorFace);
  };

  for (Index k = 0; k < nz; ++k) {
    for (Index j = 0; j < ny; ++j) {
      for (Index i = 0; i <= nx; ++i) {
        if (i == 0) {
          add_face(cell_id(0, j, k), kXMinus, kNoCell);
        } else {
          add_face(cell_id(i - 1, j, k), kXPlus,
                   i < nx ? cell_id(i, j, k) : kNoCell);
        }
      }
    }
  }
  for (Index k = 0; k < nz; ++k) {
    for (Index j = 0; j <= ny; ++j) {
      for (Index i = 0; i < nx; ++i) {
        if (j == 0) {
          add_face(cell_id(i, 0, k), kYMinus, kNoCell);
        } else {
          add_face(cell_id(i, j - 1, k), kYPlus,
                   j < ny ? cell_id(i, j, k) : kNoCell);
        }
      }
    }
  }
  for (Index k = 0; k <= nz; ++k) {
    for (Index j = 0; j < ny; ++j) {
      for (Index i = 0; i < nx; ++i) {
        if (k == 0) {
          add_face(cell_id(i, j, 0), kZMinus, kNoCell);
        } else {
          add_face(cell_id(i, j, k - 1), kZPlus,
                   k < nz ? cell_id(i, j, k) : kNoCell);
        }
      }
    }
  }

  return m;
}

std::array<Vec3, kNodesPerFace> face_area_vectors(const Mesh& m, Index f) {
  using Vec = Vec3;
  const auto& fn = m.face_nodes[f];

  std::array<Vec, kNodesPerFace> x;
  Vec center = {0.0, 0.0, 0.0};
  for (int n = 0; n < kNodesPerFace; ++n) {
    x[n] = {m.node_x[fn[n]], m.node_y[fn[n]], m.node_z[fn[n]]};
    for (int d = 0; d < 3; ++d) center[d] += x[n][d] / kNodesPerFace;
  }

  // tri[n]: area vector of triangle (p*_f, x[n], x[n+1]).
  std::array<Vec, kNodesPerFace> tri;
  Vec total = {0.0, 0.0, 0.0};
  for (int n = 0; n < kNodesPerFace; ++n) {
    const Vec& b = x[n];
    const Vec& c = x[(n + 1) % kNodesPerFace];
    const Vec u = {b[0] - center[0], b[1] - center[1], b[2] - center[2]};
    const Vec v = {c[0] - center[0], c[1] - center[1], c[2] - center[2]};
    tri[n] = {0.5 * (u[1] * v[2] - u[2] * v[1]),
              0.5 * (u[2] * v[0] - u[0] * v[2]),
              0.5 * (u[0] * v[1] - u[1] * v[0])};
    for (int d = 0; d < 3; ++d) total[d] += tri[n][d];
  }

  // Node n touches triangles n-1 and n.
  std::array<Vec, kNodesPerFace> result;
  for (int n = 0; n < kNodesPerFace; ++n) {
    const Vec& before = tri[(n + kNodesPerFace - 1) % kNodesPerFace];
    for (int d = 0; d < 3; ++d) {
      result[n][d] =
          (before[d] + tri[n][d] + total[d] / kNodesPerFace) / 3.0;
    }
  }
  return result;
}

void apply_saltzman_skew(Mesh& m, const StructuredMeshSpec& spec) {
  const double pi = std::acos(-1.0);
  const double length = spec.x_max - spec.x_min;
  for (Index p = 0; p < m.num_nodes(); ++p) {
    const double x = m.node_x[p], y = m.node_y[p];
    m.node_x[p] = x + (spec.y_max - y) * std::sin(pi * (x - spec.x_min) / length);
  }
}

}  // namespace hydro::mesh
