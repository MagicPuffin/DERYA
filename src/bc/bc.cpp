#include "bc.hpp"

#include <cmath>

namespace hydro::bc {

namespace {

// A face normal whose component orthogonal to the walls already found at a
// node is shorter than this is coplanar with them (the unit normals of faces
// on one planar wall agree to round-off).
constexpr double kCoplanarTol = 1e-8;

double dot(const Vec3& a, const Vec3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Component of v orthogonal to the orthonormal vectors in basis.
Vec3 orthogonal_part(Vec3 v, const std::vector<Vec3>& basis) {
  for (const Vec3& e : basis) {
    const double c = dot(v, e);
    for (int d = 0; d < 3; ++d) v[d] -= c * e[d];
  }
  return v;
}

// Appends the unit vector n to the orthonormal basis if it is independent of
// it, with the constraint V . n = value carried along by Gram-Schmidt:
// V . e = (value - sum_k (n . e_k) values_k) / |n - sum_k (n . e_k) e_k|.
void add_if_independent(const Vec3& n, double value, std::vector<Vec3>& basis,
                        std::vector<double>& values) {
  Vec3 w = n;
  for (std::size_t k = 0; k < basis.size(); ++k) {
    const double c = dot(n, basis[k]);
    for (int d = 0; d < 3; ++d) w[d] -= c * basis[k][d];
    value -= c * values[k];
  }
  const double norm = std::sqrt(dot(w, w));
  if (norm <= kCoplanarTol) return;
  basis.push_back({w[0] / norm, w[1] / norm, w[2] / norm});
  values.push_back(value / norm);
}

}  // namespace

void apply_pressure_bcs(const mesh::Mesh& m, const BoundarySet& bcs,
                        const std::vector<double>& cell_pressure,
                        std::vector<Vec3>& B) {
  for (mesh::Index f = 0; f < m.num_faces(); ++f) {
    const int side = m.face_boundary[f];
    if (side == mesh::kInteriorFace || bcs[side].type != BoundaryType::Outflow) continue;
    const double p_star = cell_pressure[m.face_owner[f]];
    const auto vecs = mesh::face_area_vectors(m, f);
    for (int n = 0; n < mesh::kNodesPerFace; ++n) {
      Vec3& b = B[m.face_nodes[f][n]];
      for (int d = 0; d < 3; ++d) b[d] -= p_star * vecs[n][d];
    }
  }
}

WallConstraints wall_constraints(const mesh::Mesh& m, const BoundarySet& bcs) {
  WallConstraints walls;
  walls.normals.resize(m.num_nodes());
  walls.velocity.assign(m.num_nodes(), Vec3{0.0, 0.0, 0.0});
  // Per node, V . e_k for the basis vectors e_k in walls.normals.
  std::vector<std::vector<double>> values(m.num_nodes());
  for (mesh::Index f = 0; f < m.num_faces(); ++f) {
    const int side = m.face_boundary[f];
    if (side == mesh::kInteriorFace || bcs[side].type == BoundaryType::Outflow) continue;
    Vec3 normal = {0.0, 0.0, 0.0};
    for (const auto& v : mesh::face_area_vectors(m, f)) {
      for (int d = 0; d < 3; ++d) normal[d] += v[d];
    }
    const double norm = std::sqrt(dot(normal, normal));
    for (int d = 0; d < 3; ++d) normal[d] /= norm;
    const double value =
        bcs[side].type == BoundaryType::Piston ? dot(bcs[side].velocity, normal) : 0.0;
    for (mesh::Index p : m.face_nodes[f]) {
      add_if_independent(normal, value, walls.normals[p], values[p]);
    }
  }
  for (mesh::Index p = 0; p < m.num_nodes(); ++p) {
    for (std::size_t k = 0; k < values[p].size(); ++k) {
      for (int d = 0; d < 3; ++d) walls.velocity[p][d] += values[p][k] * walls.normals[p][k][d];
    }
  }
  return walls;
}

std::vector<std::vector<Vec3>> wall_normals(const mesh::Mesh& m,
                                            const BoundarySet& bcs) {
  return wall_constraints(m, bcs).normals;
}

Vec3 solve_nodal_velocity(const Mat3& M, const Vec3& B,
                          const std::vector<Vec3>& walls) {
  // Complete the walls to an orthonormal basis of R^3, greedily taking the
  // coordinate axis with the largest part orthogonal to what is there; the
  // added vectors span the tangent space T.
  std::vector<Vec3> basis = walls;
  const auto num_walls = basis.size();
  while (basis.size() < 3) {
    Vec3 best{};
    double best_norm = -1.0;
    for (int a = 0; a < 3; ++a) {
      Vec3 axis = {0.0, 0.0, 0.0};
      axis[a] = 1.0;
      const Vec3 w = orthogonal_part(axis, basis);
      const double norm = std::sqrt(dot(w, w));
      if (norm > best_norm) {
        best = w;
        best_norm = norm;
      }
    }
    basis.push_back({best[0] / best_norm, best[1] / best_norm, best[2] / best_norm});
  }
  const std::vector<Vec3> T(basis.begin() + num_walls, basis.end());
  const int k = static_cast<int>(T.size());

  // Reduced SPD system (T^T M T) x = T^T B, solved by Gaussian elimination.
  double A[3][3];
  double r[3];
  for (int i = 0; i < k; ++i) {
    Vec3 MTi;
    for (int d = 0; d < 3; ++d) MTi[d] = dot(M[d], T[i]);
    for (int j = 0; j < k; ++j) A[j][i] = dot(T[j], MTi);
    r[i] = dot(T[i], B);
  }
  for (int i = 0; i < k; ++i) {
    for (int j = i + 1; j < k; ++j) {
      const double factor = A[j][i] / A[i][i];
      for (int l = i; l < k; ++l) A[j][l] -= factor * A[i][l];
      r[j] -= factor * r[i];
    }
  }
  double x[3];
  for (int i = k - 1; i >= 0; --i) {
    double s = r[i];
    for (int l = i + 1; l < k; ++l) s -= A[i][l] * x[l];
    x[i] = s / A[i][i];
  }

  Vec3 V = {0.0, 0.0, 0.0};
  for (int i = 0; i < k; ++i) {
    for (int d = 0; d < 3; ++d) V[d] += x[i] * T[i][d];
  }
  return V;
}

}  // namespace hydro::bc
