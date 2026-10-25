#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "eos.hpp"
#include "gradients.hpp"
#include "hydro.hpp"
#include "mesh.hpp"

using Catch::Matchers::WithinAbs;
using hydro::bc::Mat3;
using hydro::eos::IdealGasEOS;
using hydro::lagrangian::CornerValues;
using hydro::lagrangian::Vec3;
using hydro::lagrangian::cell_centers;
using hydro::lagrangian::cell_face_vectors;
using hydro::lagrangian::cell_gradients;
using hydro::lagrangian::cell_thermo;
using hydro::lagrangian::corner_impedances;
using hydro::lagrangian::make_state;
using hydro::lagrangian::nodal_pressures;
using hydro::mesh::Index;
using hydro::mesh::kFacesPerCell;
using hydro::mesh::kHexFaceNodes;
using hydro::mesh::kNodesPerFace;
using hydro::mesh::Mesh;
using hydro::mesh::StructuredMeshSpec;
using hydro::mesh::generate_structured_mesh;

namespace {

void check_vec(const Vec3& actual, const Vec3& expected, double tol) {
  for (int d = 0; d < 3; ++d) CHECK_THAT(actual[d], WithinAbs(expected[d], tol));
}

void check_mat(const Mat3& actual, const Mat3& expected, double tol) {
  for (int i = 0; i < 3; ++i) check_vec(actual[i], expected[i], tol);
}

Vec3 node_position(const Mesh& m, Index p) {
  return {m.node_x[p], m.node_y[p], m.node_z[p]};
}

// 3x2x2 mesh of [0,3]x[0,1]x[0,0.5] with every interior node displaced by a
// smooth, non-affine perturbation, so faces are non-planar and cells differ
// (as in test_hydro.cpp).
Mesh distorted_mesh() {
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 2;
  spec.nz = 2;
  spec.x_max = 3.0;
  spec.y_max = 1.0;
  spec.z_max = 0.5;
  Mesh m = generate_structured_mesh(spec);
  for (Index p = 0; p < m.num_nodes(); ++p) {
    const double x = m.node_x[p], y = m.node_y[p], z = m.node_z[p];
    const bool interior = x > 0.0 && x < 3.0 && y > 0.0 && y < 1.0 && z > 0.0 &&
                          z < 0.5;
    if (!interior) continue;
    m.node_x[p] += 0.15 * std::sin(2.0 * y + 3.0 * z);
    m.node_y[p] += 0.1 * std::cos(x + 5.0 * z);
    m.node_z[p] += 0.05 * std::sin(x * y);
  }
  return m;
}

// The single box cell [0,2]x[0,1]x[0,0.5], volume 1.
Mesh box_cell() {
  StructuredMeshSpec spec;
  spec.x_max = 2.0;
  spec.z_max = 0.5;
  return generate_structured_mesh(spec);
}

// Corner values P_cfp = f(x_p) for every cell, face and node.
template <class F>
std::vector<CornerValues> corner_values_at_nodes(const Mesh& m, F f) {
  std::vector<CornerValues> v(m.num_cells());
  for (Index c = 0; c < m.num_cells(); ++c) {
    for (int l = 0; l < kFacesPerCell; ++l) {
      for (int n = 0; n < kNodesPerFace; ++n) {
        v[c][l][n] = f(node_position(m, m.cell_nodes[c][kHexFaceNodes[l][n]]));
      }
    }
  }
  return v;
}

}  // namespace

TEST_CASE("cell_centers: mean of the cell's nodes", "[gradients]") {
  check_vec(cell_centers(box_cell())[0], {1.0, 0.5, 0.25}, 1e-15);

  const Mesh m = distorted_mesh();
  const auto x = cell_centers(m);
  REQUIRE(x.size() == static_cast<std::size_t>(m.num_cells()));
  for (Index c = 0; c < m.num_cells(); ++c) {
    Vec3 mean{};
    for (const Index p : m.cell_nodes[c]) {
      const Vec3 xp = node_position(m, p);
      for (int d = 0; d < 3; ++d) mean[d] += xp[d] / 8.0;
    }
    check_vec(x[c], mean, 1e-15);
  }
}

TEST_CASE("cell_gradients: a constant field has zero gradient (closure)",
          "[gradients]") {
  const Mesh m = distorted_mesh();
  const auto faces = cell_face_vectors(m);
  const auto P = corner_values_at_nodes(m, [](const Vec3&) { return 2.5; });
  const std::vector<Vec3> V(m.num_nodes(), Vec3{-1.0, 3.0, 0.5});
  const auto g = cell_gradients(m, faces, P, V);
  for (Index c = 0; c < m.num_cells(); ++c) {
    INFO("cell " << c);
    check_vec(g.pressure[c], {0.0, 0.0, 0.0}, 1e-13);
    check_mat(g.velocity[c], Mat3{}, 1e-13);
  }
}

TEST_CASE("cell_gradients: exact for linear fields on a distorted mesh",
          "[gradients]") {
  // Thesis Appendix A.2: 1/V_c sum_p sum_f S_pf x_p (x) n_pf = Id on any cell.
  const Mesh m = distorted_mesh();
  const auto faces = cell_face_vectors(m);
  const Vec3 a = {0.7, -1.3, 2.1};
  const auto P = corner_values_at_nodes(m, [&](const Vec3& x) {
    return 4.0 + a[0] * x[0] + a[1] * x[1] + a[2] * x[2];
  });
  const Mat3 A = {{{1.0, -0.5, 0.25}, {2.0, 0.3, -1.1}, {-0.4, 0.9, 1.7}}};
  const Vec3 b = {0.2, -0.6, 1.5};
  std::vector<Vec3> V(m.num_nodes());
  for (Index p = 0; p < m.num_nodes(); ++p) {
    const Vec3 x = node_position(m, p);
    for (int i = 0; i < 3; ++i) {
      V[p][i] = b[i] + A[i][0] * x[0] + A[i][1] * x[1] + A[i][2] * x[2];
    }
  }
  const auto g = cell_gradients(m, faces, P, V);
  for (Index c = 0; c < m.num_cells(); ++c) {
    INFO("cell " << c);
    check_vec(g.pressure[c], a, 1e-12);
    check_mat(g.velocity[c], A, 1e-12);
  }
}

TEST_CASE("cell_gradients: hand-computed values on a box cell", "[gradients]") {
  // Box [0,2]x[0,1]x[0,0.5], V_c = 1. On a planar rectangle each node gets a
  // quarter of the face area, so the face areas are 0.5 (x), 1 (y), 2 (z).
  const Mesh m = box_cell();
  const auto faces = cell_face_vectors(m);

  // P_cfp constant per face (faces -x, +x, -y, +y, -z, +z): (grad P)_j =
  // A_j (P_+j - P_-j).
  const double face_P[kFacesPerCell] = {1.0, 3.0, 2.0, 5.0, 0.0, 0.25};
  std::vector<CornerValues> P(1);
  for (int l = 0; l < kFacesPerCell; ++l) {
    for (int n = 0; n < kNodesPerFace; ++n) P[0][l][n] = face_P[l];
  }

  // V_x = x y at the nodes, bilinear: dV_x/dx = mean of y over the +x face
  // nodes times A_x = 0.5, dV_x/dy = 1 (mean of x), dV_x/dz = 0.
  std::vector<Vec3> V(m.num_nodes());
  for (Index p = 0; p < m.num_nodes(); ++p) V[p] = {m.node_x[p] * m.node_y[p], 0.0, 0.0};

  const auto g = cell_gradients(m, faces, P, V);
  check_vec(g.pressure[0], {1.0, 3.0, 0.5}, 1e-14);
  check_mat(g.velocity[0], {{{0.5, 1.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}},
            1e-14);
}

TEST_CASE("nodal_pressures: Eq. 3 at given nodal velocities", "[gradients]") {
  // Gas at rest, every node moving at +x: P_cfp = P_c + Z_cfp on the -x face,
  // P_c - Z_cfp on the +x face, P_c on the others.
  const Mesh m = box_cell();
  const auto s = make_state(m, {1.0}, {0.0}, {0.0}, {0.0}, {2.5});
  const IdealGasEOS eos(1.4);
  const auto faces = cell_face_vectors(s.mesh);
  const auto thermo = cell_thermo(s, faces, eos);
  const std::vector<Vec3> V(m.num_nodes(), Vec3{1.0, 0.0, 0.0});
  const auto Z = corner_impedances(s, faces, thermo, V);
  const auto P = nodal_pressures(s, faces, thermo, Z, V);
  const double Pc = thermo.pressure[0];
  REQUIRE_THAT(Pc, WithinAbs(1.0, 1e-14));
  const double sign[kFacesPerCell] = {1.0, -1.0, 0.0, 0.0, 0.0, 0.0};
  for (int l = 0; l < kFacesPerCell; ++l) {
    for (int n = 0; n < kNodesPerFace; ++n) {
      INFO("face " << l << ", node " << n);
      CHECK(Z[0][l][n] > 0.0);
      CHECK_THAT(P[0][l][n], WithinAbs(Pc + sign[l] * Z[0][l][n], 1e-14));
    }
  }
}
