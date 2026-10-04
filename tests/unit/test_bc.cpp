#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "bc.hpp"
#include "mesh.hpp"

using Catch::Matchers::WithinAbs;
using hydro::bc::Boundary;
using hydro::bc::BoundarySet;
using hydro::bc::BoundaryType;
using hydro::bc::Mat3;
using hydro::bc::Vec3;
using hydro::bc::apply_pressure_bcs;
using hydro::bc::solve_nodal_velocity;
using hydro::bc::wall_constraints;
using hydro::bc::wall_normals;
using hydro::mesh::Index;
using hydro::mesh::Mesh;
using hydro::mesh::StructuredMeshSpec;
using hydro::mesh::generate_structured_mesh;

namespace {

constexpr BoundaryType kSym = BoundaryType::Symmetry;
constexpr BoundaryType kOut = BoundaryType::Outflow;

void check_vec(const Vec3& actual, const Vec3& expected, double tol = 1e-13) {
  for (int d = 0; d < 3; ++d) CHECK_THAT(actual[d], WithinAbs(expected[d], tol));
}

double dot(const Vec3& a, const Vec3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 residual(const Mat3& M, const Vec3& V, const Vec3& B) {
  Vec3 r;
  for (int i = 0; i < 3; ++i) r[i] = dot(M[i], V) - B[i];
  return r;
}

Vec3 normalized(Vec3 v) {
  const double n = std::sqrt(dot(v, v));
  return {v[0] / n, v[1] / n, v[2] / n};
}

// Symmetric positive definite, with off-diagonal coupling so a wall
// constraint changes the tangential components too.
const Mat3 kM = {{{4.0, 1.0, 0.5}, {1.0, 3.0, -0.4}, {0.5, -0.4, 2.0}}};
const Vec3 kB = {1.0, -2.0, 0.7};

// The walls a node sees, as a sorted list of axes (0 = x, 1 = y, 2 = z), for
// walls that are axis-aligned.
std::vector<int> wall_axes(const std::vector<Vec3>& walls) {
  std::vector<int> axes;
  for (const Vec3& n : walls) {
    int axis = -1;
    for (int d = 0; d < 3; ++d) {
      if (std::abs(std::abs(n[d]) - 1.0) < 1e-14) axis = d;
    }
    axes.push_back(axis);
  }
  std::sort(axes.begin(), axes.end());
  return axes;
}

}  // namespace

TEST_CASE("solve_nodal_velocity: no walls solves Eq. 4 exactly", "[bc]") {
  const Vec3 V = solve_nodal_velocity(kM, kB, {});
  check_vec(residual(kM, V, kB), {0.0, 0.0, 0.0});
}

TEST_CASE("solve_nodal_velocity: one wall", "[bc]") {
  // Lagrange-multiplier conditions: V . n = 0 and M V - B = -lambda n, i.e.
  // the residual has no component along the wall.
  const Vec3 n = normalized({1.0, 1.0, 0.0});
  const Vec3 V = solve_nodal_velocity(kM, kB, {n});
  CHECK_THAT(dot(V, n), WithinAbs(0.0, 1e-14));
  const Vec3 r = residual(kM, V, kB);
  const double rn = dot(r, n);
  check_vec(r, {rn * n[0], rn * n[1], rn * n[2]});

  // Hand check with n = x: V_x = 0 and rows y, z of M V = B hold,
  // 3 V_y - 0.4 V_z = -2, -0.4 V_y + 2 V_z = 0.7 -> det = 5.84,
  // V_y = -3.72 / 5.84, V_z = 1.3 / 5.84.
  const Vec3 Vx = solve_nodal_velocity(kM, kB, {{1.0, 0.0, 0.0}});
  check_vec(Vx, {0.0, -3.72 / 5.84, 1.3 / 5.84});
}

TEST_CASE("solve_nodal_velocity: two walls leave motion along the edge", "[bc]") {
  // Walls x and z: V = t V_t with t = y, t . (M V - B) = 0 -> V_y = -2 / 3.
  const Vec3 V = solve_nodal_velocity(kM, kB, {{1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}});
  check_vec(V, {0.0, -2.0 / 3.0, 0.0});
}

TEST_CASE("solve_nodal_velocity: three walls pin the node", "[bc]") {
  const Vec3 V = solve_nodal_velocity(
      kM, kB, {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
  check_vec(V, {0.0, 0.0, 0.0});
}

TEST_CASE("wall_normals: 2x2x1 slab, symmetry on every side", "[bc]") {
  const Mesh m = generate_structured_mesh({2, 2, 1});
  const BoundarySet all_sym = {kSym, kSym, kSym, kSym, kSym, kSym};
  const auto walls = wall_normals(m, all_sym);
  REQUIRE(walls.size() == 18);

  // Nodes (i, j, 0) = j*3 + i. Every node is on a z wall; the four z faces
  // around the center node count once.
  CHECK(wall_axes(walls[0]) == std::vector<int>{0, 1, 2});  // corner
  CHECK(wall_axes(walls[1]) == std::vector<int>{1, 2});     // -y edge
  CHECK(wall_axes(walls[3]) == std::vector<int>{0, 2});     // -x edge
  CHECK(wall_axes(walls[4]) == std::vector<int>{2});        // center
  CHECK(wall_axes(walls[9 + 8]) == std::vector<int>{0, 1, 2});  // (2,2,1)
}

TEST_CASE("wall_normals: outflow sides impose no wall", "[bc]") {
  const Mesh m = generate_structured_mesh({2, 2, 1});
  const BoundarySet sod = {kOut, kOut, kSym, kSym, kSym, kSym};
  const auto walls = wall_normals(m, sod);
  CHECK(wall_axes(walls[0]) == std::vector<int>{1, 2});
  CHECK(wall_axes(walls[3]) == std::vector<int>{2});

  const BoundarySet none = {kOut, kOut, kOut, kOut, kOut, kOut};
  for (const auto& w : wall_normals(m, none)) CHECK(w.empty());
}

TEST_CASE("wall_constraints: a piston prescribes the normal velocity", "[bc]") {
  // 2x2x1 slab, piston on -x moving with (1, 0.5, 0): only its normal part
  // (1, 0, 0) is prescribed. The other sides are symmetry walls, which
  // prescribe zero.
  const Mesh m = generate_structured_mesh({2, 2, 1});
  const BoundarySet bcs = {Boundary(BoundaryType::Piston, {1.0, 0.5, 0.0}),
                           kSym, kSym, kSym, kSym, kSym};
  const auto walls = wall_constraints(m, bcs);
  CHECK(wall_axes(walls.normals[3]) == std::vector<int>{0, 2});  // -x edge
  CHECK(wall_axes(walls.normals[4]) == std::vector<int>{2});     // center
  for (Index p = 0; p < m.num_nodes(); ++p) {
    INFO("node " << p);
    const bool on_piston = m.node_x[p] == 0.0;
    check_vec(walls.velocity[p], {on_piston ? 1.0 : 0.0, 0.0, 0.0});
  }
  CHECK(wall_normals(m, bcs) == walls.normals);
}

TEST_CASE("wall_constraints: piston meeting a wall at an angle", "[bc]") {
  // Unit cube with its y = 1 nodes shifted by 0.5 in x, so the -x face (a
  // piston, velocity (2, 0, 0)) meets the -y symmetry face at 63 degrees.
  // Gram-Schmidt must carry the piston value onto the orthogonalized
  // normal: at node 0 the prescribed part u0 is in the span of the two face
  // normals n1, n2 and has u0 . n1 = (2, 0, 0) . n1, u0 . n2 = 0.
  Mesh m = generate_structured_mesh({1, 1, 1});
  for (Index p = 0; p < m.num_nodes(); ++p) {
    if (m.node_y[p] == 1.0) m.node_x[p] += 0.5;
  }
  const BoundarySet bcs = {Boundary(BoundaryType::Piston, {2.0, 0.0, 0.0}),
                           kSym, kSym, kSym, kOut, kOut};
  const auto walls = wall_constraints(m, bcs);
  REQUIRE(walls.normals[0].size() == 2);
  const Vec3 n1 = normalized({-1.0, 0.5, 0.0});  // -x face, outward
  const Vec3 n2 = {0.0, -1.0, 0.0};              // -y face
  const Vec3& u0 = walls.velocity[0];
  CHECK_THAT(dot(u0, n1), WithinAbs(2.0 * n1[0], 1e-14));
  CHECK_THAT(dot(u0, n2), WithinAbs(0.0, 1e-14));
  CHECK_THAT(u0[2], WithinAbs(0.0, 1e-14));  // no component off the normals

  // A Newton step from u0 (tangent to the walls) keeps both constraints.
  const Vec3 step = solve_nodal_velocity(kM, kB, walls.normals[0]);
  const Vec3 V = {u0[0] + step[0], u0[1] + step[1], u0[2] + step[2]};
  CHECK_THAT(dot(V, n1), WithinAbs(2.0 * n1[0], 1e-14));
  CHECK_THAT(dot(V, n2), WithinAbs(0.0, 1e-14));
}

TEST_CASE("apply_pressure_bcs: outflow force on x-boundary nodes", "[bc]") {
  // 2x1x1 cells of size 0.5 x 2 x 3: each x face has area 6, so each of its
  // nodes gets S_pf n_pf = 1.5 along the outward normal.
  StructuredMeshSpec spec;
  spec.nx = 2;
  spec.x_max = 1.0;
  spec.y_max = 2.0;
  spec.z_max = 3.0;
  const Mesh m = generate_structured_mesh(spec);
  const BoundarySet sod = {kOut, kOut, kSym, kSym, kSym, kSym};
  const std::vector<double> p = {1.0, 0.1};

  std::vector<Vec3> B(m.num_nodes(), Vec3{0.25, -0.5, 0.0});
  apply_pressure_bcs(m, sod, p, B);

  // Nodes (i, j, k) = k*6 + j*3 + i.
  for (Index n = 0; n < m.num_nodes(); ++n) {
    INFO("node " << n);
    const int i = n % 3;
    const double bx = i == 0 ? 0.25 + 1.0 * 1.5 : i == 2 ? 0.25 - 0.1 * 1.5 : 0.25;
    check_vec(B[n], {bx, -0.5, 0.0});
  }
}

TEST_CASE("outflow + symmetry: a uniform state at rest stays at rest", "[bc]") {
  // Corner node (1, 0, 0) of a single unit cell at pressure P, velocity 0:
  // Eq. 4 over the cell's three faces at the corner gives
  //   B = P (S_+x n_+x + S_-y n_-y + S_-z n_-z) = P/4 (1, -1, -1),
  //   M = Z/4 diag(1, 1, 1).
  // The outflow BC cancels the x force; the y and z walls cancel the rest.
  const Mesh m = generate_structured_mesh({1, 1, 1});
  const BoundarySet bcs = {kOut, kOut, kSym, kSym, kSym, kSym};
  const double P = 2.5, Z = 3.0;
  std::vector<Vec3> B(m.num_nodes(), Vec3{P / 4, -P / 4, -P / 4});
  apply_pressure_bcs(m, bcs, {P}, B);
  const Index corner = 1;
  check_vec(B[corner], {0.0, -P / 4, -P / 4});

  const Mat3 M = {{{Z / 4, 0.0, 0.0}, {0.0, Z / 4, 0.0}, {0.0, 0.0, Z / 4}}};
  const Vec3 V = solve_nodal_velocity(M, B[corner], wall_normals(m, bcs)[corner]);
  check_vec(V, {0.0, 0.0, 0.0});
}
