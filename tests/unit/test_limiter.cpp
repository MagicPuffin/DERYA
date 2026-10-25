#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "gradients.hpp"
#include "hydro.hpp"
#include "limiter.hpp"
#include "mesh.hpp"

using Catch::Matchers::WithinAbs;
using hydro::bc::Mat3;
using hydro::lagrangian::Vec3;
using hydro::lagrangian::barth_jespersen;
using hydro::lagrangian::cell_centers;
using hydro::lagrangian::node_neighborhoods;
using hydro::lagrangian::pressure_limiters;
using hydro::lagrangian::velocity_limiters;
using hydro::mesh::Index;
using hydro::mesh::Mesh;
using hydro::mesh::StructuredMeshSpec;
using hydro::mesh::generate_structured_mesh;

namespace {

Mesh structured(Index nx, Index ny, Index nz) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  return generate_structured_mesh(spec);
}

// Velocity components of V_c = f(x_c) per cell.
template <class F>
void cell_velocities(const Mesh& m, F f, std::vector<double>& vx,
                     std::vector<double>& vy, std::vector<double>& vz) {
  const auto xc = cell_centers(m);
  vx.resize(m.num_cells());
  vy.resize(m.num_cells());
  vz.resize(m.num_cells());
  for (Index c = 0; c < m.num_cells(); ++c) {
    const Vec3 v = f(xc[c]);
    vx[c] = v[0];
    vy[c] = v[1];
    vz[c] = v[2];
  }
}

}  // namespace

TEST_CASE("barth_jespersen: Eq. 8 by hand", "[limiter]") {
  // Controls span [1, 4] around d = 2: room 2 above, 1 below.
  const std::vector<double> a = {2.0, 1.0, 4.0};
  CHECK(barth_jespersen(a, {2.5, 1.5, 3.9, 1.1}, 2.0) == 1.0);  // inside
  CHECK_THAT(barth_jespersen(a, {5.0}, 2.0), WithinAbs(2.0 / 3.0, 1e-15));  // above
  CHECK_THAT(barth_jespersen(a, {0.0}, 2.0), WithinAbs(0.5, 1e-15));  // below
  CHECK_THAT(barth_jespersen(a, {5.0, 0.0, 2.5}, 2.0), WithinAbs(0.5, 1e-15));
  CHECK(barth_jespersen(a, {2.0, 2.0}, 2.0) == 1.0);  // b_j = d
  CHECK(barth_jespersen(a, {}, 2.0) == 1.0);
  // d at the top of the range: anything above it is cut to zero.
  CHECK(barth_jespersen({1.0, 3.0}, {3.5, 2.0}, 3.0) == 0.0);
}

TEST_CASE("node_neighborhoods: cells sharing a node, c included", "[limiter]") {
  const Mesh m = structured(3, 3, 1);
  const auto nbr = node_neighborhoods(m);
  REQUIRE(nbr.size() == 9);
  CHECK(nbr[4] == std::vector<Index>{0, 1, 2, 3, 4, 5, 6, 7, 8});
  CHECK(nbr[0] == std::vector<Index>{0, 1, 3, 4});
  CHECK(nbr[1] == std::vector<Index>{0, 1, 2, 3, 4, 5});

  // 3D: the centre cell of a 3x3x3 block touches all 27.
  const auto nbr3 = node_neighborhoods(structured(3, 3, 3));
  CHECK(nbr3[13].size() == 27);
  CHECK(nbr3[0].size() == 8);
}

TEST_CASE("limiters: a linear field away from boundaries gives chi = 1",
          "[limiter]") {
  // Uniform 5x5x5 mesh; the centre 3x3x3 cells have a full neighborhood.
  const Mesh m = structured(5, 5, 5);
  const auto xc = cell_centers(m);
  const Vec3 a = {0.7, -1.3, 2.1};
  const Mat3 A = {{{1.0, -0.5, 0.25}, {2.0, 0.3, -1.1}, {-0.4, 0.9, 1.7}}};
  const Vec3 b = {0.2, -0.6, 1.5};

  std::vector<double> P(m.num_cells());
  for (Index c = 0; c < m.num_cells(); ++c) {
    P[c] = 4.0 + a[0] * xc[c][0] + a[1] * xc[c][1] + a[2] * xc[c][2];
  }
  std::vector<double> vx, vy, vz;
  cell_velocities(m, [&](const Vec3& x) {
    Vec3 v;
    for (int i = 0; i < 3; ++i) v[i] = b[i] + A[i][0] * x[0] + A[i][1] * x[1] + A[i][2] * x[2];
    return v;
  }, vx, vy, vz);
  const std::vector<Vec3> grad_p(m.num_cells(), a);
  const std::vector<Mat3> grad_v(m.num_cells(), A);

  const auto chi_p = pressure_limiters(m, P, grad_p);
  const auto chi_v = velocity_limiters(m, vx, vy, vz, grad_v, grad_p);
  for (Index k = 1; k < 4; ++k) {
    for (Index j = 1; j < 4; ++j) {
      for (Index i = 1; i < 4; ++i) {
        const Index c = k * 25 + j * 5 + i;
        INFO("cell " << c);
        CHECK_THAT(chi_p[c], WithinAbs(1.0, 1e-12));
        CHECK_THAT(chi_v[c], WithinAbs(1.0, 1e-12));
      }
    }
  }
}

TEST_CASE("pressure_limiters: a local extremum gives chi = 0", "[limiter]") {
  const Mesh m = structured(3, 3, 1);
  std::vector<double> P(9, 1.0);
  P[4] = 2.0;  // centre cell is a maximum
  std::vector<Vec3> grad_p(9, Vec3{});
  grad_p[4] = {1.0, 0.0, 0.0};
  CHECK(pressure_limiters(m, P, grad_p)[4] == 0.0);

  // A minimum too; a zero gradient needs no limiting.
  P[4] = 0.5;
  CHECK(pressure_limiters(m, P, grad_p)[4] == 0.0);
  grad_p[4] = Vec3{};
  CHECK(pressure_limiters(m, P, grad_p)[4] == 1.0);
}

TEST_CASE("limiters: a 1D field on a slab gives one chi per column",
          "[limiter]") {
  // The property SP-V is built for: a field varying in x only is limited the
  // same in every cell of a column, so the scheme stays 1D.
  StructuredMeshSpec spec;
  spec.nx = 10;
  spec.ny = 4;
  spec.y_max = 0.4;
  spec.z_max = 0.1;
  const Mesh m = generate_structured_mesh(spec);
  const auto xc = cell_centers(m);
  std::vector<double> P(m.num_cells());
  std::vector<Vec3> grad_p(m.num_cells());
  std::vector<Mat3> grad_v(m.num_cells(), Mat3{});
  for (Index c = 0; c < m.num_cells(); ++c) {
    const double x = xc[c][0];
    P[c] = 1.0 + std::sin(5.0 * x);
    grad_p[c] = {1.3 * 5.0 * std::cos(5.0 * x), 0.0, 0.0};  // over-steep
    grad_v[c][0][0] = -1.3 * 4.0 * std::sin(4.0 * x);
  }
  std::vector<double> vx, vy, vz;
  cell_velocities(m, [](const Vec3& x) {
    return Vec3{0.5 + std::cos(4.0 * x[0]), 0.0, 0.0};
  }, vx, vy, vz);

  const auto chi_p = pressure_limiters(m, P, grad_p);
  const auto chi_v = velocity_limiters(m, vx, vy, vz, grad_v, grad_p);
  bool limited_p = false, limited_v = false;
  for (Index i = 0; i < spec.nx; ++i) {
    for (Index j = 1; j < spec.ny; ++j) {
      INFO("column " << i << ", row " << j);
      CHECK(chi_p[j * spec.nx + i] == chi_p[i]);
      CHECK(chi_v[j * spec.nx + i] == chi_v[i]);
    }
    if (chi_p[i] > 0.0 && chi_p[i] < 1.0) limited_p = true;
    if (chi_v[i] > 0.0 && chi_v[i] < 1.0) limited_v = true;
  }
  CHECK(limited_p);  // the test exercises the limiter, not just chi = 1
  CHECK(limited_v);

  // No ghost cells: at x = 0 the extrapolation to the wall has no neighbor to
  // bound it, so a monotone field is cut to first order there.
  CHECK(chi_p[0] == 0.0);
}

TEST_CASE("velocity_limiters: projection direction and the V_c = 0 fallback",
          "[limiter]") {
  // 3x3 slab of [0,1]^2, h = 1/3. V_y = 0.1 (j - 1), so the centre cell's
  // neighbors span V_y in [-0.1, 0.1]. Its dV_y/dy = 0.9 extrapolates to
  // +-0.15 at its nodes (y offsets +-1/6): limited to 0.1/0.15 = 2/3 along y.
  const Mesh m = structured(3, 3, 1);
  std::vector<double> vx(9, 0.0), vy(9), vz(9, 0.0);
  for (Index c = 0; c < 9; ++c) vy[c] = 0.1 * static_cast<double>(c / 3 - 1);
  REQUIRE(vy[4] == 0.0);
  std::vector<Mat3> grad_v(9, Mat3{});
  grad_v[4][1][1] = 0.9;
  std::vector<Vec3> grad_p(9, Vec3{});

  SECTION("V_c = 0: along the pressure gradient") {
    grad_p[4] = {0.0, 1.0, 0.0};
    CHECK_THAT(velocity_limiters(m, vx, vy, vz, grad_v, grad_p)[4],
               WithinAbs(2.0 / 3.0, 1e-14));
    grad_p[4] = {0.0, -2.0, 0.0};
    CHECK_THAT(velocity_limiters(m, vx, vy, vz, grad_v, grad_p)[4],
               WithinAbs(2.0 / 3.0, 1e-14));
  }
  SECTION("V_c = 0 and no pressure gradient: along e_x, nothing to limit") {
    CHECK(velocity_limiters(m, vx, vy, vz, grad_v, grad_p)[4] == 1.0);
  }
  SECTION("V_c != 0: along V_c, not the pressure gradient") {
    for (Index c = 0; c < 9; ++c) vx[c] = 1.0;
    grad_p[4] = {0.0, 1.0, 0.0};
    CHECK(velocity_limiters(m, vx, vy, vz, grad_v, grad_p)[4] == 1.0);
    // Turned towards y (V_c = 1e-3 e_y), the overshoot shows: the binding
    // node is the top one, room 0.1 - 1e-3 against an extrapolation of 0.15.
    for (Index c = 0; c < 9; ++c) vx[c] = 0.0;
    vy[4] = 1e-3;
    grad_p[4] = {1.0, 0.0, 0.0};
    CHECK_THAT(velocity_limiters(m, vx, vy, vz, grad_v, grad_p)[4],
               WithinAbs((0.1 - 1e-3) / 0.15, 1e-14));
  }
}
