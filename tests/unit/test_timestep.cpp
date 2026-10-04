#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "bc.hpp"
#include "eos.hpp"
#include "hydro.hpp"
#include "mesh.hpp"
#include "sync.hpp"
#include "timestep.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using hydro::eos::IdealGasEOS;
using hydro::lagrangian::HydroState;
using hydro::lagrangian::Vec3;
using hydro::lagrangian::cell_face_vectors;
using hydro::lagrangian::cell_volume;
using hydro::lagrangian::lagrangian_step;
using hydro::lagrangian::make_state;
using hydro::mesh::Index;
using hydro::mesh::Mesh;
using hydro::mesh::StructuredMeshSpec;
using hydro::mesh::generate_structured_mesh;
using hydro::sync::ReduceOp;
using hydro::sync::global_reduce;
using hydro::timestep::TimestepParams;
using hydro::timestep::characteristic_length;
using hydro::timestep::stable_dt;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Ideal gas, gamma = 1.4: a^2 = gamma (gamma - 1) eps.
const IdealGasEOS kEos(1.4);
double eps_for_sound_speed(double a) { return a * a / (1.4 * 0.4); }

// One box cell [0,2]x[0,1]x[0,0.5]: V = 1, largest face 2x1, lambda = 0.5.
Mesh box_cell() {
  StructuredMeshSpec spec;
  spec.x_max = 2.0;
  spec.z_max = 0.5;
  return generate_structured_mesh(spec);
}

HydroState state(Mesh m, const std::vector<double>& eps, const Vec3& v = {0, 0, 0}) {
  const auto n = static_cast<std::size_t>(m.num_cells());
  return make_state(std::move(m), std::vector<double>(n, 1.0),
                    std::vector<double>(n, v[0]), std::vector<double>(n, v[1]),
                    std::vector<double>(n, v[2]), eps);
}

std::vector<Vec3> node_positions(const Mesh& m) {
  std::vector<Vec3> x(m.num_nodes());
  for (Index p = 0; p < m.num_nodes(); ++p) x[p] = {m.node_x[p], m.node_y[p], m.node_z[p]};
  return x;
}

}  // namespace

TEST_CASE("global_reduce is the identity on one rank", "[timestep]") {
  CHECK(global_reduce(0.25, ReduceOp::Min) == 0.25);
  CHECK(global_reduce(-3.0, ReduceOp::Max) == -3.0);
  CHECK(global_reduce(7.5, ReduceOp::Sum) == 7.5);
}

TEST_CASE("characteristic_length: volume over the largest face area", "[timestep]") {
  const Mesh m = box_cell();
  const auto faces = cell_face_vectors(m);
  // V = 1, largest face 2 x 1: lambda = 0.5, the shortest edge.
  CHECK_THAT(characteristic_length(faces[0], cell_volume(m, 0, faces[0])),
             WithinRel(0.5, 1e-14));
}

TEST_CASE("stable_dt: acoustic CFL on one box cell", "[timestep]") {
  // a = 1, lambda = 0.5: dt = C_cfl lambda / a = 0.25 * 0.5.
  const HydroState s = state(box_cell(), {eps_for_sound_speed(1.0)});
  CHECK_THAT(stable_dt(s, kEos, {}, kInf), WithinRel(0.125, 1e-14));
  // At rest the volume criterion adds nothing.
  const std::vector<Vec3> rest(s.mesh.num_nodes(), Vec3{0, 0, 0});
  CHECK_THAT(stable_dt(s, kEos, rest, kInf), WithinRel(0.125, 1e-14));

  TimestepParams p;
  p.cfl = 0.5;
  CHECK_THAT(stable_dt(s, kEos, {}, kInf, p), WithinRel(0.25, 1e-14));
}

TEST_CASE("stable_dt: the smallest cell limit wins", "[timestep]") {
  // Two unit cubes (lambda = 1), a = 1 and a = 2: dt = 0.25 * 1 / 2.
  StructuredMeshSpec spec;
  spec.nx = 2;
  spec.x_max = 2.0;
  const HydroState s = state(generate_structured_mesh(spec),
                             {eps_for_sound_speed(1.0), eps_for_sound_speed(2.0)});
  CHECK_THAT(stable_dt(s, kEos, {}, kInf), WithinRel(0.125, 1e-14));
}

TEST_CASE("stable_dt: volume criterion bounds a cold cell", "[timestep]") {
  // eps = 0, so a = 0 and no CFL limit. With V_p = x_p the cell expands
  // homogeneously: dV/dt = sum_p n_cp . x_p = 3 V, so dt = C_V V / (3 V).
  // Compression, V_p = -x_p, gives the same limit.
  const HydroState s = state(box_cell(), {0.0});
  auto V = node_positions(s.mesh);
  CHECK_THAT(stable_dt(s, kEos, V, kInf), WithinRel(0.1 / 3.0, 1e-13));
  for (auto& v : V) v = {-v[0], -v[1], -v[2]};
  CHECK_THAT(stable_dt(s, kEos, V, kInf), WithinRel(0.1 / 3.0, 1e-13));
}

TEST_CASE("stable_dt: growth limit and dt_max", "[timestep]") {
  const HydroState s = state(box_cell(), {eps_for_sound_speed(1.0)});
  // CFL gives 0.125; growth from 0.1 caps it at 0.101.
  CHECK_THAT(stable_dt(s, kEos, {}, 0.1), WithinRel(0.101, 1e-14));
  TimestepParams p;
  p.dt_max = 0.05;
  CHECK_THAT(stable_dt(s, kEos, {}, kInf, p), WithinRel(0.05, 1e-14));
}

TEST_CASE("stable_dt: cold uniform flow has no CFL or volume limit", "[timestep]") {
  // a = 0 and the cells translate rigidly (dV/dt = 0): only the growth limit
  // and dt_max apply, and with neither there is no finite dt.
  StructuredMeshSpec spec;
  spec.nx = 3;
  const Vec3 v = {1.0, -0.5, 0.25};
  const HydroState s = state(generate_structured_mesh(spec), {0.0, 0.0, 0.0}, v);
  const std::vector<Vec3> V(s.mesh.num_nodes(), v);
  CHECK_THAT(stable_dt(s, kEos, V, 0.2), WithinRel(0.202, 1e-14));
  TimestepParams p;
  p.dt_max = 0.1;
  CHECK_THAT(stable_dt(s, kEos, V, kInf, p), WithinRel(0.1, 1e-14));
  CHECK_THROWS_AS(stable_dt(s, kEos, V, kInf), std::runtime_error);
}

TEST_CASE("lagrangian_step returns the nodal velocities it used", "[timestep]") {
  // Uniform flow: every V_p equals the flow velocity.
  const Vec3 v = {0.7, 0.0, 0.0};
  StructuredMeshSpec spec;
  spec.nx = 2;
  HydroState s = state(generate_structured_mesh(spec),
                       {eps_for_sound_speed(1.0), eps_for_sound_speed(1.0)}, v);
  using hydro::bc::BoundaryType;
  const hydro::bc::BoundarySet bcs = {BoundaryType::Outflow,  BoundaryType::Outflow,
                                      BoundaryType::Symmetry, BoundaryType::Symmetry,
                                      BoundaryType::Symmetry, BoundaryType::Symmetry};
  const auto V = lagrangian_step(s, bcs, kEos, 0.01);
  REQUIRE(V.size() == static_cast<std::size_t>(s.mesh.num_nodes()));
  for (const Vec3& Vp : V) {
    for (int d = 0; d < 3; ++d) CHECK_THAT(Vp[d], WithinAbs(v[d], 1e-13));
  }
}
