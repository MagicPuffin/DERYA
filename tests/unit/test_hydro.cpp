#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "bc.hpp"
#include "eos.hpp"
#include "gradients.hpp"
#include "hydro.hpp"
#include "limiter.hpp"
#include "mesh.hpp"
#include "scheme.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using hydro::bc::Boundary;
using hydro::bc::BoundarySet;
using hydro::bc::BoundaryType;
using hydro::bc::Mat3;
using hydro::bc::at_time;
using hydro::bc::pressure_boundary;
using hydro::eos::IdealGasEOS;
using hydro::lagrangian::CellFaceVectors;
using hydro::lagrangian::CellGradients;
using hydro::lagrangian::CellThermo;
using hydro::lagrangian::CornerState;
using hydro::lagrangian::EnergySource;
using hydro::lagrangian::HydroState;
using hydro::lagrangian::Limiter;
using hydro::lagrangian::Scheme;
using hydro::lagrangian::StageFluxes;
using hydro::lagrangian::Vec3;
using hydro::lagrangian::assemble_nodal_system;
using hydro::lagrangian::cell_face_vectors;
using hydro::lagrangian::cell_gradients;
using hydro::lagrangian::cell_thermo;
using hydro::lagrangian::cell_volume;
using hydro::lagrangian::corner_forces;
using hydro::lagrangian::corner_impedances;
using hydro::lagrangian::corner_vectors;
using hydro::lagrangian::first_order_corners;
using hydro::lagrangian::lagrangian_step;
using hydro::lagrangian::make_state;
using hydro::lagrangian::move_nodes;
using hydro::lagrangian::nodal_pressures;
using hydro::lagrangian::nodal_velocities;
using hydro::lagrangian::predictor_corrector_step;
using hydro::lagrangian::pressure_limiters;
using hydro::lagrangian::reconstructed_corners;
using hydro::lagrangian::scheme_step;
using hydro::lagrangian::second_order_fluxes;
using hydro::lagrangian::source_power;
using hydro::lagrangian::update_cells;
using hydro::lagrangian::velocity_limiters;
using hydro::mesh::Index;
using hydro::mesh::kFacesPerCell;
using hydro::mesh::kHexFaceNodes;
using hydro::mesh::kNodesPerCell;
using hydro::mesh::kNodesPerFace;
using hydro::mesh::Mesh;
using hydro::mesh::StructuredMeshSpec;
using hydro::mesh::generate_structured_mesh;

namespace {

constexpr BoundaryType kSym = BoundaryType::Symmetry;
constexpr BoundaryType kOut = BoundaryType::Outflow;
const BoundarySet kSodBcs = {kOut, kOut, kSym, kSym, kSym, kSym};
const BoundarySet kClosedBox = {kSym, kSym, kSym, kSym, kSym, kSym};

void check_vec(const Vec3& actual, const Vec3& expected, double tol = 1e-13) {
  for (int d = 0; d < 3; ++d) CHECK_THAT(actual[d], WithinAbs(expected[d], tol));
}

// 3x2x2 mesh of [0,3]x[0,1]x[0,0.5] with every interior node displaced by a
// smooth, non-affine perturbation, so faces are non-planar and cells differ.
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

// Uniform state of density rho, velocity v, specific internal energy eps.
HydroState uniform_state(Mesh m, double rho, const Vec3& v, double eps) {
  const auto n = static_cast<std::size_t>(m.num_cells());
  return make_state(std::move(m), std::vector<double>(n, rho),
                    std::vector<double>(n, v[0]), std::vector<double>(n, v[1]),
                    std::vector<double>(n, v[2]), std::vector<double>(n, eps));
}

double total_volume(const Mesh& m) {
  const auto faces = cell_face_vectors(m);
  double v = 0.0;
  for (Index c = 0; c < m.num_cells(); ++c) v += cell_volume(m, c, faces[c]);
  return v;
}

Vec3 total_momentum(const HydroState& s) {
  Vec3 p = {0.0, 0.0, 0.0};
  for (Index c = 0; c < s.mesh.num_cells(); ++c) {
    p[0] += s.mass[c] * s.vel_x[c];
    p[1] += s.mass[c] * s.vel_y[c];
    p[2] += s.mass[c] * s.vel_z[c];
  }
  return p;
}

// The cell and node data of two states are bit for bit the same.
void check_identical(const HydroState& a, const HydroState& b) {
  CHECK(a.vel_x == b.vel_x);
  CHECK(a.vel_y == b.vel_y);
  CHECK(a.vel_z == b.vel_z);
  CHECK(a.total_energy == b.total_energy);
  CHECK(a.mesh.node_x == b.mesh.node_x);
  CHECK(a.mesh.node_y == b.mesh.node_y);
  CHECK(a.mesh.node_z == b.mesh.node_z);
}

double total_energy(const HydroState& s) {
  double e = 0.0;
  for (Index c = 0; c < s.mesh.num_cells(); ++c) e += s.mass[c] * s.total_energy[c];
  return e;
}

// A smooth but non-uniform state on the distorted mesh.
HydroState wavy_state() {
  Mesh m = distorted_mesh();
  const auto n = static_cast<std::size_t>(m.num_cells());
  std::vector<double> rho(n), vx(n), vy(n), vz(n), eps(n);
  for (std::size_t c = 0; c < n; ++c) {
    rho[c] = 1.0 + 0.5 * std::sin(1.7 * c);
    vx[c] = 0.3 * std::cos(2.3 * c);
    vy[c] = -0.2 * std::sin(0.9 * c);
    vz[c] = 0.1 * std::cos(1.1 * c);
    eps[c] = 2.0 + std::cos(0.7 * c);
  }
  return make_state(std::move(m), rho, vx, vy, vz, eps);
}

// Arbitrary non-zero cell gradients, the same order of size as the state.
CellGradients wavy_gradients(Index num_cells) {
  CellGradients g;
  g.pressure.resize(num_cells);
  g.velocity.resize(num_cells);
  for (Index c = 0; c < num_cells; ++c) {
    for (int i = 0; i < 3; ++i) {
      g.pressure[c][i] = 0.4 * std::sin(1.3 * c + 2.1 * i);
      for (int j = 0; j < 3; ++j) {
        g.velocity[c][i][j] = 0.2 * std::cos(0.7 * c + 1.9 * i + 0.6 * j);
      }
    }
  }
  return g;
}

// The limited corner values of one second-order stage (decisions.md): a
// first-order solve, Green gradients from its nodal values, SP-V.
CornerState limited_corners(const HydroState& s,
                            const std::vector<CellFaceVectors>& faces,
                            const CellThermo& t, const BoundarySet& bcs) {
  const auto V = nodal_velocities(s, faces, t, bcs);
  const auto P = nodal_pressures(s, faces, t, corner_impedances(s, faces, t, V), V);
  const CellGradients g = cell_gradients(s.mesh, faces, P, V);
  return reconstructed_corners(
      s, t, g, pressure_limiters(s.mesh, t.pressure, g.pressure),
      velocity_limiters(s.mesh, s.vel_x, s.vel_y, s.vel_z, g.velocity, g.pressure));
}

// Forward-Euler step with the second-order nodal solve and corner forces on
// the given corner values.
std::vector<Vec3> corner_state_step(HydroState& s, const BoundarySet& bcs,
                                    const CornerState& corners,
                                    const std::vector<CellFaceVectors>& faces,
                                    const CellThermo& t, double dt) {
  const auto V = nodal_velocities(s, faces, t, corners, bcs);
  const auto Z = corner_impedances(s, faces, t, corners, V);
  update_cells(s, corner_forces(s, faces, corners, Z, V), V, dt);
  move_nodes(s.mesh, V, dt);
  return V;
}

}  // namespace

TEST_CASE("cell_face_vectors: owner and neighbor see opposite vectors", "[hydro]") {
  const Mesh m = distorted_mesh();
  const auto faces = cell_face_vectors(m);
  for (Index f = 0; f < m.num_faces(); ++f) {
    const auto vecs = hydro::mesh::face_area_vectors(m, f);
    for (Index c : {m.face_owner[f], m.face_neighbor[f]}) {
      if (c == hydro::mesh::kNoCell) continue;
      const double sign = c == m.face_owner[f] ? 1.0 : -1.0;
      int l = 0;
      while (m.cell_faces[c][l] != f) ++l;
      for (int n = 0; n < kNodesPerFace; ++n) {
        const Index p = m.cell_nodes[c][kHexFaceNodes[l][n]];
        int k = 0;
        while (m.face_nodes[f][k] != p) ++k;
        INFO("face " << f << ", cell " << c << ", node " << p);
        check_vec(faces[c][l][n],
                  {sign * vecs[k][0], sign * vecs[k][1], sign * vecs[k][2]}, 0.0);
      }
    }
  }
}

TEST_CASE("corner_vectors: closure and GCL (n_cp = dV_c/dx_p)", "[hydro]") {
  // V_c is cubic in the node positions; a central difference with step h has
  // error O(h^2), so h = 1e-5 resolves n_cp to ~1e-10.
  Mesh m = distorted_mesh();
  const auto faces = cell_face_vectors(m);
  const double h = 1e-5;
  for (Index c = 0; c < m.num_cells(); ++c) {
    const auto corners = corner_vectors(faces[c]);
    Vec3 closure = {0.0, 0.0, 0.0};
    for (const Vec3& n : corners) {
      for (int d = 0; d < 3; ++d) closure[d] += n[d];
    }
    check_vec(closure, {0.0, 0.0, 0.0});

    for (int k = 0; k < kNodesPerCell; ++k) {
      const Index p = m.cell_nodes[c][k];
      std::vector<double>* coord[3] = {&m.node_x, &m.node_y, &m.node_z};
      Vec3 grad;
      for (int d = 0; d < 3; ++d) {
        double& x = (*coord[d])[p];
        const double x0 = x;
        x = x0 + h;
        const double v_plus = cell_volume(m, c, cell_face_vectors(m)[c]);
        x = x0 - h;
        const double v_minus = cell_volume(m, c, cell_face_vectors(m)[c]);
        x = x0;
        grad[d] = (v_plus - v_minus) / (2.0 * h);
      }
      INFO("cell " << c << ", corner " << k);
      check_vec(corners[k], grad, 1e-9);
    }
  }
}

TEST_CASE("cell_face_vectors: triangular prisms on a polar axis", "[hydro]") {
  // Polar 3x4x2 down to r = 0: the i = 0 cells have one node twice per level
  // (corners 0 and 3, 4 and 7), so a face can list a node twice, (a, a, q, r).
  // cell_face_vectors matches corners to entries by node id, giving both
  // corners the first a's entry; that is right because the two entries are
  // equal (both (a - r) x (q - a) / 8). The corner vectors close, the corners
  // at a node sum to dV_c/dx_p, and the volumes are the triangles' times dz.
  StructuredMeshSpec spec{3, 4, 2, 0.0, 1.0, 0.0, 90.0, 0.0, 0.1};
  spec.mapping = hydro::mesh::MeshMapping::Polar;
  Mesh m = generate_structured_mesh(spec);
  const double pi = std::acos(-1.0);
  const double h = 1e-6;
  for (Index c = 0; c < m.num_cells(); ++c) {
    INFO("cell " << c);
    const auto corners = corner_vectors(cell_face_vectors(m)[c]);
    Vec3 closure = {0.0, 0.0, 0.0};
    for (const Vec3& n : corners) {
      for (int d = 0; d < 3; ++d) closure[d] += n[d];
    }
    check_vec(closure, {0.0, 0.0, 0.0}, 1e-16);
    const double r0 = (c % 3) / 3.0, r1 = (c % 3 + 1) / 3.0;
    CHECK_THAT(cell_volume(m, c, cell_face_vectors(m)[c]),
               WithinAbs(0.5 * (r1 * r1 - r0 * r0) * std::sin(pi / 8.0) * 0.05, 1e-16));
    for (int k = 0; k < kNodesPerCell; ++k) {
      const Index p = m.cell_nodes[c][k];
      Vec3 sum = {0.0, 0.0, 0.0};
      for (int k2 = 0; k2 < kNodesPerCell; ++k2) {
        if (m.cell_nodes[c][k2] != p) continue;
        for (int d = 0; d < 3; ++d) sum[d] += corners[k2][d];
      }
      std::vector<double>* coord[3] = {&m.node_x, &m.node_y, &m.node_z};
      Vec3 grad;
      for (int d = 0; d < 3; ++d) {
        double& x = (*coord[d])[p];
        const double x0 = x;
        x = x0 + h;
        const double v_plus = cell_volume(m, c, cell_face_vectors(m)[c]);
        x = x0 - h;
        const double v_minus = cell_volume(m, c, cell_face_vectors(m)[c]);
        x = x0;
        grad[d] = (v_plus - v_minus) / (2.0 * h);
      }
      INFO("corner " << k);
      check_vec(sum, grad, 1e-10);
    }
  }
}

TEST_CASE("cell_volume: box cells and a non-planar cell", "[hydro]") {
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 2;
  spec.nz = 2;
  spec.x_max = 3.0;
  spec.y_max = 1.0;
  spec.z_max = 0.5;
  const Mesh box = generate_structured_mesh(spec);
  const auto box_faces = cell_face_vectors(box);
  for (Index c = 0; c < box.num_cells(); ++c) {
    CHECK_THAT(cell_volume(box, c, box_faces[c]), WithinAbs(0.125, 1e-15));
  }

  // The prismatoid cell of test_mesh.cpp: [0,2]x[0,1]x[0,1] with bottom node 3
  // moved to (1,1,0), volume 1.75.
  StructuredMeshSpec one;
  one.x_max = 2.0;
  Mesh m = generate_structured_mesh(one);
  m.node_x[3] = 1.0;
  CHECK_THAT(cell_volume(m, 0, cell_face_vectors(m)[0]), WithinAbs(1.75, 1e-14));

  // The distorted mesh fills its box.
  CHECK_THAT(total_volume(distorted_mesh()), WithinAbs(1.5, 1e-13));
}

TEST_CASE("make_state and cell_thermo round-trip", "[hydro]") {
  const IdealGasEOS eos(1.4);
  const HydroState s = uniform_state(distorted_mesh(), 2.0, {0.3, -0.4, 1.2}, 2.5);
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  double mass = 0.0;
  for (Index c = 0; c < s.mesh.num_cells(); ++c) {
    INFO("cell " << c);
    mass += s.mass[c];
    CHECK_THAT(s.total_energy[c], WithinAbs(2.5 + 0.5 * (0.09 + 0.16 + 1.44), 1e-14));
    CHECK_THAT(t.density[c], WithinRel(2.0, 1e-13));
    CHECK_THAT(t.pressure[c], WithinRel(0.4 * 2.0 * 2.5, 1e-13));
    CHECK_THAT(t.sound_speed[c], WithinRel(std::sqrt(1.4 * 2.0 / 2.0), 1e-13));
    CHECK_THAT(t.shock_coefficient[c], WithinRel(1.2, 1e-15));
  }
  CHECK_THAT(mass, WithinRel(2.0 * 1.5, 1e-13));

  // Folding a cell inside out is caught.
  HydroState tangled = uniform_state(generate_structured_mesh({1, 1, 1}), 1.0,
                                     {0.0, 0.0, 0.0}, 1.0);
  for (Index p = 4; p < 8; ++p) tangled.mesh.node_z[p] = -1.0;
  CHECK_THROWS_AS(cell_thermo(tangled, cell_face_vectors(tangled.mesh), eos),
                  std::runtime_error);
}

TEST_CASE("cell_thermo: negative eps within round-off of the flow is clamped",
          "[hydro]") {
  // Cold gas at rest next to a shocked one (Saltzman at second order, Task
  // 15): E_c is a sum of energy fluxes of the flow's size, so a cold cell
  // holds round-off of that size. eps = -2e-14 there, with no kinetic energy
  // of its own, is clamped when another cell has E = 1, an error when no
  // cell does, and -1e-10 is always an error.
  const IdealGasEOS eos(5.0 / 3.0);
  HydroState s = uniform_state(generate_structured_mesh({3, 1, 1}), 1.0,
                               {0.0, 0.0, 0.0}, 0.0);
  const auto faces = cell_face_vectors(s.mesh);
  s.total_energy[0] = 1.0;
  s.total_energy[2] = -2e-14;
  const CellThermo t = cell_thermo(s, faces, eos);
  CHECK(t.pressure[2] == 0.0);
  CHECK(t.pressure[0] > 0.0);

  s.total_energy[0] = 0.0;
  CHECK_THROWS_AS(cell_thermo(s, faces, eos), std::runtime_error);
  s.total_energy[0] = 1.0;
  s.total_energy[2] = -1e-10;
  CHECK_THROWS_AS(cell_thermo(s, faces, eos), std::runtime_error);
}

TEST_CASE("corner_impedances and assemble_nodal_system: hand-computed corner",
          "[hydro]") {
  // Corner node 1 = (1,0,0) of a unit cell touches its +x, -y and -z faces,
  // each giving S_pf = 1/4 along the outward axis. With V_p = V_c + (dv,0,0)
  // the jump is seen only by the +x face:
  //   Z_x = rho (a + Gamma dv), Z_0 = rho a on the -y and -z faces,
  // and Eq. 4 gives
  //   M = 1/4 diag(Z_x, Z_0, Z_0),
  //   B = P/4 (1, -1, -1) + 1/4 (Z_x v_x, Z_0 v_y, Z_0 v_z).
  const IdealGasEOS eos(1.4);
  const Vec3 v = {0.3, -0.2, 0.1};
  const double dv = 0.5;
  const HydroState s = uniform_state(generate_structured_mesh({1, 1, 1}), 1.5, v, 2.0);
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  std::vector<Vec3> V(s.mesh.num_nodes(), v);
  V[1][0] += dv;
  const auto Z = corner_impedances(s, faces, t, V);
  const double rho = 1.5, a = t.sound_speed[0], P = t.pressure[0];
  const double Zx = rho * (a + 1.2 * dv), Z0 = rho * a;
  for (int l = 0; l < kFacesPerCell; ++l) {
    for (int n = 0; n < kNodesPerFace; ++n) {
      INFO("face " << l << ", node " << n);
      const bool at_1 = s.mesh.cell_nodes[0][kHexFaceNodes[l][n]] == 1;
      const bool x_face = std::abs(faces[0][l][n][0]) > 0.0;
      CHECK_THAT(Z[0][l][n], WithinRel(at_1 && x_face ? Zx : Z0, 1e-14));
    }
  }

  const auto sys = assemble_nodal_system(s, faces, t, Z);
  const Vec3 diag = {Zx / 4, Z0 / 4, Z0 / 4};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      CHECK_THAT(sys.M[1][i][j], WithinAbs(i == j ? diag[i] : 0.0, 1e-15));
    }
  }
  check_vec(sys.B[1], {P / 4 + Zx / 4 * v[0], -P / 4 + Z0 / 4 * v[1],
                       -P / 4 + Z0 / 4 * v[2]});
}

TEST_CASE("nodal_velocities: uniform flow gives V_p = V_c at every node", "[hydro]") {
  // Interior nodes: B_p = P sum_c n_cp + M_p V = M_p V, as the corner vectors
  // around the node cancel. Outflow nodes: the remaining boundary face vectors
  // are cancelled by P* = P. Symmetry walls are tangent to V.
  const IdealGasEOS eos(1.4);
  const Vec3 v = {0.7, 0.0, 0.0};
  const HydroState s = uniform_state(distorted_mesh(), 1.0, v, 2.5);
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  const auto V = nodal_velocities(s, faces, t, kSodBcs);
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    INFO("node " << p);
    check_vec(V[p], v, 1e-12);
  }

  // Same with every side outflow and an oblique velocity.
  const Vec3 w = {0.3, -0.5, 0.2};
  const HydroState s2 = uniform_state(distorted_mesh(), 1.0, w, 2.5);
  const CellThermo t2 = cell_thermo(s2, faces, eos);
  const BoundarySet all_out = {kOut, kOut, kOut, kOut, kOut, kOut};
  const auto V2 = nodal_velocities(s2, faces, t2, all_out);
  for (Index p = 0; p < s2.mesh.num_nodes(); ++p) {
    INFO("node " << p);
    check_vec(V2[p], w, 1e-12);
  }
}

TEST_CASE("nodal_velocities: uniform cold flow gives V_p = V_c", "[hydro]") {
  // eps = 0: a = 0, P = 0, and no velocity jumps, so every Z_cfp is zero and
  // M_p is singular. The forces do not depend on V_p; it stays at the mean of
  // the surrounding cell velocities, i.e. V_c.
  const IdealGasEOS eos(1.4);
  const Vec3 w = {0.3, -0.5, 0.2};
  const HydroState s = uniform_state(distorted_mesh(), 1.0, w, 0.0);
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  const BoundarySet all_out = {kOut, kOut, kOut, kOut, kOut, kOut};
  const auto V = nodal_velocities(s, faces, t, all_out);
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    INFO("node " << p);
    check_vec(V[p], w, 1e-14);
  }
}

TEST_CASE("nodal_velocities: colliding cold slabs give the two-shock velocity",
          "[hydro]") {
  // 4x1x1 unit cells, symmetry y/z walls, outflow x ends. Cold gas (eps = 0):
  // rho = 4, u = 1 in the two left cells; rho = 1, u = 0 in the two right
  // ones. At the interface nodes (x = 2) only the x faces see a jump, and
  // Eq. 4 with the two-shock Z reads
  //   rho_L Gamma (1 - V)^2 = rho_R Gamma V^2  =>  2 (1 - V) = V,  V = 2/3,
  // the strong-shock (Rankine-Hugoniot) momentum balance. Nodes inside each
  // slab move with it.
  const IdealGasEOS eos(1.4);
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.x_max = 4.0;
  const HydroState s = make_state(generate_structured_mesh(spec), {4.0, 4.0, 1.0, 1.0},
                                  {1.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
                                  {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0});
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  const auto V = nodal_velocities(s, faces, t, kSodBcs);
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    INFO("node " << p << " at x = " << s.mesh.node_x[p]);
    const double x = s.mesh.node_x[p];
    const double expected = x < 2.0 ? 1.0 : (x > 2.0 ? 0.0 : 2.0 / 3.0);
    check_vec(V[p], {expected, 0.0, 0.0}, 1e-12);
  }
}

TEST_CASE("lagrangian_step: a piston into a cold gas pushes with the shock pressure",
          "[hydro]") {
  // 2x1x1 unit cubes of cold gas at rest (rho = 1, gamma = 5/3), piston on -x
  // moving with (1, 0.5, 0), symmetry elsewhere. Every piston node is also on
  // the y and z walls, so V_p = (1, 0, 0) there. The interior nodes see only
  // cold gas at rest, so they stay at rest. On the piston face s = (V_p -
  // V_c) . n = -1 and Eq. 3 gives P = rho Gamma s^2 = 4/3, the exact
  // post-shock pressure of the piston problem; the other faces at those
  // nodes have s = 0. So cell 0 gains u = (4/3) dt (unit mass and face area)
  // and the energy the piston does as work, (4/3) dt.
  const IdealGasEOS eos(5.0 / 3.0);
  StructuredMeshSpec spec;
  spec.nx = 2;
  spec.x_max = 2.0;
  HydroState s = uniform_state(generate_structured_mesh(spec), 1.0, {0.0, 0.0, 0.0}, 0.0);
  const BoundarySet bcs = {Boundary(BoundaryType::Piston, {1.0, 0.5, 0.0}),
                           kSym, kSym, kSym, kSym, kSym};
  const double dt = 0.01;
  const auto V = lagrangian_step(s, bcs, eos, dt);
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    INFO("node " << p);
    // node_x has moved by dt V_p.
    const bool on_piston = s.mesh.node_x[p] < 0.5;
    check_vec(V[p], {on_piston ? 1.0 : 0.0, 0.0, 0.0});
  }
  CHECK_THAT(s.vel_x[0], WithinAbs(4.0 / 3.0 * dt, 1e-14));
  CHECK_THAT(s.vel_y[0], WithinAbs(0.0, 1e-14));
  CHECK_THAT(s.vel_x[1], WithinAbs(0.0, 1e-14));
  CHECK_THAT(total_energy(s), WithinAbs(4.0 / 3.0 * dt, 1e-14));
}

TEST_CASE("nodal_velocities: subnormal velocity jumps stay finite", "[hydro]") {
  // Cold gas ahead of a shock on a skewed mesh carries velocity precursors
  // that decay to subnormal values (Saltzman: 1e-315 and below). The
  // Jacobian at such a node is subnormal too, and once 1e-8 tr(J) is below
  // the smallest subnormal the regularization underflows to zero and the
  // elimination gives 0/0. 3x1x1 cold cells with
  // densities 4, 1, 1 moving at 1, 1e-320 and 0, all sides outflow (no
  // walls, so the full 3x3 system). The real jump between cells 0 and 1
  // (whose solution, 2/3, is not the initial guess 1/2) keeps Newton
  // iterating, so it also steps at the subnormal nodes at x = 2.
  const IdealGasEOS eos(5.0 / 3.0);
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.x_max = 3.0;
  const HydroState s =
      make_state(generate_structured_mesh(spec), {4.0, 1.0, 1.0}, {1.0, 1e-320, 0.0},
                 {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  const BoundarySet all_out = {kOut, kOut, kOut, kOut, kOut, kOut};
  const auto V = nodal_velocities(s, faces, t, all_out);
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    INFO("node " << p << " at x = " << s.mesh.node_x[p]);
    for (int d = 0; d < 3; ++d) CHECK(std::isfinite(V[p][d]));
    if (s.mesh.node_x[p] > 1.5) check_vec(V[p], {0.0, 0.0, 0.0}, 1e-300);
  }
}

TEST_CASE("corner_forces: P_cfp = P_c when V_p = V_c", "[hydro]") {
  // Eq. 3 with no velocity jump leaves F_cp = P_c n_cp.
  const IdealGasEOS eos(1.4);
  const Vec3 v = {0.3, -0.5, 0.2};
  const HydroState s = uniform_state(distorted_mesh(), 1.0, v, 2.5);
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  const std::vector<Vec3> V(s.mesh.num_nodes(), v);
  const auto F = corner_forces(s, faces, t, corner_impedances(s, faces, t, V), V);
  for (Index c = 0; c < s.mesh.num_cells(); ++c) {
    const auto n = corner_vectors(faces[c]);
    for (int k = 0; k < kNodesPerCell; ++k) {
      INFO("cell " << c << ", corner " << k);
      check_vec(F[c][k], {t.pressure[c] * n[k][0], t.pressure[c] * n[k][1],
                          t.pressure[c] * n[k][2]});
    }
  }
}

TEST_CASE("lagrangian_step: uniform state at rest stays at rest", "[hydro]") {
  const IdealGasEOS eos(1.4);
  HydroState s = uniform_state(distorted_mesh(), 1.0, {0.0, 0.0, 0.0}, 2.5);
  const Mesh initial = s.mesh;
  for (int step = 0; step < 10; ++step) lagrangian_step(s, kSodBcs, eos, 0.01);
  for (Index c = 0; c < s.mesh.num_cells(); ++c) {
    INFO("cell " << c);
    check_vec({s.vel_x[c], s.vel_y[c], s.vel_z[c]}, {0.0, 0.0, 0.0}, 1e-14);
    CHECK_THAT(s.total_energy[c], WithinAbs(2.5, 1e-13));
  }
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    INFO("node " << p);
    check_vec({s.mesh.node_x[p], s.mesh.node_y[p], s.mesh.node_z[p]},
              {initial.node_x[p], initial.node_y[p], initial.node_z[p]}, 1e-15);
  }
}

TEST_CASE("lagrangian_step: closed box conserves mass, volume and energy",
          "[hydro]") {
  // Non-uniform state on the distorted mesh inside symmetry walls. Energy is
  // conserved exactly by construction: at an interior node the corner forces
  // sum to zero (Eq. 4), at a wall node they sum to a wall-normal reaction
  // that does no work since V_p . n = 0. Wall nodes slide within the box
  // faces, so the total volume is fixed too.
  const IdealGasEOS eos(1.4);
  Mesh m = distorted_mesh();
  const auto n = static_cast<std::size_t>(m.num_cells());
  std::vector<double> rho(n), vx(n), vy(n), vz(n), eps(n);
  for (std::size_t c = 0; c < n; ++c) {
    rho[c] = 1.0 + 0.5 * std::sin(1.7 * c);
    vx[c] = 0.3 * std::cos(2.3 * c);
    vy[c] = -0.2 * std::sin(0.9 * c);
    vz[c] = 0.1 * std::cos(1.1 * c);
    eps[c] = 2.0 + std::cos(0.7 * c);
  }
  HydroState s = make_state(std::move(m), rho, vx, vy, vz, eps);
  const double e0 = total_energy(s);
  const double v0 = total_volume(s.mesh);
  const auto mass0 = s.mass;

  for (int step = 0; step < 50; ++step) lagrangian_step(s, kClosedBox, eos, 0.005);

  CHECK(s.mass == mass0);
  CHECK_THAT(total_energy(s), WithinRel(e0, 1e-13));
  CHECK_THAT(total_volume(s.mesh), WithinRel(v0, 1e-13));
  // The state did evolve.
  double moved = 0.0;
  for (Index c = 0; c < s.mesh.num_cells(); ++c) {
    moved = std::max(moved, std::abs(s.vel_x[c] - vx[c]));
  }
  CHECK(moved > 1e-3);
}

TEST_CASE("lagrangian_step: Sod on a thin slab", "[hydro][integration]") {
  // Sod on 100x2x2 cells of [0,1]x[0,0.02]x[0,0.02] with outflow x ends and
  // symmetry y/z walls, fixed dt = 5e-4 to t = 0.2 (acoustic CFL ~0.1).
  // Checks: the flow stays 1D (cells of a column agree, no transverse
  // velocity) and density is close to the exact solution in L1. Measured L1
  // errors with dt ~ dx/20: 0.0298, 0.0202, 0.0131, 0.0084 for nx = 50, 100,
  // 200, 400 (order ~0.6, as expected for first order with a contact and a
  // shock), with the acoustic Z; the two-shock Z gives 0.0202 at nx = 100
  // too (re-measured at that resolution only). The bound is loose regression; the Sod acceptance test (Task 6)
  // runs through the driver.
  const Index nx = 100;
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = 2;
  spec.nz = 2;
  spec.y_max = 0.02;
  spec.z_max = 0.02;
  Mesh m = generate_structured_mesh(spec);
  const auto n = static_cast<std::size_t>(m.num_cells());
  std::vector<double> rho(n), zero(n, 0.0), eps(n);
  const IdealGasEOS eos(1.4);
  for (std::size_t c = 0; c < n; ++c) {
    const bool left = static_cast<Index>(c) % nx < nx / 2;
    const double r = left ? 1.0 : 0.125;
    const double p = left ? 1.0 : 0.1;
    rho[c] = r;
    eps[c] = p / (0.4 * r);
  }
  HydroState s = make_state(std::move(m), rho, zero, zero, zero, eps);

  const double dt = 5e-4;
  for (int step = 0; step < 400; ++step) lagrangian_step(s, kSodBcs, eos, dt);

  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  for (Index c = 0; c < s.mesh.num_cells(); ++c) {
    INFO("cell " << c);
    CHECK_THAT(s.vel_y[c], WithinAbs(0.0, 1e-12));
    CHECK_THAT(s.vel_z[c], WithinAbs(0.0, 1e-12));
    CHECK_THAT(t.density[c], WithinRel(t.density[c % nx], 1e-10));
  }

  std::ifstream in(ORACLE_DATA_DIR "/sod.json");
  REQUIRE(in);
  const auto oracle = nlohmann::json::parse(in);
  REQUIRE(oracle["t"].get<double>() == 0.2);
  const auto x_ex = oracle["x"].get<std::vector<double>>();
  const auto rho_ex = oracle["density"].get<std::vector<double>>();
  auto exact = [&](double x) {
    std::size_t i = 1;
    while (i + 1 < x_ex.size() && x_ex[i] < x) ++i;
    const double w = (x - x_ex[i - 1]) / (x_ex[i] - x_ex[i - 1]);
    return (1.0 - w) * rho_ex[i - 1] + w * rho_ex[i];
  };
  double l1 = 0.0;
  for (Index i = 0; i < nx; ++i) {
    // Column i, bottom row: corners 0 and 1 give the cell's x extent.
    const auto& cn = s.mesh.cell_nodes[i];
    const double x0 = s.mesh.node_x[cn[0]], x1 = s.mesh.node_x[cn[1]];
    l1 += std::abs(t.density[i] - exact(0.5 * (x0 + x1))) * (x1 - x0);
  }
  INFO("L1 density error " << l1);
  CHECK(l1 < 0.025);
}

TEST_CASE("lagrangian_step: planar Noh with a cold gas", "[hydro][integration]") {
  // Cold (eps = 0) ideal gas, gamma = 5/3, rho = 1, u = -1 on 100x1x1 cells of
  // [0,1]x[0,0.01]x[0,0.01]: symmetry wall at x = 0, outflow (P* = 0 while
  // the gas there is cold) at x = 1, symmetry y/z walls. Fixed dt = 5e-4 to
  // t = 0.6 (CFL ~0.1 on the post-shock a + |u|). Exact solution: a shock at
  // x = D t = 0.2 (D = 1/3), gas at rest behind it with rho = 4, eps = 1/2,
  // and the untouched inflow ahead of it. The acoustic Z would make the
  // first nodal solve singular. Neither wall does work, so total energy is
  // conserved. Measured: plateau rho 3.993-4.000 on [0.04, 0.19], wall
  // heating in the 2 cells at the wall (rho 2.54, the classic Noh error),
  // shock 3 cells wide at x = 0.201, and L1 density errors 0.0376, 0.0188,
  // 0.0094, 0.0047 for nx = 50, 100, 200, 400 (dt = dx/20), first order.
  const Index nx = 100;
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.y_max = 0.01;
  spec.z_max = 0.01;
  const IdealGasEOS eos(5.0 / 3.0);
  HydroState s = uniform_state(generate_structured_mesh(spec), 1.0, {-1.0, 0.0, 0.0}, 0.0);
  const BoundarySet bcs = {kSym, kOut, kSym, kSym, kSym, kSym};
  const double e0 = total_energy(s);

  const double dt = 5e-4;
  for (int step = 0; step < 1200; ++step) lagrangian_step(s, bcs, eos, dt);

  CHECK_THAT(total_energy(s), WithinRel(e0, 1e-12));
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  double l1 = 0.0;
  for (Index c = 0; c < nx; ++c) {
    const auto& cn = s.mesh.cell_nodes[c];
    const double x0 = s.mesh.node_x[cn[0]], x1 = s.mesh.node_x[cn[1]];
    const double xc = 0.5 * (x0 + x1);
    const double eps = s.total_energy[c] - 0.5 * s.vel_x[c] * s.vel_x[c];
    INFO("cell " << c << " at x = " << xc);
    if (xc > 0.04 && xc < 0.18) {
      CHECK_THAT(t.density[c], WithinRel(4.0, 0.01));
      CHECK_THAT(eps, WithinRel(0.5, 0.01));
      CHECK_THAT(s.vel_x[c], WithinAbs(0.0, 1e-4));
    }
    if (xc > 0.25) {
      CHECK_THAT(t.density[c], WithinRel(1.0, 1e-10));
      CHECK_THAT(s.vel_x[c], WithinAbs(-1.0, 1e-10));
    }
    l1 += std::abs(t.density[c] - (xc < 0.2 ? 4.0 : 1.0)) * (x1 - x0);
  }
  INFO("L1 density error " << l1);
  CHECK(l1 < 0.025);
}

TEST_CASE("reconstructed_corners: limited linear extrapolation to the nodes",
          "[hydro]") {
  // Box cell [0,2]x[0,1]x[0,0.5]: x_c = (1, 0.5, 0.25), node 0 at the origin,
  // node 6 at (2, 1, 0.5), so x_p - x_c = -/+(1, 0.5, 0.25).
  StructuredMeshSpec spec;
  spec.x_max = 2.0;
  spec.y_max = 1.0;
  spec.z_max = 0.5;
  const IdealGasEOS eos(1.4);
  const HydroState s = uniform_state(generate_structured_mesh(spec), 1.0,
                                     {0.1, 0.2, 0.3}, 2.5);  // P_c = 1
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  CellGradients g;
  g.pressure = {{1.0, 2.0, 4.0}};  // . (x_p - x_c) = -/+3
  g.velocity = {Mat3{{{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 4.0}}}};
  const CornerState r = reconstructed_corners(s, t, g, {0.5}, {0.25});
  REQUIRE(s.mesh.cell_nodes[0][0] == 0);
  CHECK_THAT(r.pressure[0][0], WithinAbs(1.0 - 1.5, 1e-15));
  CHECK_THAT(r.pressure[0][6], WithinAbs(1.0 + 1.5, 1e-15));
  check_vec(r.velocity[0][0], {0.1 - 0.25, 0.2 - 0.25, 0.3 - 0.25}, 1e-15);
  check_vec(r.velocity[0][6], {0.1 + 0.25, 0.2 + 0.25, 0.3 + 0.25}, 1e-15);

  // Zero factors or zero gradients give the first-order corners exactly.
  const CornerState first = first_order_corners(s, t);
  const CornerState unlimited_zero = reconstructed_corners(
      s, t, CellGradients{{Vec3{}}, {Mat3{}}}, {1.0}, {1.0});
  const CornerState limited_away = reconstructed_corners(s, t, g, {0.0}, {0.0});
  for (const CornerState* c : {&unlimited_zero, &limited_away}) {
    CHECK(c->pressure == first.pressure);
    CHECK(c->velocity == first.velocity);
  }
}

TEST_CASE("second-order solve: P_cfp = P~_c(x_p) when V_p = V~_c(x_p)", "[hydro]") {
  // One cell, so each node has a single corner: with V_p equal to the corner
  // velocity there is no jump, Z_cfp = rho_c a_c, and Eq. 2.44 leaves
  // F_cp = P~_c(x_p) n_cp.
  StructuredMeshSpec spec;
  const IdealGasEOS eos(1.4);
  HydroState s = uniform_state(generate_structured_mesh(spec), 1.0, {0.1, 0.2, 0.3}, 2.5);
  s.mesh.node_x[6] += 0.2;  // a non-planar cell
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  const CornerState r = reconstructed_corners(s, t, wavy_gradients(1), {1.0}, {1.0});
  std::vector<Vec3> V(s.mesh.num_nodes());
  for (int k = 0; k < kNodesPerCell; ++k) V[s.mesh.cell_nodes[0][k]] = r.velocity[0][k];

  const auto Z = corner_impedances(s, faces, t, r, V);
  const double rho_a = t.density[0] * t.sound_speed[0];
  for (int l = 0; l < kFacesPerCell; ++l) {
    for (int n = 0; n < kNodesPerFace; ++n) CHECK_THAT(Z[0][l][n], WithinRel(rho_a, 1e-14));
  }
  const auto F = corner_forces(s, faces, r, Z, V);
  const auto n = corner_vectors(faces[0]);
  for (int k = 0; k < kNodesPerCell; ++k) {
    INFO("corner " << k);
    const double P = r.pressure[0][k];
    check_vec(F[0][k], {P * n[k][0], P * n[k][1], P * n[k][2]}, 1e-14);
  }
  // The first-order jump, against V_c, is not zero at these V_p.
  const auto Z1 = corner_impedances(s, faces, t, V);
  CHECK(Z1[0][0][0] > rho_a * (1.0 + 1e-3));
}

TEST_CASE("second-order solve: first-order corners reproduce first order exactly",
          "[hydro]") {
  // Zero gradients (or limiters) must leave the first-order scheme bit for
  // bit: the order = 1 regression check of every v0.2 task.
  const IdealGasEOS eos(1.4);
  const HydroState s = wavy_state();
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  const Index nc = s.mesh.num_cells();
  const CornerState r = reconstructed_corners(
      s, t, wavy_gradients(nc), std::vector<double>(nc, 0.0), std::vector<double>(nc, 0.0));

  const auto V1 = nodal_velocities(s, faces, t, kSodBcs);
  const auto V2 = nodal_velocities(s, faces, t, r, kSodBcs);
  CHECK(V2 == V1);
  const auto Z1 = corner_impedances(s, faces, t, V1);
  const auto Z2 = corner_impedances(s, faces, t, r, V2);
  CHECK(Z2 == Z1);
  CHECK(corner_forces(s, faces, r, Z2, V2) == corner_forces(s, faces, t, Z1, V1));
  CHECK(nodal_pressures(s, faces, r, Z2, V2) == nodal_pressures(s, faces, t, Z1, V1));

  // ... and non-zero corner values do change the solution.
  const CornerState unlimited = reconstructed_corners(
      s, t, wavy_gradients(nc), std::vector<double>(nc, 1.0), std::vector<double>(nc, 1.0));
  const auto V3 = nodal_velocities(s, faces, t, unlimited, kSodBcs);
  double change = 0.0;
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    for (int d = 0; d < 3; ++d) change = std::max(change, std::abs(V3[p][d] - V1[p][d]));
  }
  CHECK(change > 1e-3);
}

TEST_CASE("second-order solve: uniform flow gives V_p = V_c", "[hydro]") {
  // The full stage on a distorted mesh: first-order gradients vanish to
  // round-off, so the corner values are V_c and P_c to round-off, and the
  // solve keeps V_p = V_c. Warm gas only: in a cold gas that is uniform
  // everywhere, the round-off left in V~_c(x_p) (SP-V limits only its part
  // along V_c) makes the two-shock jumps ~1e-16 and the whole flux scale
  // ~1e-33, and Newton does not converge on that degenerate root (see
  // backlog, Task 12). Cold regions of real runs are judged against the
  // fluxes elsewhere.
  const IdealGasEOS eos(1.4);
  const Vec3 v = {0.3, -0.5, 0.2};
  const BoundarySet all_out = {kOut, kOut, kOut, kOut, kOut, kOut};
  const HydroState s = uniform_state(distorted_mesh(), 1.0, v, 2.5);
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo t = cell_thermo(s, faces, eos);
  const CornerState r = limited_corners(s, faces, t, all_out);
  const auto V = nodal_velocities(s, faces, t, r, all_out);
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    INFO("node " << p);
    check_vec(V[p], v);
  }
}

TEST_CASE("second-order solve: closed box conserves mass, volume and energy",
          "[hydro]") {
  // Conservation rests on the nodal balance (Eq. 4), whatever the corner
  // values: the corner forces at a node sum to zero, or to a wall-normal
  // reaction on a wall. Checked with unlimited, arbitrary corner values (the
  // hardest case) and with the limited stage.
  const IdealGasEOS eos(1.4);
  for (const bool limited : {false, true}) {
    INFO((limited ? "limited" : "unlimited"));
    HydroState s = wavy_state();
    const auto vx0 = s.vel_x;
    const double e0 = total_energy(s);
    const double v0 = total_volume(s.mesh);
    const auto mass0 = s.mass;
    const auto walls = hydro::bc::wall_normals(s.mesh, kClosedBox);
    const Index nc = s.mesh.num_cells();

    for (int step = 0; step < 50; ++step) {
      const auto faces = cell_face_vectors(s.mesh);
      const CellThermo t = cell_thermo(s, faces, eos);
      const CornerState r =
          limited ? limited_corners(s, faces, t, kClosedBox)
                  : reconstructed_corners(s, t, wavy_gradients(nc),
                                          std::vector<double>(nc, 1.0),
                                          std::vector<double>(nc, 1.0));
      if (step == 0) {
        // Nodal balance: the tangential part of sum_c F_cp vanishes.
        const auto V = nodal_velocities(s, faces, t, r, kClosedBox);
        const auto F = corner_forces(s, faces, r, corner_impedances(s, faces, t, r, V), V);
        std::vector<Vec3> sum(s.mesh.num_nodes(), Vec3{});
        double scale = 0.0;
        for (Index c = 0; c < nc; ++c) {
          for (int k = 0; k < kNodesPerCell; ++k) {
            for (int d = 0; d < 3; ++d) {
              sum[s.mesh.cell_nodes[c][k]][d] += F[c][k][d];
              scale = std::max(scale, std::abs(F[c][k][d]));
            }
          }
        }
        for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
          Vec3 r_p = sum[p];
          for (const Vec3& w : walls[p]) {
            const double rn = r_p[0] * w[0] + r_p[1] * w[1] + r_p[2] * w[2];
            for (int d = 0; d < 3; ++d) r_p[d] -= rn * w[d];
          }
          INFO("node " << p);
          check_vec(r_p, {0.0, 0.0, 0.0}, 1e-12 * scale);
        }
      }
      corner_state_step(s, kClosedBox, r, faces, t, 0.005);
    }

    CHECK(s.mass == mass0);
    CHECK_THAT(total_energy(s), WithinRel(e0, 1e-13));
    CHECK_THAT(total_volume(s.mesh), WithinRel(v0, 1e-13));
    double moved = 0.0;
    for (Index c = 0; c < nc; ++c) moved = std::max(moved, std::abs(s.vel_x[c] - vx0[c]));
    CHECK(moved > 1e-3);
  }
}

TEST_CASE("scheme_step: order 1 is lagrangian_step, bit for bit", "[hydro]") {
  const IdealGasEOS eos(1.4);
  HydroState a = wavy_state();
  HydroState b = a;
  for (int step = 0; step < 5; ++step) {
    const auto Va = lagrangian_step(a, kClosedBox, eos, 0.005);
    const auto Vb = scheme_step(b, kClosedBox, eos, 0.0, 0.005, Scheme{});
    CHECK(Va == Vb);
  }
  check_identical(a, b);
  CHECK(Scheme{}.order == 1);
  HydroState c = a;
  CHECK_THROWS_AS(scheme_step(c, kClosedBox, eos, 0.0, 0.005, Scheme{3, Limiter::SpV}),
                  std::invalid_argument);
}

TEST_CASE("predictor_corrector_step: thesis Eqs. 2.45-2.49", "[hydro]") {
  // The predictor is a forward-Euler step with the fluxes of s^n; the
  // corrector starts again from s^n and applies each stage's fluxes over
  // dt / 2: nodes move with the mean V_p, cells with the mean forces and work.
  const IdealGasEOS eos(1.4);
  const double dt = 0.01;
  for (const Limiter limiter : {Limiter::SpV, Limiter::None}) {
    INFO((limiter == Limiter::SpV ? "sp-v" : "none"));
    const HydroState s0 = wavy_state();
    HydroState s = s0;
    const auto V = predictor_corrector_step(s, kClosedBox, eos, 0.0, dt, limiter);

    const StageFluxes f0 = second_order_fluxes(s0, kClosedBox, eos, limiter);
    HydroState predicted = s0;
    update_cells(predicted, f0.forces, f0.node_velocity, dt);
    move_nodes(predicted.mesh, f0.node_velocity, dt);
    const StageFluxes f1 = second_order_fluxes(predicted, kClosedBox, eos, limiter);

    const Index nc = s0.mesh.num_cells();
    for (Index c = 0; c < nc; ++c) {
      INFO("cell " << c);
      Vec3 force = {0.0, 0.0, 0.0};
      double work = 0.0;
      for (int k = 0; k < kNodesPerCell; ++k) {
        const Index p = s0.mesh.cell_nodes[c][k];
        for (int d = 0; d < 3; ++d) {
          force[d] += 0.5 * (f0.forces[c][k][d] + f1.forces[c][k][d]);
          work += 0.5 * (f0.forces[c][k][d] * f0.node_velocity[p][d] +
                         f1.forces[c][k][d] * f1.node_velocity[p][d]);
        }
      }
      const double scale = dt / s0.mass[c];
      CHECK_THAT(s.vel_x[c], WithinAbs(s0.vel_x[c] - scale * force[0], 1e-14));
      CHECK_THAT(s.vel_y[c], WithinAbs(s0.vel_y[c] - scale * force[1], 1e-14));
      CHECK_THAT(s.vel_z[c], WithinAbs(s0.vel_z[c] - scale * force[2], 1e-14));
      CHECK_THAT(s.total_energy[c], WithinAbs(s0.total_energy[c] - scale * work, 1e-14));
    }
    for (Index p = 0; p < s0.mesh.num_nodes(); ++p) {
      INFO("node " << p);
      Vec3 mean;
      for (int d = 0; d < 3; ++d) {
        mean[d] = 0.5 * (f0.node_velocity[p][d] + f1.node_velocity[p][d]);
      }
      CHECK(V[p] == mean);
      check_vec({s.mesh.node_x[p], s.mesh.node_y[p], s.mesh.node_z[p]},
                {s0.mesh.node_x[p] + dt * mean[0], s0.mesh.node_y[p] + dt * mean[1],
                 s0.mesh.node_z[p] + dt * mean[2]},
                1e-15);
    }
    // The stages differ, so the corrector is not the predictor again.
    double differ = 0.0;
    for (Index p = 0; p < s0.mesh.num_nodes(); ++p) {
      for (int d = 0; d < 3; ++d) {
        differ = std::max(differ, std::abs(f1.node_velocity[p][d] - f0.node_velocity[p][d]));
      }
    }
    CHECK(differ > 1e-3);
  }
}

TEST_CASE("predictor_corrector_step: closed box conserves mass, momentum, energy",
          "[hydro]") {
  // Each stage's corner forces balance at every interior node and reduce to
  // a wall-normal reaction, doing no work, on the walls (Task 12), so the
  // mean of two stages does too. In a closed box the walls push, so total
  // momentum changes, but only by the wall nodes' impulse.
  const IdealGasEOS eos(1.4);
  const double dt = 0.005;
  for (const Limiter limiter : {Limiter::SpV, Limiter::None}) {
    INFO((limiter == Limiter::SpV ? "sp-v" : "none"));
    HydroState s = wavy_state();
    const auto vx0 = s.vel_x;
    const double e0 = total_energy(s);
    const double v0 = total_volume(s.mesh);
    const auto mass0 = s.mass;
    const auto walls = hydro::bc::wall_normals(s.mesh, kClosedBox);

    // Momentum over the first step: the impulse of the forces at wall nodes.
    {
      const StageFluxes f0 = second_order_fluxes(s, kClosedBox, eos, limiter);
      HydroState predicted = s;
      update_cells(predicted, f0.forces, f0.node_velocity, dt);
      move_nodes(predicted.mesh, f0.node_velocity, dt);
      const StageFluxes f1 = second_order_fluxes(predicted, kClosedBox, eos, limiter);
      Vec3 impulse = {0.0, 0.0, 0.0};
      double scale = 0.0;
      for (Index c = 0; c < s.mesh.num_cells(); ++c) {
        for (int k = 0; k < kNodesPerCell; ++k) {
          if (walls[s.mesh.cell_nodes[c][k]].empty()) continue;
          for (int d = 0; d < 3; ++d) {
            impulse[d] -= 0.5 * dt * (f0.forces[c][k][d] + f1.forces[c][k][d]);
            scale = std::max(scale, dt * std::abs(f0.forces[c][k][d]));
          }
        }
      }
      const Vec3 p0 = total_momentum(s);
      HydroState stepped = s;
      predictor_corrector_step(stepped, kClosedBox, eos, 0.0, dt, limiter);
      const Vec3 p1 = total_momentum(stepped);
      check_vec({p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]}, impulse, 1e-13 * scale);
    }

    for (int step = 0; step < 50; ++step) {
      predictor_corrector_step(s, kClosedBox, eos, 0.0, dt, limiter);
    }
    CHECK(s.mass == mass0);
    CHECK_THAT(total_energy(s), WithinRel(e0, 1e-13));
    CHECK_THAT(total_volume(s.mesh), WithinRel(v0, 1e-13));
    double moved = 0.0;
    for (Index c = 0; c < s.mesh.num_cells(); ++c) {
      moved = std::max(moved, std::abs(s.vel_x[c] - vx0[c]));
    }
    CHECK(moved > 1e-3);
  }
}

TEST_CASE("predictor_corrector_step: uniform flow stays uniform", "[hydro]") {
  // Warm uniform flow through all-outflow boundaries on the distorted mesh:
  // every stage keeps V_p = V_c, so the mesh translates rigidly.
  const IdealGasEOS eos(1.4);
  const Vec3 v = {0.3, -0.5, 0.2};
  const BoundarySet all_out = {kOut, kOut, kOut, kOut, kOut, kOut};
  HydroState s = uniform_state(distorted_mesh(), 1.0, v, 2.5);
  const Mesh initial = s.mesh;
  const double dt = 0.01;
  for (int step = 0; step < 5; ++step) {
    const auto V = predictor_corrector_step(s, all_out, eos, 0.0, dt, Limiter::SpV);
    for (Index p = 0; p < s.mesh.num_nodes(); ++p) check_vec(V[p], v);
  }
  for (Index c = 0; c < s.mesh.num_cells(); ++c) {
    INFO("cell " << c);
    check_vec({s.vel_x[c], s.vel_y[c], s.vel_z[c]}, v);
    CHECK_THAT(s.total_energy[c], WithinRel(2.5 + 0.5 * (0.09 + 0.25 + 0.04), 1e-13));
  }
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    INFO("node " << p);
    check_vec({s.mesh.node_x[p], s.mesh.node_y[p], s.mesh.node_z[p]},
              {initial.node_x[p] + 0.05 * v[0], initial.node_y[p] + 0.05 * v[1],
               initial.node_z[p] + 0.05 * v[2]});
  }
}

TEST_CASE("scheme_step: a uniform energy source heats a gas at rest uniformly",
          "[hydro]") {
  // S = 3 everywhere on a uniform gas at rest in a closed box: no pressure
  // gradient appears, so the gas stays at rest and E_c gains dt S / rho_c per
  // step, at both orders.
  const IdealGasEOS eos(1.4);
  const EnergySource source = [](const Vec3&) { return 3.0; };
  const double rho = 2.0, eps = 1.5, dt = 0.01;
  for (const int order : {1, 2}) {
    INFO("order " << order);
    HydroState s = uniform_state(distorted_mesh(), rho, {0.0, 0.0, 0.0}, eps);
    for (int step = 0; step < 5; ++step) {
      scheme_step(s, kClosedBox, eos, 0.0, dt, Scheme{order, Limiter::None}, source);
    }
    for (Index c = 0; c < s.mesh.num_cells(); ++c) {
      INFO("cell " << c);
      check_vec({s.vel_x[c], s.vel_y[c], s.vel_z[c]}, {0.0, 0.0, 0.0});
      CHECK_THAT(s.total_energy[c], WithinRel(eps + 5 * dt * 3.0 / rho, 1e-13));
    }
  }
}

TEST_CASE("scheme_step: an energy source enters each stage", "[hydro]") {
  // In a closed box the corner forces do no net work, so the total energy
  // changes only by the source: dt sum_c V_c S(x_c) at t^n's geometry at
  // order 1, and the mean of that and the predicted geometry's at order 2.
  const IdealGasEOS eos(1.4);
  const EnergySource source = [](const Vec3& x) { return 1.0 + x[0] - 2.0 * x[1] * x[2]; };
  const double dt = 0.02;
  auto sum = [](const std::vector<double>& v) {
    double t = 0.0;
    for (const double x : v) t += x;
    return t;
  };

  const HydroState s0 = wavy_state();
  const double power0 = sum(source_power(s0, cell_face_vectors(s0.mesh), source));
  REQUIRE(power0 > 0.1);

  HydroState first = s0;
  scheme_step(first, kClosedBox, eos, 0.0, dt, Scheme{1, Limiter::None}, source);
  CHECK_THAT(total_energy(first) - total_energy(s0), WithinRel(dt * power0, 1e-12));

  // The predicted geometry: the nodes of s^n moved by the first stage's V_p.
  const StageFluxes f0 = second_order_fluxes(s0, kClosedBox, eos, Limiter::SpV);
  HydroState predicted = s0;
  move_nodes(predicted.mesh, f0.node_velocity, dt);
  const double power1 =
      sum(source_power(predicted, cell_face_vectors(predicted.mesh), source));
  // Large enough against the tolerance to tell the stages apart.
  REQUIRE(std::abs(power1 - power0) > 1e-6 * power0);

  HydroState second = s0;
  scheme_step(second, kClosedBox, eos, 0.0, dt, Scheme{2, Limiter::SpV}, source);
  CHECK_THAT(total_energy(second) - total_energy(s0),
             WithinRel(0.5 * dt * (power0 + power1), 1e-10));

  // No source is the step without one, bit for bit.
  HydroState with_empty = s0, without = s0;
  scheme_step(with_empty, kClosedBox, eos, 0.0, dt, Scheme{2, Limiter::SpV}, EnergySource{});
  predictor_corrector_step(without, kClosedBox, eos, 0.0, dt, Limiter::SpV);
  check_identical(with_empty, without);
}

namespace {

// 4x1x1 cells of [0,1]x[0,0.25]x[0,0.25].
Mesh tube_mesh() {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.x_max = 1.0;
  spec.y_max = 0.25;
  spec.z_max = 0.25;
  return generate_structured_mesh(spec);
}

}  // namespace

TEST_CASE("scheme_step: a pressure boundary at the gas pressure keeps a state at rest",
          "[hydro]") {
  // P = (gamma - 1) rho eps = 1 everywhere, P* = 1 on both x sides (one a
  // constant, one a law), on the distorted mesh: nothing moves at either order.
  const IdealGasEOS eos(1.4);
  const BoundarySet bcs = {pressure_boundary(1.0),
                           pressure_boundary([](double) { return 1.0; }),
                           kSym, kSym, kSym, kSym};
  for (const int order : {1, 2}) {
    INFO("order " << order);
    HydroState s = uniform_state(distorted_mesh(), 1.0, {0.0, 0.0, 0.0}, 2.5);
    const Mesh initial = s.mesh;
    double t = 0.0;
    for (int step = 0; step < 5; ++step) {
      const auto V = scheme_step(s, bcs, eos, t, 0.01, Scheme{order, Limiter::SpV});
      t += 0.01;
      for (Index p = 0; p < s.mesh.num_nodes(); ++p) check_vec(V[p], {0.0, 0.0, 0.0}, 1e-14);
    }
    for (Index c = 0; c < s.mesh.num_cells(); ++c) {
      INFO("cell " << c);
      check_vec({s.vel_x[c], s.vel_y[c], s.vel_z[c]}, {0.0, 0.0, 0.0}, 1e-14);
      CHECK_THAT(s.total_energy[c], WithinRel(2.5, 1e-14));
    }
    for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
      INFO("node " << p);
      check_vec({s.mesh.node_x[p], s.mesh.node_y[p], s.mesh.node_z[p]},
                {initial.node_x[p], initial.node_y[p], initial.node_z[p]}, 1e-15);
    }
  }
}

TEST_CASE("lagrangian_step: a pressure above the gas's pushes the boundary in",
          "[hydro]") {
  // Gas at rest, rho = 1, P = 1, a = sqrt(1.4), Gamma = 1.2, walls but for
  // P* = P + dP on x_max. At an x_max node only the x_max face acts along x
  // (the y and z faces are walls), so Eq. 4 reads Z |V| = dP with the
  // two-shock Z = rho (a + Gamma |V|):
  //   |V| = (-a + sqrt(a^2 + 4 Gamma dP / rho)) / (2 Gamma),
  // inward (the acoustic rate dP / (rho a) for small dP). The x_max face's
  // P_cfp is then P*, so the last cell gets dV = -dt dP / (rho dx) and
  // dE = dt P* |V| / (rho dx), the boundary's work; the others stay at rest.
  const IdealGasEOS eos(1.4);
  const double rho = 1.0, P = 1.0, dP = 0.1, dt = 0.01, dx = 0.25;
  const double a = std::sqrt(1.4), gamma_s = 1.2;
  const double speed =
      (-a + std::sqrt(a * a + 4.0 * gamma_s * dP / rho)) / (2.0 * gamma_s);
  const BoundarySet bcs = {kSym, pressure_boundary(P + dP), kSym, kSym, kSym, kSym};
  HydroState s = uniform_state(tube_mesh(), rho, {0.0, 0.0, 0.0}, P / (0.4 * rho));
  const Mesh initial = s.mesh;
  const auto V = lagrangian_step(s, bcs, eos, dt);
  for (Index p = 0; p < s.mesh.num_nodes(); ++p) {
    INFO("node " << p);
    const bool pushed = initial.node_x[p] == 1.0;
    check_vec(V[p], {pushed ? -speed : 0.0, 0.0, 0.0}, 1e-12);
    CHECK_THAT(s.mesh.node_x[p], WithinAbs(initial.node_x[p] - (pushed ? dt * speed : 0.0),
                                           1e-14));
  }
  for (Index c = 0; c < 4; ++c) {
    INFO("cell " << c);
    const bool last = c == 3;
    CHECK_THAT(s.vel_x[c], WithinAbs(last ? -dt * dP / (rho * dx) : 0.0, 1e-14));
    CHECK_THAT(s.total_energy[c],
               WithinAbs(P / (0.4 * rho) + (last ? dt * (P + dP) * speed / (rho * dx) : 0.0),
                         1e-12));
  }
}

TEST_CASE("scheme_step: pressure laws are evaluated at each stage's time", "[hydro]") {
  // Thesis Eq. 2.42: the predictor uses F(t^n, U^n), the corrector stage
  // F(t^n + dt, U^(1)). P*(t) = 1 + 20 t on x_max, recording its calls.
  const IdealGasEOS eos(1.4);
  const double t = 0.3, dt = 0.01;
  std::vector<double> times;
  const auto law = [](double time) { return 1.0 + 20.0 * time; };
  const BoundarySet bcs = {kSym,
                           pressure_boundary([&times, law](double time) {
                             times.push_back(time);
                             return law(time);
                           }),
                           kSym, kSym, kSym, kSym};
  const HydroState s0 = uniform_state(tube_mesh(), 1.0, {0.0, 0.0, 0.0}, 2.5);

  HydroState first = s0;
  scheme_step(first, bcs, eos, t, dt, Scheme{1, Limiter::SpV});
  CHECK(times == std::vector<double>{t});
  HydroState by_hand = s0;
  lagrangian_step(by_hand, {kSym, pressure_boundary(law(t)), kSym, kSym, kSym, kSym}, eos,
                  dt);
  check_identical(first, by_hand);
  // The step itself takes evaluated boundaries only.
  HydroState unevaluated = s0;
  CHECK_THROWS_AS(lagrangian_step(unevaluated, bcs, eos, dt), std::logic_error);

  times.clear();
  HydroState second = s0;
  const auto V = scheme_step(second, bcs, eos, t, dt, Scheme{2, Limiter::SpV});
  CHECK(times == std::vector<double>{t, t + dt});

  // The same step by hand, each stage with its constant P*.
  const BoundarySet at_tn = {kSym, pressure_boundary(law(t)), kSym, kSym, kSym, kSym};
  const BoundarySet at_tn1 = {kSym, pressure_boundary(law(t + dt)), kSym, kSym, kSym, kSym};
  const StageFluxes f0 = second_order_fluxes(s0, at_tn, eos, Limiter::SpV);
  HydroState predicted = s0;
  update_cells(predicted, f0.forces, f0.node_velocity, dt);
  move_nodes(predicted.mesh, f0.node_velocity, dt);
  const StageFluxes f1 = second_order_fluxes(predicted, at_tn1, eos, Limiter::SpV);
  const StageFluxes f1_at_tn = second_order_fluxes(predicted, at_tn, eos, Limiter::SpV);
  double differ = 0.0;
  for (Index p = 0; p < s0.mesh.num_nodes(); ++p) {
    INFO("node " << p);
    Vec3 mean;
    for (int d = 0; d < 3; ++d) {
      mean[d] = 0.5 * (f0.node_velocity[p][d] + f1.node_velocity[p][d]);
      differ = std::max(differ, std::abs(f1.node_velocity[p][d] - f1_at_tn.node_velocity[p][d]));
    }
    CHECK(V[p] == mean);
  }
  // P* rises by 0.2 over the step, enough to tell the corrector's time apart.
  CHECK(differ > 1e-3);
}
