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
#include "hydro.hpp"
#include "mesh.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using hydro::bc::BoundarySet;
using hydro::bc::BoundaryType;
using hydro::eos::IdealGasEOS;
using hydro::lagrangian::CellThermo;
using hydro::lagrangian::HydroState;
using hydro::lagrangian::Vec3;
using hydro::lagrangian::assemble_nodal_system;
using hydro::lagrangian::cell_face_vectors;
using hydro::lagrangian::cell_thermo;
using hydro::lagrangian::cell_volume;
using hydro::lagrangian::corner_forces;
using hydro::lagrangian::corner_impedances;
using hydro::lagrangian::corner_vectors;
using hydro::lagrangian::lagrangian_step;
using hydro::lagrangian::make_state;
using hydro::lagrangian::nodal_velocities;
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

double total_energy(const HydroState& s) {
  double e = 0.0;
  for (Index c = 0; c < s.mesh.num_cells(); ++c) e += s.mass[c] * s.total_energy[c];
  return e;
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
