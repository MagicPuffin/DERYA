#include "hydro.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace hydro::lagrangian {

using mesh::Index;
using mesh::kFacesPerCell;
using mesh::kHexFaceNodes;
using mesh::kNodesPerCell;
using mesh::kNodesPerFace;

namespace {

double dot(const Vec3& a, const Vec3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 node_position(const mesh::Mesh& m, Index p) {
  return {m.node_x[p], m.node_y[p], m.node_z[p]};
}

Vec3 cell_velocity(const HydroState& s, Index c) {
  return {s.vel_x[c], s.vel_y[c], s.vel_z[c]};
}

// Relative round-off allowed in eps = E - |V|^2/2 (see cell_thermo).
constexpr double kEnergyRoundoff = 1e-13;

// Newton solve of the nodal system (see nodal_velocities). Converged when the
// tangential residual is below kNewtonTolerance times the largest nodal flux
// sum_c sum_f S_pf (|P_c| + Z_cfp |s_cfp|) over the mesh, so that cold,
// quiescent nodes are judged against the flow as a whole. Also converged when
// the step is at round-off, kStepRoundoff max_c (|V_c| + a_c): where the
// Jacobian is singular the residual and the flux scale vanish together.
constexpr double kStepRoundoff = 1e-14;
constexpr int kMaxNewtonIterations = 50;
constexpr double kNewtonTolerance = 1e-13;
// Relative regularization of the Newton Jacobian, for nodes where it is
// singular (cold cells moving alike): the step then has no component along
// the null space, so V_p keeps its initial value there.
constexpr double kNewtonRegularization = 1e-8;

// rho_c (a_c + k Gamma_c |s_cfp|), s_cfp = (V_p - V_c) . n_pf: the impedance
// Z_cfp for k = 1, and the Newton Jacobian's d(Z_cfp s_cfp)/ds_cfp for k = 2.
std::vector<CornerValues> scaled_impedances(const HydroState& s,
                                            const std::vector<CellFaceVectors>& faces,
                                            const CellThermo& thermo,
                                            const std::vector<Vec3>& node_velocity,
                                            double k) {
  const mesh::Mesh& m = s.mesh;
  std::vector<CornerValues> Z(m.num_cells());
  for (Index c = 0; c < m.num_cells(); ++c) {
    const double rho = thermo.density[c];
    const double a = thermo.sound_speed[c];
    const double gamma = thermo.shock_coefficient[c];
    const Vec3 v = cell_velocity(s, c);
    for (int l = 0; l < kFacesPerCell; ++l) {
      for (int n = 0; n < kNodesPerFace; ++n) {
        const Vec3& Vp = node_velocity[m.cell_nodes[c][kHexFaceNodes[l][n]]];
        const Vec3& Sn = faces[c][l][n];
        const double S = std::sqrt(dot(Sn, Sn));
        const Vec3 dv = {Vp[0] - v[0], Vp[1] - v[1], Vp[2] - v[2]};
        const double jump = S == 0.0 ? 0.0 : std::abs(dot(dv, Sn)) / S;
        Z[c][l][n] = rho * (a + k * gamma * jump);
      }
    }
  }
  return Z;
}

// Largest nodal flux sum_c sum_f S_pf (|P_c| + Z_cfp |s_cfp|).
double nodal_flux_scale(const HydroState& s, const std::vector<CellFaceVectors>& faces,
                        const CellThermo& thermo, const std::vector<CornerValues>& Z,
                        const std::vector<Vec3>& node_velocity) {
  const mesh::Mesh& m = s.mesh;
  std::vector<double> flux(m.num_nodes(), 0.0);
  for (Index c = 0; c < m.num_cells(); ++c) {
    const Vec3 v = cell_velocity(s, c);
    for (int l = 0; l < kFacesPerCell; ++l) {
      for (int n = 0; n < kNodesPerFace; ++n) {
        const Index p = m.cell_nodes[c][kHexFaceNodes[l][n]];
        const Vec3& Vp = node_velocity[p];
        const Vec3& Sn = faces[c][l][n];
        const double S = std::sqrt(dot(Sn, Sn));
        const Vec3 dv = {Vp[0] - v[0], Vp[1] - v[1], Vp[2] - v[2]};
        flux[p] += S * std::abs(thermo.pressure[c]) + Z[c][l][n] * std::abs(dot(dv, Sn));
      }
    }
  }
  return flux.empty() ? 0.0 : *std::max_element(flux.begin(), flux.end());
}

// Component of x tangent to the (orthonormal) wall normals.
Vec3 tangential(Vec3 x, const std::vector<Vec3>& walls) {
  for (const Vec3& w : walls) {
    const double xn = dot(x, w);
    for (int d = 0; d < 3; ++d) x[d] -= xn * w[d];
  }
  return x;
}

}  // namespace

std::vector<CellFaceVectors> cell_face_vectors(const mesh::Mesh& m) {
  std::vector<CellFaceVectors> result(m.num_cells());
  for (Index c = 0; c < m.num_cells(); ++c) {
    for (int l = 0; l < kFacesPerCell; ++l) {
      const Index f = m.cell_faces[c][l];
      const double sign = m.face_owner[f] == c ? 1.0 : -1.0;
      const auto vecs = mesh::face_area_vectors(m, f);
      // The face's node order is the owner's; match by node id.
      for (int n = 0; n < kNodesPerFace; ++n) {
        const Index p = m.cell_nodes[c][kHexFaceNodes[l][n]];
        int k = 0;
        while (m.face_nodes[f][k] != p) ++k;
        for (int d = 0; d < 3; ++d) result[c][l][n][d] = sign * vecs[k][d];
      }
    }
  }
  return result;
}

std::array<Vec3, kNodesPerCell> corner_vectors(const CellFaceVectors& s) {
  std::array<Vec3, kNodesPerCell> result{};
  for (int l = 0; l < kFacesPerCell; ++l) {
    for (int n = 0; n < kNodesPerFace; ++n) {
      Vec3& r = result[kHexFaceNodes[l][n]];
      for (int d = 0; d < 3; ++d) r[d] += s[l][n][d];
    }
  }
  return result;
}

double cell_volume(const mesh::Mesh& m, Index c, const CellFaceVectors& s) {
  // The corner vectors sum to zero, so positions may be taken relative to a
  // corner, which keeps round-off proportional to the cell size.
  const auto corners = corner_vectors(s);
  const Vec3 origin = node_position(m, m.cell_nodes[c][0]);
  double volume = 0.0;
  for (int k = 1; k < kNodesPerCell; ++k) {
    const Vec3 x = node_position(m, m.cell_nodes[c][k]);
    const Vec3 r = {x[0] - origin[0], x[1] - origin[1], x[2] - origin[2]};
    volume += dot(r, corners[k]);
  }
  return volume / 3.0;
}

HydroState make_state(mesh::Mesh m, const std::vector<double>& density,
                      const std::vector<double>& vel_x,
                      const std::vector<double>& vel_y,
                      const std::vector<double>& vel_z,
                      const std::vector<double>& internal_energy) {
  const auto faces = cell_face_vectors(m);
  HydroState s;
  s.mass.resize(m.num_cells());
  s.total_energy.resize(m.num_cells());
  for (Index c = 0; c < m.num_cells(); ++c) {
    s.mass[c] = density[c] * cell_volume(m, c, faces[c]);
    s.total_energy[c] =
        internal_energy[c] +
        0.5 * (vel_x[c] * vel_x[c] + vel_y[c] * vel_y[c] + vel_z[c] * vel_z[c]);
  }
  s.mesh = std::move(m);
  s.vel_x = vel_x;
  s.vel_y = vel_y;
  s.vel_z = vel_z;
  return s;
}

CellThermo cell_thermo(const HydroState& s,
                       const std::vector<CellFaceVectors>& faces,
                       const eos::EOS& eos) {
  const Index num_cells = s.mesh.num_cells();
  CellThermo t;
  t.density.resize(num_cells);
  t.pressure.resize(num_cells);
  t.sound_speed.resize(num_cells);
  t.shock_coefficient.resize(num_cells);
  for (Index c = 0; c < num_cells; ++c) {
    const double volume = cell_volume(s.mesh, c, faces[c]);
    if (!(volume > 0.0)) {
      throw std::runtime_error("cell " + std::to_string(c) +
                               " has non-positive volume " +
                               std::to_string(volume));
    }
    const Vec3 v = cell_velocity(s, c);
    const double rho = s.mass[c] / volume;
    const double kinetic = 0.5 * dot(v, v);
    double eps = s.total_energy[c] - kinetic;
    // eps is a difference of numbers of size `kinetic`, so in a cold, moving
    // gas it comes out a few ulp negative; that is clamped, more is an error.
    if (eps < 0.0) {
      if (eps < -kEnergyRoundoff * kinetic) {
        throw std::runtime_error("cell " + std::to_string(c) +
                                 " has negative internal energy " +
                                 std::to_string(eps));
      }
      eps = 0.0;
    }
    t.density[c] = rho;
    t.pressure[c] = eos.pressure(rho, eps);
    t.sound_speed[c] = eos.sound_speed(rho, eps);
    t.shock_coefficient[c] = eos.shock_coefficient(rho, eps);
  }
  return t;
}

std::vector<CornerValues> corner_impedances(const HydroState& s,
                                            const std::vector<CellFaceVectors>& faces,
                                            const CellThermo& thermo,
                                            const std::vector<Vec3>& node_velocity) {
  return scaled_impedances(s, faces, thermo, node_velocity, 1.0);
}

NodalSystem assemble_nodal_system(const HydroState& s,
                                  const std::vector<CellFaceVectors>& faces,
                                  const CellThermo& thermo,
                                  const std::vector<CornerValues>& impedance) {
  const mesh::Mesh& m = s.mesh;
  NodalSystem sys;
  sys.M.assign(m.num_nodes(), Mat3{});
  sys.B.assign(m.num_nodes(), Vec3{});
  for (Index c = 0; c < m.num_cells(); ++c) {
    const double P = thermo.pressure[c];
    const Vec3 v = cell_velocity(s, c);
    for (int l = 0; l < kFacesPerCell; ++l) {
      for (int n = 0; n < kNodesPerFace; ++n) {
        const Index p = m.cell_nodes[c][kHexFaceNodes[l][n]];
        const double Z = impedance[c][l][n];
        // With Sn = S_pf n_pf: S_pf (n x n) = (Sn x Sn) / S_pf.
        const Vec3& Sn = faces[c][l][n];
        const double S = std::sqrt(dot(Sn, Sn));
        if (S == 0.0) continue;
        const double vn = dot(Sn, v) / S;
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) sys.M[p][i][j] += Z * Sn[i] * Sn[j] / S;
          sys.B[p][i] += (P + Z * vn) * Sn[i];
        }
      }
    }
  }
  return sys;
}

std::vector<Vec3> nodal_velocities(const HydroState& s,
                                   const std::vector<CellFaceVectors>& faces,
                                   const CellThermo& thermo,
                                   const bc::BoundarySet& bcs) {
  // Per node, Eq. 4 reads R_p(V_p) = M_p V_p - B_p
  //   = sum_c sum_f S_pf [Z_cfp s_cfp - P_c] n_pf (+ outflow terms) = 0,
  // with s_cfp = (V_p - V_c) . n_pf, up to a wall-normal reaction. Each term is
  // the gradient of the convex rho (a s^2/2 + Gamma |s|^3/3) - P s, and the
  // Jacobian is J_p = sum S_pf rho_c (a_c + 2 Gamma_c |s_cfp|) (n_pf x n_pf).
  const mesh::Mesh& m = s.mesh;
  const auto walls = bc::wall_normals(m, bcs);

  // Initial guess: the mean velocity of the cells at each node, on its walls.
  std::vector<Vec3> V(m.num_nodes(), Vec3{});
  std::vector<int> count(m.num_nodes(), 0);
  for (Index c = 0; c < m.num_cells(); ++c) {
    const Vec3 v = cell_velocity(s, c);
    for (int k = 0; k < kNodesPerCell; ++k) {
      const Index p = m.cell_nodes[c][k];
      for (int d = 0; d < 3; ++d) V[p][d] += v[d];
      ++count[p];
    }
  }
  for (Index p = 0; p < m.num_nodes(); ++p) {
    for (int d = 0; d < 3; ++d) V[p][d] /= count[p];
    V[p] = tangential(V[p], walls[p]);
  }

  double speed = 0.0;
  for (Index c = 0; c < m.num_cells(); ++c) {
    const Vec3 v = cell_velocity(s, c);
    speed = std::max(speed, std::sqrt(dot(v, v)) + thermo.sound_speed[c]);
  }

  for (int iteration = 0; iteration < kMaxNewtonIterations; ++iteration) {
    const auto Z = scaled_impedances(s, faces, thermo, V, 1.0);
    NodalSystem system = assemble_nodal_system(s, faces, thermo, Z);
    bc::apply_pressure_bcs(m, bcs, thermo.pressure, system.B);
    std::vector<Vec3> R(m.num_nodes());
    double residual = 0.0;
    for (Index p = 0; p < m.num_nodes(); ++p) {
      Vec3 r;
      for (int d = 0; d < 3; ++d) r[d] = dot(system.M[p][d], V[p]) - system.B[p][d];
      R[p] = tangential(r, walls[p]);
      residual = std::max(residual, std::sqrt(dot(R[p], R[p])));
    }
    if (residual <= kNewtonTolerance * nodal_flux_scale(s, faces, thermo, Z, V)) {
      return V;
    }

    const auto J = assemble_nodal_system(
        s, faces, thermo, scaled_impedances(s, faces, thermo, V, 2.0)).M;
    double max_step = 0.0;
    for (Index p = 0; p < m.num_nodes(); ++p) {
      Mat3 A = J[p];
      const double trace = A[0][0] + A[1][1] + A[2][2];
      if (!(trace > 0.0)) continue;
      for (int d = 0; d < 3; ++d) A[d][d] += kNewtonRegularization * trace;
      const Vec3 step = bc::solve_nodal_velocity(A, R[p], walls[p]);
      for (int d = 0; d < 3; ++d) V[p][d] -= step[d];
      max_step = std::max(max_step, std::sqrt(dot(step, step)));
    }
    if (max_step <= kStepRoundoff * speed) return V;
  }
  throw std::runtime_error("nodal solver: Newton did not converge in " +
                           std::to_string(kMaxNewtonIterations) + " iterations");
}

std::vector<std::array<Vec3, kNodesPerCell>> corner_forces(
    const HydroState& s, const std::vector<CellFaceVectors>& faces,
    const CellThermo& thermo, const std::vector<CornerValues>& impedance,
    const std::vector<Vec3>& node_velocity) {
  const mesh::Mesh& m = s.mesh;
  std::vector<std::array<Vec3, kNodesPerCell>> F(m.num_cells());
  for (Index c = 0; c < m.num_cells(); ++c) {
    const double P = thermo.pressure[c];
    const Vec3 v = cell_velocity(s, c);
    for (int l = 0; l < kFacesPerCell; ++l) {
      for (int n = 0; n < kNodesPerFace; ++n) {
        const int k = kHexFaceNodes[l][n];
        const double Z = impedance[c][l][n];
        const Vec3& Vp = node_velocity[m.cell_nodes[c][k]];
        const Vec3& Sn = faces[c][l][n];
        const double S = std::sqrt(dot(Sn, Sn));
        if (S == 0.0) continue;
        const Vec3 dv = {v[0] - Vp[0], v[1] - Vp[1], v[2] - Vp[2]};
        const double P_cfp = P + Z * dot(dv, Sn) / S;
        for (int d = 0; d < 3; ++d) F[c][k][d] += P_cfp * Sn[d];
      }
    }
  }
  return F;
}

void update_cells(HydroState& s,
                  const std::vector<std::array<Vec3, kNodesPerCell>>& forces,
                  const std::vector<Vec3>& node_velocity, double dt) {
  const mesh::Mesh& m = s.mesh;
  for (Index c = 0; c < m.num_cells(); ++c) {
    Vec3 force = {0.0, 0.0, 0.0};
    double power = 0.0;
    for (int k = 0; k < kNodesPerCell; ++k) {
      const Vec3& F = forces[c][k];
      for (int d = 0; d < 3; ++d) force[d] += F[d];
      power += dot(F, node_velocity[m.cell_nodes[c][k]]);
    }
    const double scale = dt / s.mass[c];
    s.vel_x[c] -= scale * force[0];
    s.vel_y[c] -= scale * force[1];
    s.vel_z[c] -= scale * force[2];
    s.total_energy[c] -= scale * power;
  }
}

void move_nodes(mesh::Mesh& m, const std::vector<Vec3>& node_velocity, double dt) {
  for (Index p = 0; p < m.num_nodes(); ++p) {
    m.node_x[p] += dt * node_velocity[p][0];
    m.node_y[p] += dt * node_velocity[p][1];
    m.node_z[p] += dt * node_velocity[p][2];
  }
}

std::vector<Vec3> lagrangian_step(HydroState& s, const bc::BoundarySet& bcs,
                     const eos::EOS& eos, double dt) {
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo thermo = cell_thermo(s, faces, eos);
  const auto V = nodal_velocities(s, faces, thermo, bcs);
  // Z at the converged V_p, so that the corner forces balance at each node.
  const auto Z = corner_impedances(s, faces, thermo, V);
  const auto F = corner_forces(s, faces, thermo, Z, V);
  update_cells(s, F, V, dt);
  move_nodes(s.mesh, V, dt);
  return V;
}

}  // namespace hydro::lagrangian
