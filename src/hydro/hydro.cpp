#include "hydro.hpp"

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
  t.impedance.resize(num_cells);
  for (Index c = 0; c < num_cells; ++c) {
    const double volume = cell_volume(s.mesh, c, faces[c]);
    if (!(volume > 0.0)) {
      throw std::runtime_error("cell " + std::to_string(c) +
                               " has non-positive volume " +
                               std::to_string(volume));
    }
    const Vec3 v = cell_velocity(s, c);
    const double rho = s.mass[c] / volume;
    const double eps = s.total_energy[c] - 0.5 * dot(v, v);
    t.density[c] = rho;
    t.pressure[c] = eos.pressure(rho, eps);
    t.impedance[c] = rho * eos.sound_speed(rho, eps);
  }
  return t;
}

NodalSystem assemble_nodal_system(const HydroState& s,
                                  const std::vector<CellFaceVectors>& faces,
                                  const CellThermo& thermo) {
  const mesh::Mesh& m = s.mesh;
  NodalSystem sys;
  sys.M.assign(m.num_nodes(), Mat3{});
  sys.B.assign(m.num_nodes(), Vec3{});
  for (Index c = 0; c < m.num_cells(); ++c) {
    const double P = thermo.pressure[c];
    const double Z = thermo.impedance[c];
    const Vec3 v = cell_velocity(s, c);
    for (int l = 0; l < kFacesPerCell; ++l) {
      for (int n = 0; n < kNodesPerFace; ++n) {
        const Index p = m.cell_nodes[c][kHexFaceNodes[l][n]];
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

std::vector<Vec3> nodal_velocities(const mesh::Mesh& m, const bc::BoundarySet& bcs,
                                   const CellThermo& thermo, NodalSystem system) {
  bc::apply_pressure_bcs(m, bcs, thermo.pressure, system.B);
  const auto walls = bc::wall_normals(m, bcs);
  std::vector<Vec3> V(m.num_nodes());
  for (Index p = 0; p < m.num_nodes(); ++p) {
    V[p] = bc::solve_nodal_velocity(system.M[p], system.B[p], walls[p]);
  }
  return V;
}

std::vector<std::array<Vec3, kNodesPerCell>> corner_forces(
    const HydroState& s, const std::vector<CellFaceVectors>& faces,
    const CellThermo& thermo, const std::vector<Vec3>& node_velocity) {
  const mesh::Mesh& m = s.mesh;
  std::vector<std::array<Vec3, kNodesPerCell>> F(m.num_cells());
  for (Index c = 0; c < m.num_cells(); ++c) {
    const double P = thermo.pressure[c];
    const double Z = thermo.impedance[c];
    const Vec3 v = cell_velocity(s, c);
    for (int l = 0; l < kFacesPerCell; ++l) {
      for (int n = 0; n < kNodesPerFace; ++n) {
        const int k = kHexFaceNodes[l][n];
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

void lagrangian_step(HydroState& s, const bc::BoundarySet& bcs,
                     const eos::EOS& eos, double dt) {
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo thermo = cell_thermo(s, faces, eos);
  const auto V = nodal_velocities(s.mesh, bcs, thermo,
                                  assemble_nodal_system(s, faces, thermo));
  const auto F = corner_forces(s, faces, thermo, V);
  update_cells(s, F, V, dt);
  move_nodes(s.mesh, V, dt);
}

}  // namespace hydro::lagrangian
