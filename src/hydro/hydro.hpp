#pragma once

#include <array>
#include <vector>

#include "bc.hpp"
#include "eos.hpp"
#include "mesh.hpp"

namespace hydro::lagrangian {

using mesh::Vec3;
using bc::Mat3;

// First-order cell-centered Lagrangian step (EUCCLHYD; Colaïtis, Guisset &
// Breil, SSRN 5167115, Sec. 2.2, after Georges, Breil & Maire [29]).

// Outward nodal face area vectors S_pf n_pf (Eq. 5) of one cell, indexed
// [local face l][face node n] as in mesh::kHexFaceNodes: entry [l][n] belongs
// to local corner kHexFaceNodes[l][n], node cell_nodes[c][kHexFaceNodes[l][n]].
using CellFaceVectors =
    std::array<std::array<Vec3, mesh::kNodesPerFace>, mesh::kFacesPerCell>;

// CellFaceVectors of every cell, from the current node positions.
std::vector<CellFaceVectors> cell_face_vectors(const mesh::Mesh& m);

// Corner area vectors n_cp = sum_{f in F(c,p)} S_pf n_pf, per local corner
// (cell_nodes order). By the GCL they are dV_c/dx_p, and they sum to zero.
std::array<Vec3, mesh::kNodesPerCell> corner_vectors(const CellFaceVectors& s);

// Volume of cell c, V_c = 1/3 sum_p x_p . n_cp (exact for the face-split cell,
// V_c being cubic in the node positions with gradient n_cp).
double cell_volume(const mesh::Mesh& m, mesh::Index c, const CellFaceVectors& s);

// Hydro state: the (moving) mesh plus mass-averaged cell quantities. m_c is
// constant in the Lagrangian step; E_c is the specific total energy.
struct HydroState {
  mesh::Mesh mesh;
  std::vector<double> mass;
  std::vector<double> vel_x;
  std::vector<double> vel_y;
  std::vector<double> vel_z;
  std::vector<double> total_energy;
};

// State with the given per-cell density, velocity and specific internal
// energy: m_c = rho_c V_c, E_c = eps_c + |V_c|^2 / 2.
HydroState make_state(mesh::Mesh m, const std::vector<double>& density,
                      const std::vector<double>& vel_x,
                      const std::vector<double>& vel_y,
                      const std::vector<double>& vel_z,
                      const std::vector<double>& internal_energy);

// Cell thermodynamics at the current geometry: rho_c = m_c / V_c, and P_c,
// a_c and the shock coefficient Gamma_c from the EOS, at eps_c = E_c -
// |V_c|^2/2. Throws std::runtime_error on a cell with non-positive volume (a
// tangled mesh) or a negative eps_c beyond round-off, which is clamped to 0.
struct CellThermo {
  std::vector<double> density;
  std::vector<double> pressure;
  std::vector<double> sound_speed;
  std::vector<double> shock_coefficient;
};
CellThermo cell_thermo(const HydroState& s,
                       const std::vector<CellFaceVectors>& faces,
                       const eos::EOS& eos);

// One value per cell corner of each face, indexed like CellFaceVectors.
using CornerValues =
    std::array<std::array<double, mesh::kNodesPerFace>, mesh::kFacesPerCell>;

// Impedance of Eq. 3 per cell, face and node, in its two-shock (Dukowicz)
// form rather than the paper's acoustic Z_c = rho_c a_c (see decisions.md):
// Z_cfp = rho_c (a_c + Gamma_c |(V_p - V_c) . n_pf|). It stays positive in a
// cold gas wherever the velocity jumps.
std::vector<CornerValues> corner_impedances(const HydroState& s,
                                            const std::vector<CellFaceVectors>& faces,
                                            const CellThermo& thermo,
                                            const std::vector<Vec3>& node_velocity);

// Eq. 4, per node: M_p = sum_c sum_f S_pf Z_cfp (n_pf x n_pf) and
// B_p = sum_c sum_f S_pf [P_c n_pf + Z_cfp (n_pf x n_pf) V_c], over every face
// of every cell at p, boundary faces included (boundary conditions not yet
// applied).
struct NodalSystem {
  std::vector<Mat3> M;
  std::vector<Vec3> B;
};
NodalSystem assemble_nodal_system(const HydroState& s,
                                  const std::vector<CellFaceVectors>& faces,
                                  const CellThermo& thermo,
                                  const std::vector<CornerValues>& impedance);

// Nodal velocities V_p solving Eq. 4 with the two-shock Z_cfp, subject to the
// boundary conditions (outflow pressures in B, symmetry and piston walls on
// V_p). Z_cfp depends on V_p, so each node's system is nonlinear; it is the
// gradient of a convex function of V_p, solved by Newton's method from the
// mean of the surrounding cell velocities (moved onto the walls). Where M_p is
// singular (cold cells all moving alike) the forces do not depend on V_p and
// it keeps that initial value. Throws
// std::runtime_error if Newton does not converge.
std::vector<Vec3> nodal_velocities(const HydroState& s,
                                   const std::vector<CellFaceVectors>& faces,
                                   const CellThermo& thermo,
                                   const bc::BoundarySet& bcs);

// Corner forces F_cp = sum_{f in F(c,p)} S_pf P_cfp n_pf, per cell and local
// corner, with the nodal pressures of Eq. 3,
// P_cfp = P_c + Z_cfp (V_c - V_p) . n_pf.
std::vector<std::array<Vec3, mesh::kNodesPerCell>> corner_forces(
    const HydroState& s, const std::vector<CellFaceVectors>& faces,
    const CellThermo& thermo, const std::vector<CornerValues>& impedance,
    const std::vector<Vec3>& node_velocity);

// Eq. 1, forward Euler: m_c dV_c = -dt sum_p F_cp and
// m_c dE_c = -dt sum_p F_cp . V_p.
void update_cells(HydroState& s,
                  const std::vector<std::array<Vec3, mesh::kNodesPerCell>>& forces,
                  const std::vector<Vec3>& node_velocity, double dt);

// Eq. 2, forward Euler: x_p += dt V_p.
void move_nodes(mesh::Mesh& m, const std::vector<Vec3>& node_velocity, double dt);

// One first-order Lagrangian step of size dt: geometry, thermodynamics and
// nodal solve at t^n, then Eq. 1 and Eq. 2. Returns the nodal velocities V_p
// used (for the timestep controller's volume criterion).
std::vector<Vec3> lagrangian_step(HydroState& s, const bc::BoundarySet& bcs,
                     const eos::EOS& eos, double dt);

}  // namespace hydro::lagrangian
