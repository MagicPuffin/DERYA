#pragma once

#include <array>
#include <vector>

#include "mesh.hpp"

namespace hydro::bc {

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<Vec3, 3>;

// Boundary conditions of the EUCCLHYD nodal solver (Colaïtis, Guisset & Breil,
// SSRN 5167115, Sec. 2.2, after Maire et al. [30] and Georges, Breil & Maire
// [29]). They act on each boundary node's system M_p V_p = B (Eq. 4), assumed
// assembled over every face of every cell in C(p), boundary faces included.
//   - Symmetry: a rigid wall, V_p . n = 0 for the unit normal n of every
//     symmetry face at p.
//   - Outflow: a prescribed pressure P* = P_c of the face's owner cell, so the
//     boundary is transparent to a uniform state (zero-gradient outflow).
enum class BoundaryType { Symmetry, Outflow };

// One boundary type per domain side, indexed by mesh::Mesh::face_boundary
// (-x, +x, -y, +y, -z, +z).
using BoundarySet = std::array<BoundaryType, mesh::kFacesPerCell>;

// Adds the prescribed-pressure force of every outflow face f to its nodes:
// B_p -= P*_f S_pf n_pf, with S_pf n_pf from mesh::face_area_vectors (outward,
// since boundary faces are owned by the interior cell) and P*_f =
// cell_pressure[face_owner[f]]. B holds one entry per mesh node.
void apply_pressure_bcs(const mesh::Mesh& m, const BoundarySet& bcs,
                        const std::vector<double>& cell_pressure,
                        std::vector<Vec3>& B);

// Per node, an orthonormal basis of the normals of the symmetry faces at that
// node: empty off symmetry walls, one vector on a wall, two on an edge where
// walls meet, three at a corner. Coplanar faces count once. A face's normal is
// its total area vector, normalized.
std::vector<std::vector<Vec3>> wall_normals(const mesh::Mesh& m,
                                            const BoundarySet& bcs);

// Nodal velocity from M_p V_p = B subject to V_p . n = 0 for every n in walls
// (orthonormal, as returned by wall_normals): the Lagrange-multiplier solution,
// computed as V_p = T (T^T M_p T)^{-1} T^T B over an orthonormal basis T of the
// wall tangent space. M_p must be symmetric positive definite. With no walls
// this is the plain solve of Eq. 4; with three it is V_p = 0.
Vec3 solve_nodal_velocity(const Mat3& M, const Vec3& B,
                          const std::vector<Vec3>& walls);

}  // namespace hydro::bc
