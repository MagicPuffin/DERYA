#pragma once

#include <array>
#include <functional>
#include <vector>

#include "mesh.hpp"

namespace hydro::bc {

using Vec3 = mesh::Vec3;
using Mat3 = std::array<Vec3, 3>;

// Boundary conditions of the EUCCLHYD nodal solver (Colaïtis, Guisset & Breil,
// SSRN 5167115, Sec. 2.2, after Maire et al. [30] and Georges, Breil & Maire
// [29]). They act on each boundary node's system M_p V_p = B (Eq. 4), assumed
// assembled over every face of every cell in C(p), boundary faces included.
//   - Symmetry: a rigid wall, V_p . n = 0 for the unit normal n of every
//     symmetry face at p.
//   - Outflow: a prescribed pressure P* = P_c of the face's owner cell, so the
//     boundary is transparent to a uniform state (zero-gradient outflow).
//   - Piston: a wall moving with a prescribed velocity u_w, V_p . n = u_w . n
//     (the tangential motion stays free, as on a symmetry wall).
//   - Pressure: a prescribed pressure P*(t), given as a constant or a law of
//     time (Kidder's shell, thesis Eqs. 3.13-3.14), where Outflow has P_c.
enum class BoundaryType { Symmetry, Outflow, Piston, Pressure };

// Prescribed boundary pressure as a function of time.
using PressureLaw = std::function<double(double t)>;

// Boundary condition of one domain side. velocity is the wall velocity u_w
// of a Piston side; pressure (a constant) or pressure_law (if set) is the P*
// of a Pressure side; each is ignored on the other types. Implicitly
// constructible from a BoundaryType, so {Symmetry, Outflow, ...} initializes
// a BoundarySet.
struct Boundary {
  BoundaryType type = BoundaryType::Symmetry;
  Vec3 velocity = {0.0, 0.0, 0.0};
  double pressure = 0.0;
  PressureLaw pressure_law;

  Boundary() = default;
  Boundary(BoundaryType t, const Vec3& v = {0.0, 0.0, 0.0}) : type(t), velocity(v) {}
};

// A Pressure boundary with P* = p at all times, or P* = law(t).
Boundary pressure_boundary(double p);
Boundary pressure_boundary(PressureLaw law);

// Kidder's boundary pressure (Georges' thesis Eqs. 3.13-3.14):
// P*(t) = p0 h(t)^(-2 gamma / (gamma - 1)), h(t) = sqrt(1 - (t / t_focus)^2),
// p0 = P0(r_b) the initial pressure at the boundary radius. The law throws
// std::runtime_error at t >= t_focus (h = 0, the shell has collapsed).
PressureLaw kidder_pressure(double p0, double t_focus, double gamma);

// One boundary per domain side, indexed by mesh::Mesh::face_boundary
// (-x, +x, -y, +y, -z, +z).
using BoundarySet = std::array<Boundary, mesh::kFacesPerCell>;

// The boundaries at time t: every Pressure side with a law gets pressure =
// law(t) and no law; the other sides are copied. A step stage calls this
// once, at its own time, and hands the result to the nodal solver.
BoundarySet at_time(const BoundarySet& bcs, double t);

// Adds the prescribed-pressure force of every outflow or pressure face f to
// its nodes: B_p -= P*_f S_pf n_pf, with S_pf n_pf from
// mesh::face_area_vectors (outward, since boundary faces are owned by the
// interior cell), P*_f = cell_pressure[face_owner[f]] on an outflow face and
// the side's pressure on a pressure face. B holds one entry per mesh node.
// Throws std::logic_error on a pressure side whose law was not evaluated
// (at_time) first.
void apply_pressure_bcs(const mesh::Mesh& m, const BoundarySet& bcs,
                        const std::vector<double>& cell_pressure,
                        std::vector<Vec3>& B);

// Per node, an orthonormal basis of the normals of the wall faces (symmetry
// or piston) at that node: empty off walls, one vector on a wall, two on an
// edge where walls meet, three at a corner. Coplanar faces count once. A
// face's normal is its total area vector, normalized; a face of zero area
// (a polar mesh's axis side) is skipped.
std::vector<std::vector<Vec3>> wall_normals(const mesh::Mesh& m,
                                            const BoundarySet& bcs);

// The wall constraints of every node: normals as from wall_normals, and
// velocity, the part of V_p they prescribe, sum_k (V_p . e_k) e_k over the
// basis e_k (zero unless a piston moves along its normal). A node velocity
// V_p satisfies its walls exactly when its normal part equals velocity, so
// V_p = velocity + (any vector tangent to the walls).
struct WallConstraints {
  std::vector<std::vector<Vec3>> normals;
  std::vector<Vec3> velocity;
};
WallConstraints wall_constraints(const mesh::Mesh& m, const BoundarySet& bcs);

// Nodal velocity from M_p V_p = B subject to V_p . n = 0 for every n in walls
// (orthonormal, as returned by wall_normals): the Lagrange-multiplier solution,
// also the tangential correction (Newton step) on a moving wall,
// computed as V_p = T (T^T M_p T)^{-1} T^T B over an orthonormal basis T of the
// wall tangent space. M_p must be symmetric positive definite. With no walls
// this is the plain solve of Eq. 4; with three it is V_p = 0.
Vec3 solve_nodal_velocity(const Mat3& M, const Vec3& B,
                          const std::vector<Vec3>& walls);

}  // namespace hydro::bc
