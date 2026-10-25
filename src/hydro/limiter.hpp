#pragma once

#include <vector>

#include "hydro.hpp"

namespace hydro::lagrangian {

// SP-V slope limiter of the second-order extension (Colaïtis, Guisset &
// Breil, SSRN 5167115, Sec. 2.3, Eqs. 7-8). The limited gradients are
// chi_c (grad psi)_c, one factor chi_c in [0, 1] per cell and variable. See
// decisions.md, "Second-order extension".

// Barth-Jespersen function BJ({a_i}; {b_j}; d) of Eq. 8: the min over j of
//   gamma_j = min((max_i a_i - d) / (b_j - d), 1)  if b_j > d,
//   gamma_j = min((min_i a_i - d) / (b_j - d), 1)  if b_j < d,
//   gamma_j = 1                                     if b_j = d,
// the largest factor keeping every extrapolated value b_j within the range of
// the control values a_i. In [0, 1] when d lies in that range (the SP-V
// controls include the cell's own value); 1 for no extrapolated values. The
// preprint prints max(., 1) in the b_j < d branch, which would never limit;
// min is the standard Barth-Jespersen form and what is meant.
double barth_jespersen(const std::vector<double>& controls,
                       const std::vector<double>& extrapolated, double reference);

// Node neighborhood C(P(c)) of each cell: the cells sharing at least one node
// with c, c included, in increasing order. Boundary cells just have fewer
// (no ghost cells).
std::vector<std::vector<mesh::Index>> node_neighborhoods(const mesh::Mesh& m);

// chi_c(P) of Eq. 7: Barth-Jespersen with the P_c' of C(P(c)) as controls,
// the extrapolations P_c + (grad P)_c . (x_p - x_c) at c's own nodes, and
// reference P_c. x_c = cell_centers.
std::vector<double> pressure_limiters(const mesh::Mesh& m,
                                      const std::vector<double>& pressure,
                                      const std::vector<Vec3>& pressure_gradient);

// chi_c(V) of Eq. 7: the same on velocities projected onto u_c = V_c / |V_c|,
// i.e. controls V_c' . u_c, extrapolations (V_c + (grad V)_c (x_p - x_c)) .
// u_c, reference V_c . u_c = |V_c|. Where V_c = 0 exactly, u_c is the
// pressure-gradient direction, or e_x if that is zero too (our choice, the
// thesis SP rule). velocity_gradient[c][i][j] = dV_i/dx_j, as cell_gradients.
std::vector<double> velocity_limiters(const mesh::Mesh& m,
                                      const std::vector<double>& vel_x,
                                      const std::vector<double>& vel_y,
                                      const std::vector<double>& vel_z,
                                      const std::vector<Mat3>& velocity_gradient,
                                      const std::vector<Vec3>& pressure_gradient);

}  // namespace hydro::lagrangian
