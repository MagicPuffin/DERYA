#include "limiter.hpp"

#include <algorithm>
#include <cmath>

#include "gradients.hpp"

namespace hydro::lagrangian {

using mesh::Index;
using mesh::kNodesPerCell;

double barth_jespersen(const std::vector<double>& controls,
                       const std::vector<double>& extrapolated, double reference) {
  if (extrapolated.empty()) return 1.0;
  const auto [lo, hi] = std::minmax_element(controls.begin(), controls.end());
  double chi = 1.0;
  for (const double b : extrapolated) {
    const double db = b - reference;
    double gamma = 1.0;
    if (db > 0.0) {
      gamma = std::min((*hi - reference) / db, 1.0);
    } else if (db < 0.0) {
      gamma = std::min((*lo - reference) / db, 1.0);
    }
    chi = std::min(chi, gamma);
  }
  return chi;
}

std::vector<std::vector<Index>> node_neighborhoods(const mesh::Mesh& m) {
  std::vector<std::vector<Index>> node_cells(m.num_nodes());
  for (Index c = 0; c < m.num_cells(); ++c) {
    for (const Index p : m.cell_nodes[c]) node_cells[p].push_back(c);
  }
  std::vector<std::vector<Index>> nbr(m.num_cells());
  for (Index c = 0; c < m.num_cells(); ++c) {
    for (const Index p : m.cell_nodes[c]) {
      nbr[c].insert(nbr[c].end(), node_cells[p].begin(), node_cells[p].end());
    }
    std::sort(nbr[c].begin(), nbr[c].end());
    nbr[c].erase(std::unique(nbr[c].begin(), nbr[c].end()), nbr[c].end());
  }
  return nbr;
}

namespace {

// Offsets x_p - x_c of cell c's nodes from its reconstruction point.
std::array<Vec3, kNodesPerCell> node_offsets(const mesh::Mesh& m, Index c,
                                             const Vec3& xc) {
  std::array<Vec3, kNodesPerCell> r;
  for (int k = 0; k < kNodesPerCell; ++k) {
    const Index p = m.cell_nodes[c][k];
    r[k] = {m.node_x[p] - xc[0], m.node_y[p] - xc[1], m.node_z[p] - xc[2]};
  }
  return r;
}

double dot(const Vec3& a, const Vec3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// |V_c| at or below this fraction of the largest |V_c'| over C(P(c)) counts
// as V_c = 0 in velocity_limiters (decisions.md, "Limiter: SP-V"). Such a
// V_c is whatever reached the cell first, round-off (Sod, 1e-17) or the
// precursor of an arriving wave (Sedov, 2e-7, tangential), not the flow
// whose gradient is to be limited. Task 13 had 1e-8, enough for round-off
// only; Sedov at order 2 then fails, as at 1e-6. On Sedov at 30x30 every
// value from 1e-5 to 3e-1 runs, with the same L1 to 3e-4 relative; 1e-2 is
// inside that band, not at its edge.
constexpr double kNegligibleSpeed = 1e-2;

// a / |a|, or the zero vector for a = 0.
Vec3 unit(const Vec3& a) {
  const double n = std::sqrt(dot(a, a));
  if (n == 0.0) return Vec3{};
  return {a[0] / n, a[1] / n, a[2] / n};
}

}  // namespace

std::vector<double> pressure_limiters(const mesh::Mesh& m,
                                      const std::vector<double>& pressure,
                                      const std::vector<Vec3>& pressure_gradient) {
  const auto xc = cell_centers(m);
  const auto nbr = node_neighborhoods(m);
  std::vector<double> chi(m.num_cells());
  std::vector<double> controls;
  std::vector<double> extrapolated(kNodesPerCell);
  for (Index c = 0; c < m.num_cells(); ++c) {
    controls.clear();
    for (const Index d : nbr[c]) controls.push_back(pressure[d]);
    const auto r = node_offsets(m, c, xc[c]);
    for (int k = 0; k < kNodesPerCell; ++k) {
      extrapolated[k] = pressure[c] + dot(pressure_gradient[c], r[k]);
    }
    chi[c] = barth_jespersen(controls, extrapolated, pressure[c]);
  }
  return chi;
}

std::vector<double> velocity_limiters(const mesh::Mesh& m,
                                      const std::vector<double>& vel_x,
                                      const std::vector<double>& vel_y,
                                      const std::vector<double>& vel_z,
                                      const std::vector<Mat3>& velocity_gradient,
                                      const std::vector<Vec3>& pressure_gradient) {
  const auto xc = cell_centers(m);
  const auto nbr = node_neighborhoods(m);
  std::vector<double> chi(m.num_cells());
  std::vector<double> controls;
  std::vector<double> extrapolated(kNodesPerCell);
  for (Index c = 0; c < m.num_cells(); ++c) {
    const Vec3 Vc = {vel_x[c], vel_y[c], vel_z[c]};
    // V_c small against the surrounding speeds has no meaningful direction:
    // mirror cells of a 1D flow would project on different round-off
    // vectors, and a tangential precursor hides the incoming gradient.
    double speed = 0.0;
    for (const Index d : nbr[c]) {
      speed = std::max(speed, std::sqrt(vel_x[d] * vel_x[d] + vel_y[d] * vel_y[d] +
                                        vel_z[d] * vel_z[d]));
    }
    Vec3 u{};
    if (std::sqrt(dot(Vc, Vc)) > kNegligibleSpeed * speed) u = unit(Vc);
    if (u == Vec3{}) u = unit(pressure_gradient[c]);
    if (u == Vec3{}) u = {1.0, 0.0, 0.0};

    controls.clear();
    for (const Index d : nbr[c]) {
      controls.push_back(vel_x[d] * u[0] + vel_y[d] * u[1] + vel_z[d] * u[2]);
    }
    // (V_c + G r) . u = V_c . u + (G^T u) . r.
    const Mat3& G = velocity_gradient[c];
    Vec3 Gtu{};
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) Gtu[j] += G[i][j] * u[i];
    }
    const double reference = dot(Vc, u);
    const auto r = node_offsets(m, c, xc[c]);
    for (int k = 0; k < kNodesPerCell; ++k) {
      extrapolated[k] = reference + dot(Gtu, r[k]);
    }
    chi[c] = barth_jespersen(controls, extrapolated, reference);
  }
  return chi;
}

}  // namespace hydro::lagrangian
