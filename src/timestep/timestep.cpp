#include "timestep.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "sync.hpp"

namespace hydro::timestep {

using mesh::Index;
using mesh::kFacesPerCell;
using mesh::kNodesPerCell;
using mesh::kNodesPerFace;

namespace {

// Relative round-off below which dV_c/dt counts as zero.
constexpr double kRateRoundoff = 1e-12;

double dot(const Vec3& a, const Vec3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

}  // namespace

double characteristic_length(const lagrangian::CellFaceVectors& faces, double volume) {
  double largest = 0.0;
  for (int l = 0; l < kFacesPerCell; ++l) {
    double area = 0.0;
    for (int n = 0; n < kNodesPerFace; ++n) area += std::sqrt(dot(faces[l][n], faces[l][n]));
    largest = std::max(largest, area);
  }
  return volume / largest;
}

double stable_dt(const lagrangian::HydroState& s, const eos::EOS& eos,
                 const std::vector<Vec3>& node_velocity, double dt_prev,
                 const TimestepParams& params) {
  const mesh::Mesh& m = s.mesh;
  const auto faces = lagrangian::cell_face_vectors(m);
  const auto thermo = lagrangian::cell_thermo(s, faces, eos);
  double dt = std::min(params.max_growth * dt_prev, params.dt_max);
  for (Index c = 0; c < m.num_cells(); ++c) {
    const double volume = lagrangian::cell_volume(m, c, faces[c]);
    const double a = thermo.sound_speed[c];
    if (a > 0.0) {
      dt = std::min(dt, params.cfl * characteristic_length(faces[c], volume) / a);
    }
    if (node_velocity.empty()) continue;
    const auto n = lagrangian::corner_vectors(faces[c]);
    double rate = 0.0, magnitude = 0.0;
    for (int k = 0; k < kNodesPerCell; ++k) {
      const double term = dot(n[k], node_velocity[m.cell_nodes[c][k]]);
      rate += term;
      magnitude += std::abs(term);
    }
    // A rigid translation gives dV/dt = 0 only up to round-off of its terms.
    if (std::abs(rate) > kRateRoundoff * magnitude) {
      dt = std::min(dt, params.volume_fraction * volume / std::abs(rate));
    }
  }
  dt = sync::global_reduce(dt, sync::ReduceOp::Min);
  if (!std::isfinite(dt)) {
    throw std::runtime_error("stable_dt: no criterion limits dt (set dt_max or dt_prev)");
  }
  return dt;
}

}  // namespace hydro::timestep
