#include "driver.hpp"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "eos.hpp"
#include "mesh.hpp"
#include "scheme.hpp"
#include "timestep.hpp"

namespace hydro::app {

using lagrangian::HydroState;
using lagrangian::Vec3;
using mesh::Index;
using mesh::kNodesPerCell;

namespace {

bool inside(const io::Region& r, const Vec3& x) {
  for (int d = 0; d < 3; ++d) {
    if (x[d] < r.min[d] || x[d] >= r.max[d]) return false;
  }
  return true;
}

constexpr double kPi = 3.14159265358979323846;

}  // namespace

Vec3 taylor_green_velocity(const io::TaylorGreen& tg, const Vec3& x) {
  return {tg.c1 * std::sin(kPi * x[0]) * std::cos(kPi * x[1]),
          -tg.c1 * std::cos(kPi * x[0]) * std::sin(kPi * x[1]), 0.0};
}

double taylor_green_pressure(const io::TaylorGreen& tg, const Vec3& x) {
  return 0.25 * tg.density * tg.c1 * tg.c1 *
             (std::cos(2.0 * kPi * x[0]) + std::cos(2.0 * kPi * x[1])) +
         tg.c2;
}

double taylor_green_source(const io::TaylorGreen& tg, double gamma, const Vec3& x) {
  return kPi * tg.density * tg.c1 * tg.c1 * tg.c1 / (4.0 * (gamma - 1.0)) *
         (std::cos(3.0 * kPi * x[0]) * std::cos(kPi * x[1]) -
          std::cos(3.0 * kPi * x[1]) * std::cos(kPi * x[0]));
}

lagrangian::EnergySource energy_source(const io::Deck& deck) {
  switch (deck.source) {
    case io::EnergySourceType::TaylorGreen: {
      const io::TaylorGreen tg = *deck.taylor_green;
      const double gamma = deck.gamma;
      return [tg, gamma](const Vec3& x) { return taylor_green_source(tg, gamma, x); };
    }
    case io::EnergySourceType::None:
      break;
  }
  return {};
}

HydroState initial_state(const io::Deck& deck) {
  mesh::Mesh m = mesh::generate_structured_mesh(deck.mesh);
  if (deck.perturbation == io::MeshPerturbation::Saltzman) {
    mesh::apply_saltzman_skew(m, deck.mesh);
  }
  const auto n = static_cast<std::size_t>(m.num_cells());
  std::vector<double> rho(n), vx(n), vy(n), vz(n), eps(n);
  for (Index c = 0; c < m.num_cells(); ++c) {
    Vec3 centroid = {0.0, 0.0, 0.0};
    for (const Index p : m.cell_nodes[c]) {
      centroid[0] += m.node_x[p] / kNodesPerCell;
      centroid[1] += m.node_y[p] / kNodesPerCell;
      centroid[2] += m.node_z[p] / kNodesPerCell;
    }
    io::RegionState st = deck.background;
    if (deck.taylor_green) {
      st.density = deck.taylor_green->density;
      st.velocity = taylor_green_velocity(*deck.taylor_green, centroid);
      st.pressure = taylor_green_pressure(*deck.taylor_green, centroid);
    }
    for (const io::Region& r : deck.regions) {
      if (inside(r, centroid)) st = r.state;
    }
    rho[c] = st.density;
    vx[c] = st.velocity[0];
    vy[c] = st.velocity[1];
    vz[c] = st.velocity[2];
    eps[c] = st.pressure / ((deck.gamma - 1.0) * st.density);
  }
  return lagrangian::make_state(std::move(m), rho, vx, vy, vz, eps);
}

RunResult run(const io::Deck& deck, std::ostream& log) {
  const eos::IdealGasEOS eos(deck.gamma);
  const lagrangian::EnergySource source = energy_source(deck);
  RunResult r{initial_state(deck), 0.0, 0};
  std::vector<Vec3> node_velocity;  // of the last step; none before the first
  double dt_prev = deck.dt_initial;
  double dt = 0.0;
  while (r.t < deck.t_end) {
    if (r.steps >= deck.max_steps) {
      std::ostringstream msg;
      msg << "run: max_steps = " << deck.max_steps << " reached at t = " << r.t
          << " < t_end = " << deck.t_end;
      throw std::runtime_error(msg.str());
    }
    try {
      dt = timestep::stable_dt(r.state, eos, node_velocity, dt_prev, deck.timestep);
      // The growth limit applies to the stable dt, not to a final step cut short.
      dt_prev = dt;
      const bool last = dt >= deck.t_end - r.t;
      if (last) dt = deck.t_end - r.t;
      node_velocity = lagrangian::scheme_step(r.state, deck.boundaries, eos, dt,
                                               deck.scheme, source);
      r.t = last ? deck.t_end : r.t + dt;
    } catch (const std::runtime_error& e) {
      std::ostringstream msg;
      msg << "run: step " << r.steps + 1 << " from t = " << r.t << " failed: " << e.what();
      throw std::runtime_error(msg.str());
    }
    ++r.steps;
    if (deck.log_every > 0 && r.steps % deck.log_every == 0) {
      log << "step " << r.steps << "  t = " << r.t << "  dt = " << dt << '\n';
    }
  }
  log << "done: " << r.steps << " steps, t = " << r.t << ", last dt = " << dt << '\n';
  return r;
}

}  // namespace hydro::app
