#include "driver.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
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

// Deposits the deck's energy spike (io::EnergySpike) into s.
void deposit_energy_spike(const io::Deck& deck, HydroState& s) {
  const io::EnergySpike& spike = *deck.energy_spike;
  const mesh::Mesh& m = s.mesh;
  const double size = std::hypot(deck.mesh.x_max - deck.mesh.x_min,
                                 deck.mesh.y_max - deck.mesh.y_min,
                                 deck.mesh.z_max - deck.mesh.z_min);
  const double tol = 1e-12 * size;
  std::vector<Index> cells;
  double mass = 0.0;
  for (Index c = 0; c < m.num_cells(); ++c) {
    const bool touches = std::any_of(m.cell_nodes[c].begin(), m.cell_nodes[c].end(),
                                     [&](Index p) {
                                       return std::hypot(m.node_x[p] - spike.at[0],
                                                         m.node_y[p] - spike.at[1],
                                                         m.node_z[p] - spike.at[2]) <= tol;
                                     });
    if (!touches) continue;
    cells.push_back(c);
    mass += s.mass[c];
  }
  if (cells.empty()) {
    std::ostringstream msg;
    msg << "initial_state: no cell has a node at the energy spike's point (" << spike.at[0]
        << ", " << spike.at[1] << ", " << spike.at[2] << ")";
    throw std::runtime_error(msg.str());
  }
  const double eps = spike.energy_per_length * (deck.mesh.z_max - deck.mesh.z_min) / mass;
  for (const Index c : cells) {
    s.total_energy[c] = eps + 0.5 * (s.vel_x[c] * s.vel_x[c] + s.vel_y[c] * s.vel_y[c] +
                                     s.vel_z[c] * s.vel_z[c]);
  }
}

constexpr double kPi = 3.14159265358979323846;

void check_probe_ids(const std::vector<Index>& ids, Index count, const char* what) {
  for (const Index id : ids) {
    if (id < 0 || id >= count) {
      throw std::invalid_argument(std::string("probe: ") + what + " " + std::to_string(id) +
                                  " is not in [0, " + std::to_string(count) + ")");
    }
  }
}

// The probe lines of `probes` at this step (see Probes in driver.hpp).
void write_probes(const Probes& probes, const HydroState& s, const eos::EOS& eos,
                  const std::vector<Vec3>& node_velocity, long step, double t, double dt) {
  const mesh::Mesh& m = s.mesh;
  if (!probes.cells.empty()) {
    const auto faces = lagrangian::cell_face_vectors(m);
    const auto thermo = lagrangian::cell_thermo(s, faces, eos);
    for (const Index c : probes.cells) {
      Vec3 centroid = {0.0, 0.0, 0.0};
      for (const Index p : m.cell_nodes[c]) {
        centroid[0] += m.node_x[p] / kNodesPerCell;
        centroid[1] += m.node_y[p] / kNodesPerCell;
        centroid[2] += m.node_z[p] / kNodesPerCell;
      }
      const double eps =
          s.total_energy[c] - 0.5 * (s.vel_x[c] * s.vel_x[c] + s.vel_y[c] * s.vel_y[c] +
                                     s.vel_z[c] * s.vel_z[c]);
      const nlohmann::ordered_json line = {
          {"cell", c},
          {"step", step},
          {"t", t},
          {"dt", dt},
          {"centroid", centroid},
          {"volume", lagrangian::cell_volume(m, c, faces[c])},
          {"mass", s.mass[c]},
          {"density", thermo.density[c]},
          {"velocity", {s.vel_x[c], s.vel_y[c], s.vel_z[c]}},
          {"pressure", thermo.pressure[c]},
          {"internal_energy", eps},
          {"sound_speed", thermo.sound_speed[c]},
      };
      *probes.out << line.dump() << '\n';
    }
  }
  for (const Index p : probes.nodes) {
    nlohmann::ordered_json line = {
        {"node", p},
        {"step", step},
        {"t", t},
        {"dt", dt},
        {"position", {m.node_x[p], m.node_y[p], m.node_z[p]}},
        {"velocity", nullptr},
    };
    if (!node_velocity.empty()) line["velocity"] = node_velocity[p];
    *probes.out << line.dump() << '\n';
  }
}

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
  HydroState s = lagrangian::make_state(std::move(m), rho, vx, vy, vz, eps);
  if (deck.energy_spike) deposit_energy_spike(deck, s);
  return s;
}

RunResult run(const io::Deck& deck, std::ostream& log, const Probes& probes) {
  const eos::IdealGasEOS eos(deck.gamma);
  const lagrangian::EnergySource source = energy_source(deck);
  RunResult r{initial_state(deck), 0.0, 0};
  check_probe_ids(probes.cells, r.state.mesh.num_cells(), "cell");
  check_probe_ids(probes.nodes, r.state.mesh.num_nodes(), "node");
  const bool probing =
      probes.out != nullptr && !(probes.cells.empty() && probes.nodes.empty());
  std::vector<Vec3> node_velocity;  // of the last step; none before the first
  double dt_prev = deck.dt_initial;
  double dt = 0.0;
  if (probing) write_probes(probes, r.state, eos, node_velocity, 0, r.t, 0.0);
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
      node_velocity = lagrangian::scheme_step(r.state, deck.boundaries, eos, r.t, dt,
                                               deck.scheme, source);
      r.t = last ? deck.t_end : r.t + dt;
    } catch (const std::runtime_error& e) {
      std::ostringstream msg;
      msg << "run: step " << r.steps + 1 << " from t = " << r.t << " failed: " << e.what();
      throw std::runtime_error(msg.str());
    }
    ++r.steps;
    if (probing) write_probes(probes, r.state, eos, node_velocity, r.steps, r.t, dt);
    if (deck.log_every > 0 && r.steps % deck.log_every == 0) {
      log << "step " << r.steps << "  t = " << r.t << "  dt = " << dt << '\n';
    }
  }
  log << "done: " << r.steps << " steps, t = " << r.t << ", last dt = " << dt << '\n';
  return r;
}

}  // namespace hydro::app
