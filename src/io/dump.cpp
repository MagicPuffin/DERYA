#include "dump.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <stdexcept>
#include <vector>

namespace hydro::io {

using mesh::Index;
using mesh::kNodesPerCell;

void write_dump(const std::string& path, const lagrangian::HydroState& s,
                const eos::EOS& eos, double t, long steps) {
  const mesh::Mesh& m = s.mesh;
  const auto faces = lagrangian::cell_face_vectors(m);
  const auto thermo = lagrangian::cell_thermo(s, faces, eos);
  const auto n = static_cast<std::size_t>(m.num_cells());
  std::vector<double> cx(n), cy(n), cz(n), volume(n), eps(n);
  for (Index c = 0; c < m.num_cells(); ++c) {
    for (const Index p : m.cell_nodes[c]) {
      cx[c] += m.node_x[p] / kNodesPerCell;
      cy[c] += m.node_y[p] / kNodesPerCell;
      cz[c] += m.node_z[p] / kNodesPerCell;
    }
    volume[c] = lagrangian::cell_volume(m, c, faces[c]);
    eps[c] = s.total_energy[c] - 0.5 * (s.vel_x[c] * s.vel_x[c] + s.vel_y[c] * s.vel_y[c] +
                                        s.vel_z[c] * s.vel_z[c]);
  }
  nlohmann::json j;
  j["format"] = "hydro-dump-0";
  j["t"] = t;
  j["steps"] = steps;
  j["cells"] = {
      {"centroid_x", cx},
      {"centroid_y", cy},
      {"centroid_z", cz},
      {"volume", volume},
      {"mass", s.mass},
      {"density", thermo.density},
      {"velocity_x", s.vel_x},
      {"velocity_y", s.vel_y},
      {"velocity_z", s.vel_z},
      {"pressure", thermo.pressure},
      {"internal_energy", eps},
  };
  j["nodes"] = {{"x", m.node_x}, {"y", m.node_y}, {"z", m.node_z}};

  std::ofstream out(path);
  if (!out) throw std::runtime_error("dump: cannot open " + path + " for writing");
  out << j.dump() << '\n';
  if (!out) throw std::runtime_error("dump: failed writing " + path);
}

}  // namespace hydro::io
