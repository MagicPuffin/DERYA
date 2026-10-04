#pragma once

#include <array>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

#include "bc.hpp"
#include "mesh.hpp"
#include "timestep.hpp"

namespace hydro::io {

using mesh::Vec3;

// Uniform initial state of a region: density, velocity and pressure (the
// specific internal energy follows from the EOS).
struct RegionState {
  double density = 0.0;
  Vec3 velocity = {0.0, 0.0, 0.0};
  double pressure = 0.0;
};

// Axis-aligned box [min, max) of the initial mesh. A cell belongs to the box
// when its centroid (mean of its nodes) is inside.
struct Region {
  Vec3 min = {0.0, 0.0, 0.0};
  Vec3 max = {0.0, 0.0, 0.0};
  RegionState state;
};

// Input deck of hydro_run (JSON; schema in docs/deck_format.md). Only what
// Sod, Noh and Saltzman-type setups need: one structured mesh, an ideal gas,
// one boundary type per side, a background state overridden by boxes in
// order (the last box containing a cell wins), and time control.
struct Deck {
  mesh::StructuredMeshSpec mesh;
  double gamma = 0.0;
  bc::BoundarySet boundaries;
  RegionState background;
  std::vector<Region> regions;
  double t_end = 0.0;
  // dt^0 = min(dt_initial, CFL); infinity means CFL alone.
  double dt_initial = std::numeric_limits<double>::infinity();
  timestep::TimestepParams timestep;
  long max_steps = std::numeric_limits<long>::max();
  // Steps between progress lines in the log; 0 for none.
  long log_every = 100;
  // Final-state dump; empty when the deck names none.
  std::string output;
  // Non-fatal problems found while reading (e.g. the thin-slab guard).
  std::vector<std::string> warnings;
};

// Parses a deck. Throws std::runtime_error naming the offending key on a
// missing required key, an unknown key (to catch typos), a value of the
// wrong type, or an out-of-range value.
Deck parse_deck(const std::string& json_text);
Deck read_deck(const std::string& path);

// Thin-slab guard (decisions.md, "Lagrangian timestep"): one warning per
// single-cell direction that is thinner than the smallest cell size of the
// multi-cell directions, since lambda_c would then limit dt for no physical
// reason. Empty when every direction has several cells, or none has.
std::vector<std::string> thin_slab_warnings(const mesh::StructuredMeshSpec& spec);

}  // namespace hydro::io
