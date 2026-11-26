#pragma once

#include <array>
#include <iosfwd>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "bc.hpp"
#include "mesh.hpp"
#include "scheme.hpp"
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

// Axis-aligned box [min, max) of the initial (perturbed) mesh. A cell belongs
// to the box when its centroid (mean of its nodes) is inside.
struct Region {
  Vec3 min = {0.0, 0.0, 0.0};
  Vec3 max = {0.0, 0.0, 0.0};
  RegionState state;
};

// Taylor-Green vortex (Georges' thesis Eqs. 3.4-3.5), an initial profile in
// place of the background and regions, sampled at each cell's centroid:
//   V = c1 (sin(pi x) cos(pi y), -cos(pi x) sin(pi y), 0),
//   P = density c1^2 / 4 (cos(2 pi x) + cos(2 pi y)) + c2, rho = density.
// c2 > density c1^2 / 2 keeps P positive (decisions.md).
struct TaylorGreen {
  double density = 1.0;
  double c1 = 1.0;
  double c2 = 1.0;
};

// Prescribed energy source. TaylorGreen: thesis Eq. 3.6, with the constants
// of the deck's Taylor-Green profile, which keeps that profile stationary.
enum class EnergySourceType { None, TaylorGreen };

// Node perturbation applied to the structured mesh.
enum class MeshPerturbation { None, Saltzman };

// Input deck of hydro_run (JSON; schema in docs/deck_format.md). Only what
// Sod, Noh and Saltzman-type setups need: one structured mesh (optionally
// perturbed), an ideal gas, one boundary per side, a background state
// overridden by boxes in order (the last box containing a cell wins), and
// time control, and the scheme (first order unless the deck asks for more);
// or, instead of the background and boxes, the Taylor-Green profile, with
// its energy source.
struct Deck {
  mesh::StructuredMeshSpec mesh;
  MeshPerturbation perturbation = MeshPerturbation::None;
  double gamma = 0.0;
  bc::BoundarySet boundaries;
  RegionState background;
  std::vector<Region> regions;
  // Set: the initial state is this profile (no background, no regions).
  std::optional<TaylorGreen> taylor_green;
  EnergySourceType source = EnergySourceType::None;
  double t_end = 0.0;
  // dt^0 = min(dt_initial, CFL); infinity means CFL alone.
  double dt_initial = std::numeric_limits<double>::infinity();
  timestep::TimestepParams timestep;
  lagrangian::Scheme scheme;
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
// reason. Empty when every direction has several cells, or none has. On a
// polar mesh the y cells are arcs, measured at the middle of the innermost
// ring.
std::vector<std::string> thin_slab_warnings(const mesh::StructuredMeshSpec& spec);

}  // namespace hydro::io
