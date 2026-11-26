#include "deck.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <stdexcept>

namespace hydro::io {

using nlohmann::json;

namespace {

[[noreturn]] void fail(const std::string& path, const std::string& what) {
  throw std::runtime_error("deck: " + path + ": " + what);
}

// Object at path, with no keys outside `allowed`.
const json& object(const json& j, const std::string& path,
                   std::initializer_list<const char*> allowed) {
  if (!j.is_object()) fail(path, "expected an object");
  for (const auto& item : j.items()) {
    const bool known = std::any_of(allowed.begin(), allowed.end(),
                                   [&](const char* k) { return item.key() == k; });
    if (!known) fail(path + "." + item.key(), "unknown key");
  }
  return j;
}

const json& required(const json& j, const std::string& path, const char* key) {
  if (!j.contains(key)) fail(path + "." + key, "missing");
  return j.at(key);
}

double number(const json& j, const std::string& path) {
  if (!j.is_number()) fail(path, "expected a number");
  return j.get<double>();
}

long integer(const json& j, const std::string& path) {
  if (!j.is_number_integer()) fail(path, "expected an integer");
  return j.get<long>();
}

std::string string(const json& j, const std::string& path) {
  if (!j.is_string()) fail(path, "expected a string");
  return j.get<std::string>();
}

Vec3 vec3(const json& j, const std::string& path) {
  if (!j.is_array() || j.size() != 3) fail(path, "expected an array of 3 numbers");
  Vec3 v;
  for (int d = 0; d < 3; ++d) v[d] = number(j[d], path + "[" + std::to_string(d) + "]");
  return v;
}

// Optional number, checked with `valid` when present.
template <typename Valid>
void optional_number(const json& j, const std::string& path, const char* key,
                     double& out, Valid valid, const char* requirement) {
  if (!j.contains(key)) return;
  out = number(j.at(key), path + "." + key);
  if (!valid(out)) fail(path + "." + key, requirement);
}

// mesh.x_segments: [[x0, x1, n], ...], contiguous from mesh.min to mesh.max
// in x, n uniform cells each, adding up to mesh.cells[0]. Returns the x
// node positions.
std::vector<double> parse_x_segments(const json& j, const std::string& path,
                                     const mesh::StructuredMeshSpec& spec) {
  if (!j.is_array() || j.empty()) fail(path, "expected a non-empty array of [x0, x1, n]");
  std::vector<double> nodes = {spec.x_min};
  for (std::size_t s = 0; s < j.size(); ++s) {
    const std::string s_path = path + "[" + std::to_string(s) + "]";
    const json& seg = j[s];
    if (!seg.is_array() || seg.size() != 3) fail(s_path, "expected [x0, x1, n]");
    const double x0 = number(seg[0], s_path + "[0]");
    const double x1 = number(seg[1], s_path + "[1]");
    const long n = integer(seg[2], s_path + "[2]");
    if (n < 1) fail(s_path + "[2]", "must be >= 1");
    if (!(x1 > x0)) fail(s_path, "x1 must exceed x0");
    if (x0 != nodes.back()) {
      fail(s_path + "[0]", s == 0 ? "must equal mesh.min[0]" : "must equal the previous x1");
    }
    // As generate_structured_mesh: one segment gives the uniform nodes.
    const double dx = (x1 - x0) / n;
    for (long i = 1; i < n; ++i) nodes.push_back(x0 + i * dx);
    nodes.push_back(x1);
  }
  if (nodes.back() != spec.x_max) fail(path, "the last x1 must equal mesh.max[0]");
  if (nodes.size() != static_cast<std::size_t>(spec.nx) + 1) {
    fail(path, "the cell counts must add up to mesh.cells[0]");
  }
  return nodes;
}

mesh::StructuredMeshSpec parse_mesh(const json& j, MeshPerturbation& perturbation) {
  const std::string path = "mesh";
  object(j, path, {"cells", "min", "max", "perturbation", "x_segments", "mapping"});
  if (j.contains("perturbation")) {
    if (string(j.at("perturbation"), path + ".perturbation") != "saltzman") {
      fail(path + ".perturbation", "only \"saltzman\" is supported");
    }
    perturbation = MeshPerturbation::Saltzman;
  }
  const json& cells = required(j, path, "cells");
  if (!cells.is_array() || cells.size() != 3) fail(path + ".cells", "expected an array of 3 integers");
  std::array<long, 3> n;
  for (int d = 0; d < 3; ++d) {
    n[d] = integer(cells[d], path + ".cells[" + std::to_string(d) + "]");
    if (n[d] < 1) fail(path + ".cells", "every count must be >= 1");
  }
  const Vec3 lo = vec3(required(j, path, "min"), path + ".min");
  const Vec3 hi = vec3(required(j, path, "max"), path + ".max");
  for (int d = 0; d < 3; ++d) {
    if (!(hi[d] > lo[d])) fail(path + ".max", "must exceed mesh.min in every direction");
  }
  mesh::StructuredMeshSpec spec;
  spec.nx = static_cast<mesh::Index>(n[0]);
  spec.ny = static_cast<mesh::Index>(n[1]);
  spec.nz = static_cast<mesh::Index>(n[2]);
  spec.x_min = lo[0];
  spec.y_min = lo[1];
  spec.z_min = lo[2];
  spec.x_max = hi[0];
  spec.y_max = hi[1];
  spec.z_max = hi[2];
  if (j.contains("x_segments")) {
    if (perturbation != MeshPerturbation::None) {
      fail(path + ".x_segments", "cannot be combined with mesh.perturbation");
    }
    spec.x_nodes = parse_x_segments(j.at("x_segments"), path + ".x_segments", spec);
  }
  if (j.contains("mapping")) {
    const std::string mapping = string(j.at("mapping"), path + ".mapping");
    if (mapping == "polar") {
      if (perturbation != MeshPerturbation::None) {
        fail(path + ".mapping", "\"polar\" cannot be combined with mesh.perturbation");
      }
      if (!(spec.x_min >= 0.0)) fail(path + ".min[0]", "polar: r must be >= 0");
      const double span = spec.y_max - spec.y_min;
      if (!(span <= 180.0)) fail(path + ".max[1]", "polar: theta may span at most 180 degrees");
      if (!(span / spec.ny < 180.0)) {
        fail(path + ".cells[1]", "polar: each cell must span less than 180 degrees");
      }
      spec.mapping = mesh::MeshMapping::Polar;
    } else if (mapping != "cartesian") {
      fail(path + ".mapping", "expected \"cartesian\" or \"polar\"");
    }
  }
  return spec;
}

double parse_eos(const json& j) {
  const std::string path = "eos";
  object(j, path, {"type", "gamma"});
  if (string(required(j, path, "type"), path + ".type") != "ideal_gas") {
    fail(path + ".type", "only \"ideal_gas\" is supported");
  }
  const double gamma = number(required(j, path, "gamma"), path + ".gamma");
  if (!(gamma > 1.0)) fail(path + ".gamma", "must be > 1");
  return gamma;
}

// A pressure side: {"type": "pressure", "value": P} or {"type": "pressure",
// "law": "kidder", "p0": P0, "t_focus": t} (the law's gamma is the gas's).
bc::Boundary parse_pressure_boundary(const json& b, const std::string& path, double gamma) {
  object(b, path, {"type", "value", "law", "p0", "t_focus"});
  if (b.contains("value") == b.contains("law")) {
    fail(path, "expected exactly one of \"value\" and \"law\"");
  }
  if (b.contains("value")) {
    const double p = number(b.at("value"), path + ".value");
    if (!(p >= 0.0)) fail(path + ".value", "must be >= 0");
    if (b.contains("p0") || b.contains("t_focus")) {
      fail(path, "\"p0\" and \"t_focus\" belong to a law, not a constant value");
    }
    return bc::pressure_boundary(p);
  }
  if (string(b.at("law"), path + ".law") != "kidder") {
    fail(path + ".law", "only \"kidder\" is supported");
  }
  const double p0 = number(required(b, path, "p0"), path + ".p0");
  if (!(p0 >= 0.0)) fail(path + ".p0", "must be >= 0");
  const double t_focus = number(required(b, path, "t_focus"), path + ".t_focus");
  if (!(t_focus > 0.0)) fail(path + ".t_focus", "must be > 0");
  return bc::pressure_boundary(bc::kidder_pressure(p0, t_focus, gamma));
}

bc::BoundarySet parse_boundaries(const json& j, double gamma) {
  const std::string path = "boundaries";
  // Ordered like mesh::Mesh::face_boundary.
  static constexpr std::array<const char*, mesh::kFacesPerCell> kSides = {
      "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"};
  object(j, path, {"x_min", "x_max", "y_min", "y_max", "z_min", "z_max"});
  bc::BoundarySet bcs;
  for (int side = 0; side < mesh::kFacesPerCell; ++side) {
    const std::string key_path = path + "." + kSides[side];
    const json& b = required(j, path, kSides[side]);
    if (b.is_object()) {
      const std::string type = string(required(b, key_path, "type"), key_path + ".type");
      if (type == "pressure") {
        bcs[side] = parse_pressure_boundary(b, key_path, gamma);
        continue;
      }
      if (type != "piston") {
        fail(key_path + ".type",
             "expected \"piston\" or \"pressure\" (symmetry and outflow are plain strings)");
      }
      object(b, key_path, {"type", "velocity"});
      bcs[side] = bc::Boundary(bc::BoundaryType::Piston,
                               vec3(required(b, key_path, "velocity"), key_path + ".velocity"));
      continue;
    }
    if (!b.is_string()) fail(key_path, "expected a string or an object");
    const std::string type = b.get<std::string>();
    if (type == "symmetry") {
      bcs[side] = bc::BoundaryType::Symmetry;
    } else if (type == "outflow") {
      bcs[side] = bc::BoundaryType::Outflow;
    } else {
      fail(key_path,
           "expected \"symmetry\", \"outflow\", or a {\"type\": \"piston\" or "
           "\"pressure\", ...} object");
    }
  }
  return bcs;
}

// The state keys of a background or region object (the caller has checked
// for unknown keys).
RegionState parse_state(const json& j, const std::string& path) {
  RegionState s;
  s.density = number(required(j, path, "density"), path + ".density");
  if (!(s.density > 0.0)) fail(path + ".density", "must be > 0");
  s.pressure = number(required(j, path, "pressure"), path + ".pressure");
  if (!(s.pressure >= 0.0)) fail(path + ".pressure", "must be >= 0");
  if (j.contains("velocity")) s.velocity = vec3(j.at("velocity"), path + ".velocity");
  return s;
}

void parse_taylor_green(const json& j, Deck& deck) {
  const std::string path = "initial.taylor_green";
  object(j, path, {"density", "c1", "c2"});
  TaylorGreen tg;
  optional_number(j, path, "density", tg.density, [](double v) { return v > 0.0; },
                  "must be > 0");
  if (j.contains("c1")) tg.c1 = number(j.at("c1"), path + ".c1");
  if (j.contains("c2")) tg.c2 = number(j.at("c2"), path + ".c2");
  if (!(tg.c2 > 0.5 * tg.density * tg.c1 * tg.c1)) {
    fail(path + ".c2", "must exceed density c1^2 / 2 (the pressure would not be positive)");
  }
  deck.taylor_green = tg;
}

void parse_energy_spike(const json& j, Deck& deck) {
  const std::string path = "initial.energy_spike";
  object(j, path, {"at", "energy_per_length"});
  EnergySpike spike;
  spike.at = vec3(required(j, path, "at"), path + ".at");
  spike.energy_per_length =
      number(required(j, path, "energy_per_length"), path + ".energy_per_length");
  if (!(spike.energy_per_length > 0.0)) fail(path + ".energy_per_length", "must be > 0");
  deck.energy_spike = spike;
}

void parse_initial(const json& j, Deck& deck) {
  const std::string path = "initial";
  object(j, path, {"background", "regions", "taylor_green", "energy_spike"});
  if (j.contains("taylor_green")) {
    if (j.contains("background") || j.contains("regions") || j.contains("energy_spike")) {
      fail(path + ".taylor_green", "excludes background, regions and energy_spike");
    }
    parse_taylor_green(j.at("taylor_green"), deck);
    return;
  }
  const std::string bg_path = path + ".background";
  const json& bg = required(j, path, "background");
  object(bg, bg_path, {"density", "velocity", "pressure"});
  deck.background = parse_state(bg, bg_path);
  if (j.contains("energy_spike")) parse_energy_spike(j.at("energy_spike"), deck);
  if (!j.contains("regions")) return;
  const json& regions = j.at("regions");
  if (!regions.is_array()) fail(path + ".regions", "expected an array");
  for (std::size_t i = 0; i < regions.size(); ++i) {
    const std::string r_path = path + ".regions[" + std::to_string(i) + "]";
    const json& r = regions[i];
    object(r, r_path, {"min", "max", "density", "velocity", "pressure"});
    Region region;
    region.min = vec3(required(r, r_path, "min"), r_path + ".min");
    region.max = vec3(required(r, r_path, "max"), r_path + ".max");
    region.state = parse_state(r, r_path);
    deck.regions.push_back(region);
  }
}

void parse_time(const json& j, Deck& deck) {
  const std::string path = "time";
  object(j, path, {"t_end", "dt_initial", "cfl", "volume_fraction", "max_growth",
                   "dt_max", "max_steps"});
  deck.t_end = number(required(j, path, "t_end"), path + ".t_end");
  if (!(deck.t_end > 0.0)) fail(path + ".t_end", "must be > 0");
  auto positive = [](double v) { return v > 0.0; };
  optional_number(j, path, "dt_initial", deck.dt_initial, positive, "must be > 0");
  optional_number(j, path, "cfl", deck.timestep.cfl, positive, "must be > 0");
  optional_number(j, path, "volume_fraction", deck.timestep.volume_fraction, positive,
                  "must be > 0");
  optional_number(j, path, "max_growth", deck.timestep.max_growth,
                  [](double v) { return v >= 1.0; }, "must be >= 1");
  optional_number(j, path, "dt_max", deck.timestep.dt_max, positive, "must be > 0");
  if (j.contains("max_steps")) {
    deck.max_steps = integer(j.at("max_steps"), path + ".max_steps");
    if (deck.max_steps < 1) fail(path + ".max_steps", "must be >= 1");
  }
}

void parse_scheme(const json& j, Deck& deck) {
  const std::string path = "scheme";
  object(j, path, {"order", "limiter"});
  if (j.contains("order")) {
    const long order = integer(j.at("order"), path + ".order");
    if (order != 1 && order != 2) fail(path + ".order", "must be 1 or 2");
    deck.scheme.order = static_cast<int>(order);
  }
  if (j.contains("limiter")) {
    const std::string limiter = string(j.at("limiter"), path + ".limiter");
    if (limiter == "sp-v") {
      deck.scheme.limiter = lagrangian::Limiter::SpV;
    } else if (limiter == "none") {
      deck.scheme.limiter = lagrangian::Limiter::None;
    } else {
      fail(path + ".limiter", "expected \"sp-v\" or \"none\"");
    }
  }
}

void parse_source(const json& j, Deck& deck) {
  if (string(j, "source") != "taylor_green") fail("source", "only \"taylor_green\" is supported");
  if (!deck.taylor_green) fail("source", "\"taylor_green\" needs initial.taylor_green");
  deck.source = EnergySourceType::TaylorGreen;
}

void parse_output(const json& j, Deck& deck) {
  const std::string path = "output";
  object(j, path, {"path", "log_every"});
  if (j.contains("path")) deck.output = string(j.at("path"), path + ".path");
  if (j.contains("log_every")) {
    deck.log_every = integer(j.at("log_every"), path + ".log_every");
    if (deck.log_every < 0) fail(path + ".log_every", "must be >= 0");
  }
}

}  // namespace

Deck parse_deck(const std::string& json_text) {
  json j;
  try {
    j = json::parse(json_text);
  } catch (const json::parse_error& e) {
    throw std::runtime_error(std::string("deck: invalid JSON: ") + e.what());
  }
  object(j, "deck",
         {"mesh", "eos", "boundaries", "initial", "source", "time", "scheme", "output"});
  Deck deck;
  deck.mesh = parse_mesh(required(j, "deck", "mesh"), deck.perturbation);
  deck.gamma = parse_eos(required(j, "deck", "eos"));
  deck.boundaries = parse_boundaries(required(j, "deck", "boundaries"), deck.gamma);
  parse_initial(required(j, "deck", "initial"), deck);
  if (j.contains("source")) parse_source(j.at("source"), deck);
  parse_time(required(j, "deck", "time"), deck);
  if (j.contains("scheme")) parse_scheme(j.at("scheme"), deck);
  if (j.contains("output")) parse_output(j.at("output"), deck);
  deck.warnings = thin_slab_warnings(deck.mesh);
  return deck;
}

Deck read_deck(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("deck: cannot open " + path);
  std::ostringstream text;
  text << in.rdbuf();
  return parse_deck(text.str());
}

std::vector<std::string> thin_slab_warnings(const mesh::StructuredMeshSpec& spec) {
  const std::array<mesh::Index, 3> n = {spec.nx, spec.ny, spec.nz};
  std::array<double, 3> width = {spec.x_max - spec.x_min, spec.y_max - spec.y_min,
                                 spec.z_max - spec.z_min};
  if (spec.mapping == mesh::MeshMapping::Polar) {
    // y is theta in degrees: its cells are arcs, the shortest at the middle
    // of the innermost ring.
    const double r0 = spec.x_nodes.empty() ? spec.x_min : spec.x_nodes[0];
    const double r1 = spec.x_nodes.empty() ? spec.x_min + width[0] / spec.nx : spec.x_nodes[1];
    width[1] = 0.5 * (r0 + r1) * width[1] * 3.14159265358979323846 / 180.0;
  }
  static constexpr std::array<char, 3> kAxis = {'x', 'y', 'z'};
  double in_plane = std::numeric_limits<double>::infinity();
  for (int d = 0; d < 3; ++d) {
    if (n[d] > 1) in_plane = std::min(in_plane, width[d] / n[d]);
  }
  for (std::size_t i = 1; i < spec.x_nodes.size(); ++i) {
    in_plane = std::min(in_plane, spec.x_nodes[i] - spec.x_nodes[i - 1]);
  }
  std::vector<std::string> warnings;
  if (std::isinf(in_plane)) return warnings;  // no multi-cell direction
  for (int d = 0; d < 3; ++d) {
    if (n[d] == 1 && width[d] < in_plane) {
      std::ostringstream w;
      w << "thin slab: the single cell in " << kAxis[d] << " is " << width[d]
        << " wide, thinner than the smallest in-plane cell size " << in_plane
        << "; the CFL length lambda_c will limit dt more than the flow needs";
      warnings.push_back(w.str());
    }
  }
  return warnings;
}

}  // namespace hydro::io
