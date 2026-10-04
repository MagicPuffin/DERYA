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

mesh::StructuredMeshSpec parse_mesh(const json& j, MeshPerturbation& perturbation) {
  const std::string path = "mesh";
  object(j, path, {"cells", "min", "max", "perturbation"});
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

bc::BoundarySet parse_boundaries(const json& j) {
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
      object(b, key_path, {"type", "velocity"});
      if (string(required(b, key_path, "type"), key_path + ".type") != "piston") {
        fail(key_path + ".type", "expected \"piston\" (symmetry and outflow are plain strings)");
      }
      bcs[side] = bc::Boundary(bc::BoundaryType::Piston,
                               vec3(required(b, key_path, "velocity"), key_path + ".velocity"));
      continue;
    }
    if (!b.is_string()) fail(key_path, "expected a string or a piston object");
    const std::string type = b.get<std::string>();
    if (type == "symmetry") {
      bcs[side] = bc::BoundaryType::Symmetry;
    } else if (type == "outflow") {
      bcs[side] = bc::BoundaryType::Outflow;
    } else {
      fail(key_path, "expected \"symmetry\", \"outflow\" or {\"type\": \"piston\", ...}");
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

void parse_initial(const json& j, Deck& deck) {
  const std::string path = "initial";
  object(j, path, {"background", "regions"});
  const std::string bg_path = path + ".background";
  const json& bg = required(j, path, "background");
  object(bg, bg_path, {"density", "velocity", "pressure"});
  deck.background = parse_state(bg, bg_path);
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
  object(j, "deck", {"mesh", "eos", "boundaries", "initial", "time", "output"});
  Deck deck;
  deck.mesh = parse_mesh(required(j, "deck", "mesh"), deck.perturbation);
  deck.gamma = parse_eos(required(j, "deck", "eos"));
  deck.boundaries = parse_boundaries(required(j, "deck", "boundaries"));
  parse_initial(required(j, "deck", "initial"), deck);
  parse_time(required(j, "deck", "time"), deck);
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
  const std::array<double, 3> width = {spec.x_max - spec.x_min, spec.y_max - spec.y_min,
                                       spec.z_max - spec.z_min};
  static constexpr std::array<char, 3> kAxis = {'x', 'y', 'z'};
  double in_plane = std::numeric_limits<double>::infinity();
  for (int d = 0; d < 3; ++d) {
    if (n[d] > 1) in_plane = std::min(in_plane, width[d] / n[d]);
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
