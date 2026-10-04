#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "deck.hpp"
#include "driver.hpp"
#include "hydro.hpp"
#include "mesh.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using hydro::bc::BoundaryType;
using hydro::io::Deck;
using hydro::io::parse_deck;
using hydro::io::thin_slab_warnings;
using hydro::lagrangian::HydroState;
using hydro::mesh::Index;
using hydro::mesh::StructuredMeshSpec;
using nlohmann::json;

namespace {

// A small valid deck: 4x1x1 cells of [0,1]x[0,0.25]x[0,0.25], a two-state
// tube with walls everywhere.
json small_deck() {
  return json::parse(R"({
    "mesh": {"cells": [4, 1, 1], "min": [0, 0, 0], "max": [1, 0.25, 0.25]},
    "eos": {"type": "ideal_gas", "gamma": 1.4},
    "boundaries": {"x_min": "symmetry", "x_max": "outflow", "y_min": "symmetry",
                   "y_max": "symmetry", "z_min": "symmetry", "z_max": "symmetry"},
    "initial": {
      "background": {"density": 0.125, "pressure": 0.1},
      "regions": [{"min": [0, 0, 0], "max": [0.5, 1, 1], "density": 1,
                   "velocity": [0.5, 0, 0], "pressure": 1}]
    },
    "time": {"t_end": 0.2}
  })");
}

void check_rejected(const json& j, const std::string& message) {
  CHECK_THROWS_WITH(parse_deck(j.dump()), ContainsSubstring(message));
}

// Runs the hydro_run executable on decks/<name>.json and returns its dump.
json run_hydro(const std::string& name) {
  const auto dump_path =
      (std::filesystem::temp_directory_path() / ("hydro_test_" + name + "_dump.json")).string();
  std::filesystem::remove(dump_path);
  const std::string cmd = std::string("\"") + HYDRO_RUN_EXE + "\" \"" + DECK_DIR + "/" +
                          name + ".json\" \"" + dump_path + "\" > /dev/null";
  REQUIRE(std::system(cmd.c_str()) == 0);
  std::ifstream in(dump_path);
  REQUIRE(in);
  const json dump = json::parse(in);
  std::filesystem::remove(dump_path);
  CHECK(dump["format"] == "hydro-dump-0");
  return dump;
}

json read_oracle(const std::string& name) {
  std::ifstream in(std::string(ORACLE_DATA_DIR) + "/" + name + ".json");
  REQUIRE(in);
  return json::parse(in);
}

// Linear interpolation of an oracle field, except across the oracle's
// "jumps" (when it lists them): a point on either side of a jump between two
// samples takes the sample on its own side, since interpolating would invent
// values the exact solution never takes.
double exact_at(const json& oracle, const std::vector<double>& field, double x) {
  const auto xs = oracle["x"].get<std::vector<double>>();
  std::size_t i = 1;
  while (i + 1 < xs.size() && xs[i] < x) ++i;
  if (oracle.contains("jumps")) {
    for (const double jump : oracle["jumps"].get<std::vector<double>>()) {
      if (xs[i - 1] < jump && jump <= xs[i]) return x < jump ? field[i - 1] : field[i];
    }
  }
  const double w = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
  return (1.0 - w) * field[i - 1] + w * field[i];
}

// L1 density error of a slab dump: sum_c |rho_c - rho(x_c)| V_c / area, with
// area the fixed y-z cross-section, so V_c / area is the cell's length.
double l1_density_error(const json& dump, const json& oracle, double area) {
  const auto& cells = dump["cells"];
  const auto x = cells["centroid_x"].get<std::vector<double>>();
  const auto volume = cells["volume"].get<std::vector<double>>();
  const auto rho = cells["density"].get<std::vector<double>>();
  const auto rho_ex = oracle["density"].get<std::vector<double>>();
  double l1 = 0.0;
  for (std::size_t c = 0; c < x.size(); ++c) {
    l1 += std::abs(rho[c] - exact_at(oracle, rho_ex, x[c])) * volume[c] / area;
  }
  return l1;
}

}  // namespace

TEST_CASE("parse_deck: a valid deck, with defaults", "[io]") {
  const Deck d = parse_deck(small_deck().dump());
  CHECK(d.mesh.nx == 4);
  CHECK(d.mesh.ny == 1);
  CHECK(d.mesh.nz == 1);
  CHECK(d.mesh.x_max == 1.0);
  CHECK(d.mesh.z_max == 0.25);
  CHECK(d.gamma == 1.4);
  CHECK(d.boundaries[0] == BoundaryType::Symmetry);
  CHECK(d.boundaries[1] == BoundaryType::Outflow);
  CHECK(d.background.density == 0.125);
  CHECK(d.background.velocity == hydro::io::Vec3{0.0, 0.0, 0.0});
  REQUIRE(d.regions.size() == 1);
  CHECK(d.regions[0].max == hydro::io::Vec3{0.5, 1.0, 1.0});
  CHECK(d.regions[0].state.velocity == hydro::io::Vec3{0.5, 0.0, 0.0});
  CHECK(d.t_end == 0.2);
  CHECK(std::isinf(d.dt_initial));
  CHECK(d.timestep.cfl == 0.25);
  CHECK(d.timestep.max_growth == 1.01);
  CHECK(d.output.empty());
  CHECK(d.log_every == 100);
  CHECK(d.warnings.empty());

  json j = small_deck();
  j["time"]["dt_initial"] = 1e-4;
  j["time"]["cfl"] = 0.4;
  j["time"]["max_steps"] = 10;
  j["output"] = {{"path", "out.json"}, {"log_every", 0}};
  const Deck e = parse_deck(j.dump());
  CHECK(e.dt_initial == 1e-4);
  CHECK(e.timestep.cfl == 0.4);
  CHECK(e.max_steps == 10);
  CHECK(e.output == "out.json");
  CHECK(e.log_every == 0);
}

TEST_CASE("parse_deck: rejects bad decks, naming the key", "[io]") {
  CHECK_THROWS_WITH(parse_deck("{not json"), ContainsSubstring("invalid JSON"));
  json j = small_deck();
  j["time"]["t_ned"] = 0.2;
  check_rejected(j, "time.t_ned: unknown key");
  j = small_deck();
  j.erase("eos");
  check_rejected(j, "deck.eos: missing");
  j = small_deck();
  j["boundaries"].erase("z_max");
  check_rejected(j, "boundaries.z_max: missing");
  j = small_deck();
  j["boundaries"]["y_min"] = "periodic";
  check_rejected(j, "boundaries.y_min");
  j = small_deck();
  j["eos"]["type"] = "stiffened_gas";
  check_rejected(j, "eos.type");
  j = small_deck();
  j["mesh"]["cells"] = {4, 0, 1};
  check_rejected(j, "mesh.cells");
  j = small_deck();
  j["mesh"]["cells"] = {4.5, 1, 1};
  check_rejected(j, "mesh.cells[0]: expected an integer");
  j = small_deck();
  j["mesh"]["max"] = {1, 0, 0.25};
  check_rejected(j, "mesh.max");
  j = small_deck();
  j["initial"]["regions"][0]["density"] = -1;
  check_rejected(j, "initial.regions[0].density");
  j = small_deck();
  j["initial"]["background"]["pressure"] = "1";
  check_rejected(j, "initial.background.pressure: expected a number");
  j = small_deck();
  j["time"]["max_growth"] = 0.9;
  check_rejected(j, "time.max_growth");
}

TEST_CASE("thin_slab_warnings: a thin extrusion warns, cubic cells do not", "[io]") {
  // 100x1x1 cells of [0,1]x[0,0.01]x[0,0.005]: dx = 0.01, z only half as wide.
  StructuredMeshSpec thin{100, 1, 1, 0.0, 1.0, 0.0, 0.01, 0.0, 0.005};
  const auto w = thin_slab_warnings(thin);
  REQUIRE(w.size() == 1);
  CHECK_THAT(w[0], ContainsSubstring("in z"));

  // Cubic cells, as in the Noh test.
  CHECK(thin_slab_warnings({100, 1, 1, 0.0, 1.0, 0.0, 0.01, 0.0, 0.01}).empty());
  // No single-cell direction (the Sod slab), or no multi-cell one.
  CHECK(thin_slab_warnings({100, 2, 2, 0.0, 1.0, 0.0, 0.002, 0.0, 0.002}).empty());
  CHECK(thin_slab_warnings({1, 1, 1, 0.0, 1.0, 0.0, 0.1, 0.0, 0.5}).empty());

  // Read from a deck.
  json j = small_deck();
  j["mesh"]["max"] = {1, 0.25, 0.1};
  const Deck d = parse_deck(j.dump());
  REQUIRE(d.warnings.size() == 1);
  CHECK_THAT(d.warnings[0], ContainsSubstring("thin slab"));
}

TEST_CASE("initial_state: regions override the background", "[app]") {
  const HydroState s = hydro::app::initial_state(parse_deck(small_deck().dump()));
  REQUIRE(s.mesh.num_cells() == 4);
  const double cell_volume = 0.25 * 0.25 * 0.25;
  for (Index c = 0; c < 4; ++c) {
    INFO("cell " << c);
    const bool left = c < 2;
    const double rho = left ? 1.0 : 0.125;
    const double eps = left ? 1.0 / (0.4 * 1.0) : 0.1 / (0.4 * 0.125);
    const double u = left ? 0.5 : 0.0;
    CHECK_THAT(s.mass[c], WithinRel(rho * cell_volume, 1e-14));
    CHECK(s.vel_x[c] == u);
    CHECK_THAT(s.total_energy[c], WithinRel(eps + 0.5 * u * u, 1e-14));
  }
}

TEST_CASE("run: ends exactly at t_end; max_steps stops it", "[app]") {
  json j = small_deck();
  j["output"] = {{"log_every", 1}};
  std::ostringstream log;
  const auto r = hydro::app::run(parse_deck(j.dump()), log);
  CHECK(r.t == 0.2);
  CHECK(r.steps > 1);
  CHECK_THAT(log.str(), ContainsSubstring("step 1 "));
  CHECK_THAT(log.str(), ContainsSubstring("done: " + std::to_string(r.steps) + " steps"));

  j["time"]["max_steps"] = 2;
  std::ostringstream quiet;
  CHECK_THROWS_WITH(hydro::app::run(parse_deck(j.dump()), quiet),
                    ContainsSubstring("max_steps = 2"));
}

TEST_CASE("run: a failing step reports where it failed", "[app]") {
  // A cold gas at rest with no dt_initial: nothing bounds the first dt.
  json j = small_deck();
  j["initial"] = {{"background", {{"density", 1.0}, {"pressure", 0.0}}}};
  std::ostringstream log;
  CHECK_THROWS_WITH(hydro::app::run(parse_deck(j.dump()), log),
                    ContainsSubstring("step 1 from t = 0"));
}

TEST_CASE("hydro_run: Sod against the exact solution", "[app][integration]") {
  // First integration milestone (v0.1): decks/sod.json (100x2x2 cells, the
  // slab of the lagrangian_step Sod test, dt from the controller) run by the
  // hydro_run executable, its dump compared with tests/oracle_data/sod.json.
  // Measured: 211 steps, L1 density error 0.0195 (0.0202 with the fixed
  // dt = 5e-4 of the lagrangian_step test).
  const json dump = run_hydro("sod");
  const json oracle = read_oracle("sod");
  REQUIRE(oracle["t"].get<double>() == 0.2);
  CHECK(dump["t"].get<double>() == 0.2);
  const auto& cells = dump["cells"];
  const auto volume = cells["volume"].get<std::vector<double>>();
  const auto mass = cells["mass"].get<std::vector<double>>();
  const auto rho = cells["density"].get<std::vector<double>>();
  const auto vy = cells["velocity_y"].get<std::vector<double>>();
  const auto vz = cells["velocity_z"].get<std::vector<double>>();
  REQUIRE(rho.size() == 400);
  CHECK(dump["nodes"]["x"].size() == 101 * 3 * 3);
  for (std::size_t c = 0; c < rho.size(); ++c) {
    INFO("cell " << c);
    CHECK_THAT(vy[c], WithinAbs(0.0, 1e-12));
    CHECK_THAT(vz[c], WithinAbs(0.0, 1e-12));
    CHECK_THAT(rho[c], WithinRel(mass[c] / volume[c], 1e-12));
  }
  // The sum over the 4 columns' cells, so the mean of their L1 errors.
  const double l1 = l1_density_error(dump, oracle, 0.02 * 0.02);
  INFO("L1 density error " << l1);
  CHECK(l1 < 0.021);
}

TEST_CASE("hydro_run: planar Noh against the exact solution", "[app][integration]") {
  // decks/noh.json: cold gas (gamma = 5/3, rho = 1, u = -1) on 100x1x1 cubic
  // cells of [0,1]x[0,0.01]x[0,0.01], wall at x = 0, outflow at x = 1 (which
  // moves with the inflow, to x = 0.4 at t = 0.6), dt^0 = 1e-4, compared with
  // tests/oracle_data/noh.json: shock at x = 0.2, rho = 4, eps = 1/2 and
  // u = 0 behind it. Measured: 856 steps, plateau rho 3.990-4.001 and eps
  // 0.4999-0.5013 on [0.04, 0.18], wall heating in the cells at the wall
  // (rho 2.51, 3.64), L1 density error 0.0194 (0.0188 with the fixed
  // dt = 5e-4 of the lagrangian_step Noh test).
  const json dump = run_hydro("noh");
  const json oracle = read_oracle("noh");
  REQUIRE(oracle["t"].get<double>() == 0.6);
  CHECK(dump["t"].get<double>() == 0.6);
  const auto& cells = dump["cells"];
  const auto x = cells["centroid_x"].get<std::vector<double>>();
  const auto mass = cells["mass"].get<std::vector<double>>();
  const auto rho = cells["density"].get<std::vector<double>>();
  const auto u = cells["velocity_x"].get<std::vector<double>>();
  const auto vy = cells["velocity_y"].get<std::vector<double>>();
  const auto vz = cells["velocity_z"].get<std::vector<double>>();
  const auto eps = cells["internal_energy"].get<std::vector<double>>();
  REQUIRE(rho.size() == 100);

  // Neither the wall nor the cold outflow (P* = 0) does work: the total
  // energy stays that of the inflow, sum m |u|^2 / 2.
  double e0 = 0.0, e = 0.0;
  for (std::size_t c = 0; c < rho.size(); ++c) {
    e0 += 0.5 * mass[c];
    e += mass[c] * (eps[c] + 0.5 * (u[c] * u[c] + vy[c] * vy[c] + vz[c] * vz[c]));
  }
  CHECK_THAT(e, WithinRel(e0, 1e-12));

  for (std::size_t c = 0; c < rho.size(); ++c) {
    INFO("cell " << c << " at x = " << x[c]);
    CHECK_THAT(vy[c], WithinAbs(0.0, 1e-12));
    CHECK_THAT(vz[c], WithinAbs(0.0, 1e-12));
    if (x[c] > 0.04 && x[c] < 0.18) {
      CHECK_THAT(rho[c], WithinRel(4.0, 0.01));
      CHECK_THAT(eps[c], WithinRel(0.5, 0.01));
      CHECK_THAT(u[c], WithinAbs(0.0, 1e-4));
    }
    if (x[c] > 0.25) {
      CHECK_THAT(rho[c], WithinRel(1.0, 1e-10));
      CHECK_THAT(u[c], WithinAbs(-1.0, 1e-10));
    }
  }
  const double l1 = l1_density_error(dump, oracle, 0.01 * 0.01);
  INFO("L1 density error " << l1);
  CHECK(l1 < 0.021);
}
