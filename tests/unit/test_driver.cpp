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
  const auto dump_path =
      (std::filesystem::temp_directory_path() / "hydro_test_sod_dump.json").string();
  std::filesystem::remove(dump_path);
  const std::string cmd = std::string("\"") + HYDRO_RUN_EXE + "\" \"" + DECK_DIR +
                          "/sod.json\" \"" + dump_path + "\" > /dev/null";
  REQUIRE(std::system(cmd.c_str()) == 0);

  std::ifstream in(dump_path);
  REQUIRE(in);
  const auto dump = json::parse(in);
  CHECK(dump["format"] == "hydro-dump-0");
  CHECK(dump["t"].get<double>() == 0.2);
  const auto& cells = dump["cells"];
  const auto x = cells["centroid_x"].get<std::vector<double>>();
  const auto volume = cells["volume"].get<std::vector<double>>();
  const auto mass = cells["mass"].get<std::vector<double>>();
  const auto rho = cells["density"].get<std::vector<double>>();
  const auto vy = cells["velocity_y"].get<std::vector<double>>();
  const auto vz = cells["velocity_z"].get<std::vector<double>>();
  REQUIRE(x.size() == 400);
  CHECK(dump["nodes"]["x"].size() == 101 * 3 * 3);

  std::ifstream oin(ORACLE_DATA_DIR "/sod.json");
  REQUIRE(oin);
  const auto oracle = json::parse(oin);
  REQUIRE(oracle["t"].get<double>() == 0.2);
  const auto x_ex = oracle["x"].get<std::vector<double>>();
  const auto rho_ex = oracle["density"].get<std::vector<double>>();
  auto exact = [&](double xc) {
    std::size_t i = 1;
    while (i + 1 < x_ex.size() && x_ex[i] < xc) ++i;
    const double w = (xc - x_ex[i - 1]) / (x_ex[i] - x_ex[i - 1]);
    return (1.0 - w) * rho_ex[i - 1] + w * rho_ex[i];
  };

  // The y/z walls do not move, so V_c / (0.02 * 0.02) is the cell's length;
  // the sum is the mean L1 error of the 4 columns.
  const double area = 0.02 * 0.02;
  double l1 = 0.0;
  for (std::size_t c = 0; c < x.size(); ++c) {
    INFO("cell " << c);
    CHECK_THAT(vy[c], WithinAbs(0.0, 1e-12));
    CHECK_THAT(vz[c], WithinAbs(0.0, 1e-12));
    CHECK_THAT(rho[c], WithinRel(mass[c] / volume[c], 1e-12));
    l1 += std::abs(rho[c] - exact(x[c])) * volume[c] / area;
  }
  INFO("L1 density error " << l1);
  CHECK(l1 < 0.021);
  std::filesystem::remove(dump_path);
}
