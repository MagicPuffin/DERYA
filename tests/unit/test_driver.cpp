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
#include "eos.hpp"
#include "gradients.hpp"
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
using hydro::lagrangian::Vec3;
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

// Runs the hydro_run executable on the deck file and returns its dump,
// written to a temporary file named after `name`.
json run_hydro_deck(const std::string& deck_path, const std::string& name) {
  const auto dump_path =
      (std::filesystem::temp_directory_path() / ("hydro_test_" + name + "_dump.json")).string();
  std::filesystem::remove(dump_path);
  const std::string cmd = std::string("\"") + HYDRO_RUN_EXE + "\" \"" + deck_path +
                          "\" \"" + dump_path + "\" > /dev/null";
  REQUIRE(std::system(cmd.c_str()) == 0);
  std::ifstream in(dump_path);
  REQUIRE(in);
  const json dump = json::parse(in);
  std::filesystem::remove(dump_path);
  CHECK(dump["format"] == "hydro-dump-0");
  return dump;
}

// Runs the hydro_run executable on decks/<name>.json and returns its dump.
json run_hydro(const std::string& name) {
  return run_hydro_deck(std::string(DECK_DIR) + "/" + name + ".json", name);
}

json read_deck_json(const std::string& name) {
  std::ifstream in(std::string(DECK_DIR) + "/" + name + ".json");
  REQUIRE(in);
  return json::parse(in);
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
  CHECK(d.boundaries[0].type == BoundaryType::Symmetry);
  CHECK(d.boundaries[1].type == BoundaryType::Outflow);
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

TEST_CASE("parse_deck: piston boundaries and the Saltzman perturbation", "[io]") {
  json j = small_deck();
  j["boundaries"]["x_min"] = {{"type", "piston"}, {"velocity", {1.0, 0.0, 0.0}}};
  j["mesh"]["perturbation"] = "saltzman";
  const Deck d = parse_deck(j.dump());
  CHECK(d.boundaries[0].type == BoundaryType::Piston);
  CHECK(d.boundaries[0].velocity == hydro::io::Vec3{1.0, 0.0, 0.0});
  CHECK(d.boundaries[1].type == BoundaryType::Outflow);
  CHECK(d.perturbation == hydro::io::MeshPerturbation::Saltzman);
  CHECK(parse_deck(small_deck().dump()).perturbation == hydro::io::MeshPerturbation::None);

  // initial_state applies the skew: node (i, j) = (2, 0) of the 4x1 mesh of
  // [0,1]x[0,0.25] moves by 0.25 sin(pi/2).
  const HydroState s = hydro::app::initial_state(d);
  CHECK_THAT(s.mesh.node_x[2], WithinAbs(0.75, 1e-15));
  CHECK_THAT(s.mesh.node_x[5 + 2], WithinAbs(0.5, 1e-15));

  json bad = j;
  bad["boundaries"]["x_min"].erase("velocity");
  check_rejected(bad, "boundaries.x_min.velocity: missing");
  bad = j;
  bad["boundaries"]["x_min"]["type"] = "symmetry";
  check_rejected(bad, "boundaries.x_min.type");
  bad = j;
  bad["boundaries"]["x_min"]["speed"] = 1.0;
  check_rejected(bad, "boundaries.x_min.speed: unknown key");
  bad = j;
  bad["boundaries"]["x_max"] = 3;
  check_rejected(bad, "boundaries.x_max: expected a string or a piston object");
  bad = j;
  bad["mesh"]["perturbation"] = "random";
  check_rejected(bad, "mesh.perturbation");
}

TEST_CASE("parse_deck: the scheme block", "[io]") {
  using hydro::lagrangian::Limiter;
  const Deck d = parse_deck(small_deck().dump());
  CHECK(d.scheme.order == 1);
  CHECK(d.scheme.limiter == Limiter::SpV);

  json j = small_deck();
  j["scheme"] = {{"order", 2}, {"limiter", "none"}};
  const Deck e = parse_deck(j.dump());
  CHECK(e.scheme.order == 2);
  CHECK(e.scheme.limiter == Limiter::None);
  j["scheme"] = {{"order", 2}, {"limiter", "sp-v"}};
  CHECK(parse_deck(j.dump()).scheme.limiter == Limiter::SpV);
  j["scheme"] = json::object();
  CHECK(parse_deck(j.dump()).scheme.order == 1);

  j["scheme"] = {{"order", 3}};
  check_rejected(j, "scheme.order: must be 1 or 2");
  j["scheme"] = {{"order", 2.0}};
  check_rejected(j, "scheme.order: expected an integer");
  j["scheme"] = {{"limiter", "minmod"}};
  check_rejected(j, "scheme.limiter");
  j["scheme"] = {{"oder", 2}};
  check_rejected(j, "scheme.oder: unknown key");
  j["scheme"] = 2;
  check_rejected(j, "scheme: expected an object");
}

TEST_CASE("parse_deck: the Taylor-Green profile and source", "[io]") {
  using hydro::io::EnergySourceType;
  const json tg_initial = {{"taylor_green", {{"density", 2.0}, {"c1", 0.5}, {"c2", 3.0}}}};
  json j = small_deck();
  CHECK(!parse_deck(j.dump()).taylor_green);
  CHECK(parse_deck(j.dump()).source == EnergySourceType::None);

  j["initial"] = tg_initial;
  j["source"] = "taylor_green";
  const Deck d = parse_deck(j.dump());
  REQUIRE(d.taylor_green);
  CHECK(d.taylor_green->density == 2.0);
  CHECK(d.taylor_green->c1 == 0.5);
  CHECK(d.taylor_green->c2 == 3.0);
  CHECK(d.source == EnergySourceType::TaylorGreen);

  // Defaults (1, 1, 1), the thesis's constants; the source is optional.
  j["initial"] = {{"taylor_green", json::object()}};
  j.erase("source");
  const Deck e = parse_deck(j.dump());
  REQUIRE(e.taylor_green);
  CHECK(e.taylor_green->density == 1.0);
  CHECK(e.taylor_green->c1 == 1.0);
  CHECK(e.taylor_green->c2 == 1.0);
  CHECK(e.source == EnergySourceType::None);

  j["initial"] = {{"taylor_green", {{"c2", 0.5}}}};
  check_rejected(j, "initial.taylor_green.c2: must exceed density c1^2 / 2");
  j["initial"] = {{"taylor_green", {{"c3", 1.0}}}};
  check_rejected(j, "initial.taylor_green.c3: unknown key");
  j["initial"] = tg_initial;
  j["initial"]["background"] = {{"density", 1.0}, {"pressure", 1.0}};
  check_rejected(j, "initial.taylor_green: excludes background and regions");
  j = small_deck();
  j["source"] = "taylor_green";
  check_rejected(j, "source: \"taylor_green\" needs initial.taylor_green");
  j["initial"] = tg_initial;
  j["source"] = "gravity";
  check_rejected(j, "source: only \"taylor_green\"");
}

TEST_CASE("taylor_green: a steady solution of Euler with the energy source", "[app]") {
  // Appendix B of the thesis, checked by central differences at a few points:
  // momentum, rho (V . grad) V + grad P = 0 (V is divergence free), and total
  // energy, where the source must equal V . grad P / (gamma - 1) (Eq. B.12).
  const hydro::io::TaylorGreen tg{1.3, 0.7, 2.0};
  const double gamma = 5.0 / 3.0, h = 1e-5;
  using hydro::app::taylor_green_pressure;
  using hydro::app::taylor_green_velocity;
  for (const Vec3& x : {Vec3{0.1, 0.2, 0.0}, Vec3{0.63, 0.37, 0.5}, Vec3{0.9, 0.55, 1.0}}) {
    INFO("x = (" << x[0] << ", " << x[1] << ")");
    Vec3 grad_p{}, div_v{};
    std::array<Vec3, 2> dv{};  // dv[d] = dV/dx_d
    for (int d = 0; d < 2; ++d) {
      Vec3 a = x, b = x;
      a[d] += h;
      b[d] -= h;
      grad_p[d] = (taylor_green_pressure(tg, a) - taylor_green_pressure(tg, b)) / (2 * h);
      const Vec3 va = taylor_green_velocity(tg, a), vb = taylor_green_velocity(tg, b);
      for (int i = 0; i < 3; ++i) dv[d][i] = (va[i] - vb[i]) / (2 * h);
      div_v[0] += dv[d][d];
    }
    CHECK_THAT(div_v[0], WithinAbs(0.0, 1e-8));
    const Vec3 v = taylor_green_velocity(tg, x);
    for (int i = 0; i < 2; ++i) {
      const double advection = v[0] * dv[0][i] + v[1] * dv[1][i];
      CHECK_THAT(tg.density * advection + grad_p[i], WithinAbs(0.0, 1e-8));
    }
    const double v_grad_p = v[0] * grad_p[0] + v[1] * grad_p[1];
    CHECK_THAT(hydro::app::taylor_green_source(tg, gamma, x),
               WithinAbs(v_grad_p / (gamma - 1.0), 1e-8));
  }
  CHECK_THAT(taylor_green_pressure(tg, {0.0, 0.0, 0.0}),
             WithinRel(0.25 * 1.3 * 0.49 * 2.0 + 2.0, 1e-15));
}

TEST_CASE("initial_state: the Taylor-Green profile at the cell centroids", "[app]") {
  json j = read_deck_json("taylor_green");
  j["mesh"]["cells"] = {4, 3, 1};
  const Deck deck = parse_deck(j.dump());
  REQUIRE(deck.taylor_green);
  const HydroState s = hydro::app::initial_state(deck);
  const auto centers = hydro::lagrangian::cell_centers(s.mesh);
  REQUIRE(s.mesh.num_cells() == 12);
  for (Index c = 0; c < 12; ++c) {
    INFO("cell " << c);
    const Vec3 v = hydro::app::taylor_green_velocity(*deck.taylor_green, centers[c]);
    const double p = hydro::app::taylor_green_pressure(*deck.taylor_green, centers[c]);
    CHECK_THAT(s.mass[c], WithinRel(1.0 / 12.0, 1e-14));
    CHECK(s.vel_x[c] == v[0]);
    CHECK(s.vel_y[c] == v[1]);
    CHECK(s.vel_z[c] == 0.0);
    CHECK_THAT(s.total_energy[c], WithinRel(p / 0.4 + 0.5 * (v[0] * v[0] + v[1] * v[1]), 1e-14));
  }
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

TEST_CASE("hydro_run: Sod at second order stays 1D", "[app][integration]") {
  // decks/sod.json with scheme.order = 2 (predictor-corrector, SP-V). The
  // 100x2x2 slab's flow must stay exactly 1D: no transverse velocity and the
  // same state in the 4 cells of every column. SP-V is built for this, but
  // only if round-off V_c count as zero in its projection (limiter.hpp):
  // with exact-zero tests the columns drifted 1e-4 apart. Measured: 212
  // steps, L1 density error 0.0072 (first order: 0.0195).
  json deck = read_deck_json("sod");
  deck["scheme"] = {{"order", 2}};
  const auto deck_path =
      (std::filesystem::temp_directory_path() / "hydro_test_sod_order2.json").string();
  {
    std::ofstream out(deck_path);
    out << deck.dump();
  }
  const json dump = run_hydro_deck(deck_path, "sod_order2");
  std::filesystem::remove(deck_path);
  CHECK(dump["t"].get<double>() == 0.2);
  const auto& cells = dump["cells"];
  const auto rho = cells["density"].get<std::vector<double>>();
  const auto vx = cells["velocity_x"].get<std::vector<double>>();
  const auto vy = cells["velocity_y"].get<std::vector<double>>();
  const auto vz = cells["velocity_z"].get<std::vector<double>>();
  REQUIRE(rho.size() == 400);
  for (std::size_t c = 0; c < rho.size(); ++c) {
    INFO("cell " << c);
    CHECK_THAT(vy[c], WithinAbs(0.0, 1e-12));
    CHECK_THAT(vz[c], WithinAbs(0.0, 1e-12));
    // Cell (i, j, k) is c = i + 100 (j + 2 k); compare with (i, 0, 0).
    CHECK_THAT(rho[c], WithinRel(rho[c % 100], 1e-12));
    CHECK_THAT(vx[c], WithinAbs(vx[c % 100], 1e-12));
  }
  const double l1 = l1_density_error(dump, read_oracle("sod"), 0.02 * 0.02);
  INFO("L1 density error " << l1);
  CHECK(l1 < 0.008);
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

TEST_CASE("hydro_run: Saltzman piston against the exact solution", "[app][integration]") {
  // decks/saltzman.json: cold gas (gamma = 5/3, rho = 1) at rest on Saltzman's
  // skewed 100x10x1 mesh of [0,1]x[0,0.1]x[0,0.01], pushed by a piston at
  // u = 1 from x = 0, walls elsewhere, dt^0 = 1e-4, compared with
  // tests/oracle_data/saltzman.json: piston at x = 0.6, shock at x = 0.8,
  // rho = 4, u = 1, eps = 1/2 between them. First order on this mesh gives a
  // shock front tilted along the initial skew and a density overshoot where
  // it meets the top wall, so the checks are on averages and bounds.
  // Measured: 1052 steps (2 min in a Debug build, 11 s optimized); shock in
  // rows 0-9 between x = 0.789 and 0.811; on [0.65, 0.75] rho 3.87-4.22 per
  // cell, mass-weighted means rho 4.009, u 0.998, eps 0.504; total energy
  // 7.907e-4 vs the exact piston work 8e-4 (-1.2%, start-up heating at the
  // piston); gas ahead of x = 0.85 untouched (|u| < 5e-12); L1 density
  // error 0.0454.
  const json dump = run_hydro("saltzman");
  const json oracle = read_oracle("saltzman");
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
  const std::size_t nx = 100, ny = 10;
  REQUIRE(rho.size() == nx * ny);

  // The piston does work P u A t = (4/3) t (0.1 x 0.01) on gas that starts
  // with none.
  double e = 0.0;
  for (std::size_t c = 0; c < rho.size(); ++c) {
    e += mass[c] * (eps[c] + 0.5 * (u[c] * u[c] + vy[c] * vy[c] + vz[c] * vz[c]));
  }
  CHECK_THAT(e, WithinRel(4.0 / 3.0 * 0.6 * 0.1 * 0.01, 0.02));

  double m_sum = 0.0, rho_sum = 0.0, u_sum = 0.0, eps_sum = 0.0;
  for (std::size_t c = 0; c < rho.size(); ++c) {
    INFO("cell " << c << " at x = " << x[c]);
    CHECK_THAT(vz[c], WithinAbs(0.0, 1e-12));
    if (x[c] > 0.65 && x[c] < 0.75) {
      CHECK_THAT(rho[c], WithinRel(4.0, 0.07));
      m_sum += mass[c];
      rho_sum += mass[c] * rho[c];
      u_sum += mass[c] * u[c];
      eps_sum += mass[c] * eps[c];
    }
    if (x[c] > 0.85) {
      CHECK_THAT(rho[c], WithinRel(1.0, 1e-5));
      CHECK_THAT(u[c], WithinAbs(0.0, 1e-10));
      CHECK_THAT(vy[c], WithinAbs(0.0, 1e-10));
    }
  }
  CHECK_THAT(rho_sum / m_sum, WithinRel(4.0, 0.01));
  CHECK_THAT(u_sum / m_sum, WithinRel(1.0, 0.01));
  CHECK_THAT(eps_sum / m_sum, WithinRel(0.5, 0.02));

  // Shock position in each row (cells j*nx + i): the last cell with rho > 2.5
  // and the next one bracket it, within 2 initial cells of x = 0.8.
  for (std::size_t j = 0; j < ny; ++j) {
    std::size_t last = 0;
    for (std::size_t i = 0; i < nx; ++i) {
      if (rho[j * nx + i] > 2.5) last = i;
    }
    INFO("row " << j);
    REQUIRE(last + 1 < nx);
    CHECK(x[j * nx + last] > 0.78);
    CHECK(x[j * nx + last + 1] < 0.82);
  }

  const double l1 = l1_density_error(dump, oracle, 0.1 * 0.01);
  INFO("L1 density error " << l1);
  CHECK(l1 < 0.05);
}

namespace {

// Pressure errors of a Taylor-Green run against the exact (initial)
// pressure at each cell's centroid, thesis Eq. 3.7: L1, L2, Linf.
std::array<double, 3> taylor_green_errors(const Deck& deck, const HydroState& s) {
  using namespace hydro::lagrangian;
  const auto faces = cell_face_vectors(s.mesh);
  const CellThermo thermo = cell_thermo(s, faces, hydro::eos::IdealGasEOS(deck.gamma));
  const auto centers = cell_centers(s.mesh);
  double volume = 0.0, l1 = 0.0, l2 = 0.0, linf = 0.0;
  for (Index c = 0; c < s.mesh.num_cells(); ++c) {
    const double v = cell_volume(s.mesh, c, faces[c]);
    const double dp = std::abs(thermo.pressure[c] -
                               hydro::app::taylor_green_pressure(*deck.taylor_green, centers[c]));
    volume += v;
    l1 += v * dp;
    l2 += v * dp * dp;
    linf = std::max(linf, dp);
  }
  return {l1 / volume, std::sqrt(l2 / volume), linf};
}

}  // namespace

TEST_CASE("hydro_run: Taylor-Green vortex converges at order ~2", "[app][integration]") {
  // decks/taylor_green.json (thesis Sec. 3.2) at 10x10, 20x20 and 40x40, to
  // t = 0.5. The thesis's Table 3.1(b) (predictor-corrector, Green
  // gradients, unlimited), L1 / L2 pressure errors: 3.72e-2 / 4.96e-2,
  // 2.11e-2 / 2.37e-2 (its "2.11e-3" is a typo: the order it gives, 0.82,
  // needs 2.11e-2), 8.29e-3 / 9.61e-3, then 2.39e-3 / 2.83e-3 at 80x80.
  // Measured at order 2: 2.81e-2 / 3.98e-2, 2.03e-2 / 2.28e-2, 8.08e-3 /
  // 9.44e-3 (2.31e-3 / 2.76e-3 at 80x80, L1 order 1.81 vs the thesis's
  // 1.79; 160x160 and up via scripts/taylor_green_convergence.py). At order
  // 1: L1 1.37e-1, 9.09e-2, 5.29e-2 (orders 0.59, 0.78; 0.89 at 80x80).
  // About 20 s at order 2, 6 s at order 1, optimized.
  const json base = read_deck_json("taylor_green");
  const std::array<int, 3> sizes = {10, 20, 40};
  std::array<std::array<std::array<double, 3>, 3>, 2> errors{};  // [order][mesh][norm]
  for (const int order : {1, 2}) {
    for (std::size_t m = 0; m < sizes.size(); ++m) {
      json j = base;
      j["mesh"]["cells"] = {sizes[m], sizes[m], 1};
      j["scheme"]["order"] = order;
      const Deck deck = parse_deck(j.dump());
      std::ostringstream log;
      const auto r = hydro::app::run(deck, log);
      CHECK(r.t == 0.5);
      errors[order - 1][m] = taylor_green_errors(deck, r.state);
    }
  }
  auto rate = [&](int order, std::size_t m, int norm) {
    return std::log2(errors[order - 1][m - 1][norm] / errors[order - 1][m][norm]);
  };
  for (int norm = 0; norm < 3; ++norm) {
    INFO("norm " << norm << ": order 2 errors " << errors[1][0][norm] << ", "
         << errors[1][1][norm] << ", " << errors[1][2][norm] << "; order 1 "
         << errors[0][0][norm] << ", " << errors[0][1][norm] << ", " << errors[0][2][norm]);
    for (std::size_t m = 0; m < sizes.size(); ++m) {
      CHECK(errors[1][m][norm] < errors[0][m][norm]);
      if (m > 0) {
        CHECK(errors[1][m][norm] < errors[1][m - 1][norm]);
        CHECK(errors[0][m][norm] < errors[0][m - 1][norm]);
      }
    }
  }
  // L1 and L2: the order-2 rate rises towards 2 and passes the order-1 one;
  // within a factor 1.5 of the thesis's errors.
  for (int norm = 0; norm < 2; ++norm) {
    INFO("norm " << norm << ": order 2 rates " << rate(2, 1, norm) << ", "
         << rate(2, 2, norm) << "; order 1 " << rate(1, 1, norm) << ", " << rate(1, 2, norm));
    CHECK(rate(2, 2, norm) > rate(2, 1, norm));
    CHECK(rate(2, 2, norm) > 1.2);
    CHECK(rate(2, 2, norm) > rate(1, 2, norm) + 0.3);
    CHECK(rate(1, 2, norm) > 0.6);
    CHECK(rate(1, 2, norm) < 1.2);
  }
  const std::array<std::array<double, 3>, 2> thesis = {{{3.72e-2, 2.11e-2, 8.29e-3},
                                                         {4.96e-2, 2.37e-2, 9.61e-3}}};
  for (int norm = 0; norm < 2; ++norm) {
    for (std::size_t m = 0; m < sizes.size(); ++m) {
      INFO("norm " << norm << ", " << sizes[m] << "x" << sizes[m]);
      CHECK(errors[1][m][norm] < 1.5 * thesis[norm][m]);
      CHECK(errors[1][m][norm] > thesis[norm][m] / 1.5);
    }
  }
}
