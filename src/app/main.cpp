// hydro_run <deck.json> [dump.json] [--probe cell_id]... [--probe-node node_id]...
//           [--probe-out path]: runs a deck to its t_end and writes the final
// state to the dump path (the argument, else the deck's output.path). Probed
// cells and nodes are traced as JSON lines (app::Probes) to --probe-out, else
// to stderr.
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "deck.hpp"
#include "driver.hpp"
#include "dump.hpp"
#include "eos.hpp"

namespace {

constexpr const char* kUsage =
    "usage: hydro_run <deck.json> [dump.json] [--probe cell_id]... "
    "[--probe-node node_id]... [--probe-out path]\n";

// A non-negative id, the whole string; -1 if it is not one.
hydro::mesh::Index parse_id(const std::string& s) {
  if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos) return -1;
  try {
    return static_cast<hydro::mesh::Index>(std::stoll(s));
  } catch (const std::exception&) {
    return -1;
  }
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> positional;
  hydro::app::Probes probes;
  std::string probe_out;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--probe" || arg == "--probe-node" || arg == "--probe-out") {
      if (i + 1 == argc) {
        std::cerr << "hydro_run: " << arg << " needs a value\n" << kUsage;
        return 2;
      }
      const std::string value = argv[++i];
      if (arg == "--probe-out") {
        probe_out = value;
        continue;
      }
      const hydro::mesh::Index id = parse_id(value);
      if (id < 0) {
        std::cerr << "hydro_run: " << arg << " " << value << ": not an id\n" << kUsage;
        return 2;
      }
      (arg == "--probe" ? probes.cells : probes.nodes).push_back(id);
    } else if (arg.rfind("--", 0) == 0) {
      std::cerr << "hydro_run: unknown option " << arg << '\n' << kUsage;
      return 2;
    } else {
      positional.push_back(arg);
    }
  }
  if (positional.empty() || positional.size() > 2) {
    std::cerr << kUsage;
    return 2;
  }
  try {
    const hydro::io::Deck deck = hydro::io::read_deck(positional[0]);
    for (const auto& w : deck.warnings) std::cerr << "warning: " << w << '\n';
    const std::string output = positional.size() == 2 ? positional[1] : deck.output;
    if (output.empty()) {
      std::cerr << "hydro_run: no dump path (give one, or set output.path in the deck)\n";
      return 2;
    }
    std::ofstream probe_file;
    probes.out = &std::cerr;
    if (!probe_out.empty()) {
      probe_file.open(probe_out);
      if (!probe_file) {
        std::cerr << "hydro_run: cannot open " << probe_out << " for writing\n";
        return 1;
      }
      probes.out = &probe_file;
    }
    const auto result = hydro::app::run(deck, std::cout, probes);
    hydro::io::write_dump(output, result.state, hydro::eos::IdealGasEOS(deck.gamma),
                          result.t, result.steps);
    std::cout << "wrote " << output << '\n';
  } catch (const std::exception& e) {
    std::cerr << "hydro_run: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
