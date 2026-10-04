// hydro_run <deck.json> [dump.json]: runs a deck to its t_end and writes the
// final state to the dump path (the argument, else the deck's output.path).
#include <exception>
#include <iostream>
#include <string>

#include "deck.hpp"
#include "driver.hpp"
#include "dump.hpp"
#include "eos.hpp"

int main(int argc, char** argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "usage: hydro_run <deck.json> [dump.json]\n";
    return 2;
  }
  try {
    const hydro::io::Deck deck = hydro::io::read_deck(argv[1]);
    for (const auto& w : deck.warnings) std::cerr << "warning: " << w << '\n';
    const std::string output = argc == 3 ? argv[2] : deck.output;
    if (output.empty()) {
      std::cerr << "hydro_run: no dump path (give one, or set output.path in the deck)\n";
      return 2;
    }
    const auto result = hydro::app::run(deck, std::cout);
    hydro::io::write_dump(output, result.state, hydro::eos::IdealGasEOS(deck.gamma),
                          result.t, result.steps);
    std::cout << "wrote " << output << '\n';
  } catch (const std::exception& e) {
    std::cerr << "hydro_run: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
