#pragma once

#include <string>

#include "eos.hpp"
#include "hydro.hpp"

namespace hydro::io {

// Writes the state at time t, after `steps` steps, as JSON (format
// "hydro-dump-0", described in docs/deck_format.md): per cell the centroid
// (mean of its nodes), volume, mass, density, velocity, pressure and specific
// internal energy, plus the node coordinates. A minimal snapshot for
// plotting and comparison; the restart schema is deferred to v1.0. Throws
// std::runtime_error if the file cannot be written.
void write_dump(const std::string& path, const lagrangian::HydroState& s,
                const eos::EOS& eos, double t, long steps);

}  // namespace hydro::io
