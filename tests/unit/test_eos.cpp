#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

#include "eos.hpp"

using Catch::Matchers::WithinRel;
using hydro::eos::IdealGasEOS;

TEST_CASE("IdealGasEOS matches hand-computed pressure and sound speed", "[eos]") {
  const IdealGasEOS eos(1.4);
  const double rho = 1.0;
  const double eps = 2.5;

  // p = (gamma - 1) * rho * eps = 0.4 * 1.0 * 2.5
  const double p_expected = 1.0;
  CHECK_THAT(eos.pressure(rho, eps), WithinRel(p_expected, 1e-14));

  // a = sqrt(gamma * p / rho) = sqrt(1.4 * 1.0 / 1.0) = sqrt(1.4)
  const double a_expected = std::sqrt(1.4);
  CHECK_THAT(eos.sound_speed(rho, eps), WithinRel(a_expected, 1e-14));
}
