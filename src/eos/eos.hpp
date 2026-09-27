#pragma once

#include <cmath>

namespace hydro::eos {

// Abstract equation of state: closes the Euler system given density and
// specific internal energy.
class EOS {
public:
  virtual ~EOS() = default;

  virtual double pressure(double rho, double eps) const = 0;
  virtual double sound_speed(double rho, double eps) const = 0;
};

// Ideal gas: p = (gamma - 1) * rho * eps, a = sqrt(gamma * p / rho).
class IdealGasEOS final : public EOS {
public:
  explicit IdealGasEOS(double gamma) : gamma_(gamma) {}

  double gamma() const { return gamma_; }

  double pressure(double rho, double eps) const override {
    return (gamma_ - 1.0) * rho * eps;
  }

  double sound_speed(double rho, double eps) const override {
    return std::sqrt(gamma_ * pressure(rho, eps) / rho);
  }

private:
  double gamma_;
};

}  // namespace hydro::eos
