#pragma once

#include "constants.hpp"
#include "dopr5.hpp"
#include <array>
#include <cmath>

struct Geodesic {
  using YVector = std::array<double, 3>;

  double rs, R_star;
  StepperDopr5<3, Geodesic> stepper;

  Geodesic(double M_star, double R_star, double orbit_atol = 1e-10,
           double orbit_rtol = 1e-10)
      : rs(M_star * schwarzschild_radius_of_sun_in_km), R_star(R_star),
        stepper(*this, orbit_atol, orbit_rtol) {
    stepper.add_event(
        [R_star](double, YVector const& y) { return y[0] - R_star; });
  }

  void
  init(double r0, double psi0, double alpha0) {
    stepper.init(0.0, 1e-3 * R_star, {r0, psi0, alpha0});
  }

  double
  l_old() const {
    return stepper.x_old;
  }
  double
  l_new() const {
    return stepper.x_old + stepper.h_old;
  }

  bool
  step() {
    stepper.do_step(0.1 * stepper.y_old[0]);
    return stepper.detect_event() == 0;
  }

  YVector
  interp(double l) const {
    return stepper.dense_out(l);
  }

  void
  operator()(double, YVector const& y, YVector& dydx) const {
    auto [r, psi, alpha] = y;
    double L = std::sqrt(1 - rs / r);
    dydx[0] = L * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * L) * (1 - 3 * rs / (2 * r));
  }
};
