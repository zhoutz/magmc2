#pragma once

#include "constants.hpp"
#include "doubles.hpp"
#include "polarization.hpp"
#include "ran.hpp"
#include <cmath>
#include <concepts>
#include <stdexcept>

struct InitializedPhoton {
  double3 n, e1, e2; // fixed, right-handed orbital basis: e1 x e2 = n
  double energy_inf; // photon energy hbar * omega at infinity, in keV
  // Radius in km; orbital and propagation angles in radians.
  double r0, psi0, alpha0;
  AdiabaticPolarization pol;
};

template <typename T>
concept PhotonInitializer = requires(T const& t, Ran& ran) {
  { t.sample_initial_photon(ran) } -> std::same_as<InitializedPhoton>;
};

struct InitFD11 {
  struct Options {
    // Equal-radius caps about the magnetic +/-z axes, uniform per unit area.
    // pi/2 with north_fraction = 0.5 gives uniform full-surface emission.
    // FD11's 5-degree caps use 5*pi/180; north_fraction = 0 selects the south.
    double cap_half_angle = pi / 2; // radians, in (0, pi/2]
    double north_fraction = 0.5;  // probability of emission from the north cap
    double e_mode_fraction = 1.0; // probability of E mode; 0.5 is unpolarized
  };

  double R_star, T_inf;
  Options options;

  InitFD11(double R_star, double T_inf) : InitFD11(R_star, T_inf, Options{}) {}

  InitFD11(double R_star, double T_inf, Options options)
      : R_star(R_star), T_inf(T_inf), options(options) {
    if (!std::isfinite(R_star) || R_star <= 0)
      throw std::invalid_argument("R_star must be positive and finite");
    if (!std::isfinite(T_inf) || T_inf <= 0)
      throw std::invalid_argument("T_inf must be positive and finite");
    if (!std::isfinite(options.cap_half_angle) || options.cap_half_angle <= 0 ||
        options.cap_half_angle > pi / 2)
      throw std::invalid_argument(
          "cap_half_angle must be finite and in (0, pi/2]");
    if (!std::isfinite(options.north_fraction) || options.north_fraction < 0 ||
        options.north_fraction > 1)
      throw std::invalid_argument(
          "north_fraction must be finite and in [0, 1]");
    if (!std::isfinite(options.e_mode_fraction) ||
        options.e_mode_fraction < 0 || options.e_mode_fraction > 1)
      throw std::invalid_argument(
          "e_mode_fraction must be finite and in [0, 1]");
  }

  double
  sample_planck_number(Ran& ran) const {
    // FD11 Sections 2.2 and 3.4: sample photon NUMBER, not energy, so
    // p(x) = x^2 / [2*zeta(3)*(exp(x)-1)], x = energy_inf / T_inf.
    // Expanding 1/(exp(x)-1) gives Gamma(shape=3, rate=m) components
    // with probabilities 1/[zeta(3)*m^3], m = 1, 2, ... .
    double m, acceptance;
    do {
      m = std::floor(1 / std::sqrt(ran.U()));
      // Proposal q(m) = 1/m^2 - 1/(m+1)^2; q(m)*acceptance = 3/(4*m^3).
      // Use a double for m: Ran's smallest positive U is 2^-64, so m <= 2^32.
      acceptance = 3 * (m + 1) * (m + 1) / (4 * m * (2 * m + 1));
    } while (ran.U() >= acceptance);

    double x = 0;
    for (int i = 0; i < 3; ++i)
      x -= std::log(ran.U());
    return x / m;
  }

  // One equal-photon-number sample. Both caps share T_inf (k_B*T in keV)
  // and cap_half_angle; north_fraction sets their relative photon luminosity.
  // Defaults reproduce FD11's uniform-temperature, full-surface, pure E model.
  // rs is the Schwarzschild radius and R_star the stellar radius, both in km.
  InitializedPhoton
  sample_initial_photon(Ran& ran) const {
    // The Fortran seed routine blueshifts its LOCAL energy. Here we store the
    // conserved energy at infinity; propagation uses energy_inf/sqrt(1-rs/r).
    double energy_inf = T_inf * sample_planck_number(ran);
    if (!std::isfinite(energy_inf) || energy_inf <= 0)
      throw std::runtime_error(
          "Initial photon energy is not finite and positive");

    double pole = ran.U() < options.north_fraction ? 1 : -1;
    double sin_half_cap = std::sin(options.cap_half_angle / 2);
    // d = 1 - |cos(theta)|, avoiding subtraction at small cap radii.
    double d = 2 * sin_half_cap * sin_half_cap * ran.U();
    double rho = std::sqrt(d * (2 - d));
    double c, s;
    ran.point_on_unit_circle(c, s);
    double3 e1{rho * c, rho * s, pole * (1 - d)};
    double3 n = ran.unit_perp_to(e1);
    double3 e2 = cross(n, e1);

    // Isotropic surface intensity gives outgoing flux p(mu) = 2*mu, mu in
    // [0,1]. This is mustar = sqrt(RAND) in magmc_2019/src/physics.F90, seed().
    // Choosing n independently keeps the basis regular even for radial rays.
    double alpha0 = std::acos(std::sqrt(ran.U()));
    AdiabaticPolarization pol = ran.U() < options.e_mode_fraction
                                    ? AdiabaticPolarization::E
                                    : AdiabaticPolarization::O;
    return {n, e1, e2, energy_inf, R_star, 0, alpha0, pol};
  }
};
