#pragma once

#include "constants.hpp"
#include "doubles.hpp"
#include "ran.hpp"
#include <cmath>
#include <stdexcept>

enum class Polarization { O, E };

struct InitializedPhoton {
  double3 n, e1, e2; // fixed, right-handed orbital basis: e1 x e2 = n
  double omega_inf;  // photon energy hbar * omega at infinity, in keV
  // Radius in km; orbital and propagation angles in radians.
  double r0, psi0, alpha0;
  Polarization pol;
};

struct InitialPhotonOptions {
  // Equal-radius caps about the magnetic +/-z axes, uniform per unit area.
  // pi/2 with north_fraction = 0.5 gives uniform full-surface emission.
  // FD11's 5-degree caps use 5*pi/180; north_fraction = 0 selects the south.
  double cap_half_angle = pi / 2; // radians, in (0, pi/2]
  double north_fraction = 0.5;    // probability of emission from the north cap
  double e_mode_fraction = 1.0;   // probability of E mode; 0.5 is unpolarized
};

namespace initialize_fd11_detail {

inline
double
open_unit(Ran& ran) {
  // Ran::U() can return zero, or round its integer conversion to one.
  double u;
  do {
    u = ran.U();
  } while (u <= 0 || u >= 1);
  return u;
}

inline
double
sample_planck_number(Ran& ran) {
  // FD11 Sections 2.2 and 3.4: sample photon NUMBER, not energy, so
  // p(x) = x^2 / [2*zeta(3)*(exp(x)-1)], x = omega_inf / T_inf.
  // Expanding 1/(exp(x)-1) gives Gamma(shape=3, rate=m) components
  // with probabilities 1/[zeta(3)*m^3], m = 1, 2, ... .
  double m, acceptance;
  do {
    m = std::floor(1 / std::sqrt(open_unit(ran)));
    // Proposal q(m) = 1/m^2 - 1/(m+1)^2; q(m)*acceptance = 3/(4*m^3).
    // Use a double for m: Ran's smallest positive U is 2^-64, so m <= 2^32.
    acceptance = 3 * (m + 1) * (m + 1) / (4 * m * (2 * m + 1));
  } while (open_unit(ran) >= acceptance);

  double x = 0;
  for (int i = 0; i < 3; ++i)
    x -= std::log(open_unit(ran));
  return x / m;
}

} // namespace initialize_fd11_detail

// One equal-photon-number sample. Both caps share T_inf (k_B*T in keV)
// and cap_half_angle; north_fraction sets their relative photon luminosity.
// Defaults reproduce FD11's uniform-temperature, full-surface, pure E model.
// rs is the Schwarzschild radius and R_star the stellar radius, both in km.
inline
InitializedPhoton
sample_initial_photon_fd11(Ran& ran, double rs, double R_star, double T_inf,
                           InitialPhotonOptions const& options = {}) {
  if (!std::isfinite(R_star) || R_star <= 0)
    throw std::invalid_argument("R_star must be positive and finite");
  if (!std::isfinite(rs) || rs < 0 || rs >= R_star)
    throw std::invalid_argument(
        "rs must be finite and satisfy 0 <= rs < R_star");
  if (!std::isfinite(T_inf) || T_inf <= 0)
    throw std::invalid_argument("T_inf must be positive and finite");
  if (!std::isfinite(options.cap_half_angle) || options.cap_half_angle <= 0 ||
      options.cap_half_angle > pi / 2)
    throw std::invalid_argument(
        "cap_half_angle must be finite and in (0, pi/2]");
  if (!std::isfinite(options.north_fraction) || options.north_fraction < 0 ||
      options.north_fraction > 1)
    throw std::invalid_argument("north_fraction must be finite and in [0, 1]");
  if (!std::isfinite(options.e_mode_fraction) || options.e_mode_fraction < 0 ||
      options.e_mode_fraction > 1)
    throw std::invalid_argument("e_mode_fraction must be finite and in [0, 1]");

  using initialize_fd11_detail::open_unit;
  // The Fortran seed routine blueshifts its LOCAL energy. Here we store the
  // conserved energy at infinity; propagation uses omega_inf/sqrt(1-rs/r).
  double omega_inf = T_inf * initialize_fd11_detail::sample_planck_number(ran);
  if (!std::isfinite(omega_inf) || omega_inf <= 0)
    throw std::runtime_error(
        "Initial photon energy is not finite and positive");

  double pole = open_unit(ran) < options.north_fraction ? 1 : -1;
  double sin_half_cap = std::sin(options.cap_half_angle / 2);
  // d = 1 - |cos(theta)|, avoiding subtraction at small cap radii.
  double d = 2 * sin_half_cap * sin_half_cap * open_unit(ran);
  double rho = std::sqrt(d * (2 - d));
  double c, s;
  ran.point_on_unit_circle(c, s);
  double3 e1{rho * c, rho * s, pole * (1 - d)};
  double3 n = ran.unit_perp_to(e1);
  double3 e2 = cross(n, e1);

  // Isotropic surface intensity gives outgoing flux p(mu) = 2*mu, mu in [0,1].
  // This is mustar = sqrt(RAND) in magmc_2019/src/physics.F90, seed().
  // Choosing n independently keeps the basis regular even for radial rays.
  double alpha0 = std::acos(std::sqrt(open_unit(ran)));
  Polarization pol = open_unit(ran) < options.e_mode_fraction ? Polarization::E
                                                              : Polarization::O;
  return {n, e1, e2, omega_inf, R_star, 0, alpha0, pol};
}
