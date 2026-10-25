#pragma once

#include "bfield.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "magnus.hpp"
#include "polarization.hpp"
#include "quad.hpp"
#include "ran.hpp"
#include <algorithm>
#include <cmath>

template <PhotonInitializer Init, typename Scat>
struct PhotonEvolution {
  struct Options {
    double eta_couple = 1e-3, escape_tolerance = 1e-3;
  } options;

  BField const& bfield;
  Init const& initializer;
  Scat const& resonant; // The scattering model must outlive this evolution.
  Ran ran;
  Geodesic geodesic;
  Quad quad;

  double3 n, e1, e2;
  double omega_inf, tau;
  Polarization pol;

  PhotonEvolution(double M_star, double R_star, int seed, BField const& bfield,
                  Init const& init, Scat const& scat,
                  Options options = Options{})
      : options(options), bfield(bfield), initializer(init), resonant(scat),
        ran(seed), geodesic(M_star, R_star) {}

  struct Basis {
    double3 r_hat, psi_hat, k_hat, s1, s2, B, omega;
    double Bperp2, r, c, s, kappa, Bmag, mu, twist, J;
    Basis(double radius, double psi, double alpha, double3 n, double3 e1,
          double3 e2, BField const& field, double energy, double rs) {
      r = radius;
      r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
      psi_hat = cross(n, r_hat);
      k_hat = std::cos(alpha) * r_hat + std::sin(alpha) * psi_hat;
      s1 = cross(n, k_hat);
      s2 = n;
      double muz = std::clamp(r_hat.z, -1.0, 1.0);
      double rho = std::hypot(r_hat.x, r_hat.y);
      double3 theta_hat =
          rho > 0 ? double3{r_hat.x * muz / rho, r_hat.y * muz / rho, -rho}
                  : double3{std::copysign(1.0, muz), 0, 0};
      double3 phi_hat = rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0}
                                : double3{0, 1, 0};
      auto field_result = field.evaluate(r, muz);
      double3 B_vec = field_result.B;
      twist = field_result.twist;
      J = field_result.J;
      Bmag = field_result.Bmag;
      B = B_vec.x * r_hat + B_vec.y * theta_hat + B_vec.z * phi_hat;
      mu = Bmag > 0 ? std::clamp(dot(k_hat, B) / Bmag, -1.0, 1.0) : 0;
      double L = std::sqrt(1 - rs / r);
      double B1 = dot(B, s1), B2 = dot(B, s2);
      Bperp2 = B1 * B1 + B2 * B2;
      double diff = B1 * B1 - B2 * B2, twice = 2 * B1 * B2;
      c = Bperp2 > 0 ? diff / Bperp2 : 1;
      s = Bperp2 > 0 ? twice / Bperp2 : 0;
      double K = 2.013919e-19 * energy / L;
      kappa = K * Bperp2;
      omega = {K * diff, K * twice, 0};
    }
  };

  Basis
  basis(Geodesic::YVector y) const {
    auto [r, psi, alpha] = y;
    return Basis(r, psi, alpha, n, e1, e2, bfield, omega_inf, geodesic.rs);
  }

  Basis
  get_basis(double l) const {
    auto y = geodesic.interp(l);
    return basis(y);
  }

  void
  init() {
    init(initializer.sample_initial_photon(ran));
  }
  void
  init(InitializedPhoton const& photon) {
    geodesic.init(photon.r0, photon.psi0, photon.alpha0);
    n = photon.n;
    e1 = photon.e1;
    e2 = photon.e2;
    omega_inf = photon.omega_inf;
    tau = std::log(ran.U());
    pol.set_adiabatic(photon.pol);
  }

  template <class PolarizationAt>
  bool
  integrate_optical_depth(double left, double right,
                          PolarizationAt const& pol_at) {
    auto event = resonant.integrate_optical_depth(
        left, right, [&](double l) { return get_basis(l); }, pol_at, omega_inf,
        geodesic.rs, tau, ran, quad);
    if (!event)
      return false;
    scattered_data = {geodesic.interp(event->l), event->beta, event->icomp,
                      event->iroot};
    pol = event->incoming;
    return true;
  }

  void
  record_escape(double l, Geodesic::YVector y) {
    auto [r, psi, alpha] = y;
    double impact = r * std::sin(alpha) / std::sqrt(1 - geodesic.rs / r);
    // u = r/rho maps [r,infinity) to [1,0], with a regular integrand at 0.
    auto integrand = [&](double u) {
      double t = impact / r;
      return t / std::sqrt(1 - t * t * u * u * (1 - geodesic.rs * u / r));
    };
    double psi_inf = psi + quad.qags(integrand, 0, 1, 1e-11, 1e-11);
    double3 k = to_unit(std::cos(psi_inf) * e1 + std::sin(psi_inf) * e2);
    Stokes sky = pol.observer_stokes(k, n);
    escaped_data = {.omega_inf = omega_inf,
                    .muk = std::clamp(k.z, -1.0, 1.0),
                    .Q = sky.x,
                    .U = sky.y,
                    .V = sky.z};
  }

  enum class EvolveResult { Absorbed, Escaped, Scattered };

  struct EscapedData {
    double omega_inf, muk;
    double Q, U, V;
  } escaped_data;

  struct ScatteredData {
    Geodesic::YVector r_psi_alpha;
    double beta;
    std::size_t icomp;
    int iroot;
  } scattered_data;

  EvolveResult
  evolve() {
    while (true) {
      bool absorbed = geodesic.step();
      if (pol.is_adiabatic) {
        auto f_couple = [&](double l) {
          Basis basis = get_basis(l);
          return options.eta_couple * basis.kappa * basis.r - 1;
        };
        if (f_couple(geodesic.l_old()) < 0) {
          Basis b = get_basis(geodesic.l_old());
          pol.set_stokes(b.c, b.s);
        }
      }
      double left = geodesic.l_old(), right = geodesic.l_new();
      if (pol.is_adiabatic) {
        if (integrate_optical_depth(left, right, [&](double) { return pol; }))
          return EvolveResult::Scattered;
      } else {
        StepperMagnus stepper;
        stepper.init(geodesic.l_old(), geodesic.l_new(),
                     0.1 * (geodesic.l_new() - geodesic.l_old()), pol.stokes);
        auto omega = [&](double l) { return get_basis(l).omega; };
        while (stepper.do_step(omega)) {
          if (integrate_optical_depth(
                  stepper.x_old, stepper.x_new,
                  [&](double l) { return Polarization(stepper.dense_out(l)); }))
            return EvolveResult::Scattered;

          pol.set_stokes(stepper.y_new);
          auto y = geodesic.interp(stepper.x_new);
          double change = (stepper.y_new + (-1.0) * stepper.y_old).length() *
                          y[0] / (stepper.x_new - stepper.x_old);
          if (change < options.escape_tolerance &&
                  resonant.resonance_escape(
                      stepper.x_new, geodesic,
                      [&](Geodesic::YVector y) { return basis(y); },
                      omega_inf) ||
              (y[0] > 1e4 * geodesic.R_star && std::cos(y[2]) > 0)) {
            record_escape(stepper.x_new, y);
            return EvolveResult::Escaped;
          }
          stepper.update_old();
        }
      }
      if (absorbed || geodesic.stepper.y_new[0] <= geodesic.R_star)
        return EvolveResult::Absorbed;
      geodesic.stepper.update_old();
    }
  }

  void
  perform_scattering() {
    auto b = basis(scattered_data.r_psi_alpha);
    init(resonant.perform_scattering(b, scattered_data.beta, omega_inf, ran));
  }
};
