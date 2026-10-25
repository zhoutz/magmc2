#pragma once

#include "bfield.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "magnus.hpp"
#include "polarization.hpp"
#include "quad.hpp"
#include "ran.hpp"
#include "rotation_cache.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>

struct Basis {
  double3 r_hat, psi_hat, k_hat, s1, s2, B, omega;
  double Bperp2, r, c, s, kappa, Bmag, mu, twist, J, energy_local;

  Basis() = default;

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
    double3 phi_hat =
        rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
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
    energy_local = energy / L;
    double K = 2.013919e-19 * energy_local;
    kappa = K * Bperp2;
    omega = {K * diff, K * twice, 0};
  }
};

template <PhotonInitializer Init, typename Scat>
struct PhotonEvolution {
  struct Options {
    double eta_couple = 1e-3, delta_stokes = 1e-3;
    bool thinning = false;
    double polarization_tolerance = 1e-6;
    double cache_tolerance = 1e-7;
    bool profile = false;
  } options;

  BField const& bfield;
  Init const& initializer;
  Scat const& resonant;
  Ran ran;
  Geodesic geodesic;
  Quad quad;

  double3 n, e1, e2;
  double energy_inf, tau;
  Polarization pol;

  PhotonEvolution(double M_star, double R_star, int seed, BField const& bfield,
                  Init const& init, Scat const& scat, Options options = {})
      : options(options), bfield(bfield), initializer(init), resonant(scat),
        ran(seed), geodesic(M_star, R_star) {}

  Basis
  basis(Geodesic::State y) const {
    auto [r, psi, alpha] = y;
    return Basis(r, psi, alpha, n, e1, e2, bfield, energy_inf, geodesic.rs);
  }

  Basis
  basis(double l) const {
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
    energy_inf = photon.energy_inf;
    tau = std::log(ran.U());
    pol.set_adiabatic(photon.pol);
  }

  template <class PolarizationAt>
  bool
  integrate_optical_depth(double left, double right,
                          PolarizationAt const& pol_at) {
    auto event = resonant.integrate_optical_depth(
        left, right, [&](double l) { return basis(l); }, pol_at, tau, ran,
        quad);
    if (!event) {
      return false;
    } else {
      scatter_data = *event;
      return true;
    }
  }

  void
  record_escape(Geodesic::State y) {
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
    escape_data = {.energy_inf = energy_inf,
                   .muk = std::clamp(k.z, -1.0, 1.0),
                   .Q = sky.x,
                   .U = sky.y,
                   .V = sky.z};
  }

  enum class EvolveResult { Absorbed, Escaped, Scattered };

  struct EscapeData {
    double energy_inf, muk;
    double Q, U, V;
  } escape_data;

  Scat::ScatterData scatter_data;

  EvolveResult
  evolve() {
    return options.thinning ? evolve_thinning() : evolve_reference();
  }

  EvolveResult
  evolve_reference() {
    while (true) {
      int geodesic_event = geodesic.step();
      if (pol.is_adiabatic) {
        auto f_couple = [&](double l) {
          Basis b = basis(l);
          return options.eta_couple * b.kappa * b.r - 1;
        };
        if (f_couple(geodesic.l_old()) < 0) {
          Basis b = basis(geodesic.l_old());
          pol.set_stokes(b.c, b.s);
        }
      }
      if (pol.is_adiabatic) {
        if (integrate_optical_depth(geodesic.l_old(), geodesic.l_new(),
                                    [&](double) { return pol; }))
          return EvolveResult::Scattered;
      } else {
        StepperMagnus stepper(options.polarization_tolerance,
                              options.polarization_tolerance);
        stepper.init(geodesic.l_old(), geodesic.l_new(),
                     0.1 * (geodesic.l_new() - geodesic.l_old()), pol.stokes);
        auto stokes_rotation = [&](double l) { return basis(l).omega; };
        while (stepper.do_step(stokes_rotation)) {
          if (integrate_optical_depth(
                  stepper.x_old, stepper.x_new,
                  [&](double l) { return Polarization(stepper.interp(l)); }))
            return EvolveResult::Scattered;

          pol.set_stokes(stepper.y_new);
          Geodesic::State y = geodesic.interp(stepper.x_new);
          bool outward_and_far_enough =
              std::cos(y[2]) > 0 && y[0] > 1e3 * geodesic.R_star;
          double stokes_change = (stepper.y_new - stepper.y_old).length() *
                                 y[0] / (stepper.x_new - stepper.x_old);
          bool is_resonant_escaped = resonant.resonance_escape(
              1e-3 * stepper.h_old, basis(stepper.x_old), basis(stepper.x_new));
          if (outward_and_far_enough && stokes_change < options.delta_stokes &&
              is_resonant_escaped && stepper.phase_advance < pi / 2) {
            record_escape(y);
            return EvolveResult::Escaped;
          }
          stepper.update_old();
        }
      }
      if (geodesic_event == 0 && std::cos(geodesic.stepper.y_new[2]) < 0)
        return EvolveResult::Absorbed;
      if (geodesic_event == 1 && std::cos(geodesic.stepper.y_new[2]) > 0) {
        record_escape(geodesic.stepper.y_new);
        return EvolveResult::Escaped;
      }

      geodesic.stepper.update_old();
    }
  }

  unsigned long long candidates = 0, rejections = 0, magnus_steps = 0;
  double envelope_seconds = 0, polarization_seconds = 0;

  bool
  possible_resonance() const {
    // Bound the dense quartic by its power coefficients on t in [0,1].
    auto const& st = geodesic.stepper;
    double a = st.rcont1[0];
    double c[] = {st.rcont2[0] + st.rcont3[0],
                  st.rcont4[0] + st.rcont5[0] - st.rcont3[0],
                  -st.rcont4[0] - 2 * st.rcont5[0], st.rcont5[0]};
    double low = a, high = a;
    for (double x : c) {
      low += std::min(x, 0.0);
      high += std::max(x, 0.0);
    }
    low -= 1e-11 * (1 + std::abs(low));
    high += 1e-11 * (1 + std::abs(high));
    if (!(low > geodesic.rs && bfield.angular_lower > 0 &&
          std::isfinite(bfield.angular_upper)))
      return true;
    double min_x = B_to_energy * 0.5 * bfield.B_pole * bfield.angular_lower *
                   std::pow(bfield.R_star / high, 2 + bfield.p) *
                   std::sqrt(1 - geodesic.rs / low) / energy_inf;
    double max_x = B_to_energy * 0.5 * bfield.B_pole * bfield.angular_upper *
                   std::pow(bfield.R_star / low, 2 + bfield.p) / energy_inf;
    bool possible = false;
    resonant.for_each_component([&](std::size_t, auto const& comp) {
      if (comp.current_fraction == 0)
        return;
      double beta =
          std::max(std::abs(comp.dist.b_min), std::abs(comp.dist.b_max));
      if (!(beta < 1)) {
        possible = true;
        return;
      }
      double xmax = std::sqrt((1 + beta) / (1 - beta)) * (1 + 1e-12);
      possible |= !(comp.charge_mass_ratio * max_x < 1 / xmax ||
                    comp.charge_mass_ratio * min_x > xmax);
    });
    return possible;
  }

  EvolveResult
  evolve_thinning() {
    while (true) {
      int event = geodesic.step();
      double left = geodesic.l_old(), right = geodesic.l_new();
      auto get_basis = [&](double l) { return basis(l); };
      if (pol.is_adiabatic) {
        auto b = basis(left);
        if (options.eta_couple * b.kappa * b.r < 1)
          pol.set_stokes(b.c, b.s);
      }
      auto clock_start = options.profile
                             ? std::chrono::steady_clock::now()
                             : std::chrono::steady_clock::time_point{};
      auto cuts = possible_resonance()
                      ? resonant.integration_cuts(left, right, get_basis, true)
                      : std::vector<double>{};
      if (options.profile)
        envelope_seconds += std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - clock_start)
                                .count();
      RotationCache cache;
      cache.tolerance = options.cache_tolerance;
      int coefficient_calls = 0;
      auto exact_rotation = [&](double l) { return basis(l).omega; };
      auto rotation = [&](double l) {
        if (++coefficient_calls <= 180 || options.cache_tolerance == 0)
          return exact_rotation(l);
        if (cache.leaves.empty())
          cache.build(left, right, exact_rotation);
        return cache.at(l, exact_rotation);
      };
      double cursor = left;
      while (cursor < right) {
        clock_start = options.profile ? std::chrono::steady_clock::now()
                                      : std::chrono::steady_clock::time_point{};
        auto candidate =
            resonant.next_candidate(cursor, right, cuts, get_basis, tau, quad);
        if (options.profile)
          envelope_seconds +=
              std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                            clock_start)
                  .count();
        double end = candidate.value_or(right);
        clock_start = options.profile ? std::chrono::steady_clock::now()
                                      : std::chrono::steady_clock::time_point{};
        if (!pol.is_adiabatic) {
          StepperMagnus stepper(options.polarization_tolerance,
                                options.polarization_tolerance);
          stepper.init(cursor, end, end - cursor, pol.stokes);
          while (stepper.do_step(rotation, false)) {
            ++magnus_steps;
            pol.set_stokes(stepper.y_new);
            auto y = geodesic.interp(stepper.x_new);
            double change =
                (stepper.y_new - stepper.y_old).length() * y[0] / stepper.h_old;
            // Candidates at this endpoint take precedence over escape.
            if (!(candidate && stepper.x_new == end) && std::cos(y[2]) > 0 &&
                y[0] > 1e3 * geodesic.R_star && change < options.delta_stokes &&
                stepper.phase_advance < pi / 2 &&
                resonant.resonance_escape(stepper.h_old, basis(stepper.x_old),
                                          basis(stepper.x_new))) {
              record_escape(y);
              return EvolveResult::Escaped;
            }
            stepper.update_old();
          }
        }
        if (options.profile)
          polarization_seconds +=
              std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                            clock_start)
                  .count();
        if (!candidate)
          break;
        ++candidates;
        auto scatter = resonant.accept_candidate(basis(end), pol, ran);
        if (scatter) {
          scatter_data = *scatter;
          return EvolveResult::Scattered;
        }
        ++rejections;
        tau = std::log(ran.U());
        cursor = end;
      }
      if (event == 0 && std::cos(geodesic.stepper.y_new[2]) < 0)
        return EvolveResult::Absorbed;
      if (event == 1 && std::cos(geodesic.stepper.y_new[2]) > 0) {
        record_escape(geodesic.stepper.y_new);
        return EvolveResult::Escaped;
      }
      geodesic.stepper.update_old();
    }
  }

  void
  perform_scattering() {
    init(resonant.perform_scattering(scatter_data, energy_inf, ran));
  }
};
