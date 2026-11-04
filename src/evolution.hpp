#pragma once

#include "bfield.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "magnus.hpp"
#include "polarization.hpp"
#include "quad.hpp"
#include "ran.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

struct Basis {
  double3 r_hat, B, omega;
  double3 k_hat, screen_1, screen_2;
  double r, c, s, kappa, Bmag, mu, J, energy_local;

  Basis() = default;

  Basis(double radius, double psi, double alpha, double3 n, double3 e1,
        double3 e2, BField const& field, double energy, double rs) {
    r = radius;
    r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    double3 psi_hat = cross(n, r_hat);
    k_hat = std::cos(alpha) * r_hat + std::sin(alpha) * psi_hat;
    double3 s1 = cross(n, k_hat);
    screen_1 = s1;
    screen_2 = n;
    double muz = std::clamp(r_hat.z, -1.0, 1.0);
    double rho = std::sqrt(r_hat.x * r_hat.x + r_hat.y * r_hat.y);
    double3 theta_hat =
        rho > 0 ? double3{r_hat.x * muz / rho, r_hat.y * muz / rho, -rho}
                : double3{std::copysign(1.0, muz), 0, 0};
    double3 phi_hat =
        rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    auto field_result = field.evaluate(r, muz);
    double3 B_vec = field_result.B;
    J = field_result.J;
    Bmag = field_result.Bmag;
    B = B_vec.x * r_hat + B_vec.y * theta_hat + B_vec.z * phi_hat;
    mu = Bmag > 0 ? std::clamp(dot(k_hat, B) / Bmag, -1.0, 1.0) : 0;
    double L = std::sqrt(1 - rs / r);
    double B1 = dot(B, s1), B2 = dot(B, n);
    double Bperp2 = B1 * B1 + B2 * B2;
    double diff = B1 * B1 - B2 * B2, twice = 2 * B1 * B2;
    c = Bperp2 > 0 ? diff / Bperp2 : 1;
    s = Bperp2 > 0 ? twice / Bperp2 : 0;
    if (Bmag > 0 && cross(k_hat, B / Bmag).length() <= 1e-14) {
      c = 1;
      s = 0;
    }
    energy_local = energy / L;
    double K = 2.013919e-19 * energy_local;
    kappa = K * Bperp2;
    omega = {K * diff, K * twice, 0};
  }
};

struct TransportCoefficients {
  double Bmag, energy_local, J, mu;
  double3 omega;
};

class TransportCache {
  using Values = std::array<double, 7>;
  struct Segment {
    double left, right;
    std::array<Values, 5> differences;
    bool direct = false;
    Values
    at(double l) const {
      double u = 4 * ((l - left) / (right - left));
      Values out;
      for (int j = 0; j < 7; ++j)
        out[j] = differences[0][j] +
                 u * (differences[1][j] +
                      (u - 1) * (differences[2][j] / 2 +
                                 (u - 2) * (differences[3][j] / 6 +
                                            (u - 3) * differences[4][j] / 24)));
      return out;
    }
  };
  std::vector<Segment> segments;
  double tolerance;

  template <class F>
  static
  Values
  sample(F const& exact, double l) {
    auto b = exact(l);
    return {b.Bmag, b.energy_local, b.J, b.mu, b.omega.x, b.omega.y, b.omega.z};
  }

  template <class F>
  void
  add(double left, double right, F const& exact, int depth) {
    Segment s{left, right, {}};
    std::array<Values, 9> values;
    for (int k = 0; k < 9; ++k)
      values[k] = sample(exact, std::lerp(left, right, k / 8.));
    auto d = std::array<Values, 5>{values[0], values[2], values[4], values[6],
                                   values[8]};
    for (int k = 0; k < 5; ++k) {
      s.differences[k] = d[0];
      for (int i = 0; i < 4 - k; ++i)
        for (int j = 0; j < 7; ++j)
          d[i][j] = d[i + 1][j] - d[i][j];
    }
    Values scale{};
    for (auto const& v : values)
      for (int j = 0; j < 7; ++j)
        scale[j] = std::max(scale[j], std::abs(v[j]));
    scale[3] = 1;
    scale[4] = scale[5] = scale[6] = std::sqrt(
        scale[4] * scale[4] + scale[5] * scale[5] + scale[6] * scale[6]);
    bool accurate = true;
    for (int k : {1, 3, 5, 7}) {
      auto estimate = s.at(std::lerp(left, right, k / 8.));
      for (int j = 0; j < 7; ++j)
        accurate &= std::abs(estimate[j] - values[k][j]) <=
                    tolerance * std::max(scale[j], 1e-100);
    }
    if (!accurate && depth < 4) {
      double mid = std::midpoint(left, right);
      add(left, mid, exact, depth + 1);
      add(mid, right, exact, depth + 1);
    } else {
      s.direct = !accurate;
      segments.push_back(s);
    }
  }

public:
  explicit TransportCache(double tolerance = 1e-3) : tolerance(tolerance) {}
  template <class F>
  void
  reset(double left, double right, F const& exact) {
    segments.clear();
    if (right > left)
      add(left, right, exact, 0);
  }
  template <class F>
  TransportCoefficients
  at(double l, F const& exact) const {
    auto it = std::lower_bound(
        segments.begin(), segments.end(), l,
        [](Segment const& s, double l) { return s.right < l; });
    if (it == segments.end() || l < it->left || it->direct) {
      auto b = exact(l);
      return {b.Bmag, b.energy_local, b.J, b.mu, b.omega};
    }
    auto v = it->at(l);
    // Positive endpoint samples do not guarantee a positive polynomial.
    if (!(v[0] > 0 && v[1] > 0 && v[2] >= 0 && std::abs(v[3]) <= 1)) {
      auto b = exact(l);
      return {b.Bmag, b.energy_local, b.J, b.mu, b.omega};
    }
    return {v[0], v[1], v[2], v[3], {v[4], v[5], v[6]}};
  }
};

enum class EvolveResult { Absorbed, Escaped, Scattered };

template <PhotonInitializer Init, typename Scat>
struct PhotonEvolution {
  struct Options {
    double eta_couple = 1e-3, delta_stokes = 1e-3;
    Geodesic::Options geodesic;
    // Optional monitoring of the direct reference solver; no physics changes.
    std::function<void(double, double)> progress;
  } options;

  BField const& bfield;
  Init const& initializer;
  Scat resonant;
  Ran ran;
  Geodesic geodesic;
  Quad quad;
  TransportCache transport_cache;

  double3 n, e1, e2;
  double energy_inf, tau;
  Polarization pol;

  PhotonEvolution(double M_star, double R_star, int seed, BField const& bfield,
                  Init const& init, Scat const& scat, Options options = {})
      : options(options), bfield(bfield), initializer(init), resonant(scat),
        ran(seed), geodesic(M_star, R_star, options.geodesic) {}

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
    if constexpr (requires { resonant.reset(); }) resonant.reset();
    geodesic.init(photon.r0, photon.psi0, photon.alpha0);
    n = photon.n;
    e1 = photon.e1;
    e2 = photon.e2;
    energy_inf = photon.energy_inf;
    tau = std::log(ran.U());
    pol = photon.pol;
  }

  void
  record_escape(Geodesic::State y) {
    if (pol.is_adiabatic) {
      auto b = basis(y);
      pol.set_stokes(b.c, b.s);
    }
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

  struct EscapeData {
    double energy_inf, muk;
    double Q, U, V;
  } escape_data;

  Scat::ScatterData scatter_data;

  EvolveResult
  evolve() {
    while (true) {
      if (options.progress)
        options.progress(geodesic.l_old(), geodesic.stepper.y_old[0]);
      int geodesic_event = geodesic.step();
      if (pol.is_adiabatic) {
        Basis b = basis(geodesic.l_old());
        if (options.eta_couple * b.kappa * b.r < 1)
          pol.set_stokes(b.c, b.s);
      }

      double cursor = geodesic.l_old();
      auto exact = [&](double l) { return basis(l); };
      transport_cache.reset(cursor, geodesic.l_new(), exact);
      auto coefficients = [&](double l) {
        return transport_cache.at(l, exact);
      };
      while (cursor < geodesic.l_new()) {
        auto candidate = [&] {
          if constexpr (requires { Scat::finite_width; })
            return resonant.next_candidate(cursor, geodesic.l_new(), exact, tau, quad);
          else
            return resonant.next_candidate(cursor, geodesic.l_new(), coefficients, tau, quad);
        }();
        double end = candidate.value_or(geodesic.l_new());

        if (!pol.is_adiabatic) {
          StepperMagnusWithoutDense stepper;
          stepper.init(cursor, end, end - cursor, pol.stokes);
          auto stokes_rotation = [&](double l) {
            return coefficients(l).omega;
          };
          while (stepper.do_step(stokes_rotation)) {
            pol.set_stokes(stepper.y_new);
            Geodesic::State y = geodesic.interp(stepper.x_new);
            if constexpr (!requires { Scat::finite_width; }) {
              bool outward_and_far_enough =
                  std::cos(y[2]) > 0 && y[0] > 1e3 * geodesic.R_star;
              double stokes_change = (stepper.y_new - stepper.y_old).length() *
                                     y[0] / (stepper.x_new - stepper.x_old);
              if (outward_and_far_enough &&
                  stokes_change < options.delta_stokes &&
                  resonant.resonance_escape(stepper.h_old, basis(stepper.x_old),
                                            basis(y))) {
                record_escape(y);
                return EvolveResult::Escaped;
              }
            }
            stepper.update_old();
          }
        }

        if (!candidate)
          break;

        auto scatter = resonant.accept_candidate(basis(end), pol, ran);
        if (scatter) {
          scatter_data = *scatter;
          return EvolveResult::Scattered;
        }
        tau = std::log(ran.U());
        cursor = end;
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

  void
  perform_scattering() {
    init(resonant.perform_scattering(scatter_data, energy_inf, ran));
  }
};
