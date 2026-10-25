#pragma once

#include "bfield.hpp"
#include "constants.hpp"
#include "doubles.hpp"
#include "find_root.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "polarization.hpp"
#include "quad.hpp"
#include "ran.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

// Scattering physics uses local geometry supplied by the transport solver.
// Random draws and quadrature share the caller's generator and workspace.
template <typename Comp>
struct ResonantScattering {
  struct Options {
    double quad_atol = 1e-8, quad_rtol = 1e-8;
    // n_knots sets the constructor's velocity grid; n_scan is used per
    // interval.
    int n_knots = 8, n_scan = 8;
  };

  BField const& bfield;
  Comp comp;
  std::vector<double> knots;

  ResonantScattering(BField const& bfield, Comp const& comp,
                     Options const& options = Options{})
      : bfield(bfield), comp(comp) {
    for (int i = 0; i < options.n_knots; ++i)
      knots.push_back(comp.charge_sign *
                      std::lerp(comp.dist.b_min, comp.dist.b_max,
                                double(i) / (options.n_knots - 1)));
  }

  struct Event {
    double l, beta;
    Polarization incoming;
  };

  template <class Basis>
  static
  double
  x_e(Basis const& b, double omega_inf, double rs) {
    double L = std::sqrt(1 - rs / b.r);
    return B_to_omega * b.Bmag * L / omega_inf;
  }

  static
  double
  discriminant(double x, double mu) {
    return std::fma(x, x, (mu - 1) * (mu + 1));
  }

  // Sorted roots of (x^2+mu^2) beta^2 - 2 mu beta + (1-x^2)=0.
  // Compute one root without cancellation, the other by their product.
  static
  std::array<double, 2>
  resonance_roots(double x, double mu, double D) {
    double a = x * x + mu * mu;
    double q = mu + std::copysign(x * std::sqrt(D), mu);
    std::array<double, 2> roots{q / a, ((1 - x) * (1 + x)) / q};
    if (roots[1] < roots[0])
      std::swap(roots[0], roots[1]);
    return roots;
  }

  static
  double
  overlap(double qb, double v, double mup, int charge_sign) {
    double w = 0.25 * ((1 + mup * mup) + (mup * mup - 1) * qb) +
               0.5 * charge_sign * mup * v;
    if (!std::isfinite(w) || w < -1e-10 || w > 1 + 1e-10)
      throw std::runtime_error("Nonphysical polarization overlap");
    return std::clamp(w, 0.0, 1.0);
  }

  struct Rate {
    int iroot = 0; // 0 = lower beta root, 1 = upper beta root
    double beta = 0, weight = 0;
  };
  using Rates = std::array<Rate, 2>;

  template <class Basis>
  Rates
  rates(Basis const& b, Polarization const& pol, double omega_inf,
        double rs) const {
    Rates result{};
    Stokes sb = pol.in_B(b.c, b.s);
    double qb = sb.x, v = sb.z;

    if (comp.current_fraction == 0 || b.twist == 0)
      return result;
    double x = comp.charge_mass_ratio * x_e(b, omega_inf, rs);
    double D = discriminant(x, b.mu);
    // Exact coalescence is an integration boundary, not a sample.
    if (!(D > 0) || !(x > 0))
      return result;
    auto roots = resonance_roots(x, b.mu, D);
    for (int j = 0; j < 2; ++j) {
      double beta = roots[j];
      double f = comp.dist.f(comp.charge_sign * beta);
      if (!(beta > -1 && beta < 1) || f == 0)
        continue;
      double t = (1 - beta) * (1 + beta);
      double mup = std::clamp((b.mu - beta) / (1 - beta * b.mu), -1.0, 1.0);
      // Equivalent to (1-beta*mu)^2(1-beta^2)/abs(mu-beta),
      // but avoids cancellation at coalescing roots.
      double jacobian = x * x * t * std::sqrt(t) / std::sqrt(D);
      double w = comp.current_fraction * (bfield.p + 1) * pi * b.twist *
                 comp.dist.inv_abs_b_mean / b.r * f * jacobian *
                 overlap(qb, v, mup, comp.charge_sign);
      if (!std::isfinite(w) || w < 0)
        throw std::runtime_error("Invalid differential optical depth");
      result[j] = {j, beta, w};
    }
    return result;
  }
  static
  double
  total_rate(Rates const& rates) {
    double sum = 0;
    for (auto const& r : rates)
      sum += r.weight;
    return sum;
  }
  static
  Rate
  select_rate(Rates const& rates, Ran& ran) {
    double sum = total_rate(rates);
    if (!(sum > 0 && std::isfinite(sum)))
      throw std::runtime_error("No positive rate at scattering event");
    double target = ran.U() * sum, cumulative = 0;
    Rate last{};
    for (auto const& r : rates) {
      if (r.weight == 0)
        continue;
      last = r;
      cumulative += r.weight;
      if (target < cumulative)
        return r;
    }
    return last; // Only roundoff in the final cumulative sum.
  }

  template <class BasisAt>
  std::vector<double>
  integration_cuts(double left, double right, BasisAt const& get_basis,
                   double omega_inf, double rs, Options const& options) const {
    std::vector<double> cuts{left, std::midpoint(left, right), right};
    struct ScanPoint {
      double l, x, mu;
    };
    auto sample = [&](double l) {
      auto b = get_basis(l);
      return ScanPoint{l, comp.charge_mass_ratio * x_e(b, omega_inf, rs), b.mu};
    };
    // All boundary functions share this grid. Cache only the geometry they
    // need; root-finding trial points still evaluate the actual trajectory.
    std::vector<ScanPoint> grid;
    grid.reserve(options.n_scan + 1);
    grid.push_back(sample(left));
    for (int j = 1; j <= options.n_scan; ++j)
      grid.push_back(
          sample(std::lerp(left, right, double(j) / options.n_scan)));

    auto scan = [&](auto const& value) {
      auto f = [&](double l) {
        auto point = sample(l);
        return value(point.x, point.mu);
      };
      double a = grid.front().l, fa = value(grid.front().x, grid.front().mu);
      for (size_t j = 1; j < grid.size(); ++j) {
        double b = grid[j].l, fb = value(grid[j].x, grid[j].mu);
        if (fa == 0)
          cuts.push_back(a);
        if (fb == 0)
          cuts.push_back(b);
        if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0))
          cuts.push_back(ridders_root(f, a, b, fa, fb));
        a = b;
        fa = fb;
      }
    };

    scan([](double x, double mu) { return discriminant(x, mu); });
    for (double beta : knots) {
      double inv_gamma = std::sqrt((1 - beta) * (1 + beta));
      scan([&](double x, double mu) { return x * inv_gamma + beta * mu - 1; });
    }

    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    return cuts;
  }

  // tau is the remaining (negative) optical-depth threshold. At an event it
  // retains its value at the beginning of the containing integration segment.
  // The caller applies the returned event and resets tau after scattering.
  template <class BasisAt, class PolarizationAt>
  std::optional<Event>
  integrate_optical_depth(double left, double right, BasisAt const& get_basis,
                          PolarizationAt const& pol_at, double omega_inf,
                          double rs, double& tau, Ran& ran, Quad& quad,
                          Options const& options) const {
    auto cuts =
        integration_cuts(left, right, get_basis, omega_inf, rs, options);
    auto integrand = [&](double l) {
      return total_rate(rates(get_basis(l), pol_at(l), omega_inf, rs));
    };
    for (size_t j = 0; j + 1 < cuts.size(); ++j) {
      double a = cuts[j], b = cuts[j + 1];
      if (!(a < b))
        continue;
      auto integral = [&](double end) {
        return quad.qags(integrand, a, end, options.quad_atol,
                         options.quad_rtol);
      };
      double dtau = integral(b);
      if (dtau < 0)
        throw std::runtime_error("Negative integrated optical depth");
      if (dtau > 0 && tau + dtau > 0) {
        // Find the point where the optical depth reaches zero.
        auto residual = [&](double l) { return tau + integral(l); };
        double l = ridders_root(residual, a, b, tau, tau + dtau);
        Polarization incoming = pol_at(l);
        Rate event =
            select_rate(rates(get_basis(l), incoming, omega_inf, rs), ran);
        return Event{l, event.beta, incoming};
      }
      tau += dtau;
    }
    return std::nullopt;
  }

  template <class BasisForState>
  bool
  resonance_escape(double l, Geodesic const& geodesic,
                   BasisForState const& basis, double omega_inf) const {
    auto y = geodesic.interp(l);
    if (std::cos(y[2]) <= 0)
      return false; // Infinity reconstruction is only for the outgoing branch.
    auto b = basis(y);
    double rs = geodesic.rs;
    double h = 1e-3 * geodesic.stepper.h_old;
    auto minus = basis(geodesic.interp(l - h));
    double xe = x_e(b, omega_inf, rs);
    double dx = (xe - x_e(minus, omega_inf, rs)) / (h);
    double dmu = (b.mu - minus.mu) / (h);
    double x = comp.charge_mass_ratio * xe;
    double tx = comp.charge_mass_ratio * dx;
    return discriminant(x, b.mu) < 0 && tx < 0 &&
           std::abs(b.mu - x * dmu / tx) < 1;
  }

  static
  double
  sample_mup(Ran& ran) {
    double const u = ran.U();
    if (u < 0.75)
      return (8.0 / 3.0) * u - 1.0;
    double const y = 8.0 * u - 7.0;
    double m = std::max(std::abs(y), ran.U());
    m = std::max(m, ran.U());
    return std::copysign(m, y);
  }

  template <class Basis>
  static
  InitializedPhoton
  perform_scattering(Basis b, double beta, double omega_inf, Ran& ran) {
    double mup = sample_mup(ran);
    double mu_out = std::clamp((mup + beta) / (1 + beta * mup), -1.0, 1.0);
    double3 bhat = b.B / b.Bmag;
    double3 k = mu_out * bhat +
                std::sqrt((1 - mu_out) * (1 + mu_out)) * ran.unit_perp_to(bhat);
    k = to_unit(k);
    double3 normal = cross(b.r_hat, k);
    double cross_len = normal.length();
    normal = cross_len > 1e-14 ? normal / cross_len : ran.unit_perp_to(b.r_hat);
    double energy = omega_inf * (1 - beta * b.mu) / (1 - beta * mu_out);
    auto mode = ran.U() < 1 / (1 + mup * mup) ? AdiabaticPolarization::E
                                              : AdiabaticPolarization::O;
    return InitializedPhoton{.n = normal,
                             .e1 = b.r_hat,
                             .e2 = cross(normal, b.r_hat),
                             .omega_inf = energy,
                             .r0 = b.r,
                             .psi0 = 0,
                             .alpha0 = std::atan2(cross_len, dot(k, b.r_hat)),
                             .pol = mode};
  }
};
