#pragma once

#include "constants.hpp"
#include "doubles.hpp"
#include "evolution.hpp"
#include "find_root.hpp"
#include "initialize.hpp"
#include "polarization.hpp"
#include "quad.hpp"
#include "ran.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

// Scattering physics uses local geometry supplied by the transport solver.
// Random draws and quadrature share the caller's generator and workspace.
template <typename... Comps>
struct ResonantScattering {
  static constexpr std::size_t n_components = sizeof...(Comps);
  static_assert(n_components > 0,
                "At least one scattering component is required");
  struct Options {
    double quad_atol = 1e-6, quad_rtol = 1e-6;
    // n_knots sets the constructor's velocity grid;
    // n_scan is used per interval.
    // n_od_sample is used to skip intervals with no resonance.
    int n_knots = 4, n_scan = 4, n_od_sample = 4;
  } options;

  // Distributions and their bounds use signed velocities along B.
  std::tuple<Comps...> components;
  std::array<std::vector<double>, n_components> knots;

  template <class F>
  void
  for_each_component(F const& f) const {
    [&]<std::size_t... I>(std::index_sequence<I...>) {
      (f(I, std::get<I>(components)), ...);
    }(std::index_sequence_for<Comps...>{});
  }

  ResonantScattering(std::tuple<Comps...> components, Options options = {})
      : options(options), components(std::move(components)) {
    if (options.n_knots < 2 || options.n_scan < 2 || options.n_od_sample < 2)
      throw std::invalid_argument("Invalid resonance integration grid");
    for_each_component([&](std::size_t i, auto const& comp) {
      if (!std::isfinite(comp.current_fraction) || comp.current_fraction < 0 ||
          !std::isfinite(comp.charge_mass_ratio) || comp.charge_mass_ratio <= 0)
        throw std::invalid_argument(
            "Invalid current fraction or charge/mass ratio");
      for (int j = 0; j < options.n_knots; ++j)
        knots[i].push_back(std::lerp(comp.dist.b_min, comp.dist.b_max,
                                     double(j) / (options.n_knots - 1)));
    });
  }

  struct ScatterData {
    Basis basis;
    double beta;
  };

  static
  double
  x_e(double Bmag, double energy_local) {
    return B_to_energy * Bmag / energy_local;
  }

  static
  double
  discriminant(double x, double mu) {
    return std::fma(x, x, (mu - 1) * (mu + 1));
  }

  // Sorted roots of (x^2+mu^2) beta^2 - 2 mu beta + (1-x^2)=0.
  // Compute one root without cancellation, the other by their product.
  struct ResonantRoot {
    double beta, sqrtD, mup;
  };

  static
  std::array<ResonantRoot, 2>
  resonance_roots(double x, double mu, double D) {
    double a = x * x + mu * mu;
    double sqrtD = std::sqrt(D);
    double q = mu + std::copysign(x * sqrtD, mu);
    std::array<ResonantRoot, 2> roots{{
        {.beta = q / a, .sqrtD = sqrtD, .mup = sqrtD / x},
        {.beta = ((1 - x) * (1 + x)) / q, .sqrtD = sqrtD, .mup = -sqrtD / x},
    }};
    if (roots[1].beta < roots[0].beta) {
      // mup for the lower root is positive, for the upper root is negative.
      // so do not swap the mup values when swapping the roots.
      std::swap(roots[0].beta, roots[1].beta);
    }
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
    std::size_t icomp = 0;
    int jroot = 0; // 0 = lower beta root, 1 = upper beta root
    double beta = 0, weight = 0;
  };
  using Rates = std::array<Rate, 2 * n_components>;

  Rates
  rates(Basis const& b, Polarization const& pol) const {
    Stokes sb = pol.in_B(b.c, b.s);
    return rates_impl(b, [qb = sb.x, v = sb.z](double mup, int charge_sign) {
      return overlap(qb, v, mup, charge_sign);
    });
  }

  template <class Geometry>
  Rates
  rates_enveloped(Geometry const& b) const {
    return rates_impl(b, [](double mup, int) { return (1 + mup * mup) / 2; });
  }

  template <class Geometry, class OverlapFunc>
  Rates
  rates_impl(Geometry const& b, OverlapFunc const& get_overlap) const {
    Rates result{};

    if (!std::isfinite(b.J) || b.J < 0)
      throw std::runtime_error("Invalid local current density");
    for (std::size_t i = 0; i < n_components; ++i)
      for (int j = 0; j < 2; ++j)
        result[2 * i + j] = {i, j, 0, 0};
    if (b.J == 0)
      return result;
    if (!std::isfinite(b.Bmag) || b.Bmag <= 0)
      throw std::runtime_error(
          "Positive current requires positive magnetic strength");
    // J is J_physical/c, expressed with lengths in km, so rates are per km.
    double current_factor = 4 * pi * pi * b.J / b.Bmag;
    double xe = x_e(b.Bmag, b.energy_local);
    for_each_component([&](std::size_t i, auto const& comp) {
      if (comp.current_fraction == 0)
        return;
      double x = comp.charge_mass_ratio * xe;
      double D = discriminant(x, b.mu);
      // Exact coalescence is an integration boundary, not a sample.
      if (!(D > 0) || !(x > 0))
        return;
      auto roots = resonance_roots(x, b.mu, D);
      for (int j = 0; j < 2; ++j) {
        double beta = roots[j].beta;
        double sqrtD = roots[j].sqrtD;
        double mup = roots[j].mup;
        double f = comp.dist.f(beta);
        if (!(beta > -1 && beta < 1) || f == 0)
          continue;
        double t = (1 - beta) * (1 + beta);
        // Equivalent to (1-beta*mu)^2(1-beta^2)/abs(mu-beta),
        // but avoids cancellation at coalescing roots.
        double jacobian = x * x * t * std::sqrt(t) / sqrtD;
        double w = current_factor * comp.current_fraction *
                   comp.dist.inv_abs_b_mean * f * jacobian *
                   get_overlap(mup, comp.charge_sign);
        if (!std::isfinite(w) || w < 0)
          throw std::runtime_error("Invalid differential optical depth");
        result[2 * i + j] = {i, j, beta, w};
      }
    });
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
  integration_cuts(double left, double right, BasisAt const& get_basis) const {
    std::vector<double> cuts{left, right};
    struct ScanPoint {
      double l, xe, mu;
    };
    auto sample = [&](double l) {
      auto b = get_basis(l);
      return ScanPoint{l, x_e(b.Bmag, b.energy_local), b.mu};
    };
    // All boundary functions share this grid. Cache only the geometry they
    // need; root-finding trial points still evaluate the actual trajectory.
    std::vector<ScanPoint> grid;
    grid.reserve(options.n_scan);
    for (int j = 0; j < options.n_scan; ++j)
      grid.push_back(
          sample(std::lerp(left, right, double(j) / (options.n_scan - 1))));

    auto scan = [&](auto const& value) {
      auto f = [&](double l) {
        auto point = sample(l);
        return value(point.xe, point.mu);
      };
      double a = grid.front().l, fa = value(grid.front().xe, grid.front().mu);
      for (int j = 1; j < (int)grid.size(); ++j) {
        double b = grid[j].l, fb = value(grid[j].xe, grid[j].mu);
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

    for_each_component([&](std::size_t i, auto const& comp) {
      if (comp.current_fraction == 0)
        return;
      double ratio = comp.charge_mass_ratio;
      scan([&](double xe, double mu) { return discriminant(ratio * xe, mu); });
      for (double beta : knots[i]) {
        double inv_gamma = std::sqrt((1 - beta) * (1 + beta));
        scan([&](double xe, double mu) {
          return ratio * xe * inv_gamma + beta * mu - 1;
        });
      }
    });

    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    return cuts;
  }

  template <class BasisAt>
  std::optional<double>
  next_candidate(double left, double right, BasisAt const& get_basis,
                 double& tau, Quad& quad) const {
    auto envelope = [&](double l) {
      return total_rate(rates_enveloped(get_basis(l)));
    };
    auto cuts = integration_cuts(left, right, get_basis);
    for (size_t j = 0; j + 1 < cuts.size(); ++j) {
      double a = cuts[j], b = cuts[j + 1];
      {
        bool skip = true;
        for (int i_sample = 0; i_sample < options.n_od_sample; ++i_sample) {
          double l_sample =
              std::lerp(a, b, double(i_sample) / (options.n_od_sample - 1));
          if (envelope(l_sample) > 0) {
            skip = false;
            break;
          }
        }
        if (skip)
          continue;
      }
      auto integral = [&](double end) {
        return quad.qags(envelope, a, end, options.quad_atol,
                         options.quad_rtol);
      };
      double dtau = integral(b);
      if (tau + dtau > 0) {
        auto residual = [&](double l) { return tau + integral(l); };
        double l =
            ridders_root(residual, a, b, tau, tau + dtau, 1e-5 * (b - a));
        return l;
      }
      tau += dtau;
    }
    return std::nullopt;
  }

  std::optional<ScatterData>
  accept_candidate(Basis const& b, Polarization const& pol, Ran& ran) const {
    auto actual = rates(b, pol);
    double rate = total_rate(actual), bound = total_rate(rates_enveloped(b));
    if (ran.U() * bound < rate)
      return ScatterData{b, select_rate(actual, ran).beta};
    return std::nullopt;
  }

  bool
  resonance_escape(double h, Basis const& b_old, Basis const& b_new) const {
    double xe = x_e(b_new.Bmag, b_new.energy_local);
    double mu = b_new.mu;
    double dx = (xe - x_e(b_old.Bmag, b_old.energy_local)) / h;
    double dmu = (b_new.mu - b_old.mu) / h;
    bool escaped = true;
    for_each_component([&](std::size_t, auto const& comp) {
      if (comp.current_fraction == 0)
        return;
      double x = comp.charge_mass_ratio * xe;
      double tx = comp.charge_mass_ratio * dx;
      escaped = escaped && discriminant(x, mu) < 0 && tx < 0 &&
                std::abs(mu - x * dmu / tx) < 1;
    });
    return escaped;
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

  static
  InitializedPhoton
  perform_scattering(ScatterData scattered_data, double energy_inf, Ran& ran) {
    auto const& b = scattered_data.basis;
    double beta = scattered_data.beta;
    double mup = sample_mup(ran);
    double mu_out = std::clamp((mup + beta) / (1 + beta * mup), -1.0, 1.0);
    double3 bhat = b.B / b.Bmag;
    double3 k = mu_out * bhat +
                std::sqrt((1 - mu_out) * (1 + mu_out)) * ran.unit_perp_to(bhat);
    k = to_unit(k);
    double3 normal = cross(b.r_hat, k);
    double cross_len = normal.length();
    normal = cross_len > 1e-14 ? normal / cross_len : ran.unit_perp_to(b.r_hat);
    double energy_inf_new =
        energy_inf * (1 - beta * b.mu) / (1 - beta * mu_out);
    auto mode = ran.U() < 1 / (1 + mup * mup) ? AdiabaticPolarization::E
                                              : AdiabaticPolarization::O;
    return InitializedPhoton{.n = normal,
                             .e1 = b.r_hat,
                             .e2 = cross(normal, b.r_hat),
                             .energy_inf = energy_inf_new,
                             .r0 = b.r,
                             .psi0 = 0,
                             .alpha0 = std::atan2(cross_len, dot(k, b.r_hat)),
                             .pol = mode};
  }
};
