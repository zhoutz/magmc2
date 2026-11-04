#pragma once

#include "component.hpp"
#include "evolution.hpp"
#include "qed_response.hpp"
#include <map>
#include <tuple>

namespace qed {
struct Boost {
  double x, mu, gamma;
};
inline
Boost
boost(double epsilon, double mu, double beta) {
  double gamma = 1 / std::sqrt((1 - beta) * (1 + beta));
  return {gamma * epsilon * (1 - beta * mu), (mu - beta) / (1 - beta * mu),
          gamma};
}
template <class D>
auto
support(D const& d) {
  if constexpr (requires { d.support_intervals(); })
    return d.support_intervals();
  else
    return std::vector<std::pair<double, double>>{{d.b_min, d.b_max}};
}
inline
std::vector<double>
speed_roots(double b, double epsilon, double mu) {
  double h = b - epsilon * epsilon * (1 - mu) * (1 + mu) / 2;
  if (h <= 0)
    return {};
  double R = h / epsilon, A = R * R + mu * mu,
         disc = std::fma(R, R, mu * mu - 1);
  if (disc < -64 * std::numeric_limits<double>::epsilon() * std::max(1., A))
    return {};
  disc = std::max(0., disc);
  double q = mu + std::copysign(R * std::sqrt(disc), mu);
  std::vector<double> out;
  if (q == 0)
    out.push_back(0);
  else {
    out.push_back(q / A);
    out.push_back((1 - R) * (1 + R) / q);
  }
  sort_cuts(out);
  out.erase(std::remove_if(out.begin(), out.end(),
                           [&](double v) {
                             return std::abs(v) >= 1 ||
                                    std::abs(boost(epsilon, mu, v).x - h) >
                                        1e-9 * std::max(h, epsilon);
                           }),
            out.end());
  return out;
}
inline
std::vector<double>
speed_cuts(double b, double epsilon, double mu, RestWidths widths, double left,
           double right) {
  std::vector<double> cuts{left, right};
  auto add = [&](double v) {
    if (v > left && v < right)
      cuts.push_back(v);
  };
  add(mu);
  auto roots = speed_roots(b, epsilon, mu);
  // Also resolve a nearby peak just outside the support and stationary minima.
  double stationary = std::clamp(mu, left, right);
  roots.push_back(stationary);
  for (double r : roots) {
    add(r);
    auto v = boost(epsilon, mu, r);
    double E = std::sqrt(1 + 2 * b + v.x * v.x * v.mu * v.mu);
    double derivative = epsilon * std::pow(v.gamma, 3) * (r - mu) / (1 + v.x);
    double second = epsilon * std::pow(v.gamma, 3) / (1 + v.x);
    for (int spin : {-1, 1}) {
      double width = widths.at(spin) * std::sqrt(1 + 2 * b) / E;
      // Quadratic estimate stays finite when two roots coalesce.
      double detuning =
          (2 * v.x + epsilon * epsilon * (1 - mu) * (1 + mu) - 2 * b) /
          (1 + v.x + E);
      double resolved =
          r == stationary ? std::max(width, std::abs(detuning)) : width;
      double scale =
          resolved / (std::abs(derivative) +
                      std::sqrt(derivative * derivative + second * resolved));
      for (double d = scale; d < 4; d *= 4) {
        require(d > 0, Failure::Convergence, "QED speed linewidth underflow");
        add(r - d);
        add(r + d);
      }
    }
  }
  sort_cuts(cuts);
  return cuts;
}
template <class F>
Integral<4>
integrate_speed(F const& f, double b, double epsilon, double mu,
                RestWidths widths, double left, double right, Options options) {
  // Exact coordinate change only: retain the complete numerator, spin widths,
  // interference and nonresonant wings at every node. A logarithmic distance
  // coordinate resolves both linewidths without hundreds of tiny beta panels.
  auto roots = speed_roots(b, epsilon, mu);
  std::vector<double> branches{left, right};
  if (mu > left && mu < right)
    branches.push_back(mu);
  sort_cuts(branches);
  Integral<4> total;
  for (size_t j = 1; j < branches.size(); ++j) {
    double lo = branches[j - 1], hi = branches[j], mid = std::midpoint(lo, hi);
    double center = std::clamp(mu, lo, hi),
           distance = std::numeric_limits<double>::infinity();
    for (double root : roots) {
      if ((mid < mu) != (root < mu) && root != mu)
        continue;
      double d = std::abs(root - std::clamp(root, lo, hi));
      if (d < distance) {
        distance = d;
        center = root;
      }
    }
    auto v = boost(epsilon, mu, center);
    if (!std::isfinite(v.gamma)) {
      center = mid;
      v = boost(epsilon, mu, center);
    }
    double E = std::sqrt(1 + 2 * b + v.x * v.x * v.mu * v.mu);
    double delta = (2 * v.x + epsilon * epsilon * (1 - mu) * (1 + mu) - 2 * b) /
                   (1 + v.x + E);
    double derivative =
        epsilon * std::pow(v.gamma, 3) * (center - mu) / (1 + v.x);
    double second = epsilon * std::pow(v.gamma, 3) / (1 + v.x);
    double width = widths.plus * std::sqrt(1 + 2 * b) / E;
    double resolved = std::max(width, std::abs(delta));
    double scale =
        resolved / (std::abs(derivative) +
                    std::sqrt(derivative * derivative + second * resolved));
    if (!(scale > 0 && scale < hi - lo)) {
      auto part =
          integrate<4>(f, {lo, hi}, options.speed_rtol, 0, options.max_panels);
      total.panels.insert(total.panels.end(), part.panels.begin(),
                          part.panels.end());
      for (int k = 0; k < 4; ++k) {
        total.value[k] += part.value[k];
        total.error[k] += part.error[k];
      }
      continue;
    }
    auto t = [&](double beta) { return std::asinh((beta - center) / scale); };
    auto beta = [&](double t) { return center + scale * std::sinh(t); };
    auto transformed = [&](double t) {
      return scaled<4>(f(beta(t)), scale * std::cosh(t));
    };
    std::vector<double> cuts{t(lo), t(hi)};
    if (center > lo && center < hi)
      cuts.push_back(0);
    double broad = widths.minus / widths.plus * scale;
    for (double offset : {-broad, -scale, scale, broad})
      if (center + offset > lo && center + offset < hi)
        cuts.push_back(t(center + offset));
    auto part = integrate<4>(transformed, cuts, options.speed_rtol, 0,
                             options.max_panels);
    for (auto p : part.panels) {
      p.left = std::clamp(beta(p.left), lo, hi);
      p.right = std::clamp(beta(p.right), lo, hi);
      p.asinh_center = center;
      p.asinh_scale = scale;
      total.panels.push_back(p);
    }
    for (int k = 0; k < 4; ++k) {
      total.value[k] += part.value[k];
      total.error[k] += part.error[k];
    }
  }
  return total;
}
inline
double3
perpendicular(double3 v) {
  double3 axis =
      std::abs(v.x) <= std::abs(v.y) && std::abs(v.x) <= std::abs(v.z)
          ? double3{1, 0, 0}
          : (std::abs(v.y) <= std::abs(v.z) ? double3{0, 1, 0}
                                            : double3{0, 0, 1});
  return to_unit(axis - dot(axis, v) * v);
}
} // namespace qed

template <class... Comps>
struct QEDScattering {
  static constexpr size_t n_components = sizeof...(Comps);
  static constexpr bool finite_width = true;
  using Options = qed::Options;
  Options options;
  std::tuple<Comps...> components;
  // The transport owns a copy of this object. These bounded exact-key caches
  // are never shared by worker threads and do not interpolate the response.
  mutable std::map<double, qed::RestWidths> width_cache;
  struct LocalRateResponse {
    qed::Response K;
    double Lambda = 0;
    double negligible_bound = 0;
    std::array<qed::Integral<4>, n_components> populations;
  };
  using Key = std::array<double, 4>;
  mutable std::map<Key, LocalRateResponse> rate_cache;
  struct ScatterData {
    Basis basis;
    double beta;
    size_t population;
    int charge_sign;
    Stokes incoming_B;
  };
  static
  size_t
  select_population(std::array<double, n_components> const& weights, double u) {
    double sum = 0;
    for (double w : weights) {
      qed::require(w >= 0 && std::isfinite(w), qed::Failure::Convergence,
                   "QED invalid population weight");
      sum += w;
    }
    qed::require(sum > 0 && u > 0 && u < 1, qed::Failure::Convergence,
                 "QED empty population selection");
    double target = u * sum;
    for (size_t i = 0; i < n_components; ++i) {
      if (target < weights[i])
        return i;
      target -= weights[i];
    }
    throw qed::Error(qed::Failure::Convergence,
                     "QED population normalization mismatch");
  }
  template <class F>
  void
  for_each(F const& f) const {
    [&]<size_t... J>(std::index_sequence<J...>) {
      (f(J, std::get<J>(components)), ...);
    }(std::index_sequence_for<Comps...>{});
  }
  explicit QEDScattering(std::tuple<Comps...> c, Options opts = {})
      : options(opts), components(std::move(c)) {
    options.validate();
    for_each([&](size_t, auto const& comp) {
      qed::require(
          std::isfinite(comp.current_fraction) && comp.current_fraction >= 0 &&
              comp.charge_mass_ratio == 1 && std::abs(comp.charge_sign) == 1,
          qed::Failure::InvalidInput,
          "QED kernel supports unit-charge leptons only");
      qed::require(std::isfinite(comp.dist.inv_abs_b_mean) &&
                       comp.dist.inv_abs_b_mean > 0,
                   qed::Failure::DensityClosure,
                   "QED missing nonzero current-density closure");
      double previous = -1;
      for (auto [a, b] : qed::support(comp.dist)) {
        qed::require(std::isfinite(a) && std::isfinite(b) && a >= previous &&
                         b > a && a >= -1 && b <= 1,
                     qed::Failure::InvalidInput, "Invalid velocity support");
        qed::require(comp.charge_sign == -1 ? b <= 0 : a >= 0,
                     qed::Failure::DensityClosure,
                     "QED current closure requires charge-aligned "
                     "unidirectional streams");
        previous = b;
      }
    });
  }
  void
  reset() {
    rate_cache.clear();
    width_cache.clear();
  }
  qed::RestWidths
  widths(double b) const {
    qed::require(b > 0 && std::isfinite(b), qed::Failure::InvalidInput,
                 "QED natural-width cache requires positive finite field");
    if (auto i = width_cache.find(b); i != width_cache.end())
      return i->second;
    auto w = qed::rest_widths(b, options.width_rtol);
    if (width_cache.size() >= 128)
      width_cache.clear();
    width_cache.emplace(b, w);
    return w;
  }
  template <class C>
  void
  check_support(Basis const& local, C const& comp) const {
    double epsilon = local.energy_local / electron_rest_keV;
    for (auto [lo, hi] : qed::support(comp.dist)) {
      for (double v : {lo, hi}) {
        // The co-propagating exactly axial tail has x -> 0; every other
        // |beta| -> 1 tail has unbounded x and fails the full-domain check.
        double x =
            std::abs(v) == 1
                ? (v * local.mu == 1 ? 0
                                     : std::numeric_limits<double>::infinity())
                : qed::boost(epsilon, local.mu, v).x;
        double invariant =
            2 * x + epsilon * epsilon * (1 - local.mu) * (1 + local.mu);
        qed::require(
            invariant < 8, qed::Failure::PairDomain,
            std::format("QED pair domain: r={}, B={}, energy={}, mu={}, "
                        "charge={}, beta=[{},{}], endpoint={}, C={}",
                        local.r, local.Bmag, local.energy_local, local.mu,
                        comp.charge_sign, lo, hi, v, invariant));
      }
    }
  }
  template <class C>
  qed::Response
  speed_response(Basis const& local, C const& comp, double beta) const {
    double b = local.Bmag / critical_field_gauss;
    auto incident =
        qed::boost(local.energy_local / electron_rest_keV, local.mu, beta);
    qed::AngularResponse angular{comp.charge_sign, b,         incident.x,
                                 incident.mu,      widths(b), options};
    double f = comp.dist.f(beta);
    qed::require(std::isfinite(f) && f >= 0, qed::Failure::DensityClosure,
                 "Invalid velocity number PDF");
    double prefactor = comp.current_fraction * local.J *
                       comp.dist.inv_abs_b_mean / elementary_charge_esu *
                       classical_electron_radius_cm *
                       classical_electron_radius_cm;
    return qed::scaled<4>(angular.reduced_integral(),
                          prefactor * f * (1 - beta * local.mu));
  }
  template <class C, class F>
  qed::Integral<4>
  integrate_population(Basis const& local, C const& comp, F const& f) const {
    double b = local.Bmag / critical_field_gauss,
           epsilon = local.energy_local / electron_rest_keV;
    auto w = widths(b);
    qed::Integral<4> total;
    for (auto [a, z] : qed::support(comp.dist)) {
      std::vector<double> knots{a, z};
      if constexpr (requires { comp.dist.integration_knots(); })
        for (double beta : comp.dist.integration_knots())
          if (beta > a && beta < z)
            knots.push_back(beta);
      qed::sort_cuts(knots);
      for (size_t segment = 1; segment < knots.size(); ++segment) {
        auto part =
            qed::integrate_speed(f, b, epsilon, local.mu, w, knots[segment - 1],
                                 knots[segment], options);
        total.panels.insert(total.panels.end(), part.panels.begin(),
                            part.panels.end());
        for (int j = 0; j < 4; ++j) {
          total.value[j] += part.value[j];
          total.error[j] += part.error[j];
        }
      }
    }
    return total;
  }
  LocalRateResponse
  local_response(Basis const& local) const {
    qed::require(std::isfinite(local.J) && local.J >= 0,
                 qed::Failure::DensityClosure, "QED invalid local current");
    LocalRateResponse out;
    bool active = false;
    for_each([&](size_t, auto const& c) { active |= c.current_fraction > 0; });
    if (local.J == 0 || !active)
      return out;
    qed::require(local.Bmag > 0 && std::isfinite(local.Bmag) &&
                     local.energy_local > 0 &&
                     std::isfinite(local.energy_local) &&
                     std::isfinite(local.mu) && std::abs(local.mu) <= 1,
                 qed::Failure::InvalidInput,
                 "QED invalid local geometry or photon energy");
    Key key{local.Bmag, local.energy_local, local.J, local.mu};
    if (auto i = rate_cache.find(key); i != rate_cache.end())
      return i->second;
    // Validate every population before beginning any normalization integral.
    for_each([&](size_t, auto const& c) {
      if (c.current_fraction > 0)
        check_support(local, c);
    });
    double b = local.Bmag / critical_field_gauss,
           epsilon = local.energy_local / electron_rest_keV;
    auto w = widths(b);
    if (options.rate_atol > 0 &&
        epsilon * epsilon * (1 - local.mu) * (1 + local.mu) / (2 * b) > 20) {
      double bound = 0;
      for_each([&](size_t, auto const& comp) {
        if (comp.current_fraction == 0)
          return;
        auto upper = [&](double beta) {
          auto v = qed::boost(epsilon, local.mu, beta);
          qed::AngularResponse angular{comp.charge_sign, b, v.x, v.mu, w,
                                       options};
          double f = comp.dist.f(beta);
          qed::require(f >= 0 && std::isfinite(f), qed::Failure::DensityClosure,
                       "Invalid velocity PDF");
          double prefactor = comp.current_fraction * local.J *
                             comp.dist.inv_abs_b_mean / elementary_charge_esu *
                             classical_electron_radius_cm *
                             classical_electron_radius_cm;
          double value =
              f == 0 ? 0
                     : std::exp(angular.log_bound() + std::log(prefactor) +
                                std::log(f) + std::log1p(-beta * local.mu));
          return qed::Response{{value, value, 0, 0}, {}};
        };
        auto integral = integrate_population(local, comp, upper);
        bound += (integral.value[0] + integral.error[0]) * (1 + 1e-8);
      });
      if (bound <= options.rate_atol) {
        // A numerical zero with a nonzero certified angular bound and a
        // resolved velocity-integration error, never a claim of zero opacity
        // or an escape condition. The path integrator carries this error.
        out.negligible_bound = bound;
        out.K.error = {bound, bound, bound / 2, bound / 2};
        if (rate_cache.size() >= 8)
          rate_cache.clear();
        rate_cache.emplace(key, out);
        return out;
      }
    }
    for_each([&](size_t index, auto const& comp) {
      if (comp.current_fraction == 0)
        return;
      auto f = [&](double beta) { return speed_response(local, comp, beta); };
      auto& total = out.populations[index];
      total = integrate_population(local, comp, f);
      for (int j = 0; j < 4; ++j) {
        out.K.value[j] += total.value[j];
        out.K.error[j] += total.error[j];
      }
    });
    out.Lambda = qed::envelope(out.K);
    if (rate_cache.size() >= 8)
      rate_cache.clear();
    rate_cache.emplace(key, out);
    return out;
  }

  // Split the path using the finite-field resonance conditions at support
  // endpoints and at the speed minimum beta=mu (entry/exit and tangency).
  template <class BasisAt>
  std::vector<double>
  path_cuts(double left, double right, BasisAt const& exact) const {
    std::vector<double> cuts{left, right};
    auto add = [&](double l) {
      if (l > left && l < right)
        cuts.push_back(l);
    };
    std::vector<double> markers;
    for_each([&](size_t, auto const& c) {
      if (c.current_fraction > 0)
        for (auto [lo, hi] : qed::support(c.dist)) {
          if (lo > -1)
            markers.push_back(lo);
          if (hi < 1)
            markers.push_back(hi);
        }
    });
    markers.push_back(2); // stationary beta=mu when it lies inside a support
    for (double marker : markers) {
      auto detune = [&](double l) {
        auto g = exact(l);
        double beta = marker == 2 ? g.mu : marker;
        if (std::abs(beta) >= 1)
          return 1.;
        double e = g.energy_local / electron_rest_keV;
        return 2 * qed::boost(e, g.mu, beta).x +
               e * e * (1 - g.mu) * (1 + g.mu) -
               2 * g.Bmag / critical_field_gauss;
      };
      auto scan = [&](auto&& self, double a, double z, double fa, double fz,
                      int depth) -> void {
        double m = std::midpoint(a, z), fm = detune(m);
        double curvature = std::abs(fm - (fa + fz) / 2);
        double closest = std::min({std::abs(fa), std::abs(fm), std::abs(fz)});
        if (depth == 4 && (fm - fa) * (fz - fm) < 0) {
          // Locate a turning point even when the layer only touches the
          // support boundary and there is no sign-changing root.
          double sign = fm < fa ? 1. : -1., lo = a, hi = z;
          constexpr double golden = .3819660112501051;
          double u = lo + golden * (hi - lo), v = hi - golden * (hi - lo);
          double fu = sign * detune(u), fv = sign * detune(v);
          for (int k = 0; k < 64; ++k) {
            if (fu < fv) {
              hi = v;
              v = u;
              fv = fu;
              u = lo + golden * (hi - lo);
              fu = sign * detune(u);
            } else {
              lo = u;
              u = v;
              fu = fv;
              v = hi - golden * (hi - lo);
              fv = sign * detune(v);
            }
          }
          double extremum = std::midpoint(lo, hi);
          add(extremum);
          auto g = exact(extremum);
          double field = g.Bmag / critical_field_gauss;
          if (field > 0 && curvature > 0) {
            auto w = widths(field);
            double second = 8 * curvature / ((z - a) * (z - a));
            for (double gamma : {w.minus, w.plus}) {
              double invariant_width = std::sqrt(1 + 2 * field) * gamma;
              double scale = std::sqrt(
                  2 * std::max(invariant_width, std::abs(detune(extremum))) /
                  second);
              for (double d = scale; d < right - left; d *= 4) {
                add(extremum - d);
                add(extremum + d);
              }
            }
          }
        }
        if (depth < 4 ||
            (depth < 24 && curvature > 0.02 * std::max(closest, 1e-300))) {
          self(self, a, m, fa, fm, depth + 1);
          self(self, m, z, fm, fz, depth + 1);
          return;
        }
        if ((fa > 0 && fz > 0) || (fa < 0 && fz < 0))
          return;
        double r = ridders_root(detune, a, z, fa, fz);
        add(r);
        auto g = exact(r);
        double b = g.Bmag / critical_field_gauss;
        double slope = std::abs((fz - fa) / (z - a));
        if (!(slope > 0 && b > 0))
          return;
        auto w = widths(b);
        // Invariant detuning has half-width E*Gamma = a*Gamma_rest.
        for (double gamma : {w.minus, w.plus}) {
          double scale = std::sqrt(1 + 2 * b) * gamma / slope;
          for (double d = scale; d < right - left; d *= 4) {
            qed::require(d > 0, qed::Failure::Convergence,
                         "Path linewidth underflow");
            add(r - d);
            add(r + d);
          }
        }
      };
      scan(scan, left, right, detune(left), detune(right), 0);
    }
    qed::sort_cuts(cuts);
    return cuts;
  }
  template <class BasisAt>
  std::optional<double>
  next_candidate(double left, double right, BasisAt const& exact, double& tau,
                 Quad&) const {
    auto rate = [&](double l) {
      auto r = local_response(exact(l));
      return qed::Estimate<1>{{r.Lambda},
                              {r.negligible_bound > 0
                                   ? r.negligible_bound
                                   : qed::response_error(r.K)}};
    };
    auto cuts = path_cuts(left, right, exact);
    return qed::candidate_from_rate(rate, cuts, tau, options.path_rtol,
                                    options.rate_atol, options.max_panels);
  }
  std::optional<ScatterData>
  accept_candidate(Basis const& local, Polarization const& pol,
                   Ran& ran) const {
    auto r = local_response(local);
    auto input = pol.in_B(local.c, local.s);
    qed::check_stokes(input);
    double rate = qed::contract(r.K.value, input), bound = r.Lambda;
    double tolerance = qed::response_error(r.K) + 1e-13 * bound;
    qed::require(rate >= -tolerance && rate <= bound + tolerance,
                 qed::Failure::Convergence, "QED candidate envelope violated");
    if (bound == 0 || ran.U() * bound >= rate)
      return std::nullopt;
    std::array<double, n_components> weights{};
    for (size_t i = 0; i < n_components; ++i)
      weights[i] = qed::polarized(r.populations[i], input).value[0];
    size_t selected = select_population(weights, ran.U());
    std::optional<ScatterData> result;
    for_each([&](size_t index, auto const& comp) {
      if (index != selected)
        return;
      auto cdf = qed::polarized(r.populations[index], input);
      auto f = [&](double beta) {
        return qed::polarized(speed_response(local, comp, beta), input);
      };
      double beta = qed::inverse_cdf(f, cdf, ran.U(), options.speed_rtol * 4);
      result = ScatterData{local, beta, index, comp.charge_sign, input};
    });
    qed::require(result.has_value(), qed::Failure::Convergence,
                 "QED population selection failed");
    return result;
  }
  InitializedPhoton
  perform_scattering(ScatterData const& event, double energy_inf,
                     Ran& ran) const {
    auto const& g = event.basis;
    double beta = event.beta, b = g.Bmag / critical_field_gauss;
    auto incident = qed::boost(g.energy_local / electron_rest_keV, g.mu, beta);
    qed::AngularResponse angular{event.charge_sign, b,         incident.x,
                                 incident.mu,       widths(b), options};
    double u_mu = ran.U(), u_phi = ran.U();
    auto sample = qed::sample_outgoing(angular, event.incoming_B, u_mu, u_phi);
    double mu = (sample.muf + beta) / (1 + beta * sample.muf);
    mu = std::clamp(mu, -1., 1.);
    double3 bhat = g.B / g.Bmag;
    double3 transverse = g.k_hat - g.mu * bhat;
    // At exact axial incidence choose O along screen_1. Since eO=-mu*t1,
    // t1=-mu*screen_1 gives the same incoming Jones coefficients.
    double3 t1 =
        transverse.length() > 1e-14 ? to_unit(transverse) : -g.mu * g.screen_1;
    double3 t2 = cross(bhat, t1);
    double3 radial = std::cos(sample.phi) * t1 + std::sin(sample.phi) * t2;
    double3 k = to_unit(mu * bhat + std::sqrt((1 - mu) * (1 + mu)) * radial);
    double3 normal = cross(g.r_hat, k);
    double norm = normal.length();
    normal = norm > 1e-14 ? normal / norm : qed::perpendicular(g.r_hat);
    double3 s1 = cross(normal, k);
    double3 O = -mu * radial + std::sqrt((1 - mu) * (1 + mu)) * bhat;
    double c = dot(O, s1), s = dot(O, normal), normcs = std::hypot(c, s);
    c /= normcs;
    s /= normcs;
    Stokes out =
        Polarization::from_B(sample.stokes_B, c * c - s * s, 2 * c * s);
    double energy = energy_inf / g.energy_local * electron_rest_keV *
                    incident.gamma * sample.kernel.kinematics.y *
                    (1 + beta * sample.muf);
    return {.n = normal,
            .e1 = g.r_hat,
            .e2 = cross(normal, g.r_hat),
            .energy_inf = energy,
            .r0 = g.r,
            .psi0 = 0,
            .alpha0 = std::atan2(norm, dot(k, g.r_hat)),
            .pol = Polarization(out)};
  }
};
