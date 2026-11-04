#pragma once

#include "qed_kernel.hpp"

namespace qed {
struct Options {
  double width_rtol = 1e-9, angular_rtol = 1e-7, speed_rtol = 1e-6,
         path_rtol = 1e-6;
  // Absolute opacity error budget in km^-1, used only with a complete
  // angular/velocity upper bound. Set zero to disable this numerical bound.
  double rate_atol = 1e-30;
  size_t max_panels = 8192;
  bool bound_azimuth_tails = true;
  void
  validate() const {
    for (double t : {width_rtol, angular_rtol, speed_rtol, path_rtol})
      require(std::isfinite(t) && t >= 1e-12 && t < 0.1, Failure::InvalidInput,
              "QED tolerances must be in [1e-12,0.1)");
    require(max_panels >= 16, Failure::InvalidInput,
            "QED panel limit too small");
    require(rate_atol >= 0 && std::isfinite(rate_atol), Failure::InvalidInput,
            "QED absolute rate tolerance must be finite and nonnegative");
  }
};

inline
std::vector<double>
angular_cuts(double b, double x, double mui) {
  std::vector<double> cuts{-1, 1};
  // Resolve the outgoing Gaussian boundary layers, including axial incidence.
  for (double end : {-1., 1.}) {
    double y = recoil(x, mui, end).y;
    double scale = b / (y * y);
    for (double d = scale; d < 1; d *= 4) {
      require(d > 0, Failure::Convergence, "Angular boundary layer underflow");
      double point = end * (1 - d);
      if (point > -1 && point < 1)
        cuts.push_back(point);
    }
  }
  sort_cuts(cuts);
  return cuts;
}

struct AngularResponse {
  int charge;
  double b, x, mui;
  RestWidths widths;
  Options options;

  // Conservative all-angle upper bound on the largest eigenvalue of
  // integral (M/r_e)^dagger(M/r_e). Used ONLY for representational underflow,
  // never as a phenomenological opacity cutoff. |V_out|<=1; each V_in is
  // bounded by 4 exp(-ki_perp^2/4b) max(1,ki_perp/sqrt(2b)).
  double
  log_bound() const {
    double C = 2 * x + x * x * (1 - mui) * (1 + mui);
    double ymaxperp = C / (std::sqrt(1 + C) + 1);
    double ymin = x / (1 + 2 * x), sum = 0;
    for (int n : {0, 1}) {
      double E = std::sqrt(1 + 2 * n * b + x * x * mui * mui);
      for (int s : {-1, 1}) {
        if (n == 0 && s == 1)
          continue;
        double width = n ? widths.at(s) * std::sqrt(1 + 2 * b) / E : 0;
        double ds = (C - 2 * n * b) / (1 + x + E);
        sum += 1 / std::hypot(ds, width / 2) +
               1 / std::hypot(1 + x + E, width / 2);
        // Both u-channel energy signs have this common conservative bound.
        sum += 4 * (1 + x) / (ymin * (2 - ymaxperp) + 2 * n * b);
      }
    }
    double h2 = x * x * (1 - mui) * (1 + mui) / (2 * b);
    double logA = std::log(4.) - h2 / 2 + 0.5 * std::log(std::max(1., h2)) +
                  std::log(sum);
    return std::log(32 * pi) + 2 * std::log1p(x) + 2 * logA;
  }
  bool
  underflows() const {
    return log_bound() < std::log(std::numeric_limits<double>::denorm_min());
  }
  bool
  polar_underflows(double muf) const {
    double y = recoil(x, mui, muf).y;
    double h2 = y * y * (1 - muf) * (1 + muf) / (2 * b);
    // Bound the outgoing vertex as well as the incoming one. This covers
    // every azimuth, including nodes not yet visited by quadrature.
    double log_out =
        std::min(0., std::log(4.) - h2 / 2 + 0.5 * std::log(std::max(1., h2)));
    return log_bound() + 2 * log_out <
           std::log(std::numeric_limits<double>::denorm_min());
  }
  std::vector<double>
  phi_cuts(double muf, double end = 2 * pi) const {
    double y = recoil(x, mui, muf).y;
    double phase =
        x * y * std::sqrt((1 - mui) * (1 + mui) * (1 - muf) * (1 + muf)) / b;
    double count = std::max(4., std::ceil(phase));
    require(count <= options.max_panels, Failure::Convergence,
            "QED azimuth phase exceeds integration budget");
    std::vector<double> cuts{0};
    for (size_t j = 1; j <= static_cast<size_t>(count); ++j)
      cuts.push_back(end * j / count);
    return cuts;
  }
  Response
  at(double muf, double phi) const {
    return response(differential(charge, b, x, mui, muf, phi, widths).reduced);
  }
  Response
  polar(double muf, std::array<double, 2> absolute_budget = {}) const {
    if (polar_underflows(muf))
      return {{},
              {std::numeric_limits<double>::denorm_min(),
               std::numeric_limits<double>::denorm_min(), 0, 0}};
    // Positive periodic trapezoidal quadrature; seed from the largest
    // interference phase frequency, then double and reuse all old nodes.
    size_t n = 8;
    AzimuthKernel line(charge, b, x, mui, muf, widths);
    if (absolute_budget[0] > 0 && absolute_budget[1] > 0) {
      auto bound = line.integral_bound();
      if (bound.error[0] <= absolute_budget[0] &&
          bound.error[1] <= absolute_budget[1])
        return bound;
    }
    double y = recoil(x, mui, muf).y;
    double phase =
        x * y * std::sqrt((1 - mui) * (1 + mui) * (1 - muf) * (1 + muf)) / b;
    while (n < 4 * phase + 8 && n < options.max_panels)
      n *= 2;
    Response sum;
    for (size_t i = 0; i < n; ++i) {
      auto v = response(line.reduced(2 * pi * i / n));
      for (int j = 0; j < 4; ++j)
        sum.value[j] += v.value[j];
    }
    while (n * 2 <= options.max_panels) {
      auto old = qed::scaled<4>(sum, 2 * pi / n);
      for (size_t i = 0; i < n; ++i) {
        auto v = response(line.reduced(2 * pi * (i + 0.5) / n));
        for (int j = 0; j < 4; ++j)
          sum.value[j] += v.value[j];
      }
      n *= 2;
      auto result = qed::scaled<4>(sum, 2 * pi / n);
      double scale =
          std::max(std::abs(result.value[0]), std::abs(result.value[1]));
      bool good = true;
      for (int j = 0; j < 4; ++j) {
        result.error[j] = std::max(4 * std::abs(result.value[j] - old.value[j]),
                                   64 * std::numeric_limits<double>::epsilon() *
                                       std::abs(result.value[j]));
        good &= result.error[j] <=
                std::max(options.angular_rtol * 0.1 * scale,
                         64 * std::numeric_limits<double>::denorm_min());
      }
      if (good)
        return result;
    }
    throw Error(Failure::Convergence,
                std::format("QED periodic azimuth quadrature exhausted: b={}, "
                            "x={}, mui={}, muf={}, phase={}, n={}",
                            b, x, mui, muf, phase, n));
  }
  Integral<4>
  reduced_integral() const {
    require(b > 0 && std::isfinite(b) && std::abs(charge) == 1,
            Failure::InvalidInput,
            "QED angular response requires a positive field and unit charge");
    check_domain(x, mui);
    if (underflows()) {
      Integral<4> out;
      out.error[0] = std::numeric_limits<double>::denorm_min();
      out.error[1] = std::numeric_limits<double>::denorm_min();
      out.panels.push_back({-1, 1, static_cast<Response const&>(out)});
      return out;
    }
    auto cuts = angular_cuts(b, x, mui);
    std::array<double, 2> absolute_budget{};
    if (options.bound_azimuth_tails && x * x / b > 4) {
      // A positive coarse rule sets a provisional scale, NEVER certifies
      // accuracy. The subsequent adaptive integral includes all omitted-tail
      // bounds and must independently satisfy its global error budget.
      auto pilot = [&](double mu) {
        Response sum;
        if (polar_underflows(mu))
          return sum;
        AzimuthKernel line(charge, b, x, mui, mu, widths);
        for (int i = 0; i < 8; ++i) {
          auto v = response(line.reduced(2 * pi * i / 8));
          for (int j = 0; j < 4; ++j)
            sum.value[j] += 2 * pi / 8 * v.value[j];
        }
        return sum;
      };
      for (size_t i = 1; i < cuts.size(); ++i) {
        auto p = rule<4>(pilot, cuts[i - 1], cuts[i]);
        for (int j = 0; j < 2; ++j)
          absolute_budget[j] += p.value[j] * options.angular_rtol / 4096;
      }
    }
    auto f = [&](double mu) { return polar(mu, absolute_budget); };
    auto out =
        integrate<4>(f, cuts, options.angular_rtol, 0, options.max_panels);
    envelope(out);
    return out;
  }
  Response
  sigma() const {
    return scaled<4>(reduced_integral(), classical_electron_radius_cm *
                                             classical_electron_radius_cm);
  }
};

inline
Estimate<1>
polarized(Response const& r, Stokes s) {
  double value = contract(r.value, s);
  double err = std::abs(1 + s.x) * r.error[0] / 2 +
               std::abs(1 - s.x) * r.error[1] / 2 + std::abs(s.y) * r.error[2] +
               std::abs(s.z) * r.error[3];
  require(value >= -err - 1e-13 * std::abs(r.value[0]), Failure::Convergence,
          "QED negative polarized response");
  // Roundoff alone may make an exactly dark state infinitesimally negative.
  return {{std::max(0., value)}, {err}};
}
inline
Integral<1>
polarized(Integral<4> const& r, Stokes s) {
  Integral<1> out;
  auto v = polarized(static_cast<Response const&>(r), s);
  out.value = v.value;
  out.error = v.error;
  for (auto const& p : r.panels)
    out.panels.push_back({p.left, p.right, polarized(p.integral, s),
                          p.asinh_center, p.asinh_scale});
  return out;
}

struct OutgoingSample {
  double muf, phi;
  DifferentialKernel kernel;
  Stokes stokes_B;
};
inline
OutgoingSample
sample_outgoing(AngularResponse const& angular, Stokes input, double u_mu,
                double u_phi) {
  auto total = polarized(angular.reduced_integral(), input);
  auto polar = [&](double mu) { return polarized(angular.polar(mu), input); };
  double mu = inverse_cdf(polar, total, u_mu, angular.options.angular_rtol * 4);
  AzimuthKernel line(angular.charge, angular.b, angular.x, angular.mui, mu,
                     angular.widths);
  auto azimuth = [&](double phi) {
    return polarized(response(line.reduced(phi)), input);
  };
  auto conditional =
      integrate<1>(azimuth, angular.phi_cuts(mu), angular.options.angular_rtol,
                   0, angular.options.max_panels);
  double phi = inverse_cdf(azimuth, conditional, u_phi,
                           angular.options.angular_rtol * 4);
  auto kernel = differential(angular.charge, angular.b, angular.x, angular.mui,
                             mu, phi, angular.widths);
  return {mu, phi, kernel, outgoing_stokes(kernel.reduced, input)};
}
} // namespace qed
