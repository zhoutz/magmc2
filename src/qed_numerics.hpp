#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace qed {
enum class Failure { InvalidInput, PairDomain, DensityClosure, Convergence };
inline
char const*
failure_name(Failure kind) {
  switch (kind) {
  case Failure::InvalidInput:
    return "invalid-input";
  case Failure::PairDomain:
    return "pair-domain";
  case Failure::DensityClosure:
    return "density-closure";
  case Failure::Convergence:
    return "numerical-convergence";
  }
  return "unknown";
}
struct Error : std::runtime_error {
  Failure kind;
  Error(Failure kind, std::string const& message)
      : std::runtime_error(message), kind(kind) {}
};
inline
void
require(bool ok, Failure kind, std::string const& message) {
  if (!ok)
    throw Error(kind, message);
}

// Error is an absolute bound/estimate in the same units as value. Inner
// quadrature errors are integrated with positive weights at the next layer.
template <size_t N>
struct Estimate {
  std::array<double, N> value{}, error{};
};
template <size_t N>
struct Panel {
  double left, right;
  Estimate<N> integral;
  // Optional exact beta=center+scale*sinh(t) integration coordinate.
  double asinh_center = 0, asinh_scale = 0;
};
template <size_t N>
struct Integral : Estimate<N> {
  std::vector<Panel<N>> panels;
};
template <size_t N>
Estimate<N>
scaled(Estimate<N> a, double scale) {
  for (size_t j = 0; j < N; ++j) {
    a.value[j] *= scale;
    a.error[j] *= std::abs(scale);
  }
  return a;
}

// Positive 15-point Kronrod rule, embedded 7-point Gauss rule. No shared
// workspace, so arbitrary nesting and concurrent photons are safe.
template <size_t N, class F>
Estimate<N>
rule(F const& f, double left, double right) {
  constexpr double x[] = {0.9914553711208126, 0.9491079123427585,
                          0.8648644233597691, 0.7415311855993945,
                          0.5860872354676911, 0.4058451513773972,
                          0.2077849550078985, 0};
  constexpr double wk[] = {0.02293532201052922, 0.06309209262997855,
                           0.1047900103222502,  0.1406532597155259,
                           0.1690047266392679,  0.1903505780647854,
                           0.2044329400752989,  0.2094821410847278};
  constexpr double wg[] = {0.1294849661688697, 0.2797053914892767,
                           0.3818300505051189, 0.4179591836734694};
  double mid = std::midpoint(left, right), half = (right - left) / 2;
  Estimate<N> out;
  std::array<double, N> low{}, magnitude{};
  for (int i = 0; i < 8; ++i) {
    for (int sign : {-1, 1}) {
      if (i == 7 && sign == 1)
        continue;
      auto v = f(mid + sign * half * x[i]);
      for (size_t j = 0; j < N; ++j) {
        require(std::isfinite(v.value[j]) && std::isfinite(v.error[j]) &&
                    v.error[j] >= 0,
                Failure::Convergence,
                "QED quadrature: nonfinite integrand or invalid error");
        out.value[j] += wk[i] * v.value[j];
        out.error[j] += wk[i] * v.error[j];
        magnitude[j] += wk[i] * std::abs(v.value[j]);
        if (i % 2 == 1)
          low[j] += wg[i / 2] * v.value[j];
      }
    }
  }
  for (size_t j = 0; j < N; ++j) {
    out.error[j] =
        half *
        (out.error[j] +
         std::max(4 * std::abs(out.value[j] - low[j]),
                  64 * std::numeric_limits<double>::epsilon() * magnitude[j]));
    out.value[j] *= half;
  }
  return out;
}

inline
void
sort_cuts(std::vector<double>& cuts) {
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
}

template <size_t N, class F>
Integral<N>
integrate(F const& f, std::vector<double> cuts, double rtol, double atol = 0,
          size_t limit = 8192) {
  require(rtol > 0 && std::isfinite(rtol) && atol >= 0 && std::isfinite(atol) &&
              cuts.size() >= 2,
          Failure::InvalidInput, "QED quadrature: invalid options");
  sort_cuts(cuts);
  Integral<N> out;
  if (cuts.size() < 2)
    return out;
  size_t count = 0;
  auto evaluate = [&](double a, double b) {
    require(++count <= limit, Failure::Convergence,
            "QED quadrature: panel budget exhausted");
    return Panel<N>{a, b, rule<N>(f, a, b)};
  };
  for (size_t i = 1; i < cuts.size(); ++i)
    out.panels.push_back(evaluate(cuts[i - 1], cuts[i]));
  for (;;) {
    out.value.fill(0);
    out.error.fill(0);
    size_t worst = 0;
    double worst_error = -1;
    for (size_t i = 0; i < out.panels.size(); ++i) {
      auto const& p = out.panels[i];
      double local_error = 0;
      for (size_t j = 0; j < N; ++j) {
        out.value[j] += p.integral.value[j];
        out.error[j] += p.integral.error[j];
        local_error = std::max(local_error, p.integral.error[j]);
      }
      if (local_error > worst_error) {
        worst = i;
        worst_error = local_error;
      }
    }
    // Control the error of the WHOLE integral, as in globally adaptive GSL
    // quadrature. Requiring a relative tolerance independently in every
    // exponentially small angular tail wastes work without improving the
    // normalization. All positive panels and their errors remain in the sum.
    double scale = 0;
    for (double z : out.value)
      scale = std::max(scale, std::abs(z));
    double budget = std::max(rtol * scale + atol,
                             64 * std::numeric_limits<double>::denorm_min());
    bool good = true;
    for (double e : out.error)
      good &= e <= budget;
    if (good) {
      std::sort(out.panels.begin(), out.panels.end(),
                [](auto const& a, auto const& b) { return a.left < b.left; });
      return out;
    }
    double a = out.panels[worst].left, b = out.panels[worst].right;
    double m = std::midpoint(a, b);
    require(m > a && m < b, Failure::Convergence,
            std::format("QED quadrature: unresolved interval [{}, {}]", a, b));
    auto first = evaluate(a, m), second = evaluate(m, b);
    out.panels[worst] = first;
    out.panels.push_back(second);
  }
}

// Invert a scalar positive measure assembled from the *same* quadrature
// panels as its normalization. Plateaus have zero probability.
template <class F>
double
inverse_cdf(F const& f, Integral<1> const& total, double u, double rtol) {
  require(u > 0 && u < 1 && total.value[0] > 0, Failure::InvalidInput,
          "QED CDF: invalid probability or zero mass");
  double target = u * total.value[0];
  for (auto const& p : total.panels) {
    double mass = p.integral.value[0];
    require(mass >= 0, Failure::Convergence, "QED CDF: negative panel mass");
    if (target >= mass) {
      target -= mass;
      continue;
    }
    auto coordinate = [&](double beta) {
      return p.asinh_scale > 0
                 ? std::asinh((beta - p.asinh_center) / p.asinh_scale)
                 : beta;
    };
    auto physical = [&](double t) {
      return p.asinh_scale > 0 ? p.asinh_center + p.asinh_scale * std::sinh(t)
                               : t;
    };
    auto integrand = [&](double t) {
      return scaled<1>(f(physical(t)),
                       p.asinh_scale > 0 ? p.asinh_scale * std::cosh(t) : 1);
    };
    double start = coordinate(p.left), lo = start, hi = coordinate(p.right),
           flo = 0, fhi = mass;
    for (int k = 0; k < 64; ++k) {
      double fraction = std::clamp((target - flo) / (fhi - flo), 0.1, 0.9);
      double x = std::lerp(lo, hi, fraction);
      if (!(x > lo && x < hi))
        return physical(std::midpoint(lo, hi));
      auto part = integrate<1>(integrand, {start, x}, rtol, rtol * mass * 0.05);
      double val = part.value[0];
      require(val >= -part.error[0] &&
                  val <= mass + part.error[0] + p.integral.error[0],
              Failure::Convergence, "QED CDF: inconsistent partial integral");
      if (std::abs(val - target) <= rtol * mass + p.integral.error[0])
        return physical(x);
      if (val < target) {
        lo = x;
        flo = val;
      } else {
        hi = x;
        fhi = val;
      }
    }
    throw Error(Failure::Convergence, "QED CDF inversion did not converge");
  }
  throw Error(Failure::Convergence, "QED CDF: inconsistent total mass");
}

template <class F>
std::optional<double>
candidate_from_rate(F const& rate, std::vector<double> const& cuts, double& tau,
                    double rtol, double rate_atol, size_t limit) {
  for (size_t j = 1; j < cuts.size(); ++j) {
    auto integral =
        integrate<1>(rate, {cuts[j - 1], cuts[j]}, rtol,
                     2 * rate_atol * (cuts[j] - cuts[j - 1]), limit);
    if (tau + integral.value[0] > 0)
      return inverse_cdf(rate, integral, -tau / integral.value[0], rtol * 4);
    tau += integral.value[0];
  }
  return std::nullopt;
}
} // namespace qed
