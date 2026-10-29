#pragma once
#include "diagnostics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>

struct QuadratureEstimate {
  double integral, error;
};

// Embedded Gauss-7/Kronrod-15 rule. The error scaling follows QUADPACK's
// resasc/resabs safeguards. This is only a fast first attempt: difficult
// intervals still use the existing adaptive QAGS integrator.
template <class F>
QuadratureEstimate
gauss_kronrod15(F const& fn, double a, double b) {
  if (a == b)
    return {0, 0};
  constexpr std::array<double, 8> x{.9914553711208126392, .9491079123427585245,
                                    .8648644233597690728, .7415311855993944399,
                                    .5860872354676911303, .4058451513773971669,
                                    .2077849550078984676, 0};
  constexpr std::array<double, 8> wk{
      .02293532201052922496, .06309209262997855329, .1047900103222501838,
      .1406532597155259187,  .1690047266392679028,  .1903505780647854099,
      .2044329400752988924,  .2094821410847278280};
  constexpr std::array<double, 4> wg{.1294849661688696933, .2797053914892766679,
                                     .3818300505051189449,
                                     .4179591836734693878};
  double mid = std::midpoint(a, b), half = .5 * (b - a);
  MC_COUNT(quad_samples);
  double center = fn(mid), kronrod = wk[7] * center, gauss = wg[3] * center;
  double abs_sum = wk[7] * std::abs(center);
  std::array<double, 7> low, high;
  for (int i = 0; i < 7; ++i) {
    MC_COUNT(quad_samples);
    MC_COUNT(quad_samples);
    low[i] = fn(mid - half * x[i]);
    high[i] = fn(mid + half * x[i]);
    double sum = low[i] + high[i];
    kronrod += wk[i] * sum;
    abs_sum += wk[i] * (std::abs(low[i]) + std::abs(high[i]));
    if (i % 2 == 1)
      gauss += wg[i / 2] * sum;
  }
  double average = .5 * kronrod;
  double asc = wk[7] * std::abs(center - average);
  for (int i = 0; i < 7; ++i)
    asc += wk[i] * (std::abs(low[i] - average) + std::abs(high[i] - average));
  asc *= std::abs(half);
  abs_sum *= std::abs(half);
  double error = std::abs(half * (kronrod - gauss));
  if (asc > 0 && error > 0)
    error = asc * std::min(1., std::pow(200 * error / asc, 1.5));
  error =
      std::max(error, 50 * std::numeric_limits<double>::epsilon() * abs_sum);
  return {half * kronrod, error};
}
