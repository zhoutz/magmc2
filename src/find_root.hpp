#pragma once

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

template <class T>
double
bisect_root(T const& func, double x1, double x2, double xacc, double f1,
            double f2) {
  if (xacc < 0)
    throw std::invalid_argument("xacc must be non-negative");
  if (std::isnan(f1) || std::isnan(f2))
    throw std::runtime_error("NaN endpoint function value");
  if ((f1 < 0 && f2 < 0) || (f1 > 0 && f2 > 0))
    throw std::runtime_error("Root must be bracketed for bisection");
  if (f1 == 0)
    return x1;
  if (f2 == 0)
    return x2;
  if (f1 > 0) {
    std::swap(x1, x2);
  }
  int i = 64;
  while (i--) {
    if (std::abs(x1 - x2) < xacc)
      return std::midpoint(x1, x2);
    double xmid = std::midpoint(x1, x2);
    if (xmid == x1 || xmid == x2)
      return xmid;

    double fmid = func(xmid);
    if (std::isnan(fmid))
      throw std::runtime_error("NaN function value");

    if (fmid == 0)
      return xmid;
    if (fmid < 0)
      x1 = xmid;
    else
      x2 = xmid;
  }
  throw std::runtime_error("Too many bisections in bisect_root");
}
