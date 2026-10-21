#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

template <class F>
double
bisect_root(F const& func, double x1, double x2, double xacc, double f1,
            double f2) {
  if (xacc < 0)
    throw std::invalid_argument("xacc must be non-negative");
  if (std::isnan(f1) || std::isnan(f2))
    throw std::runtime_error("NaN endpoint function value");
  if ((f1 < 0 && f2 < 0) || (f1 > 0 && f2 > 0))
    throw std::runtime_error("Root must be bracketed");
  if (f1 == 0)
    return x1;
  if (f2 == 0)
    return x2;
  if (f1 > 0) {
    std::swap(x1, x2);
  }
  int i = 64;
  while (i--) {
    if (std::abs(x1 - x2) <= xacc)
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

template <class F>
double
ridders_root(F const& func, double x1, double x2, double xacc, double f1,
             double f2) {
  if (xacc < 0)
    throw std::invalid_argument("xacc must be non-negative");
  if (std::isnan(f1) || std::isnan(f2))
    throw std::runtime_error("NaN endpoint function value");
  if ((f1 < 0 && f2 < 0) || (f1 > 0 && f2 > 0))
    throw std::runtime_error("Root must be bracketed");
  if (f1 == 0)
    return x1;
  if (f2 == 0)
    return x2;

  double xl, xh, fl, fh;
  if (x1 < x2) {
    xl = x1;
    xh = x2;
    fl = f1;
    fh = f2;
  } else {
    xl = x2;
    xh = x1;
    fl = f2;
    fh = f1;
  }

  auto update_bracket = [&](double x, double fx) {
    if (std::signbit(fx) == std::signbit(fl)) {
      xl = x;
      fl = fx;
    } else {
      xh = x;
      fh = fx;
    }
  };
  int i = 64;
  while (i--) {
    if ((xh - xl) <= xacc)
      return std::midpoint(xl, xh);
    double xm = std::midpoint(xl, xh);
    if (xm == xl || xm == xh)
      return xm;

    double fm = func(xm);
    if (std::isnan(fm))
      throw std::runtime_error("NaN function value");
    if (fm == 0.0)
      return xm;

    double scale = std::max({std::abs(fl), std::abs(fh), std::abs(fm)});
    double a = fl / scale;
    double b = fh / scale;
    double c = fm / scale;
    double q = std::sqrt(std::abs(a)) * std::sqrt(std::abs(b));
    double s = std::hypot(c, q);

    double xn = xm;
    if (s > 0.0) {
      if (std::signbit(fl) != std::signbit(fm)) {
        xn = xl + (xm - xl) * (q / s) * (q / (s + std::abs(c)));
      } else {
        xn = xh - (xh - xm) * (q / s) * (q / (s + std::abs(c)));
      }
    }

    update_bracket(xm, fm);

    if (xl < xn && xn < xh) {
      double fn = func(xn);
      if (std::isnan(fn))
        throw std::runtime_error("NaN function value");
      if (fn == 0.0)
        return xn;
      update_bracket(xn, fn);
    }
  }
  throw std::runtime_error("ridders_root exceeded maximum iterations");
}
