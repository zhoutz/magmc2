#pragma once

#include <cmath>
#include <cstdlib>
#include <format>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

struct Quad {
  static constexpr size_t limit = 8192;
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);
  Quad() { gsl_set_error_handler_off(); }
  ~Quad() { gsl_integration_workspace_free(workspace); }

  template <class F>
  double
  smooth(F& fn, double a, double b, double epsabs, double epsrel) {
    gsl_function gf;
    gf.function = +[](double x, void* p) { return (*static_cast<F*>(p))(x); };
    gf.params = &fn;
    double result, error, absolute, asc;
    gsl_integration_qk15(&gf, a, b, &result, &error, &absolute, &asc);
    if (std::isfinite(result) && std::isfinite(error) &&
        error <= std::max(epsabs, epsrel * std::abs(result)) &&
        (error == 0 || error != asc))
      return result;
    return qags(fn, a, b, epsabs, epsrel);
  }

  template <class F>
  double
  qags(F& fn, double a, double b, double epsabs, double epsrel) {
    gsl_function gf;
    gf.function = +[](double x, void* p) { return (*static_cast<F*>(p))(x); };
    gf.params = &fn;
    double result, abserr;
    int status = gsl_integration_qags(&gf, a, b, epsabs, epsrel, limit,
                                      workspace, &result, &abserr);

    if (!std::isfinite(result) || !std::isfinite(abserr)) {
      throw std::runtime_error(std::format(
          "qags failed: {} (status={}, result={}, abserr={})\n",
          status == GSL_SUCCESS ? "non-finite output" : gsl_strerror(status),
          status, result, abserr));
    }
    return result;
  }

  Quad(Quad const&) = delete;
  Quad& operator=(Quad const&) = delete;
};
