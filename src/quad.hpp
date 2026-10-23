#pragma once

#include <cmath>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <print>

struct Quad {
  static constexpr size_t limit = 8192;
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);
  Quad() { gsl_set_error_handler_off(); }
  ~Quad() { gsl_integration_workspace_free(workspace); }
  Quad(Quad const&) = delete;

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
      std::println(
          stderr, "qags failed: {} (status={}, result={}, abserr={})\n",
          status == GSL_SUCCESS ? "non-finite output" : gsl_strerror(status),
          status, result, abserr);
      std::abort();
    }
    return result;
  }
};
