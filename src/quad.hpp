#pragma once

#include <algorithm>
#include <cmath>
#include <exception>
#include <format>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <limits>
#include <stdexcept>

struct Quad {
  static constexpr size_t limit = 8192;
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);
  Quad() {
    if (!workspace)
      throw std::bad_alloc();
    gsl_set_error_handler_off();
  }
  Quad(Quad const&) = delete;
  Quad& operator=(Quad const&) = delete;
  ~Quad() { gsl_integration_workspace_free(workspace); }

  template <class F>
  double
  qags(F& fn, double a, double b, double epsabs, double epsrel) {
    if (a == b)
      return 0;
    // Never unwind a C++ exception through the GSL callback.
    struct Context {
      F& fn;
      std::exception_ptr error;
    } context{fn, {}};
    gsl_function gf;
    gf.function = +[](double x, void* p) {
      auto& ctx = *static_cast<Context*>(p);
      if (!ctx.error) {
        try {
          return ctx.fn(x);
        } catch (...) {
          ctx.error = std::current_exception();
        }
      }
      return std::numeric_limits<double>::quiet_NaN();
    };
    gf.params = &context;
    double result, abserr;
    int status = gsl_integration_qags(&gf, a, b, epsabs, epsrel, limit,
                                      workspace, &result, &abserr);
    if (context.error)
      std::rethrow_exception(context.error);
    double tolerance = std::max(epsabs, epsrel * std::abs(result));
    if (!std::isfinite(result) || !std::isfinite(abserr) ||
        status != GSL_SUCCESS || abserr > tolerance)
      throw std::runtime_error(std::format(
          "qags failed on [{},{}]: {} (result={}, error={}, tolerance={})", a,
          b, gsl_strerror(status), result, abserr, tolerance));
    return result;
  }
};
