#pragma once

#ifdef MAGMC_DIAGNOSTICS
#include <chrono>
#include <cstdio>
namespace diagnostic {
inline unsigned long long basis = 0, field = 0, scans = 0, roots = 0,
                          quad_samples = 0, magnus_accept = 0,
                          magnus_reject = 0, orbit_steps = 0;
inline double optical_seconds = 0, polarization_seconds = 0;
struct Timer {
  double& total;
  std::chrono::steady_clock::time_point start =
      std::chrono::steady_clock::now();
  explicit Timer(double& total) : total(total) {}
  ~Timer() {
    total +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
  }
};
inline
void
report() {
  std::fprintf(
      stderr,
      "{\"basis\":%llu,\"field\":%llu,\"scans\":%llu,\"roots\":%llu,\"quad_"
      "samples\":%llu,\"magnus_accept\":%llu,\"magnus_reject\":%llu,\"orbit_"
      "steps\":%llu,\"optical_seconds\":%.9g,\"polarization_seconds\":%.9g}\n",
      basis, field, scans, roots, quad_samples, magnus_accept, magnus_reject,
      orbit_steps, optical_seconds, polarization_seconds);
}
} // namespace diagnostic
#define MC_COUNT(name) (++diagnostic::name)
#define MC_TIMER(name) diagnostic::Timer timer_##name(diagnostic::name)
#else
#define MC_COUNT(name) ((void)0)
#define MC_TIMER(name) ((void)0)
#endif
