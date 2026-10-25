#pragma once
#include <cstdint>

struct TransportStats {
  std::uint64_t basis = 0, cuts = 0, roots = 0, quadrature = 0, integrand = 0,
                polarization = 0, excluded = 0, magnus = 0;
};
#ifdef TRANSPORT_STATS
inline TransportStats transport_stats;
#define TRANSPORT_COUNT(name) (++transport_stats.name)
#else
#define TRANSPORT_COUNT(name) ((void)0)
#endif
