#include "bfield.hpp"
#include "component.hpp"
#include "evolution.hpp"
#include "resonant.hpp"
#include <charconv>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <tuple>

int
main(int argc, char** argv) {
  auto integer = [](char const* text) {
    int value;
    auto end = text + std::string_view(text).size();
    auto [ptr, err] = std::from_chars(text, end, value);
    if (err != std::errc{} || ptr != end || value < 0)
      throw std::invalid_argument("Expected a nonnegative integer");
    return value;
  };
  int first = argc > 1 ? integer(argv[1]) : 42;
  int count = argc > 2 ? integer(argv[2]) : 2000;
  if (count > std::numeric_limits<int>::max() - first)
    throw std::invalid_argument("Seed range overflows int");
  BField field("table/bfield_t10.txt", 1e14, 10);
  InitFD11 init(10, 0.4);
  auto electron = CurrentComponent{NegateDist(PowerLaw2(0.2, 0.8)), -1, 0.5, 1};
  auto positron = CurrentComponent{PowerLaw2(0.2, 0.8), 1, 0.5, 1};
  ResonantScattering scattering(std::make_tuple(electron, positron));
#ifdef TRANSPORT_CURRENT
  using Evolution = PhotonEvolution<InitFD11, decltype(scattering)>;
  Evolution::Options options;
  options.fast_transport =
      !(argc > 3 && std::string_view(argv[3]) == "--reference");
#endif
  for (int i = 0; i < count; ++i) {
    auto start = std::chrono::steady_clock::now();
    PhotonEvolution photon(1.4, 10, first + i, field, init, scattering
#ifdef TRANSPORT_CURRENT
                           ,
                           options
#endif
    );
    photon.init();
    int scatters = 0;
    for (;;) {
      auto result = photon.evolve();
      if (result == decltype(photon)::EvolveResult::Scattered) {
        ++scatters;
        photon.perform_scattering();
      } else {
        if (result == decltype(photon)::EvolveResult::Escaped) {
          auto const& e = photon.escape_data;
          std::printf("%d %d E %.17g %.17g %.17g %.17g %.17g\n", first + i,
                      scatters, e.energy_inf, e.muk, e.Q, e.U, e.V);
        } else
          std::printf("%d %d A\n", first + i, scatters);
        break;
      }
    }
#ifdef TRANSPORT_STATS
    auto t =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::fprintf(stderr, "%d %.9g %d\n", first + i, t, scatters);
#else
    (void)start;
#endif
  }
#ifdef TRANSPORT_STATS
  auto const& s = transport_stats;
  std::fprintf(stderr,
               "basis=%llu cuts=%llu roots=%llu quadrature=%llu integrand=%llu "
               "polarization=%llu excluded=%llu magnus=%llu\n",
               (unsigned long long)s.basis, (unsigned long long)s.cuts,
               (unsigned long long)s.roots, (unsigned long long)s.quadrature,
               (unsigned long long)s.integrand,
               (unsigned long long)s.polarization,
               (unsigned long long)s.excluded, (unsigned long long)s.magnus);
#endif
}
