#pragma once

#include <cstdint>

// A bijective 64-bit mixer assigns distinct initial seeds to photon IDs
// within a run. It does not consume a worker's transport random stream.
inline std::uint64_t mix_photon_seed(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

inline std::uint64_t photon_seed(std::uint64_t seed, std::uint64_t id) {
  return mix_photon_seed(mix_photon_seed(seed) ^ id);
}
