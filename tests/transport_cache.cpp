#include "transport_cache.hpp"
#include "photon_stream.hpp"
#include <iostream>
#include <set>
#include <stdexcept>

void require(bool ok, char const* message) {
  if (!ok) throw std::runtime_error(message);
}

int main() {
  TransportCache cache;
  cache.tolerance = 1e-10;
  int calls = 0;
  auto polynomial = [&](double x) {
    ++calls;
    double p = 1 + x * (2 + x * (3 + x * (4 + x * 5)));
    return TransportCache::Values{p, 0.1 * p, p, p, x - 0.5, p, p};
  };
  cache.build(0, 1, polynomial);
  require(cache.leaves.size() == 1, "Quartic must fit without subdivision");
  for (int i = 0; i <= 1000; ++i) {
    double x = i / 1000.0;
    int before = calls;
    auto y = cache.at(x, polynomial);
    require(calls == before, "Accepted cache evaluated expensive callback");
    auto exact = polynomial(x);
    for (int k = 0; k < 7; ++k)
      require(std::abs(y[k] - exact[k]) < 1e-12, "Quartic interpolation error");
  }
  auto smooth = [](double x) {
    double v = std::exp(x) + std::sin(3 * x);
    return TransportCache::Values{v, v, v, v, v, v, v};
  };
  cache.leaves.clear();
  cache.build(0, 1, smooth);
  require(cache.leaves.size() > 1, "Non-polynomial needs adaptive subdivision");
  for (int i = 0; i <= 10000; ++i) {
    double x = i / 10000.0;
    require(std::abs(cache.at(x, smooth)[0] - smooth(x)[0]) < 1e-9,
            "Smooth interpolation exceeds sampled error budget");
  }
  auto difficult = [](double x) {
    double v = std::sin(2000000 * x);
    return TransportCache::Values{v, v, v, v, v, v, v};
  };
  cache.leaves.clear();
  cache.build(0, 1, difficult);
  require(std::any_of(cache.leaves.begin(), cache.leaves.end(),
                     [](auto const& leaf) { return leaf.exact; }), "Missing fallback");
  for (auto const& leaf : cache.leaves)
    if (leaf.exact) {
      double x = std::lerp(leaf.a, leaf.b, 0.37);
      require(cache.at(x, difficult) == difficult(x), "Fallback is not exact");
    }
  require(cache.at(-1, smooth) == smooth(-1), "Outside cache must evaluate exactly");
  std::set<std::uint64_t> seeds;
  for (int id = 0; id < 10000; ++id) seeds.insert(photon_seed(110701, id));
  require(seeds.size() == 10000, "Repeated per-photon seed");
  std::cout << "Transport interpolation, fallback and photon seeds passed\n";
}
