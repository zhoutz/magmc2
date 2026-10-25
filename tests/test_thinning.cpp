#include "component.hpp"
#include "resonant.hpp"
#include "rotation_cache.hpp"
#include <cassert>
#include <iostream>

int
main() {
  auto comp = CurrentComponent{.dist = PowerLaw2(0.2, 0.8),
                               .charge_sign = 1,
                               .current_fraction = 1,
                               .charge_mass_ratio = 1};
  ResonantScattering scattering(std::make_tuple(comp));
  Ran ran(98765);
  Basis b{};
  b.Bmag = 1 / B_to_energy;
  b.energy_local = 1;
  b.J = b.Bmag / (4 * pi * pi);
  b.c = 1;
  b.s = 0;
  b.mu = 0.6;
  for (int i = 0; i < 10000; ++i) {
    b.mu = ran.U(-1, 1);
    b.Bmag = ran.U(0.1, 3) / B_to_energy;
    auto pol = Polarization(ran.point_on_unit_sphere());
    double actual = scattering.total_rate(scattering.rates(b, pol));
    double upper = scattering.envelope_rate(b);
    assert(actual >= 0 && actual <= upper * (1 + 1e-12));
  }
  b.mu = 0.6;
  b.Bmag = 0.85 / B_to_energy;
  Polarization pol(Stokes{0, 1, 0});
  double lambda = scattering.total_rate(scattering.rates(b, pol));
  double bound = scattering.envelope_rate(b);
  assert(lambda > 0 && bound > lambda);
  Quad quad;
  double sum = 0;
  int survived = 0;
  constexpr int n = 20000;
  for (int i = 0; i < n; ++i) {
    double x = 0;
    for (;;) {
      double target = std::log(ran.U());
      double right = x + 50 / bound;
      std::vector<double> cuts{x, right};
      auto get_basis = [&](double) { return b; };
      auto candidate =
          scattering.next_candidate(x, right, cuts, get_basis, target, quad);
      assert(candidate);
      assert(std::abs((*candidate - x) * bound + target) < 1e-5);
      x = *candidate;
      if (scattering.accept_candidate(b, pol, ran))
        break;
    }
    sum += x * lambda;
    survived += x * lambda > 1;
  }
  assert(std::abs(sum / n - 1) < 0.035);
  assert(std::abs(double(survived) / n - std::exp(-1.0)) < 0.02);
  b.J = 0;
  double target = -1;
  std::vector<double> cuts{0, 1};
  assert(!scattering.next_candidate(
      0, 1, cuts, [&](double) { return b; }, target, quad));
  assert(target == -1);
  // Identical populations must be selected in their current-fraction ratio.
  auto first = comp, second = comp;
  first.current_fraction = 0.25;
  second.current_fraction = 0.75;
  ResonantScattering mixed(std::make_tuple(first, second));
  b.J = b.Bmag / (4 * pi * pi);
  auto weights = mixed.rates(b, pol);
  double total = mixed.total_rate(weights);
  std::array<int, 4> observed{};
  for (int i = 0; i < 100000; ++i) {
    auto r = mixed.select_rate(weights, ran);
    ++observed[2 * r.icomp + r.jroot];
  }
  for (int i = 0; i < 4; ++i)
    assert(std::abs(observed[i] / 100000.0 - weights[i].weight / total) < 0.01);
  // Full-support distributions cannot use the finite-support exclusion.
  auto thermal = CurrentComponent{.dist = Boltzmann(0.3),
                                  .charge_sign = 1,
                                  .current_fraction = 1,
                                  .charge_mass_ratio = 1};
  ResonantScattering thermal_scattering(std::make_tuple(thermal));
  assert(thermal_scattering.envelope_rate(b) > 0);
  // Root coalescence and the two sides of a velocity cutoff.
  b.mu = 0.6;
  b.Bmag = 0.8 / B_to_energy;
  auto merged = scattering.resonance_roots(0.8, 0.6, 0);
  assert(std::abs(merged[0].beta - 0.6) < 1e-14);
  assert(std::abs(merged[1].beta - 0.6) < 1e-14);
  assert(comp.dist.f(0.2) == 0 && comp.dist.f(std::nextafter(0.2, 1.0)) > 0);
  assert(comp.dist.f(0.8) == 0 && comp.dist.f(std::nextafter(0.8, 0.0)) > 0);
  // Resolve an integrable coalescence singularity, including candidates
  // near its endpoint, against independent tight integration.
  auto varying = [&](double l) {
    Basis local = b;
    local.Bmag = std::sqrt(0.64 + l) / B_to_energy;
    return local;
  };
  auto envelope = [&](double l) {
    return scattering.envelope_rate(varying(l));
  };
  auto singular_cuts = scattering.integration_cuts(0, 0.04, varying, true);
  double mass = quad.qags(envelope, 0, 0.04, 1e-10, 1e-10);
  assert(mass > 0);
  for (double fraction : {0.01, 0.2, 0.8}) {
    double threshold = -fraction * mass;
    auto x = scattering.next_candidate(0, 0.04, singular_cuts, varying,
                                       threshold, quad);
    assert(x && *x > 0 && *x < 0.04);
    double measured = quad.qags(envelope, 0, *x, 1e-10, 1e-10);
    assert(std::abs(measured - fraction * mass) < 1e-4 * mass);
  }
  StepperMagnus mag;
  mag.init(0, 10, 10, Stokes{0, 1, 0});
  auto omega = [](double) { return double3{100, 0, 0}; };
  while (mag.do_step(omega, false))
    mag.update_old();
  assert(std::abs(mag.y_old.x) < 1e-12);
  assert(std::abs(mag.y_old.y - std::cos(1000.0)) < 1e-11);
  assert(std::abs(mag.y_old.z - std::sin(1000.0)) < 1e-11);
  assert(std::abs(mag.y_old.length() - 1) < 1e-12);
  bool threw = false;
  try {
    mag.interp(5);
  } catch (std::logic_error const&) {
    threw = true;
  }
  assert(threw);
  RotationCache cache;
  auto polynomial = [](double x) {
    return double3{1 + x * x, x * x * x * x, 0};
  };
  cache.build(0, 1, polynomial);
  for (int i = 0; i < 100; ++i) {
    double x = ran.U();
    assert((cache.at(x, polynomial) - polynomial(x)).length() < 1e-12);
  }
  Geodesic inward(1.4, 10);
  inward.init(100, 0, 2.0);
  bool turned = false;
  for (int i = 0; i < 150; ++i) {
    int event = inward.step();
    assert(event != 0);
    turned |= std::cos(inward.stepper.y_new[2]) > 0;
    inward.stepper.update_old();
  }
  assert(turned);
  BField field("table/bfield_t10.txt", 1e14, 10);
  assert(field.angular_lower > 0 && std::isfinite(field.angular_upper));
  InitFD11 initializer(10, 0.4);
  for (int seed = 0; seed < 30; ++seed) {
    PhotonEvolution ray(1.4, 10, seed, field, initializer, scattering);
    ray.init();
    for (int step = 0; step < 100; ++step) {
      ray.geodesic.step();
      if (!ray.possible_resonance()) {
        for (int j = 0; j <= 16; ++j) {
          double l =
              std::lerp(ray.geodesic.l_old(), ray.geodesic.l_new(), j / 16.0);
          assert(scattering.envelope_rate(ray.basis(l)) == 0);
        }
      }
      ray.geodesic.stepper.update_old();
    }
  }
  for (int i = 0; i < 10000; ++i) {
    double mu = double(i) / 9999;
    auto f = field.evaluate(10, mu);
    assert(f.Bmag >= 0.5e14 * field.angular_lower);
    assert(f.Bmag <= 0.5e14 * field.angular_upper);
  }
  Geodesic geo(1.4, 10);
  geo.init(10, 0, 0.9);
  double impact = 10 * std::sin(0.9) / std::sqrt(1 - geo.rs / 10);
  for (int i = 0; i < 100; ++i) {
    geo.step();
    auto y = geo.stepper.y_new;
    double actual = y[0] * std::sin(y[2]) / std::sqrt(1 - geo.rs / y[0]);
    assert(std::abs(actual / impact - 1) < 1e-7);
    geo.stepper.update_old();
  }
  std::cout
      << "thinning, polarization, cache, field bounds, geodesic: passed\n";
}
