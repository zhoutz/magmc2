#include "bfield.hpp"
#include "component.hpp"
#include "evolution.hpp"
#include "resonant.hpp"
#include <cassert>
#include <cstdio>
#include <random>
#include <tuple>

int
main() {
  BField field("table/bfield_t10.txt", 1e14, 10);
  std::mt19937_64 rng(83126);
  std::uniform_real_distribution<double> uniform(0, 1);
  // Enclosures must hold over cells, both hemispheres, poles and wide ranges.
  for (int i = 0; i < 3000; ++i) {
    double a = 2 * uniform(rng) - 1, b = 2 * uniform(rng) - 1;
    if (a > b)
      std::swap(a, b);
    if (i % 3 == 0)
      b = std::min(1.0, a + 1e-5);
    if (i % 50 == 0) {
      a = -1;
      b = 1;
    }
    double r = 10 + 1e5 * uniform(rng), end = r * (1 + uniform(rng));
    auto bounds = field.strength_bounds({r, end}, {a, b});
    auto global = field.strength_bounds({r, end});
    for (int j = 0; j <= 16; ++j) {
      double mu = std::lerp(a, b, j / 16.0),
             radius = std::lerp(r, end, uniform(rng));
      double strength = field.evaluate(radius, mu).Bmag;
      assert(bounds.lo <= strength && strength <= bounds.hi);
      assert(global.lo <= strength && strength <= global.hi);
    }
  }
  // A vanishing current must not request expensive polarization interpolation.
  auto electron = CurrentComponent{NegateDist(PowerLaw2(0.2, 0.8)), -1, 0.5, 1};
  auto positron = CurrentComponent{PowerLaw2(0.2, 0.8), 1, 0.5, 1};
  ResonantScattering scattering(std::make_tuple(electron, positron));
  Basis b{};
  b.J = 0;
  bool called = false;
  auto pol = [&] {
    called = true;
    return Polarization();
  };
  assert(scattering.total_rate(scattering.rates_at(b, pol)) == 0 && !called);
  b.J = 1;
  b.Bmag = 1e14;
  b.energy_local = 1;
  b.mu = 0;
  b.c = 1;
  b.s = 0;
  assert(scattering.total_rate(scattering.rates_at(b, pol)) == 0 && !called);
  assert(!scattering.possible_energy({0, 0.1}));
  assert(!scattering.possible_energy({4, 100}));
  assert(scattering.possible_energy({0.5, 2}));
  // Excluded intervals include both inward and outward propagation. Probe the
  // exact trajectory to check that an exclusion never hides an active rate.
  InitFD11 init(10, 0.4);
  for (int i = 0; i < 40; ++i) {
    PhotonEvolution photon(1.4, 10, 42 + i, field, init, scattering);
    photon.init();
    if (i % 2)
      photon.geodesic.init(100, 0, 2.3);
    for (int k = 0; k < 100; ++k) {
      int event = photon.geodesic.step();
      if (!photon.scattering_possible()) {
        for (int j = 0; j <= 16; ++j) {
          auto basis = photon.basis(std::lerp(
              photon.geodesic.l_old(), photon.geodesic.l_new(), j / 16.0));
          assert(scattering.total_rate(
                     scattering.rates(basis, Polarization())) == 0);
        }
      }
      if (event >= 0)
        break;
      photon.geodesic.stepper.update_old();
    }
  }
  // Endpoint evolution must be bitwise identical with dense output disabled.
  StepperMagnus dense, sparse;
  dense.init(0, 1, 0.01, {0, 0, 1});
  sparse.init(0, 1, 0.01, {0, 0, 1});
  auto omega = [](double x) { return double3{100 + 3 * x, 2 * x, 0}; };
  while (dense.do_step(omega)) {
    assert(sparse.do_step(omega, false));
    assert(dense.x_new == sparse.x_new);
    assert(dense.y_new.x == sparse.y_new.x && dense.y_new.y == sparse.y_new.y &&
           dense.y_new.z == sparse.y_new.z);
    bool threw = false;
    try {
      sparse.interp(sparse.x_new);
    } catch (std::logic_error const&) {
      threw = true;
    }
    assert(threw);
    dense.update_old();
    sparse.update_old();
  }
  assert(!sparse.do_step(omega, false));
  // Coalescence and distribution cutoffs stay on the reference integration
  // path. Exercise a synthetic crossing and optical-depth event selection.
  auto basis_at = [&](double l) {
    Basis v{};
    v.Bmag = (0.5 + l) / B_to_energy;
    v.energy_local = 1;
    v.J = v.Bmag / (4 * pi * pi);
    v.mu = 0.6;
    v.c = 1;
    v.s = 0;
    return v;
  };
  auto cuts = scattering.integration_cuts(0, 1, basis_at);
  assert(std::any_of(cuts.begin(), cuts.end(),
                     [](double l) { return std::abs(l - 0.3) < 1e-12; }));
  Quad quad;
  Ran ran(123);
  double tau = -1e-3;
  auto result = scattering.integrate_optical_depth(
      0, 1, basis_at, [](double) { return Polarization(); }, tau, ran, quad);
  assert(result && std::isfinite(result->beta));
  std::puts("Transport bounds, exclusions, dense output and resonance "
            "boundaries passed");
}
