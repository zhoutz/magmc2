#include "resonant.hpp"
#include "component.hpp"
#include <cassert>
#include <iostream>
#include <limits>

namespace {
// A symmetric distribution makes both resonance branches observable.
struct UniformDistribution {
  double b_min = -0.9, b_max = 0.9, inv_abs_b_mean = 1 / 0.45;
  double f(double beta) const {
    return beta > b_min && beta < b_max ? 1 / 1.8 : 0;
  }
};

// Local inputs deliberately independent of PhotonEvolution.
struct Basis {
  double r = 20, Bmag = 1, mu = 0, J = 0.5, c = 1, s = 0;
  double3 r_hat{1, 0, 0}, B{0, 0, 1};
};

bool close(double a, double b) {
  return std::abs(a - b) <= 1e-11 * std::max({1.0, std::abs(a), std::abs(b)});
}

void test_multiple_components() {
  auto electron = CurrentComponent{NegateDist(PowerLaw2(0.2, 0.9)), -1, 0.25, 1.0};
  auto positron = CurrentComponent{Boltzmann(0.4), 1, 0.75, 0.75};
  ResonantScattering mixture(std::make_tuple(electron, positron));
  ResonantScattering electrons(std::make_tuple(electron));
  ResonantScattering positrons(std::make_tuple(positron));
  using Mixture = decltype(mixture);
  static_assert(std::tuple_size_v<Mixture::Rates> == 4);
  Mixture::Options options;
  options.n_knots = 3;
  options.n_scan = 16;
  ResonantScattering configured(std::make_tuple(electron, positron), options);
  options.n_scan = 4;
  assert(configured.options.n_scan == 16); // The model owns its configuration.
  assert(configured.knots[0].size() == 3 && configured.knots[1].size() == 3);
  Basis b;
  Polarization pol(AdiabaticPolarization::E);
  double energy = B_to_omega / 2;
  auto rates = mixture.rates(b, pol, energy, 0);
  auto er = electrons.rates(b, pol, energy, 0);
  auto pr = positrons.rates(b, pol, energy, 0);
  assert(rates[0].weight > 0 && rates[0].beta < 0);
  assert(rates[1].weight == 0 && rates[2].weight == 0);
  assert(rates[3].weight > 0 && rates[3].beta > 0);
  assert(!close(std::abs(rates[0].beta), rates[3].beta));
  for (std::size_t i = 0; i < rates.size(); ++i) {
    assert(rates[i].icomp == i / 2 && rates[i].iroot == int(i % 2));
    assert(close(rates[i].weight, i < 2 ? er[i].weight : pr[i - 2].weight));
  }
  double total = Mixture::total_rate(rates);
  assert(close(total, electrons.total_rate(er) + positrons.total_rate(pr)));
  std::get<0>(electrons.components).current_fraction *= 2;
  assert(close(electrons.rates(b, pol, energy, 0)[0].weight, 2 * er[0].weight));
  std::get<0>(electrons.components).current_fraction /= 2;
  assert(mixture.knots[0].front() == -0.9 && mixture.knots[0].back() == -0.2);
  assert(mixture.knots[1].front() == 0 && mixture.knots[1].back() == 1);

  // Recover the original self-similar prefactor without a field model.
  double p = 0.8, twist = 0.5;
  b.J = (p + 1) * twist * b.Bmag / (4 * pi * b.r);
  auto self_similar = mixture.rates(b, pol, energy, 0);
  double beta = self_similar[0].beta;
  double jacobian = std::pow(1 - beta * b.mu, 2) * (1 - beta * beta) /
                    std::abs(b.mu - beta);
  double old_rate = electron.current_fraction * (p + 1) * pi * twist /
                    b.r * electron.dist.inv_abs_b_mean * electron.dist.f(beta) *
                    jacobian * 0.5;
  assert(close(self_similar[0].weight, old_rate));
  b.J = 0.5;

  // The charge sign changes only the circular overlap, not velocity support.
  auto reversed_charge = electron;
  reversed_charge.charge_sign = 1;
  ResonantScattering reversed(std::make_tuple(reversed_charge));
  Polarization circular(Stokes{0, 0, 1});
  double wm = electrons.rates(b, circular, energy, 0)[0].weight;
  double wp = reversed.rates(b, circular, energy, 0)[0].weight;
  assert(wp > wm);
  assert(close(reversed.rates(b, pol, energy, 0)[0].weight, er[0].weight));
  circular.set_stokes({0, 0, -1});
  assert(close(electrons.rates(b, circular, energy, 0)[0].weight, wp));

  // The joint draw consumes one uniform variate, including across species.
  Ran ran(110701), reference(110701);
  int n_electron = 0;
  constexpr int n_samples = 20000;
  for (int i = 0; i < n_samples; ++i) {
    auto chosen = Mixture::select_rate(rates, ran);
    bool choose_electron = reference.U() * total < rates[0].weight;
    assert(chosen.icomp == (choose_electron ? 0u : 1u));
    assert(chosen.iroot == (choose_electron ? 0 : 1));
    n_electron += choose_electron;
  }
  assert(ran.uint64() == reference.uint64());
  assert(std::abs(double(n_electron) / n_samples - rates[0].weight / total) < 0.02);

  Quad quad;
  double tau = -0.75 * total;
  auto event = mixture.integrate_optical_depth(
      0, 1, [&](double) { return b; }, [&](double) { return pol; },
      energy, 0, tau, ran, quad);
  assert(event && close(event->l, 0.75));
  auto selected = rates[2 * event->icomp + event->iroot];
  assert(selected.weight > 0 && selected.beta == event->beta);

  // Each species contributes its own resonance and support crossings.
  auto varying_basis = [&](double l) {
    Basis point = b;
    point.Bmag = 0.5 + 2 * l;
    return point;
  };
  auto cuts = mixture.integration_cuts(0, 1, varying_basis, B_to_omega, 0);
  auto expected = electrons.integration_cuts(0, 1, varying_basis, B_to_omega, 0);
  auto pcuts = positrons.integration_cuts(0, 1, varying_basis, B_to_omega, 0);
  expected.insert(expected.end(), pcuts.begin(), pcuts.end());
  std::sort(expected.begin(), expected.end());
  expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
  assert(cuts == expected);
  auto contains = [&](double l) {
    return std::any_of(cuts.begin(), cuts.end(), [&](double c) { return close(c, l); });
  };
  assert(contains(0.25)); // x_e=1
  assert(contains((1 / positron.charge_mass_ratio - 0.5) / 2));
  for (double v : mixture.knots[0])
    assert(contains((1 / std::sqrt(1 - v * v) - 0.5) / 2));
  std::get<1>(mixture.components).current_fraction = 0;
  assert(mixture.integration_cuts(0, 1, varying_basis, B_to_omega, 0) ==
         electrons.integration_cuts(0, 1, varying_basis, B_to_omega, 0));
  assert(close(Mixture::total_rate(mixture.rates(b, pol, energy, 0)),
               electrons.total_rate(er)));

  // Outward ray: species 0 has escaped but species 1 still has resonances.
  Geodesic geodesic(0, 1);
  geodesic.init(20, 0, 0);
  geodesic.step();
  double l = geodesic.l_new();
  double escape_energy = 2 * B_to_omega / geodesic.interp(l)[0];
  auto escape_basis = [&](Geodesic::YVector y) {
    Basis point;
    point.r = y[0];
    point.Bmag = 1 / y[0];
    point.J = 0; // Local zero current alone must not establish escape.
    return point;
  };
  assert(electrons.resonance_escape(l, geodesic, escape_basis, escape_energy));
  std::get<1>(mixture.components).charge_mass_ratio = 3;
  std::get<1>(mixture.components).current_fraction = 0.75;
  assert(!mixture.resonance_escape(l, geodesic, escape_basis, escape_energy));
  std::get<1>(mixture.components).current_fraction = 0;
  assert(mixture.resonance_escape(l, geodesic, escape_basis, escape_energy));

  b.J = 0;
  assert(Mixture::total_rate(mixture.rates(b, pol, energy, 0)) == 0);
  for (double invalid : {-1.0, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
    b.J = invalid;
    bool rejected = false;
    try { mixture.rates(b, pol, energy, 0); }
    catch (std::runtime_error const&) { rejected = true; }
    assert(rejected);
  }
  b.J = 0.5;
  b.Bmag = 0;
  bool rejected = false;
  try { mixture.rates(b, pol, energy, 0); }
  catch (std::runtime_error const&) { rejected = true; }
  assert(rejected);
}
} // namespace

int main() {
  CurrentComponent component{UniformDistribution{}, -1, 1.0, 1.0};
  ResonantScattering scattering(std::make_tuple(component));
  using Scattering = decltype(scattering);
  Basis b;
  Polarization pol(AdiabaticPolarization::E);
  // With rs=0 and Bmag=1, these energies give x=2, 1/2 and 1 exactly.
  double energy = B_to_omega / 2;
  auto rates = scattering.rates(b, pol, energy, 0);
  assert(rates[0].weight > 0 && rates[1].weight > 0);
  assert(close(rates[0].beta, -std::sqrt(0.75)));
  assert(close(rates[1].beta, std::sqrt(0.75)));
  assert(close(rates[0].weight, rates[1].weight));

  std::get<0>(scattering.components).current_fraction = 0;
  assert(Scattering::total_rate(scattering.rates(b, pol, energy, 0)) == 0);
  std::get<0>(scattering.components).current_fraction = 1;
  b.J = 0;
  assert(Scattering::total_rate(scattering.rates(b, pol, energy, 0)) == 0);
  b.J = 0.5;
  assert(Scattering::discriminant(0.5, 0) < 0);
  assert(Scattering::total_rate(scattering.rates(b, pol, 2 * B_to_omega, 0)) == 0);
  assert(Scattering::discriminant(1, 0) == 0);
  assert(Scattering::total_rate(scattering.rates(b, pol, B_to_omega, 0)) == 0);
  auto near = scattering.rates(b, pol, B_to_omega / (1 + 1e-8), 0);
  assert(std::isfinite(Scattering::total_rate(near)));
  assert(Scattering::total_rate(near) > 0);

  bool rejected = false;
  Ran ran(42), untouched(42);
  try {
    Scattering::select_rate({}, ran);
  } catch (std::runtime_error const&) {
    rejected = true;
  }
  assert(rejected);
  assert(ran.uint64() == untouched.uint64());

  Quad quad;
  auto basis_at = [&](double) { return b; };
  auto pol_at = [&](double) { return pol; };
  double rate = Scattering::total_rate(rates);
  double tau = -2 * rate;
  auto event = scattering.integrate_optical_depth(
      0, 1, basis_at, pol_at, energy, 0, tau, ran, quad);
  assert(!event && close(tau, -rate));
  assert(ran.uint64() == untouched.uint64());

  // Constant rate: the event must occur at l=3/4, in the second segment.
  // Integration preserves that segment's starting threshold until reinit.
  pol.set_stokes({0, 0, 1});
  rate = Scattering::total_rate(scattering.rates(b, pol, energy, 0));
  tau = -0.75 * rate;
  event = scattering.integrate_optical_depth(
      0, 1, basis_at, pol_at, energy, 0, tau, ran, quad);
  assert(event && close(event->l, 0.75));
  assert(event->icomp == 0 && event->iroot >= 0 && event->iroot < 2);
  assert(close(tau, -0.25 * rate));
  assert(!event->incoming.is_adiabatic && event->incoming.stokes.z == 1);
  assert(close(std::abs(event->beta), std::sqrt(0.75)));

  for (double beta : {-0.8, 0.0, 0.8}) {
    for (int i = 0; i < 100; ++i) {
      auto photon = scattering.perform_scattering(b, beta, energy, ran);
      assert(std::isfinite(photon.omega_inf) && photon.omega_inf > 0);
      assert(photon.r0 == b.r && photon.psi0 == 0);
      assert(std::isfinite(photon.alpha0));
      assert(close(photon.n.length(), 1) && close(photon.e1.length(), 1));
      assert(close(photon.e2.length(), 1));
      assert(close(dot(photon.n, photon.e1), 0));
      assert(close(dot(cross(photon.e1, photon.e2), photon.n), 1));
      double3 k = std::cos(photon.alpha0) * photon.e1 +
                  std::sin(photon.alpha0) * photon.e2;
      assert(close(k.length(), 1));
      // Elastic scattering in the particle rest frame.
      assert(close(photon.omega_inf * (1 - beta * dot(k, b.B)),
                   energy * (1 - beta * b.mu)));
      assert(photon.pol == AdiabaticPolarization::E ||
             photon.pol == AdiabaticPolarization::O);
    }
  }
  test_multiple_components();
  std::cout << "Resonant scattering checks passed\n";
}
