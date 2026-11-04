#include "qed_scattering.hpp"
#include <chrono>
#include <iostream>
#include <random>

using namespace qed;
static int checks = 0;
void
check(bool ok, char const* message) {
  ++checks;
  if (!ok)
    throw std::runtime_error(message);
}
void
near(double a, double b, double rel = 1e-10, double abs = 1e-14) {
  if (std::abs(a - b) > abs + rel * std::max(std::abs(a), std::abs(b)))
    throw std::runtime_error(std::format("Mismatch: {} vs {}", a, b));
  ++checks;
}
Jones
positive_charge_reference(double b, double x, double mui, double muf,
                          double phi, RestWidths widths) {
  // Direct positive-charge Hamiltonian states and swapped oscillator orders;
  // deliberately does not call the production reversed-field mapping.
  auto positive = [](int n, int s, int eta, double b, double q) {
    auto c = spinor(n, s, -eta, b, -q);
    return Spinor{std::conj(c[3]), -std::conj(c[2]), -std::conj(c[1]),
                  std::conj(c[0])};
  };
  auto vtx = [&](Spinor l, int nl, double al, Spinor r, int nr, double ar,
                 double3 k, double3 e) {
    auto A = coupling(e);
    Complex out = 0;
    for (int j = 0; j < 4; ++j)
      for (int h = 0; h < 4; ++h) {
        if (A[j][h] == Complex{} || l[j] == Complex{} || r[h] == Complex{})
          continue;
        out += std::conj(l[j]) * A[j][h] * r[h] *
               overlap(nl - j % 2, nr - h % 2, b, k.x, al, ar);
      }
    return out;
  };
  auto kin = recoil(x, mui, muf);
  double y = kin.y;
  double3 ki{x * std::sqrt(1 - mui * mui), 0, x * mui};
  double3 kf{y * std::sqrt(1 - muf * muf) * std::cos(phi),
             y * std::sqrt(1 - muf * muf) * std::sin(phi), y * muf};
  auto ei = real_basis(mui, 0), ef = real_basis(muf, phi);
  auto ci = positive(0, -1, 1, b, 0), cf = positive(0, -1, 1, b, kin.qf);
  double af = (ki.y - kf.y) / b;
  Jones out{};
  for (int diagram = 0; diagram < 2; ++diagram) {
    double V = diagram ? 1 - y : 1 + x, q = diagram ? -kf.z : ki.z,
           ac = diagram ? -kf.y / b : ki.y / b;
    for (int n : {0, 1})
      for (int s : {-1, 1})
        for (int eta : {-1, 1}) {
          if (n == 0 && s == 1)
            continue;
          auto c = positive(n, s, eta, b, q);
          double E = std::sqrt(1 + 2 * n * b + q * q);
          Complex D{V - eta * E,
                    n ? eta * widths.at(s) * std::sqrt(1 + 2 * b) / (2 * E)
                      : 0};
          for (int row = 0; row < 2; ++row)
            for (int col = 0; col < 2; ++col) {
              if (diagram == 0)
                out[row][col] += vtx(cf, 0, af, c, n, ac, -1. * kf, ef[row]) *
                                 vtx(c, n, ac, ci, 0, 0, ki, ei[col]) / D;
              else
                out[row][col] += vtx(cf, 0, af, c, n, ac, ki, ei[col]) *
                                 vtx(c, n, ac, ci, 0, 0, -1. * kf, ef[row]) / D;
            }
        }
  }
  return out;
}
void
kernel_tests() {
  for (double b : {1e-4, 0.01, 0.1, 1., 10.})
    for (double q : {-2., 0., 0.37})
      for (int n : {0, 1}) {
        double p = std::sqrt(2 * n * b), E = std::sqrt(1 + 2 * n * b + q * q);
        Matrix4 H{{{1, 0, q, -I * p},
                   {0, 1, I * p, -q},
                   {q, -I * p, -1, 0},
                   {I * p, -q, 0, -1}}};
        std::vector<Spinor> columns;
        for (int s : {-1, 1})
          for (int eta : {-1, 1}) {
            if (n == 0 && s == 1)
              continue;
            auto c = spinor(n, s, eta, b, q);
            columns.push_back(c);
            for (int j = 0; j < 4; ++j) {
              Complex v = 0;
              for (int k = 0; k < 4; ++k)
                v += H[j][k] * c[k];
              near(std::abs(v - eta * E * c[j]), 0, 0, 2e-14);
            }
          }
        for (size_t i = 0; i < columns.size(); ++i)
          for (size_t j = 0; j < columns.size(); ++j) {
            Complex v = 0;
            for (int k = 0; k < 4; ++k)
              v += std::conj(columns[i][k]) * columns[j][k];
            near(std::abs(v - Complex(i == j ? 1 : 0)), 0, 0, 2e-14);
          }
        double V = .17;
        Matrix4 resolvent{};
        for (int s : {-1, 1})
          for (int eta : {-1, 1}) {
            if (n == 0 && s == 1)
              continue;
            auto c = spinor(n, s, eta, b, q);
            for (int i = 0; i < 4; ++i)
              for (int j = 0; j < 4; ++j)
                resolvent[i][j] += c[i] * std::conj(c[j]) / (V - eta * E);
          }
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            Complex ref = (V * double(i == j) + H[i][j]) / (V * V - E * E);
            if (n == 0 && (i % 2 == 0 || j % 2 == 0))
              ref = 0;
            near(std::abs(ref - resolvent[i][j]), 0, 0, 2e-14);
          }
      }
  std::array<double, 4> bs{0.01, 0.1, 1, 10};
  std::array<double, 4> gm{9.4636115403e-7, 7.5278889649e-5, 2.0013591911e-3,
                           1.0684385576e-2};
  std::array<double, 4> gp{4.6850710298e-9, 3.4288668450e-6, 5.3626257901e-4,
                           6.8566235677e-3};
  for (int k = 0; k < 4; ++k) {
    auto w = rest_widths(bs[k]);
    near(w.minus, gm[k], 1e-10, 0);
    near(w.plus, gp[k], 1e-10, 0);
    double a = std::sqrt(1 + 2 * bs[k]);
    near(w.plus / w.minus, 2 * bs[k] / ((a + 1) * (a + 1)), 1e-10);
  }
  auto k = recoil(.5, 1, 0);
  near(k.y, (3 - std::sqrt(5.)) / 2);
  near(k.y + k.Ef, 1.5);
  for (double b : {.03, .4, 1., 3.})
    for (double x : {.01, .1, .7, 2.5})
      for (double mu : {-.8, -.2, .4, .9}) {
        auto kin = recoil(x, 1, mu);
        double y = kin.y, z = x * y * (1 - mu);
        auto t = amplitude(b, x, 1, mu, 0, {});
        double F = pi * y * y * std::exp(-y * y * (1 - mu * mu) / (2 * b)) /
                   (x * (2 * x - y - z));
        double terms = 1 / ((x - b) * (x - b)) + 1 / std::pow(x + b - z, 2);
        for (int row = 0; row < 2; ++row) {
          double T = row == 0 ? (2 + x) * (x - z) - 2 * y : x * (x - z);
          double actual = pi * y / (x * kin.J) *
                          (std::norm(t[row][0]) + std::norm(t[row][1]));
          near(actual, F / 2 * T * terms, 2e-9, 1e-12);
        }
      }
  for (double b : {.01, .4, 3.})
    for (double mui : {-.8, .3, 1.}) {
      double x = .2, muf = -.2, phi = .7;
      auto kin = recoil(x, mui, muf);
      double y = kin.y, E = kin.Ef, q = kin.qf;
      double G = std::exp(-(x * x * (1 - mui * mui) + y * y * (1 - muf * muf)) /
                          (4 * b));
      double zeta = x * std::sqrt(1 - mui * mui) * y *
                    std::sqrt(1 - muf * muf) * std::sin(phi) / (2 * b);
      double Fs = (E + 1 + q * mui) /
                  (std::sqrt(2 * E * (E + 1)) * (2 + x * (1 - mui * mui)));
      double Fu = (E + 1 + q * muf) /
                  (std::sqrt(2 * E * (E + 1)) * (2 - y * (1 - muf * muf)));
      double p0 = (x * std::sqrt(1 - mui * mui) -
                   y * std::sqrt(1 - muf * muf) * std::cos(phi)) *
                  y * std::sqrt(1 - muf * muf) * std::sin(phi) / (2 * b);
      Complex ref = G * std::exp(I * p0) *
                    std::sqrt((1 - mui * mui) * (1 - muf * muf)) *
                    (Fs * std::exp(-I * zeta) + Fu * std::exp(I * zeta));
      auto t = amplitude(b, x, mui, muf, phi, {}, 0);
      near(std::abs(t[0][0] - ref), 0, 0, 2e-13);
      near(std::abs(t[0][1]) + std::abs(t[1][0]) + std::abs(t[1][1]), 0);
      auto w = rest_widths(b);
      auto a = amplitude(b, x, mui, muf, phi, w);
      auto shifted = amplitude(b, x, mui, muf, phi, w, 1, .37);
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
          near(std::norm(a[i][j]), std::norm(shifted[i][j]), 1e-11);
    }
  Jones phase{{{1, 0}, {0, I}}};
  auto s = outgoing_stokes(phase, {0, 1, 0});
  near(s.x, 0);
  near(s.y, 0);
  near(s.z, -1);
  for (Stokes input : {Stokes{1, 0, 0},
                       {-1, 0, 0},
                       {0, 1, 0},
                       {0, -1, 0},
                       {0, 0, 1},
                       {0, 0, -1},
                       {.2, .1, .3}}) {
    auto d = differential(-1, .4, .2, .3, -.2, .7, rest_widths(.4));
    auto output = outgoing_stokes(d.reduced, input);
    check(output.length() <= 1 + 1e-12, "output positivity");
    if (input.length() == 1)
      near(output.length(), 1);
    check(contract(response(d.reduced).value, input) >= 0,
          "negative polarized weight");
  }
  bool failed = false;
  try {
    check_domain(4, 1);
  } catch (Error const& e) {
    failed = e.kind == Failure::PairDomain;
  }
  check(failed, "domain not enforced");
  for (double b : {.01, .1, 1., 3.})
    for (double mui : {-.7, .2, .8})
      for (double phi : {.3, 2.}) {
        auto w = rest_widths(b);
        auto pos = positive_charge_reference(b, .3, mui, -.4, phi, w);
        auto map = amplitude(b, .3, -mui, .4, -phi, w);
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 2; ++j)
            near(std::abs(pos[i][j] - map[i][j]), 0, 0, 3e-13);
      }
  for (double b : {1e-5, 1e-7, 1e-10}) {
    auto w = rest_widths(b);
    near(w.minus / (fine_structure * b * b), 4. / 3, 10 * b, 1e-12);
    near(w.plus / (fine_structure * b * b * b), 2. / 3, 10 * b, 1e-12);
  }
  for (double b : {.001, .03, .4, 3.})
    for (double mui : {-.7, .2, 1.})
      for (double muf : {-1., -.4, .9})
        for (int charge : {-1, 1}) {
          auto w = rest_widths(b);
          AzimuthKernel line(charge, b, .2, mui, muf, w);
          for (double phi : {.123, 1.7, 3.8, 6.}) {
            auto direct = differential(charge, b, .2, mui, muf, phi, w).reduced;
            auto fast = line.reduced(phi);
            auto rd = response(direct), rf = response(fast);
            for (int j = 0; j < 4; ++j)
              near(rd.value[j], rf.value[j], 1e-10, 1e-12);
            auto sd = outgoing_stokes(direct, {.2, .1, .3}),
                 sf = outgoing_stokes(fast, {.2, .1, .3});
            near(sd.x, sf.x, 1e-10, 1e-12);
            near(sd.y, sf.y, 1e-10, 1e-12);
            near(sd.z, sf.z, 1e-10, 1e-12);
          }
        }
}
void
response_tests() {
  Options opts;
  AngularResponse a{-1, .4, .2, .3, rest_widths(.4), opts};
  auto start = std::chrono::steady_clock::now();
  auto total = a.reduced_integral();
  std::cout << "Angular response " << total.value[0] << " +/- "
            << response_error(total) << " in "
            << std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                             start)
                   .count()
            << " s\n";
  check(envelope(total) > 0, "empty angular response");
  for (Stokes input : {Stokes{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {.2, .1, .3}}) {
    // Independent composite Simpson + periodic sum, no response contraction
    // inside the matrix update used to evaluate the polarized intensity.
    double direct = 0;
    int nm = 80, np = 64;
    auto rho = density(input);
    for (int i = 0; i <= nm; ++i) {
      double mu = -1 + 2. * i / nm,
             weight = i == 0 || i == nm ? 1 : (i % 2 ? 4 : 2);
      for (int j = 0; j < np; ++j) {
        auto m =
            differential(-1, .4, .2, .3, mu, 2 * pi * j / np, a.widths).reduced;
        Complex t = 0;
        for (int row = 0; row < 2; ++row)
          for (int k = 0; k < 2; ++k)
            for (int l = 0; l < 2; ++l)
              t += m[row][k] * rho[k][l] * std::conj(m[row][l]);
        direct += weight * t.real() * 2 / (3 * nm) * 2 * pi / np;
      }
    }
    near(contract(total.value, input), direct, 2e-7);
    auto sample = sample_outgoing(a, input, .37, .63);
    near(sample.kernel.kinematics.y + sample.kernel.kinematics.Ef, 1.2);
    check(sample.stokes_B.length() <= 1 + 1e-12, "sampled state invalid");
  }
  auto e = CurrentComponent{.dist = NegateDist(PowerLaw2(.2, .8)),
                            .charge_sign = -1,
                            .current_fraction = .5,
                            .charge_mass_ratio = 1.};
  auto p = CurrentComponent{.dist = PowerLaw2(.2, .8),
                            .charge_sign = 1,
                            .current_fraction = .5,
                            .charge_mass_ratio = 1.};
  QEDScattering scattering(std::make_tuple(e, p));
  Basis g{};
  g.r = 10;
  g.Bmag = .4 * critical_field_gauss;
  g.energy_local = .2 * electron_rest_keV;
  g.J = 1e10;
  g.mu = .3;
  start = std::chrono::steady_clock::now();
  auto local = scattering.local_response(g);
  std::cout << "Local response " << local.Lambda << " in "
            << std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                             start)
                   .count()
            << " s\n";
  check(local.Lambda > 0, "empty local response");
  auto invalid_geometry = g;
  invalid_geometry.Bmag = std::numeric_limits<double>::quiet_NaN();
  bool rejected_cache = false;
  try {
    scattering.local_response(invalid_geometry);
  } catch (Error const& e) {
    rejected_cache = e.kind == Failure::InvalidInput;
  }
  check(rejected_cache, "nonfinite cache key bypassed input validation");
  auto refined = scattering;
  refined.options.angular_rtol *= .1;
  refined.options.speed_rtol *= .1;
  refined.options.width_rtol *= .1;
  refined.reset();
  auto local_refined = refined.local_response(g);
  near(local.Lambda, local_refined.Lambda, 2e-6, 0);
  std::cout << "Speed refinement relative rate change="
            << std::abs(local.Lambda / local_refined.Lambda - 1) << '\n';
  Ran species_random(6127);
  int selected_first = 0;
  for (int i = 0; i < 100000; ++i)
    selected_first += decltype(scattering)::select_population(
                          {3., 1.}, species_random.U()) == 0;
  near(selected_first / 100000., .75, 0, .008);
  // A resolved real event uses the actual collision-weighted speed PDF.
  g.B = {0, 0, g.Bmag};
  g.r_hat = {1, 0, 0};
  g.k_hat = {std::sqrt(1 - g.mu * g.mu), 0, g.mu};
  g.screen_1 = {-g.mu, 0, std::sqrt(1 - g.mu * g.mu)};
  g.screen_2 = {0, -1, 0};
  g.c = 1;
  g.s = 0;
  Ran random(90210);
  Polarization incoming(Stokes{.2, .1, .3});
  std::optional<decltype(scattering)::ScatterData> event;
  for (int i = 0; i < 100 && !event; ++i)
    event = scattering.accept_candidate(g, incoming, random);
  check(event.has_value(), "candidate never accepted");
  std::cout << "Accepted QED target beta=" << event->beta << std::endl;
  check(event->charge_sign * event->beta > 0,
        "selected velocity reflected twice");
  double einf = .9 * g.energy_local;
  auto out = scattering.perform_scattering(*event, einf, random);
  check(!out.pol.is_adiabatic, "coherent output projected to a mode");
  auto kout = std::cos(out.alpha0) * out.e1 + std::sin(out.alpha0) * out.e2;
  double beta = event->beta, eps = g.energy_local / electron_rest_keV;
  double epsout = out.energy_inf / .9 / electron_rest_keV;
  auto star = boost(eps, g.mu, beta), outstar = boost(epsout, kout.z, beta);
  double qf = star.x * star.mu - outstar.x * outstar.mu,
         Ef = std::hypot(1., qf);
  near(eps + star.gamma, epsout + star.gamma * (Ef + beta * qf), 1e-11);
  near(eps * g.mu + star.gamma * beta,
       epsout * kout.z + star.gamma * (qf + beta * Ef), 1e-11);
  near(out.r0, g.r);
  near(dot(kout, out.n), 0, 0, 1e-12);
  Quad quad;
  double tau = -.1;
  auto point = [&](double) { return g; };
  auto candidate = scattering.next_candidate(0, 1000, point, tau, quad);
  check(candidate.has_value(), "constant-medium candidate not found");
  near(*candidate, .1 / local.Lambda, 1e-5);
  double split_tau = -.1, mid = .05 / local.Lambda;
  check(!scattering.next_candidate(0, mid, point, split_tau, quad),
        "event placed after interval boundary");
  auto second = scattering.next_candidate(mid, 1000, point, split_tau, quad);
  check(second.has_value(), "residual candidate clock lost");
  near(*second, *candidate, 1e-5);
  g.J = 0;
  near(scattering.local_response(g).Lambda, 0);
  g.J = 1;
  g.energy_local = 4 * electron_rest_keV;
  bool failed = false;
  try {
    scattering.local_response(g);
  } catch (Error const& e) {
    failed = e.kind == Failure::PairDomain;
  }
  check(failed, "population support domain not checked");
  auto thermal = CurrentComponent{.dist = Boltzmann(.2),
                                  .charge_sign = 1,
                                  .current_fraction = 1.,
                                  .charge_mass_ratio = 1.};
  QEDScattering thermal_scattering(std::make_tuple(thermal));
  g.energy_local = .2 * electron_rest_keV;
  g.mu = .3;
  failed = false;
  try {
    thermal_scattering.local_response(g);
  } catch (Error const& e) {
    failed = e.kind == Failure::PairDomain;
  }
  check(failed, "unsupported thermal tail silently truncated");
  g.mu = 1;
  check(thermal_scattering.local_response(g).Lambda > 0,
        "co-propagating axial thermal tail incorrectly rejected");
  auto tail = g;
  tail.Bmag = 1e-8 * critical_field_gauss;
  tail.energy_local = .003 * electron_rest_keV;
  tail.J = 1;
  tail.mu = .9;
  auto bounded = scattering.local_response(tail);
  check(bounded.negligible_bound > 0 &&
            bounded.negligible_bound <= scattering.options.rate_atol,
        "absolute opacity bound was not retained");
  auto exact_tail = scattering;
  exact_tail.options.rate_atol = 0;
  exact_tail.reset();
  auto reference = exact_tail.local_response(tail);
  check(reference.Lambda > 0 && reference.Lambda <= bounded.negligible_bound,
        "complete angular/speed opacity bound violated");
}
void
numerical_tests() {
  Options opts;
  for (double offset : {0., 1e-5, .003}) {
    double epsilon = .02, mu = .4;
    double b = epsilon * std::sqrt(1 - mu * mu) +
               epsilon * epsilon * (1 - mu * mu) / 2 + offset;
    auto w = rest_widths(b);
    double lo = .1, hi = .8;
    auto density = [&](double beta) {
      auto v = boost(epsilon, mu, beta);
      double E = std::sqrt(1 + 2 * b + v.x * v.x * v.mu * v.mu);
      double delta =
          (2 * v.x + epsilon * epsilon * (1 - mu * mu) - 2 * b) / (1 + v.x + E);
      double minus = w.minus * std::sqrt(1 + 2 * b) / (2 * E),
             plus = w.plus * std::sqrt(1 + 2 * b) / (2 * E);
      return Response{{minus / (delta * delta + minus * minus) +
                           plus / (delta * delta + plus * plus),
                       0, 0, 0},
                      {}};
    };
    auto direct =
        integrate<4>(density, speed_cuts(b, epsilon, mu, w, lo, hi), 1e-7);
    auto mapped = integrate_speed(density, b, epsilon, mu, w, lo, hi, opts);
    near(direct.value[0], mapped.value[0], 3e-6, 0);
    auto cdf = polarized(mapped, {0, 0, 0});
    auto f = [&](double beta) { return polarized(density(beta), {0, 0, 0}); };
    for (double u : {.01, .2, .5, .9, .99}) {
      double beta = inverse_cdf(f, cdf, u, 4e-6);
      auto part = integrate_speed(density, b, epsilon, mu, w, lo, beta, opts);
      near(part.value[0] / mapped.value[0], u, 0, 1e-5);
    }
  }
  // Endpoint root and two roots are physically guided rather than detected
  // from a sparse scan of a narrow positive function.
  for (double beta : {.2, .8}) {
    double e = .02, mu = .4, x = boost(e, mu, beta).x;
    double b = x + e * e * (1 - mu * mu) / 2;
    auto roots = speed_roots(b, e, mu);
    check(std::any_of(roots.begin(), roots.end(),
                      [&](double r) { return std::abs(r - beta) < 1e-12; }),
          "missing endpoint resonance");
  }
  auto flat = [](double) { return Estimate<1>{{1}, {}}; };
  auto uniform = integrate<1>(flat, {0, 1}, 1e-8);
  near(inverse_cdf(flat, uniform, .37, 1e-8), .37);
  auto gap = [](double x) {
    return Estimate<1>{{x < .4 || x > .6 ? 1. : 0.}, {}};
  };
  auto disconnected = integrate<1>(gap, {0, .4, .6, 1}, 1e-8);
  near(disconnected.value[0], .8);
  check(inverse_cdf(gap, disconnected, .75, 1e-8) > .6, "CDF plateau sampled");
  bool failed = false;
  try {
    integrate<1>([](double x) { return Estimate<1>{{std::exp(20 * x)}, {}}; },
                 {0, 1}, 1e-12, 0, 1);
  } catch (Error const& e) {
    failed = e.kind == Failure::Convergence;
  }
  check(failed, "panel exhaustion not reported");
  failed = false;
  try {
    envelope(Response{{1, -1, 0, 0}, {}});
  } catch (Error const& e) {
    failed = e.kind == Failure::Convergence;
  }
  check(failed, "negative response not reported");
  // Constant response thinning and true depth inversion, independent random
  // trajectories. This tests the no-selective-extinction survival law.
  Ran ran(20261103);
  double sum = 0, survive = 0;
  int N = 100000;
  for (int i = 0; i < N; ++i) {
    double path = 0;
    do {
      path -= std::log(ran.U()) / 1.1;
    } while (ran.U() >= .7 / 1.1);
    sum += path;
    survive += path > 2;
  }
  near(sum / N, 1 / .7, 0, .02);
  near(survive / N, std::exp(-1.4), 0, .006);
  std::cout << "Thinning mean=" << sum / N << ", P(L>2)=" << survive / N
            << '\n';
  // Variable rate with vacuum rotation S=(cos(2l),sin(2l),0). The envelope
  // integral l+0.1 l^2 and the physical depth have independent closed forms.
  auto physical_depth = [](double l) {
    return .75 * l + .075 * l * l + .125 * std::sin(2 * l) +
           .05 * (l * std::sin(2 * l) / 2 + (std::cos(2 * l) - 1) / 4);
  };
  double depth_sum = 0, tail = 0;
  for (int i = 0; i < 20000; ++i) {
    double l = 0;
    for (;;) {
      double tau = std::log(ran.U());
      auto envelope = [](double position) {
        return Estimate<1>{{1 + .2 * position}, {}};
      };
      auto candidate =
          candidate_from_rate(envelope, {l, 100}, tau, 1e-9, 0, 8192);
      check(candidate.has_value(), "variable-rate candidate clock failed");
      l = *candidate;
      double rate = .75 + .15 * l + (.25 + .05 * l) * std::cos(2 * l);
      if (ran.U() * (1 + .2 * l) < rate)
        break;
    }
    double depth = physical_depth(l);
    depth_sum += depth;
    tail += depth > 2;
  }
  near(depth_sum / 20000, 1, 0, .025);
  near(tail / 20000, std::exp(-2.), 0, .012);
  std::cout << "Variable-ray physical-depth mean=" << depth_sum / 20000 << '\n';
}

void
convergence_tests() {
  Options opts;
  AngularResponse a{-1, .4, .2, .3, rest_widths(.4), opts};
  auto coarse = a.reduced_integral();
  a.options.angular_rtol = 1e-9;
  auto fine = a.reduced_integral();
  for (int j = 0; j < 4; ++j)
    near(coarse.value[j], fine.value[j], 1e-7, 1e-10);
  std::cout << "Angular refinement relative trace change="
            << std::abs(coarse.value[0] / fine.value[0] - 1) << '\n';
  AngularResponse low{-1, .1, 1e-5, .3, {}, opts};
  auto total = low.reduced_integral();
  near(contract(total.value, {1, 0, 0}), 8 * pi / 3 * (1 - .3 * .3), 3e-5);
  near(contract(total.value, {-1, 0, 0}), 0, 0, 1e-6);
  double n = 1e15;
  near(1e5 * n * thomson_cross_section_cm2, 6.65245873e-5, 1e-8);
  // Complete angular bound checks and a provably below-double-range response.
  check(std::log(envelope(fine)) < a.log_bound(), "angular bound violated");
  AngularResponse tails{-1,  2.66958e-8, .0035, .985, rest_widths(2.66958e-8),
                        opts};
  auto bounded = tails.reduced_integral();
  tails.options.bound_azimuth_tails = false;
  auto explicit_tails = tails.reduced_integral();
  for (int j = 0; j < 4; ++j)
    near(bounded.value[j], explicit_tails.value[j], 1e-6,
         bounded.error[j] + explicit_tails.error[j]);
  AngularResponse tiny{-1, 1e-14, .001, .3, rest_widths(1e-14), opts};
  check(tiny.underflows(), "missing certified underflow bound");
  near(tiny.reduced_integral().value[0], 0);
  // Independent angular sampling check in the controlled longitudinal limit.
  // Analytic CDF is (2+3mu-mu^3)/4; azimuth is uniform.
  Ran rng(8147);
  std::vector<double> us, vs;
  for (int i = 0; i < 256; ++i) {
    auto s = sample_outgoing(low, {1, 0, 0}, rng.U(), rng.U());
    double mu = s.muf;
    us.push_back(.5 + .75 * mu - .25 * mu * mu * mu);
    vs.push_back(s.phi / (2 * pi));
    near(s.stokes_B.x, 1, 0, 1e-6);
  }
  for (auto samples : {us, vs}) {
    std::sort(samples.begin(), samples.end());
    double ks = 0;
    for (size_t i = 0; i < samples.size(); ++i)
      ks = std::max({ks, std::abs(samples[i] - double(i) / samples.size()),
                     std::abs(samples[i] - double(i + 1) / samples.size())});
    check(ks < .11,
          "sampled angular CDF disagrees with analytic low-energy limit");
    std::cout << "Low-energy angular KS=" << ks << '\n';
  }
}

void
transport_tests() {
  BField field("table/bfield_t10.txt", 1e8, 10);
  InitFD11 init(10, .4);
  auto stream = CurrentComponent{.dist = PowerLaw2(.2, .8),
                                 .charge_sign = 1,
                                 .current_fraction = 0.,
                                 .charge_mass_ratio = 1.};
  QEDScattering scat(std::make_tuple(stream));
  std::array<PhotonEvolution<InitFD11, decltype(scat)>::EscapeData, 2> result;
  for (int i = 0; i < 2; ++i) {
    PhotonEvolution<InitFD11, decltype(scat)>::Options options;
    options.geodesic.escape_factor = i ? 2e5 : 1e5;
    PhotonEvolution photon(1.4, 10, 42, field, init, scat, options);
    photon.init();
    check(photon.evolve() == EvolveResult::Escaped,
          "zero density escape failed");
    result[i] = photon.escape_data;
    near(std::hypot(result[i].Q, result[i].U, result[i].V), 1, 1e-8);
  }
  near(result[0].energy_inf, result[1].energy_inf);
  near(result[0].muk, result[1].muk, 0, 2e-5);
  near(result[0].Q, result[1].Q, 0, 2e-5);
  near(result[0].U, result[1].U, 0, 2e-5);
  near(result[0].V, result[1].V, 0, 2e-5);
  PhotonEvolution photon(1.4, 10, 43, field, init, scat);
  InitializedPhoton inward{.n = {0, 0, 1},
                           .e1 = {1, 0, 0},
                           .e2 = {0, 1, 0},
                           .energy_inf = .4,
                           .r0 = 11,
                           .psi0 = 0,
                           .alpha0 = pi,
                           .pol = Polarization(Stokes{0, 0, 1})};
  photon.init(inward);
  check(photon.evolve() == EvolveResult::Absorbed, "inward absorption failed");
  photon.init(inward);
  photon.pol.set_adiabatic(AdiabaticPolarization::E);
  photon.record_escape({11, 0, 0});
  near(std::hypot(photon.escape_data.Q, photon.escape_data.U,
                  photon.escape_data.V),
       1);
  std::cout
      << "Zero-density boundary convergence and inward absorption passed\n";
}

void
smoke_tests() {
  // Deliberately weak but nonzero magnetosphere: resolves the full kernel
  // while allowing a short reproducible end-to-end test of both boundaries.
  BField field("table/bfield_t10.txt", 1e8, 10);
  InitFD11 init(10, .4);
  auto e = CurrentComponent{.dist = NegateDist(PowerLaw2(.2, .8)),
                            .charge_sign = -1,
                            .current_fraction = .5,
                            .charge_mass_ratio = 1.};
  auto p = CurrentComponent{.dist = PowerLaw2(.2, .8),
                            .charge_sign = 1,
                            .current_fraction = .5,
                            .charge_mass_ratio = 1.};
  QEDScattering scat(std::make_tuple(e, p));
  std::array<PhotonEvolution<InitFD11, decltype(scat)>::EscapeData, 2> escaped;
  for (int outer = 0; outer < 2; ++outer) {
    PhotonEvolution<InitFD11, decltype(scat)>::Options options;
    options.geodesic.escape_factor = outer ? 2e5 : 1e5;
    if (std::getenv("MAGMC_QED_TRACE"))
      options.progress = [](double l, double r) {
        std::cout << "Smoke flight l=" << l << " r=" << r << '\n';
      };
    auto start = std::chrono::steady_clock::now();
    PhotonEvolution photon(1.4, 10, 42, field, init, scat, options);
    photon.init();
    int count = 0;
    for (;;) {
      auto result = photon.evolve();
      if (result == EvolveResult::Scattered) {
        ++count;
        photon.perform_scattering();
        continue;
      }
      check(result == EvolveResult::Escaped,
            "weak-field smoke photon absorbed");
      escaped[outer] = photon.escape_data;
      break;
    }
    std::cout << "QED smoke outer=" << options.geodesic.escape_factor
              << " scatters=" << count
              << " energy=" << escaped[outer].energy_inf << " seconds="
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - start)
                     .count()
              << std::endl;
  }
  near(escaped[0].energy_inf, escaped[1].energy_inf);
  near(escaped[0].muk, escaped[1].muk, 0, 2e-5);
  near(escaped[0].Q, escaped[1].Q, 0, 2e-5);
  near(escaped[0].U, escaped[1].U, 0, 2e-5);
  near(escaped[0].V, escaped[1].V, 0, 2e-5);
}
int
main(int argc, char** argv) try {
  std::cout << std::unitbuf;
  auto start = std::chrono::steady_clock::now();
  kernel_tests();
  numerical_tests();
  if (argc > 1) {
    response_tests();
    convergence_tests();
    transport_tests();
  }
  if (argc > 1 && std::string(argv[1]) == "smoke")
    smoke_tests();
  std::cout << "Passed " << checks << " checks in "
            << std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                             start)
                   .count()
            << " s\n";
} catch (std::exception const& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
