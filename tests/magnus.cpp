// clang-format off
// From the repository root (create build/ if necessary):
//   g++-16 -std=c++23 -O2 -Wall -Wextra -pedantic -I src tests/magnus.cpp -o build/test_magnus
//   build/test_magnus
// clang-format on
#include "magnus.hpp" // First include: the header must be self-contained.
#include "geodesic.hpp"
#include <iostream>
#include <string>
#include <type_traits>

namespace {
void
require(bool value, char const* message) {
  if (!value)
    throw std::runtime_error(message);
}

double
distance(Stokes a, Stokes b) {
  auto d = a + (-1.0) * b;
  return std::hypot(d.x, d.y, d.z);
}

void
near(double actual, double expected, double tolerance, char const* message) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
    throw std::runtime_error(std::string(message) +
                             ": actual=" + std::to_string(actual) +
                             ", expected=" + std::to_string(expected));
}

void
same(Stokes a, Stokes b, char const* message) {
  require(a.x == b.x && a.y == b.y && a.z == b.z, message);
}

template <class Exception, class F>
void
throws(F const& f, char const* message) {
  try {
    f();
  } catch (Exception const&) {
    return;
  }
  throw std::runtime_error(message);
}

// Independent Rodrigues formula, without using the implementation's quaternion.
Stokes
rotate(Stokes s, double3 axis, double angle) {
  return std::cos(angle) * s + std::sin(angle) * cross(axis, s) +
         ((1 - std::cos(angle)) * dot(axis, s)) * axis;
}

// omega(t) = Rz(rate*t) (kappa, 0, z).
// In the corotating frame the generator is constant: (kappa, 0, z-rate).
struct TurningAxis {
  double kappa = 2.1, rate = 1.3, z = 0.4;
  double3
  operator()(double t) const {
    return {kappa * std::cos(rate * t), kappa * std::sin(rate * t), z};
  }
  Stokes
  exact(double t, Stokes initial) const {
    double3 w{kappa, 0, z - rate};
    double length = w.length();
    return rotate(rotate(initial, w / length, t * length), {0, 0, 1}, rate * t);
  }
};

void
test_constant_and_cache() {
  for (double frequency : {0.0, 1e-12, 2.0, 1000.0}) {
    for (Stokes initial : {Stokes{0, 1, 0}, Stokes{0.2, -0.3, 0.1}, Stokes{}}) {
      int calls = 0;
      bool callable = true;
      auto omega = [&](double x) -> double3 {
        require(callable, "Dense output must not evaluate the generator");
        require(x >= 3 && x <= 5, "Generator outside initialized interval");
        ++calls;
        return {frequency, 0, 0};
      };
      using Stepper = StepperMagnus;
      static_assert(
          !std::is_constructible_v<Stepper, decltype(omega)&&, double, double>);
      Stepper stepper(1e-10, 1e-10);
      stepper.init(3, 5, 2, initial);
      require(stepper.do_step(omega), "Constant step accepted");
      require(calls == 9,
              "Six trial Gauss samples and three accepted endpoints");
      callable = false;
      for (double t : {0., 0.07, 0.31, 0.5, 0.67, 0.91, 1.}) {
        double x = 3 + 2 * t;
        auto s = stepper.dense_out(x);
        near(distance(s, rotate(initial, {1, 0, 0}, frequency * (x - 3))), 0,
             3e-12, "Constant generator analytic dense solution");
        near(s.length(), initial.length(), 3e-14,
             "Preserve polarization degree");
      }
      same(stepper.dense_out(3), initial, "Exact initial endpoint");
      same(stepper.dense_out(5), stepper.y_new, "Exact final endpoint");
      auto midpoint = stepper.dense_out(4);
      near(distance(stepper.dense_out(std::nextafter(4.0, 3.0)), midpoint), 0,
           2e-12, "Continuous left midpoint limit");
      near(distance(stepper.dense_out(std::nextafter(4.0, 5.0)), midpoint), 0,
           2e-12, "Continuous right midpoint limit");
      require(calls == 9, "Queries reuse cached coefficients");
      throws<std::logic_error>([&] { stepper.do_step(omega); },
                               "Uncommitted step guard");
      auto final = stepper.y_new;
      stepper.update_old();
      same(stepper.y_old, final, "Commit accepted state");
      require(stepper.x_old == 5, "Commit exact interval endpoint");
      throws<std::logic_error>([&] { stepper.dense_out(4); },
                               "Invalidate committed cache");
      require(!stepper.do_step(omega), "Completed interval returns false");
    }
  }
  auto x = Rotation::exponential({0.3, 0, 0});
  auto y = Rotation::exponential({0, 0.7, 0});
  near(distance(Rotation::compose(y, x).apply({0, 0, 1}),
                y.apply(x.apply({0, 0, 1}))),
       0, 5e-16, "Composition order");
}

void
test_convergence() {
  TurningAxis omega;
  Stokes initial{0.2, -0.4, 0.7};
  double previous = 0;
  for (int n : {4, 8, 16, 32}) {
    // Loose tolerances and an explicit cap isolate the discretization order.
    StepperMagnus stepper(1., 0.);
    double h = 1.0 / n;
    stepper.init(0, 1, h, initial);
    int steps = 0;
    while (stepper.do_step(omega, h)) {
      require(stepper.h_old == h, "Fixed convergence step");
      stepper.update_old();
      ++steps;
    }
    require(steps == n, "Fixed step count");
    double error = distance(stepper.y_old, omega.exact(1, initial));
    std::cout << "global n=" << n << " error=" << error
              << " ratio=" << (previous ? previous / error : 0) << '\n';
    if (previous)
      require(previous / error > 14 && previous / error < 18,
              "Global fourth order");
    previous = error;
  }
  previous = 0;
  for (double h : {0.4, 0.2, 0.1, 0.05}) {
    double start = 0.2;
    StepperMagnus stepper(1., 0.);
    stepper.init(start, start + h, h, omega.exact(start, initial));
    require(stepper.do_step(omega), "Local convergence step");
    double error = 0;
    for (double t : {0.13, 0.37, 0.63, 0.87}) {
      double x = std::lerp(stepper.x_old, stepper.x_new, t);
      auto s = stepper.dense_out(x);
      error = std::max(error, distance(s, omega.exact(x, initial)));
      near(s.length(), initial.length(), 2e-14,
           "Dense norm for varying generator");
    }
    std::cout << "dense h=" << h << " error=" << error
              << " ratio=" << (previous ? previous / error : 0) << '\n';
    if (previous)
      require(previous / error > 25 && previous / error < 40,
              "Dense local fifth order");
    previous = error;
  }
}

void
test_adaptation_and_intervals() {
  TurningAxis model;
  Stokes initial{0.2, -0.4, 0.7};
  double previous = 0;
  int previous_steps = 0;
  for (double tolerance : {1e-5, 1e-8, 1e-11}) {
    int calls = 0;
    auto omega = [&](double x) {
      ++calls;
      require(x >= 0 && x <= 2, "Rejected trials stay in interval");
      return model(x);
    };
    StepperMagnus stepper(tolerance, tolerance);
    stepper.init(0, 2, 2, initial);
    require(stepper.do_step(omega), "Adaptive step accepted");
    require(calls > 9 && (calls - 3) % 6 == 0,
            "Rejected trials skip dense endpoints");
    require(stepper.h_old < 2, "Rejection reduces step");
    int steps = 1;
    auto check_error = [&] {
      double scale = tolerance * (1 + std::max(stepper.y_old.length(),
                                               stepper.y_new.length()));
      require(stepper.y_err.length() <= scale,
              "Accepted error within tolerance");
    };
    check_error();
    stepper.update_old();
    while (stepper.do_step(omega)) {
      check_error();
      stepper.update_old();
      require(++steps < 10000, "Adaptive termination");
    }
    double error = distance(stepper.y_old, model.exact(2, initial));
    std::cout << "tol=" << tolerance << " steps=" << steps << " error=" << error
              << '\n';
    if (previous) {
      require(error < previous / 30, "Tighter tolerance improves solution");
      require(steps > previous_steps, "Tighter tolerance refines steps");
    }
    previous = error;
    previous_steps = steps;
  }

  double left = 0, right = 0.7;
  auto omega = [&](double x) -> double3 {
    require(x >= left && x <= right, "Callback outside current outer interval");
    return {0, 0, 2};
  };
  StepperMagnus stepper(1e-10, 1e-10);
  Stokes state{1, 0, 0};
  double suggested = 10;
  for (double end : {0.7, 1.03, 1.87}) {
    right = end;
    stepper.init(left, right, suggested, state);
    int count = 0;
    while (stepper.do_step(omega, 0.2)) {
      require(stepper.h_old <= 0.2 && stepper.x_new <= right,
              "Step capped by both bounds");
      stepper.update_old();
      require(++count < 100, "Bounded interval terminates");
    }
    require(stepper.x_old == right, "Exact end after truncation");
    state = stepper.y_old;
    suggested = stepper.h_new;
    near(distance(state, {std::cos(2 * right), std::sin(2 * right), 0}), 0,
         1e-13, "Continuous state across outer intervals");
    left = right;
  }

  // No accidental normalization, and zero state works with relative-only error.
  StepperMagnus zero(0., 1e-8);
  left = 0;
  right = 1;
  zero.init(0, 1, 1, {});
  require(zero.do_step(omega), "Zero state with relative-only tolerance");
  same(zero.y_new, {}, "Zero state remains zero");
}

void
test_invalid_and_lifecycle() {
  auto omega = [](double) -> double3 { return {0, 0, 1}; };
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  for (auto tolerances :
       {std::array{0., 0.}, std::array{-1., 1.}, std::array{1., -1.},
        std::array{nan, 1.}, std::array{1., nan}, std::array{inf, 1.},
        std::array{1., inf}})
    throws<std::invalid_argument>(
        [&] { StepperMagnus s(tolerances[0], tolerances[1]); },
        "Reject invalid tolerances");
  StepperMagnus stepper(1e-8, 0.);
  throws<std::logic_error>([&] { stepper.do_step(omega); },
                           "Require initialization");
  throws<std::logic_error>([&] { stepper.dense_out(0); },
                           "No cache before initialization");
  throws<std::logic_error>([&] { stepper.update_old(); }, "No state to commit");
  for (double h : {0., -1., nan, inf})
    throws<std::invalid_argument>([&] { stepper.init(0, 1, h, {}); },
                                  "Reject invalid initial step");
  for (auto bounds :
       {std::array{1., 0.}, std::array{nan, 1.}, std::array{0., nan},
        std::array{-inf, 0.}, std::array{0., inf}})
    throws<std::invalid_argument>(
        [&] { stepper.init(bounds[0], bounds[1], 1, {}); },
        "Reject invalid interval");
  for (Stokes s : {Stokes{nan, 0, 0}, Stokes{0, inf, 0}, Stokes{0, 0, -inf}})
    throws<std::invalid_argument>([&] { stepper.init(0, 1, 1, s); },
                                  "Reject invalid Stokes");
  stepper.init(0, 1, 1, {1, 0, 0});
  throws<std::logic_error>([&] { stepper.dense_out(0); },
                           "No cache before acceptance");
  for (double h : {0., -1., nan})
    throws<std::invalid_argument>([&] { stepper.do_step(omega, h); },
                                  "Reject invalid maximum step");
  require(stepper.do_step(omega), "Valid state after invalid maximum step");
  for (double x : {-1., std::nextafter(1., 2.), nan, inf})
    throws<std::out_of_range>([&] { stepper.dense_out(x); },
                              "Reject dense extrapolation");
  stepper.init(2, 2, 1, {0.2, 0.3, 0.4});
  require(!stepper.do_step(omega), "Empty interval");
  same(stepper.y_old, {0.2, 0.3, 0.4},
       "Empty interval preserves initial state");
  throws<std::logic_error>([&] { stepper.dense_out(0.5); },
                           "Reinit discards previous cache");
  stepper.init(1, std::nextafter(1., 2.), 1, {});
  throws<std::runtime_error>([&] { stepper.do_step(omega); },
                             "Unrepresentable half steps");
  stepper.init(1, 2, 1e-30, {});
  throws<std::runtime_error>([&] { stepper.do_step(omega); }, "Step underflow");
  for (bool endpoints_only : {false, true}) {
    auto bad = [&](double x) -> double3 {
      return (!endpoints_only || x == 0) ? double3{nan, 0, 0}
                                         : double3{0, 0, 1};
    };
    StepperMagnus broken(1e-8, 0.);
    broken.init(0, 1, 1, {1, 0, 0});
    throws<std::runtime_error>([&] { broken.do_step(bad); },
                               "Reject nonfinite generator");
    throws<std::logic_error>([&] { broken.dense_out(0.5); },
                             "Failed trial exposes no cache");
    same(broken.y_old, {1, 0, 0}, "Failed trial preserves committed state");
  }
}

// The quadrature and root finder are independent test fixtures, not a physical
// scattering model. Opacity is 2 + Q/2, S=(cos(l),sin(l),0), hence
// tau(l)=2*l+sin(l)/2. This exercises the intended three-level nesting.
template <class F>
double
integrate(F const& f, double a, double b) {
  if (a == b)
    return 0;
  constexpr int n = 64;
  double sum = f(a) + f(b);
  for (int i = 1; i < n; ++i)
    sum += (i % 2 ? 4 : 2) * f(std::lerp(a, b, double(i) / n));
  return sum * (b - a) / (3 * n);
}

void
test_nested_optical_depth() {
  // Flat-space radial ray: geodesic r=10+l is also analytic.
  Geodesic geodesic(0, 10);
  geodesic.init(10, 0, 0);
  double left = 0, right = 0;
  int coefficient_calls = 0, orbit_steps = 0, pol_steps = 0;
  auto omega = [&](double l) -> double3 {
    require(l >= left && l <= right,
            "Magnus callback within live geodesic cache");
    near(geodesic.interp(l)[0], 10 + l, 1e-12,
         "Geometry at coefficient samples");
    ++coefficient_calls;
    return {0, 0, 1};
  };
  StepperMagnus polarization(1e-10, 1e-10);
  Stokes state{1, 0, 0};
  double h_pol = 0.003, accumulated = 0;
  constexpr double event = 0.37;
  double target = 2 * event + std::sin(event) / 2;
  while (++orbit_steps < 100) {
    bool absorbed = geodesic.step();
    require(!absorbed, "Outward radial test ray is not absorbed");
    left = geodesic.l_old();
    right = geodesic.l_new();
    polarization.init(left, right, h_pol, state);
    while (polarization.do_step(omega, 0.04)) {
      ++pol_steps;
      int cached_calls = coefficient_calls;
      auto opacity = [&](double l) {
        require(l >= left && l <= right, "Opacity within live geodesic cache");
        near(geodesic.interp(l)[0], 10 + l, 1e-12,
             "Geometry during quadrature");
        return 2 + polarization.dense_out(l).x / 2;
      };
      double a = polarization.x_old, b = polarization.x_new;
      // Splitting at the cached half-step join also supports general
      // generators.
      double mid = std::midpoint(a, b);
      auto depth_to = [&](double x) {
        return integrate(opacity, a, std::min(x, mid)) +
               (x > mid ? integrate(opacity, mid, x) : 0);
      };
      double dtau = depth_to(b);
      near(dtau, 2 * (b - a) + (std::sin(b) - std::sin(a)) / 2, 2e-13,
           "Analytic optical depth per polarization step");
      if (accumulated + dtau >= target) {
        double lo = a, hi = b;
        for (int i = 0; i < 48; ++i) {
          double x = std::midpoint(lo, hi);
          if (accumulated + depth_to(x) < target)
            lo = x;
          else
            hi = x;
        }
        double scattered = std::midpoint(lo, hi);
        near(scattered, event, 2e-12, "Locate scattering event");
        near(distance(polarization.dense_out(scattered),
                      {std::cos(event), std::sin(event), 0}),
             0, 3e-12, "Scattering Stokes before committing step");
        require(coefficient_calls == cached_calls,
                "Optical depth and root search use only cache");
        require(orbit_steps > 1 && pol_steps > orbit_steps,
                "Exercise all nested levels");
        return; // Leave the enclosing geodesic step uncommitted at scattering.
      }
      require(coefficient_calls == cached_calls,
              "Quadrature does not resample generator");
      accumulated += dtau;
      polarization.update_old();
    }
    state = polarization.y_old;
    h_pol = polarization.h_new;
    near(accumulated, 2 * right + std::sin(right) / 2, 2e-12,
         "Optical depth across geodesic intervals");
    geodesic.stepper.update_old();
  }
  throw std::runtime_error("Scattering event not found");
}
} // namespace

int
main() try {
  test_constant_and_cache();
  test_convergence();
  test_adaptation_and_intervals();
  test_invalid_and_lifecycle();
  test_nested_optical_depth();
  std::cout << "All Magnus checks passed.\n";
} catch (std::exception const& e) {
  std::cerr << "FAIL: " << e.what() << '\n';
  return 1;
}
