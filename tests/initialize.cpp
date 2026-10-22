// From the repository root:
//   mkdir -p build
//   g++-16 -std=c++23 -O2 -I src tests/initialize.cpp -o build/test_initialize
//   build/test_initialize
#include "initialize.hpp"
#include "geodesic.hpp"
#include <array>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>

namespace {

constexpr double R_star = 10;
constexpr double M_star = 1.4;
constexpr double rs = M_star * schwarzschild_radius_of_sun_in_km;
constexpr double T_inf = 0.4;
constexpr double zeta3 = 1.2020569031595942854;
constexpr double mean_x = pi * pi * pi * pi / (30 * zeta3);
constexpr double mean_x2 = 12 * 1.0369277551433699263 / zeta3;
constexpr double mean_x4 = 360 * 1.0083492773819228268 / zeta3;

void
require(bool condition, std::string_view message) {
  if (!condition)
    throw std::runtime_error(std::string(message));
}

void
near(double value, double expected, double tolerance,
     std::string_view message) {
  if (!std::isfinite(value) || std::abs(value - expected) > tolerance) {
    std::ostringstream out;
    out << std::setprecision(17) << message << ": got " << value
        << ", expected " << expected << " +/- " << tolerance;
    throw std::runtime_error(out.str());
  }
}

void
mean_near(double value, double expected, double variance, int count,
          std::string_view message) {
  near(value, expected, 6 * std::sqrt(variance / count) + 1e-12, message);
}

void
probability_near(int hits, int count, double probability,
                 std::string_view message) {
  mean_near(double(hits) / count, probability, probability * (1 - probability),
            count, message);
}

void
same_vector(double3 a, double3 b) {
  require(a.x == b.x && a.y == b.y && a.z == b.z, "Reproducible vector");
}

void
same_geometry_and_mode(InitializedPhoton const& a, InitializedPhoton const& b) {
  same_vector(a.n, b.n);
  same_vector(a.e1, b.e1);
  same_vector(a.e2, b.e2);
  require(a.r0 == b.r0 && a.psi0 == b.psi0 && a.alpha0 == b.alpha0 &&
              a.pol == b.pol,
          "Reproducible geometry and mode");
}

void
check_geometry(InitializedPhoton const& p) {
  require(std::isfinite(p.omega_inf) && p.omega_inf > 0,
          "Positive finite energy");
  require(p.r0 == R_star && p.psi0 == 0, "Starts at stellar surface, psi = 0");
  require(std::isfinite(p.alpha0) && p.alpha0 >= 0 && p.alpha0 <= pi / 2,
          "Outward emission angle");
  require(p.pol == Polarization::E || p.pol == Polarization::O, "Valid mode");
  near(p.e1.length(), 1, 1e-12, "Unit e1");
  near(p.e2.length(), 1, 1e-12, "Unit e2");
  near(p.n.length(), 1, 1e-12, "Unit plane normal");
  near(dot(p.e1, p.e2), 0, 1e-12, "e1 perpendicular to e2");
  near(dot(p.e1, p.n), 0, 1e-12, "e1 perpendicular to n");
  near(dot(p.e2, p.n), 0, 1e-12, "e2 perpendicular to n");
  near((cross(p.e1, p.e2) + (-1.0) * p.n).length(), 0, 1e-12,
       "Right-handed orbital basis");
  auto k = std::cos(p.alpha0) * p.e1 + std::sin(p.alpha0) * p.e2;
  near(k.length(), 1, 1e-12, "Unit initial photon direction");
  require(dot(k, p.e1) >= -1e-15, "Photon points outward");
}

double
planck_cdf(double upper) {
  // Independent integration of the photon-number density, without the mixture.
  auto density = [](double x) {
    return x == 0 ? 0 : x * x / (2 * zeta3 * std::expm1(x));
  };
  constexpr int intervals = 8192;
  double h = upper / intervals;
  double sum = density(0) + density(upper);
  for (int i = 1; i < intervals; ++i)
    sum += (i % 2 ? 4 : 2) * density(i * h);
  return sum * h / 3;
}

void
test_distribution(InitialPhotonOptions const& options, unsigned long long seed,
                  std::string_view name) {
  constexpr int count = 200000;
  constexpr std::array energy_edges{0.5, 1.0, 3.0, 10.0};
  constexpr std::array area_edges{0.25, 0.5, 0.75};
  std::array<int, energy_edges.size()> energy_hits{};
  std::array<int, area_edges.size()> area_hits{}, angle_hits{};
  int north = 0, e_mode = 0;
  double sum_x = 0, sum_x2 = 0, sum_mu = 0, sum_mu2 = 0;
  double sum_area = 0, sum_cos_phi = 0, sum_sin_phi = 0;
  double sum_cos_azimuth = 0, sum_cos_azimuth2 = 0, sum_sin_azimuth = 0;
  double cap_area = 1 - std::cos(options.cap_half_angle);
  Ran ran(seed);
  for (int i = 0; i < count; ++i) {
    auto p = sample_initial_photon_fd11(ran, rs, R_star, T_inf, options);
    check_geometry(p);
    double area = (1 - std::abs(p.e1.z)) / cap_area;
    require(area >= -1e-12 && area <= 1 + 1e-12,
            "Position inside selected cap");
    double mu = std::cos(p.alpha0), x = p.omega_inf / T_inf;
    sum_area += area;
    sum_mu += mu;
    sum_mu2 += mu * mu;
    sum_x += x;
    sum_x2 += x * x;
    north += p.e1.z > 0;
    e_mode += p.pol == Polarization::E;
    for (std::size_t j = 0; j < energy_edges.size(); ++j)
      energy_hits[j] += x <= energy_edges[j];
    for (std::size_t j = 0; j < area_edges.size(); ++j) {
      area_hits[j] += area <= area_edges[j];
      angle_hits[j] += mu <= area_edges[j];
    }
    double rho = std::hypot(p.e1.x, p.e1.y);
    require(rho > 0, "Nonpolar position in statistical sample");
    sum_cos_phi += p.e1.x / rho;
    sum_sin_phi += p.e1.y / rho;
    double3 phi_hat{-p.e1.y / rho, p.e1.x / rho, 0};
    double cos_azimuth = dot(p.e2, phi_hat);
    sum_cos_azimuth += cos_azimuth;
    sum_cos_azimuth2 += cos_azimuth * cos_azimuth;
    sum_sin_azimuth += dot(p.e2, cross(p.e1, phi_hat));
  }

  probability_near(north, count, options.north_fraction, "North cap fraction");
  probability_near(e_mode, count, options.e_mode_fraction, "E mode fraction");
  mean_near(sum_area / count, 0.5, 1.0 / 12, count, "Uniform surface area");
  mean_near(sum_cos_phi / count, 0, 0.5, count,
            "Uniform position azimuth (cos)");
  mean_near(sum_sin_phi / count, 0, 0.5, count,
            "Uniform position azimuth (sin)");
  mean_near(sum_cos_azimuth / count, 0, 0.5, count, "Emission azimuth (cos)");
  mean_near(sum_sin_azimuth / count, 0, 0.5, count, "Emission azimuth (sin)");
  mean_near(sum_cos_azimuth2 / count, 0.5, 1.0 / 8, count,
            "Emission azimuth second moment");
  mean_near(sum_mu / count, 2.0 / 3, 1.0 / 18, count, "Mean emission cosine");
  mean_near(sum_mu2 / count, 0.5, 1.0 / 12, count,
            "Mean squared emission cosine");
  mean_near(sum_x / count, mean_x, mean_x2 - mean_x * mean_x, count,
            "Mean blackbody photon energy / T_inf");
  mean_near(sum_x2 / count, mean_x2, mean_x4 - mean_x2 * mean_x2, count,
            "Second blackbody energy moment");
  for (std::size_t j = 0; j < energy_edges.size(); ++j)
    probability_near(energy_hits[j], count, planck_cdf(energy_edges[j]),
                     "Planck photon-number CDF");
  for (std::size_t j = 0; j < area_edges.size(); ++j) {
    probability_near(area_hits[j], count, area_edges[j], "Uniform area CDF");
    probability_near(angle_hits[j], count, area_edges[j] * area_edges[j],
                     "Cosine-weighted emission CDF");
  }
  std::cout << name << ": " << count
            << " photons, mean energy/T_inf = " << std::setprecision(10)
            << sum_x / count << " (theory " << mean_x << ")\n";
}

void
test_boundaries_and_geodesic() {
  std::array<InitialPhotonOptions, 5> configurations{{
      {},
      {5 * pi / 180, 0, 0},
      {5 * pi / 180, 1, 1},
      {1e-12, 0.5,
       0.5}, // z rounds to +/-1; the transverse position is nonzero.
      {std::numeric_limits<double>::min(), 0.5, 0.5}, // Exact poles in double.
  }};
  Ran ran(1033);
  Geodesic geodesic(M_star, R_star);
  for (auto const& options : configurations) {
    for (int i = 0; i < 128; ++i) {
      auto p = sample_initial_photon_fd11(ran, rs, R_star, T_inf, options);
      check_geometry(p);
      require(options.north_fraction != 0 || p.e1.z < 0, "South cap endpoint");
      require(options.north_fraction != 1 || p.e1.z > 0, "North cap endpoint");
      require(options.e_mode_fraction != 0 || p.pol == Polarization::O,
              "Pure O endpoint");
      require(options.e_mode_fraction != 1 || p.pol == Polarization::E,
              "Pure E endpoint");
      geodesic.init(p.r0, p.psi0, p.alpha0);
      require(!geodesic.step(), "No spurious surface absorption on emission");
      auto [r, psi, alpha] = geodesic.interp(geodesic.l_new());
      require(std::isfinite(r) && r > R_star && std::isfinite(psi) &&
                  psi >= 0 && std::isfinite(alpha) && alpha >= 0 &&
                  alpha <= p.alpha0,
              "First geodesic step propagates outward");
      double initial_b =
          R_star * std::sin(p.alpha0) / std::sqrt(1 - rs / R_star);
      double final_b = r * std::sin(alpha) / std::sqrt(1 - rs / r);
      near(final_b, initial_b, 1e-8,
           "Conserved Schwarzschild impact parameter");
    }
  }
}

void
test_reproducibility_and_scaling() {
  Ran base(182), repeat(182), hot(182), flat(182), compact(182);
  for (int i = 0; i < 128; ++i) {
    auto a = sample_initial_photon_fd11(base, rs, R_star, T_inf);
    for (auto b : {sample_initial_photon_fd11(repeat, rs, R_star, T_inf),
                   sample_initial_photon_fd11(flat, 0, R_star, T_inf),
                   sample_initial_photon_fd11(
                       compact, std::nextafter(R_star, 0), R_star, T_inf)}) {
      same_geometry_and_mode(a, b);
      require(a.omega_inf == b.omega_inf,
              "Reproducible energy independent of rs");
    }
    auto b = sample_initial_photon_fd11(hot, rs, R_star, 2 * T_inf);
    same_geometry_and_mode(a, b);
    require(b.omega_inf == 2 * a.omega_inf,
            "Energy scales linearly with T_inf");
  }
}

template <class Exception, class F>
void
expect_throw(F const& f, std::string_view message) {
  try {
    f();
  } catch (Exception const&) {
    return;
  }
  throw std::runtime_error(std::string(message));
}

void
test_invalid_parameters() {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  Ran ran(771), original = ran;
  auto invalid = [&](double radius_s, double radius, double temperature,
                     InitialPhotonOptions options = {}) {
    expect_throw<std::invalid_argument>(
        [&] {
          sample_initial_photon_fd11(ran, radius_s, radius, temperature,
                                     options);
        },
        "Invalid photon parameters must throw invalid_argument");
    require(ran.u == original.u && ran.v == original.v && ran.w == original.w,
            "Invalid parameters must not consume random numbers");
  };
  for (double value : {-1.0, 0.0, inf, -inf, nan}) {
    invalid(rs, value, T_inf);
    invalid(rs, R_star, value);
  }
  for (double value : {-1.0, R_star, 2 * R_star, inf, -inf, nan})
    invalid(value, R_star, T_inf);
  for (double value : {-1.0, 0.0, std::nextafter(pi / 2, inf), inf, -inf, nan})
    invalid(rs, R_star, T_inf, {.cap_half_angle = value});
  for (double value : {-0.1, std::nextafter(1.0, inf), inf, -inf, nan}) {
    invalid(rs, R_star, T_inf, {.north_fraction = value});
    invalid(rs, R_star, T_inf, {.e_mode_fraction = value});
  }
}

void
test_random_endpoints_and_energy_range() {
  // These public Ran states make the next int64() exactly 0 or UINT64_MAX:
  // the u state inverts its LCG to zero, and the v state inverts its xorshifts.
  for (bool upper : {false, true}) {
    Ran ran(1);
    ran.u = 13090081054908427669ULL;
    ran.v = upper ? 9259612767413907632ULL : 0;
    ran.w = 0;
    Ran probe = ran, photon_rng = ran;
    require(probe.U() == (upper ? 1.0 : 0.0), "Random endpoint fixture");
    double u = initialize_fd11_detail::open_unit(ran);
    require(u > 0 && u < 1, "Redraw random endpoints");
    check_geometry(sample_initial_photon_fd11(photon_rng, rs, R_star, T_inf));
  }

  for (double temperature : {std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::denorm_min()}) {
    Ran ran(7032);
    int unrepresentable = 0, valid = 0;
    for (int i = 0; i < 2000; ++i) {
      InitializedPhoton p;
      try {
        p = sample_initial_photon_fd11(ran, rs, R_star, temperature);
      } catch (std::runtime_error const&) {
        ++unrepresentable;
        continue;
      }
      check_geometry(p);
      ++valid;
    }
    require(unrepresentable > 0 && valid > 0,
            "Energy range test must exercise both rejection and valid samples");
  }
}

} // namespace

int
main() {
  try {
    test_invalid_parameters();
    test_random_endpoints_and_energy_range();
    test_reproducibility_and_scaling();
    test_boundaries_and_geodesic();
    test_distribution({}, 7774, "Full surface, pure E");
    test_distribution({5 * pi / 180, 0.5, 0.5}, 10345,
                      "Equal 5-degree caps, unpolarized");
    test_distribution({5 * pi / 180, 0.3, 0.25}, 9081,
                      "Unequal caps and E/O mixture");
    std::cout << "All FD11 initialization checks passed.\n";
    return 0;
  } catch (std::exception const& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
