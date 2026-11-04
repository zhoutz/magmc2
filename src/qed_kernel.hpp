#pragma once

#include "constants.hpp"
#include "polarization.hpp"
#include "qed_numerics.hpp"
#include <complex>

namespace qed {
using Complex = std::complex<double>;
constexpr Complex I{0, 1};
using Spinor = std::array<Complex, 4>;
using Matrix4 = std::array<Spinor, 4>;
using Jones = std::array<std::array<Complex, 2>, 2>;
// Hermitian entries (OO, EE, Re(OE), Im(OE)). Keeping the two positive
// diagonals separately avoids subtracting s0-sQ for a strongly suppressed E
// mode. The equivalent Pauli coefficients are ((OO+EE)/2,(OO-EE)/2,Re,Im).
using Response = Estimate<4>;
inline constexpr char kernel_version[] = "lll-n01-st-width-v1";

inline
double
contract(std::array<double, 4> const& a, Stokes s) {
  return (1 + s.x) * a[0] / 2 + (1 - s.x) * a[1] / 2 + a[2] * s.y + a[3] * s.z;
}
inline
double
response_error(Response const& r) {
  return std::max(r.error[0], r.error[1]) + std::hypot(r.error[2], r.error[3]);
}
inline
double
envelope(Response const& r) {
  double mean = (r.value[0] + r.value[1]) / 2;
  double length =
      std::hypot((r.value[0] - r.value[1]) / 2, r.value[2], r.value[3]);
  require(std::isfinite(mean + length) &&
              mean - length >= -response_error(r) - 1e-13 * std::abs(mean),
          Failure::Convergence, "QED response is not positive semidefinite");
  return mean + length;
}
inline
void
check_stokes(Stokes s) {
  require(std::isfinite(s.length()) && s.length() <= 1 + 1e-10,
          Failure::InvalidInput, "QED: nonphysical incident Stokes vector");
}
inline
Jones
density(Stokes s) {
  check_stokes(s);
  return {{{(1 + s.x) / 2, Complex(s.y, s.z) / 2.},
           {Complex(s.y, -s.z) / 2., (1 - s.x) / 2}}};
}
inline
Response
response(Jones const& m) {
  double a = 0, d = 0;
  Complex c = 0;
  for (int k = 0; k < 2; ++k) {
    a += std::norm(m[k][0]);
    d += std::norm(m[k][1]);
    c += std::conj(m[k][0]) * m[k][1];
  }
  return {{a, d, c.real(), c.imag()}, {}};
}
inline
Stokes
outgoing_stokes(Jones const& m, Stokes input) {
  auto rho = density(input);
  Jones a{};
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          a[i][j] += m[i][k] * rho[k][l] * std::conj(m[j][l]);
  double trace = (a[0][0] + a[1][1]).real();
  require(trace > 0 && std::isfinite(trace), Failure::Convergence,
          "QED: selected zero-weight outgoing state");
  Stokes s{(a[0][0] - a[1][1]).real() / trace, 2 * a[0][1].real() / trace,
           2 * a[0][1].imag() / trace};
  check_stokes(s);
  return s; // In particular, do not normalize a mixed state's Stokes length.
}

inline
Spinor
spinor(int n, int spin, int eta, double b, double q) {
  double p = std::sqrt(2 * n * b), a = std::sqrt(1 + 2 * n * b),
         E = std::hypot(a, q);
  double N = 2 * std::sqrt(E * a * (a + 1) * (E + a));
  Spinor c;
  if (eta == 1 && spin == -1)
    c = {-I * q * p, (E + a) * (a + 1), -I * p * (E + a), -q * (a + 1)};
  else if (eta == 1)
    c = {(E + a) * (a + 1), -I * q * p, q * (a + 1), I * p * (E + a)};
  else if (spin == -1)
    c = {-p * (E + a), I * q * (a + 1), q * p, I * (E + a) * (a + 1)};
  else
    c = {I * q * (a + 1), -p * (E + a), -I * (E + a) * (a + 1), -q * p};
  for (auto& z : c)
    z /= N;
  return c;
}
inline
Matrix4
coupling(double3 e) {
  Complex minus = e.x - I * e.y, plus = e.x + I * e.y;
  return {{{0, 0, e.z, minus},
           {0, 0, plus, -e.z},
           {e.z, minus, 0, 0},
           {plus, -e.z, 0, 0}}};
}
inline
Complex
overlap(int l, int m, double b, double kx, double aL, double aR) {
  if (l < 0 || m < 0)
    return 0;
  require(l + m <= 1, Failure::InvalidInput, "Unsupported oscillator overlap");
  double d = std::sqrt(b) * (aL - aR), v = kx / std::sqrt(b);
  Complex g = std::exp(Complex(-(d * d + v * v) / 4, kx * (aL + aR) / 2));
  if (l == 1)
    return g * (-d + I * v) / std::sqrt(2.);
  if (m == 1)
    return g * (d + I * v) / std::sqrt(2.);
  return g;
}
inline
int
order(int n, int j) {
  return n - (j % 2 == 0 ? 1 : 0);
}
inline
Complex
vertex(Spinor const& l, int nl, double al, Spinor const& r, int nr, double ar,
       double b, double3 k, double3 e) {
  auto A = coupling(e);
  Complex result = 0;
  for (int j = 0; j < 4; ++j)
    for (int h = 0; h < 4; ++h) {
      if (A[j][h] == Complex{} || l[j] == Complex{} || r[h] == Complex{})
        continue;
      result += std::conj(l[j]) * A[j][h] * r[h] *
                overlap(order(nl, j), order(nr, h), b, k.x, al, ar);
    }
  return result;
}

struct RestWidths {
  double minus = 0, plus = 0, minus_error = 0, plus_error = 0;
  double
  at(int spin) const {
    return spin == -1 ? minus : plus;
  }
};
inline
RestWidths
rest_widths(double b, double rtol = 1e-9) {
  require(b > 0 && std::isfinite(b), Failure::InvalidInput,
          "QED widths require b>0");
  double a = std::sqrt(1 + 2 * b);
  double Ad = std::sqrt((a + 1) / (2 * a)), Bd = std::sqrt(b / (a * (a + 1)));
  // Scale the two integrands separately to retain the b^3 spin-plus width.
  auto f = [&](double nu) {
    double w = 2 * b / (a + std::sqrt(1 + 2 * b * nu * nu));
    double E = std::hypot(1., w * nu), J = 1 + w * nu * nu / E;
    double F = std::sqrt((E + 1) / (2 * E)),
           G = w * nu / std::sqrt(2 * E * (E + 1));
    double sn = std::sqrt((1 - nu) * (1 + nu));
    double h = w * sn / std::sqrt(2 * b),
           g = std::exp(-w * w * sn * sn / (4 * b));
    double om = g * (F * Bd * nu + G * Ad * h * sn), em = -g * F * Bd;
    double op = -g * (F * Bd * h * sn + G * Ad * nu), ep = -g * G * Ad;
    return Estimate<2>{
        {fine_structure * w / J * (om * om + em * em) / (b * b),
         fine_structure * w / J * (op * op + ep * ep) / (b * b * b)},
        {}};
  };
  auto r = integrate<2>(f, {-1, 0, 1}, rtol);
  RestWidths result{r.value[0] * b * b, r.value[1] * b * b * b,
                    r.error[0] * b * b, r.error[1] * b * b * b};
  require(result.minus > 0 && result.plus > 0, Failure::Convergence,
          "QED natural width underflow");
  return result;
}

struct Kinematics {
  double y, qf, Ef, J;
};
inline
void
check_domain(double x, double mui) {
  require(x > 0 && std::isfinite(x) && std::isfinite(mui) && std::abs(mui) <= 1,
          Failure::InvalidInput, "QED requires x>0 and |mu_i|<=1");
  double C = 2 * x + x * x * (1 - mui) * (1 + mui);
  require(C < 8, Failure::PairDomain,
          std::format("QED full-angle domain: x={}, mu_i={}, C={} >= 8", x, mui,
                      C));
}
inline
Kinematics
recoil(double x, double mui, double muf) {
  require(std::isfinite(muf) && std::abs(muf) <= 1, Failure::InvalidInput,
          "QED requires |mu_f|<=1");
  double C = 2 * x + x * x * (1 - mui) * (1 + mui), D = 1 + x * (1 - mui * muf);
  double delta = (1 - muf) * (1 + muf) + std::pow(muf + x * (muf - mui), 2);
  double y = C / (D + std::sqrt(delta)), q = x * mui - y * muf,
         E = std::hypot(1., q);
  return {y, q, E, 1 - q * muf / E};
}
inline
std::array<double3, 2>
real_basis(double mu, double phi) {
  return {{{-mu * std::cos(phi), -mu * std::sin(phi),
            std::sqrt((1 - mu) * (1 + mu))},
           {std::sin(phi), -std::cos(phi), 0}}};
}

// Exposed spectral amplitude is also the independent reference for tests.
// initial_center changes only a common phase. max_level=0 and zero widths
// are algebraic test controls, never production resonance prescriptions.
inline
std::array<Jones, 2>
amplitude_terms(double b, double x, double mui, double muf, double phi,
                RestWidths widths, int max_level = 1,
                double initial_center = 0) {
  auto kin = recoil(x, mui, muf);
  double y = kin.y;
  double3 ki{x * std::sqrt((1 - mui) * (1 + mui)), 0, x * mui};
  double3 kf{y * std::sqrt((1 - muf) * (1 + muf)) * std::cos(phi),
             y * std::sqrt((1 - muf) * (1 + muf)) * std::sin(phi), y * muf};
  auto ei = real_basis(mui, 0), ef = real_basis(muf, phi);
  auto ci = spinor(0, -1, 1, b, 0), cf = spinor(0, -1, 1, b, kin.qf);
  double af = initial_center + (kf.y - ki.y) / b;
  std::array<Jones, 2> out{};
  for (int diagram = 0; diagram < 2; ++diagram) {
    double V = diagram == 0 ? 1 + x : 1 - y, q = diagram == 0 ? ki.z : -kf.z;
    double ac = initial_center + (diagram == 0 ? -ki.y : kf.y) / b;
    for (int n = 0; n <= max_level; ++n) {
      double a = std::sqrt(1 + 2 * n * b), E = std::hypot(a, q);
      double invariant =
          diagram == 0 ? 2 * x + x * x * (1 - mui) * (1 + mui) - 2 * n * b
                       : -2 * y + y * y * (1 - muf) * (1 + muf) - 2 * n * b;
      for (int spin : {-1, 1}) {
        if (n == 0 && spin == 1)
          continue;
        double width = n == 0 ? 0 : widths.at(spin) * a / E;
        for (int eta : {-1, 1}) {
          double detuning = V - eta * E;
          if (std::abs(detuning) < 0.5 * E)
            detuning = invariant / (V + eta * E);
          Complex D{detuning, eta * width / 2};
          require(std::norm(D) > 0, Failure::Convergence,
                  "QED unresolved zero-width pole");
          auto c = spinor(n, spin, eta, b, q);
          std::array<Complex, 2> left, right;
          for (int p = 0; p < 2; ++p) {
            if (diagram == 0) {
              left[p] = vertex(cf, 0, af, c, n, ac, b, -1. * kf, ef[p]);
              right[p] = vertex(c, n, ac, ci, 0, initial_center, b, ki, ei[p]);
            } else {
              left[p] = vertex(cf, 0, af, c, n, ac, b, ki, ei[p]);
              right[p] =
                  vertex(c, n, ac, ci, 0, initial_center, b, -1. * kf, ef[p]);
            }
          }
          for (int row = 0; row < 2; ++row)
            for (int col = 0; col < 2; ++col)
              out[diagram][row][col] +=
                  (diagram == 0 ? left[row] * right[col]
                                : left[col] * right[row]) /
                  D;
        }
      }
    }
  }
  return out;
}
inline
Jones
amplitude(double b, double x, double mui, double muf, double phi,
          RestWidths widths, int max_level = 1, double initial_center = 0) {
  auto terms =
      amplitude_terms(b, x, mui, muf, phi, widths, max_level, initial_center);
  auto out = terms[0];
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      out[i][j] += terms[1][i][j];
  return out;
}

// Exact finite Fourier representation of the vertices at fixed outgoing mu.
// After removing each guiding-centre phase, n=0 contributes harmonic 0 and
// n=1 contributes harmonics +/-1: e_z I_01 and e_+/- I_00 are linear in
// exp(+/-i phi). Three direct spectral sums therefore recover ALL coefficients,
// not an interpolated physical response. Quadrature still evaluates positive
// M^dagger M at every azimuth and resolves the full exp(+/-i zeta)
// interference.
struct AzimuthKernel {
  std::array<std::array<Jones, 3>, 2> coefficients{};
  double phase_scale, factor;
  int charge;
  AzimuthKernel(int charge, double b, double x, double mui, double muf,
                RestWidths widths)
      : charge(charge) {
    if (charge == 1) {
      mui = -mui;
      muf = -muf;
    }
    auto k = recoil(x, mui, muf);
    double ki = x * std::sqrt((1 - mui) * (1 + mui));
    double kf = k.y * std::sqrt((1 - muf) * (1 + muf));
    phase_scale = ki * kf / (2 * b);
    factor = std::sqrt(k.y / (x * k.J));
    for (int node = 0; node < 3; ++node) {
      double phi = 2 * pi * node / 3;
      double common = (ki - kf * std::cos(phi)) * kf * std::sin(phi) / (2 * b);
      auto terms = amplitude_terms(b, x, mui, muf, phi, widths);
      for (int d = 0; d < 2; ++d)
        for (int h = -1; h <= 1; ++h) {
          Complex weight =
              std::exp(-I * (common +
                             (d == 0 ? -1. : 1.) * phase_scale * std::sin(phi) +
                             h * phi)) /
              3.;
          for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
              coefficients[d][h + 1][i][j] += weight * terms[d][i][j];
        }
    }
  }
  Jones
  reduced(double phi) const {
    if (charge == 1)
      phi = -phi;
    double cp = std::cos(phi), sp = std::sin(phi);
    double cz = std::cos(phase_scale * sp), sz = std::sin(phase_scale * sp);
    Jones out{};
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j) {
        std::array<Complex, 2> term;
        for (int d = 0; d < 2; ++d)
          term[d] =
              coefficients[d][1][i][j] +
              cp * (coefficients[d][0][i][j] + coefficients[d][2][i][j]) +
              I * sp * (coefficients[d][2][i][j] - coefficients[d][0][i][j]);
        out[i][j] =
            factor * (cz * (term[0] + term[1]) + I * sz * (term[1] - term[0]));
      }
    return out; // Common phase omitted; all density-matrix operations
                // invariant.
  }
  Response
  integral_bound() const {
    // Triangle inequality for the exact six Fourier coefficients, followed
    // by Cauchy--Schwarz for the off-diagonal entry. Bounds the ENTIRE
    // azimuthal integral, not just values at a few sampling points.
    Response bound;
    for (int col = 0; col < 2; ++col) {
      double norm = 0;
      for (int row = 0; row < 2; ++row) {
        double sum = 0;
        for (auto const& diagram : coefficients)
          for (auto const& harmonic : diagram)
            sum += std::abs(harmonic[row][col]);
        norm += sum * sum;
      }
      bound.error[col] = std::max(2 * pi * factor * factor * norm,
                                  std::numeric_limits<double>::denorm_min());
    }
    bound.error[2] = bound.error[3] =
        std::sqrt(bound.error[0]) * std::sqrt(bound.error[1]);
    return bound;
  }
};

struct DifferentialKernel {
  Kinematics kinematics;
  Jones M;       // cm
  Jones reduced; // M/r_e, for well-scaled quadrature and state conditioning
};
inline
DifferentialKernel
differential(int charge, double b, double x, double mui, double muf, double phi,
             RestWidths widths) {
  require((charge == -1 || charge == 1) && b > 0 && std::isfinite(b) &&
              std::isfinite(phi),
          Failure::InvalidInput, "QED: invalid charge, field or azimuth");
  check_domain(x, mui);
  auto kin = recoil(x, mui, muf);
  auto t = charge == -1 ? amplitude(b, x, mui, muf, phi, widths)
                        : amplitude(b, x, -mui, -muf, -phi, widths);
  DifferentialKernel result{kin, {}, t};
  double factor = std::sqrt(kin.y / (x * kin.J));
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      result.reduced[i][j] *= factor;
      result.M[i][j] = classical_electron_radius_cm * result.reduced[i][j];
      require(std::isfinite(std::abs(result.reduced[i][j])),
              Failure::Convergence, "QED nonfinite differential amplitude");
    }
  return result;
}
} // namespace qed
