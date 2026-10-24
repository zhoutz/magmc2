#pragma once

#include "doubles.hpp"
// #include <variant>

enum class AdiabaticPolarization { O, E };
using Stokes = double3;

struct Polarization {
  bool is_adiabatic;
  AdiabaticPolarization apol;
  Stokes stokes;

  Polarization(AdiabaticPolarization pol = AdiabaticPolarization::E)
      : is_adiabatic(true), apol(pol), stokes(0, 0, 0) {}
  Polarization(Stokes stokes)
      : is_adiabatic(false), apol(AdiabaticPolarization::O), stokes(stokes) {}

  static
  Stokes
  eigenmode_B(AdiabaticPolarization apol) {
    double sign = (apol == AdiabaticPolarization::O) ? 1 : -1;
    return {sign, 0, 0};
  }

  static
  Stokes
  from_B(Stokes S, double c2, double s2) {
    return {c2 * S.x - s2 * S.y, s2 * S.x + c2 * S.y, S.z};
  }

  static
  Stokes
  to_B(Stokes S, double c2, double s2) {
    return {c2 * S.x + s2 * S.y, -s2 * S.x + c2 * S.y, S.z};
  }

  static
  Stokes
  eigenmode(AdiabaticPolarization apol, double cos2chi, double sin2chi) {
    return from_B(eigenmode_B(apol), cos2chi, sin2chi);
  }

  void
  set_adiabatic(AdiabaticPolarization pol) {
    is_adiabatic = true;
    apol = pol;
  }

  void
  set_stokes(Stokes s) {
    is_adiabatic = false;
    stokes = s;
  }

  void
  set_stokes(double cos2chi, double sin2chi) {
    set_stokes(eigenmode(apol, cos2chi, sin2chi));
  }

  Stokes
  in_B(double cos2chi, double sin2chi) const {
    if (is_adiabatic)
      return eigenmode_B(apol);
    else
      return to_B(stokes, cos2chi, sin2chi);
  }

  Stokes
  observer_stokes(double3 k, double3 normal) const {
    double3 axis{0, 0, 1};
    double3 u = axis + (-dot(axis, k)) * k;
    if (u.length() < 1e-12) {
      axis = {1, 0, 0};
      u = axis + (-dot(axis, k)) * k;
    }
    u = to_unit(u);
    double c = dot(u, cross(normal, k)), s = dot(u, normal);
    double norm = std::hypot(c, s);
    c /= norm;
    s /= norm;
    double c2 = c * c - s * s, s2 = 2 * c * s;
    return {c2 * stokes.x + s2 * stokes.y, -s2 * stokes.x + c2 * stokes.y,
            stokes.z};
  }
};
