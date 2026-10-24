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
  eigenmode(AdiabaticPolarization apol, double cos2chi, double sin2chi) {
    double sign = (apol == AdiabaticPolarization::O) ? 1 : -1;
    return {sign * cos2chi, sign * sin2chi, 0};
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
};
