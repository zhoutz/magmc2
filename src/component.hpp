#pragma once

#include "distribution.hpp"

template <ParticleDistribution Dist>
struct CurrentComponent {
  Dist dist;
  int charge_sign;
  double current_fraction;
  double charge_mass_ratio;
};
