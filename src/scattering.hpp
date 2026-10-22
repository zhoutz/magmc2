#pragma once

#include "distribution.hpp"

template <ParticleDistribution Dist>
struct ResonantScattering {
  Dist dist;
  int charge_sign;
  double current_fraction;
};
