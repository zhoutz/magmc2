#pragma once

#include "distribution.hpp"

template <ParticleDistribution Dist>
struct CurrentComponent {
  Dist dist; // Signed velocity along B: evaluate directly at beta.
  int charge_sign;
  double current_fraction; // Nonnegative current share; not auto-normalized.
  double charge_mass_ratio;
};
