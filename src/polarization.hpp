#pragma once

#include "doubles.hpp"
#include <variant>

enum class AdiabaticPolarization { O, E };
using Stokes = double3;
using Polarization = std::variant<AdiabaticPolarization, Stokes>;
