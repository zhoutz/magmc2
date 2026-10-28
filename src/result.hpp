#pragma once

#include "highfive/highfive.hpp" // IWYU pragma: export
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

struct ResultPhoton {
  int id, n_scattered;
  double energy_inf, muk;
  double Q, U, V;
};

struct ResultPhotonGroup {
  std::vector<ResultPhoton> photons;
  int n_absorbed, n_escaped;
  int seed;
  ResultPhotonGroup(int seed) : n_absorbed(0), n_escaped(0), seed(seed) {}
  void
  add_absorbed() {
    ++n_absorbed;
  }
  void
  add_escaped(ResultPhoton const& photon) {
    ++n_escaped;
    photons.push_back(photon);
  }
};

inline
HighFive::CompoundType
create_result_photon_type() {
  static_assert(std::is_standard_layout_v<ResultPhoton>);
  static_assert(std::is_trivially_copyable_v<ResultPhoton>);
  // Match the C++ memory layout, including any padding between records.
  return HighFive::CompoundType(
      {{"id", HighFive::create_datatype<int>(), offsetof(ResultPhoton, id)},
       {"n_scattered", HighFive::create_datatype<int>(),
        offsetof(ResultPhoton, n_scattered)},
       {"energy_inf", HighFive::create_datatype<double>(),
        offsetof(ResultPhoton, energy_inf)},
       {"muk", HighFive::create_datatype<double>(),
        offsetof(ResultPhoton, muk)},
       {"Q", HighFive::create_datatype<double>(), offsetof(ResultPhoton, Q)},
       {"U", HighFive::create_datatype<double>(), offsetof(ResultPhoton, U)},
       {"V", HighFive::create_datatype<double>(), offsetof(ResultPhoton, V)}},
      sizeof(ResultPhoton));
}

HIGHFIVE_REGISTER_TYPE(ResultPhoton, create_result_photon_type)

struct OutputHDF5 {
  HighFive::File file;
  OutputHDF5(std::string const& filename)
      : file(filename, HighFive::File::OpenOrCreate) {}

  void
  write_photon_group(std::string const& group_name,
                     ResultPhotonGroup const& group) {
    auto group_h5 = file.createGroup(group_name);
    group_h5.createAttribute("n_absorbed", group.n_absorbed);
    group_h5.createAttribute("n_escaped", group.n_escaped);
    group_h5.createAttribute("seed", group.seed);
    // A one-dimensional compound dataset: one ResultPhoton per element.
    group_h5.createDataSet("photons", group.photons);
  }
};
