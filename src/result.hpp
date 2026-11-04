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
  int n_absorbed, n_escaped, n_failed, n_domain_failed;
  int seed;
  ResultPhotonGroup(int seed)
      : n_absorbed(0), n_escaped(0), n_failed(0), n_domain_failed(0), seed(seed) {}
  void
  add_failed(bool domain_failure = false) {
    ++n_failed;
    if (domain_failure) ++n_domain_failed;
  }
  void
  add_absorbed() {
    ++n_absorbed;
  }
  void
  add_escaped(ResultPhoton const& photon) {
    photons.push_back(photon);
    ++n_escaped;
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
  write_photon_group(std::string const& dset_name,
                     ResultPhotonGroup const& group) {
    // A one-dimensional compound dataset: one ResultPhoton per element.
    auto dset = file.createDataSet(dset_name, group.photons);
    dset.createAttribute("n_absorbed", group.n_absorbed);
    dset.createAttribute("n_escaped", group.n_escaped);
    dset.createAttribute("n_failed", group.n_failed);
    dset.createAttribute("n_domain_failed", group.n_domain_failed);
    dset.createAttribute("seed", group.seed);
  }
};
