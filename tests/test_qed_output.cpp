#include "qed_numerics.hpp"
#include "result.hpp"
#include <filesystem>
#include <iostream>

int
main() try {
  std::filesystem::create_directories("build/test-fixtures");
  std::string path = "build/test-fixtures/qed-output.h5";
  std::filesystem::remove(path);
  {
    OutputHDF5 output(path);
    output.file.createAttribute("scattering_kernel",
                                std::string("lll-n01-st-width-v1"));
    ResultPhotonGroup group(42);
    // A domain failure is one failed photon, followed by a successful photon.
    group.add_failed(true);
    group.add_failed(false);
    group.add_absorbed();
    group.add_escaped({.id = 3,
                       .n_scattered = 1,
                       .energy_inf = .7,
                       .muk = .2,
                       .Q = .3,
                       .U = .4,
                       .V = .5});
    output.write_photon_group("seed_42", group);
    output.file.flush();
  }
  {
    HighFive::File file(path, HighFive::File::ReadOnly);
    auto dataset = file.getDataSet("seed_42");
    int failed = 0, domain = 0, absorbed = 0, escaped = 0;
    dataset.getAttribute("n_failed").read(failed);
    dataset.getAttribute("n_domain_failed").read(domain);
    dataset.getAttribute("n_absorbed").read(absorbed);
    dataset.getAttribute("n_escaped").read(escaped);
    if (failed != 2 || domain != 1 || absorbed != 1 || escaped != 1)
      throw std::runtime_error("Incorrect failure accounting");
    std::vector<ResultPhoton> rows;
    dataset.read(rows);
    if (rows.size() != 1 || rows[0].id != 3 || rows[0].V != .5)
      throw std::runtime_error("Photon output schema did not round-trip");
  }
  std::cout << "QED HDF5 failure accounting and photon fields passed\n";
} catch (std::exception const& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
