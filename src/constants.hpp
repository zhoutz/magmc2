#pragma once

constexpr double pi = 3.14159265358979323846;

// https://pdg.lbl.gov/2024/reviews/rpp2024-rev-astrophysical-constants.pdf
constexpr double schwarzschild_radius_of_sun_in_km = 2.9532501; // km

// compute electron's cyclotron frequency, B in G(Gauss), energy in keV.
constexpr double B_to_energy = 1.157676359638892e-11;

// Gaussian cgs charge; microscopic QED energies are in electron rest masses.
constexpr double electron_rest_keV = 510.99895;
constexpr double critical_field_gauss = 4.414e13;
constexpr double fine_structure = 7.2973525693e-3;
constexpr double classical_electron_radius_cm = 2.81794032e-13;
constexpr double elementary_charge_esu = 4.80320471257e-10;
constexpr double thomson_cross_section_cm2 =
    (8 * pi / 3) * classical_electron_radius_cm * classical_electron_radius_cm;
