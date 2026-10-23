# Photon evolution

Run from the repository root:

```sh
make main
make test
```

The build uses C++23, GCC 16 and GSL (Homebrew paths by default). The magnetic
field table `table/bfield_t10.txt` must already exist. `make main` initializes
one photon, calls `evolve()` once, and performs one collision if the returned
event is `Scattered`. It does not run a photon ensemble.

## State and event interface

`PhotonEvolution` in `src/main.cpp` owns its tuple of `CurrentComponent`s.
The field and initializer remain caller-owned references. `init()` samples a
source photon; `init(InitializedPhoton const&)` starts an explicitly specified
free flight. Each initialization sets `tau=0` and draws `tau_target=-log(U)`.

`evolve()` returns the first `Absorbed`, `Escaped`, or `Scattered` event. For a
scattering event, `scattered_data` records the local geodesic state, path length,
component index, signed resonant velocity, root index (0: lower beta; 1: upper
beta), and incident polarization. Call `perform_scattering()` before evolving
again: it changes the direction, energy and normal mode, rebuilds the orbital
basis, and draws a new target optical depth without moving the interaction
point. A terminal event requires reinitialization before further evolution.

Lengths are in km, fields in gauss, and `omega_inf` is photon energy at infinity
in keV. `V=+2 Im(A1 conj(A2))`. Between collisions, the transported screen is
`s1=n cross k`, `s2=n`. Electrons use negative beta, protons positive beta,
with the positive-speed distribution evaluated at `charge_sign*beta`.

## Integration and numerical controls

The local opacity sums complete component/root weights, including current
fractions, mean-speed normalization, species-specific cyclotron energies,
resonance Jacobians, and the instantaneous polarization overlap. Component
selection uses these same weights evaluated at the located event, rather than
integrated interval weights or distribution values alone.

Both adiabatic intervals and individual accepted Magnus intervals are split at
the union of every active component's velocity knots and `D=0` crossings. The
Magnus dense-output midpoint join is also included. Defaults are eight speed
knots including support endpoints and eight scan intervals. Bisection locates
boundary crossings; direct GSL QAGS integration resolves each resulting piece.
At an optical-depth crossing, another bisection locates the scattering position
using the integral from that piece's left endpoint. Polarization dense output
is consumed before its accepted step is committed.

`Options` exposes quadrature tolerances (both `1e-8`), scattering/coupling root
position tolerance (`1e-10 km`), Stokes tolerances (both `1e-6`), coupling and
escape thresholds (both `1e-3`), knot/scan counts and maximum sampled phase
increment (`0.25 rad`). Specify these when constructing the evolution object;
knots are generated at construction. Refine scans/knots, quadrature tolerances,
and polarization steps separately when checking convergence. Boundary detection
uses sign changes on the scan mesh, as in the reference implementation; finite
scanning is not a proof of detecting arbitrarily narrow multiple crossings.

QAGS failures and nonphysical opacity/overlap values raise exceptions. The
integrator does not silently drop failed integrals or renormalize Stokes vectors.
The Magnus error estimate controls endpoint error, while the additional phase
cap resolves oscillations even when the input is an exact eigenmode. Existing
one-argument `do_step(omega)` calls retain their previous uncapped behavior.

## Escape and observer convention

Escape is checked at the end of a resolved Stokes propagation step. The ray
must be outgoing; every active species must satisfy the resonance-plane tangent
criterion, and `norm(delta S)*r/delta l < 1e-3`. Tangents use symmetric local RK4
trial propagation, including gravitational bending and redshift. There is no
separate frozen polarization state.

At escape, integrate the remaining Schwarzschild deflection to obtain `k_inf`.
The Stokes components at termination are transported in the orbital screen to
infinity, then rotated into the sky basis defined by the fixed global +z axis
projected perpendicular to `k_inf`. If its projection length is below `1e-12`,
use global +x. The output rotation does not consult the magnetic field.
`escaped_data` contains energy, `muk=k_inf.z`, sky Q/U/V, the escape position and
path length, and `k_inf`.

Both escape tests are local prescriptions from the note, not global error
bounds. Tests include continuation to four times the escape radius for a bent
ray, independent transformed quadrature near coalescing resonant roots,
refinement of scattering locations, outgoing angular/mode statistics, and
repeated electron/proton collisions through absorption or escape.
