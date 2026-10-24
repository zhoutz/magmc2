# FD11 Figure 7

This calculation uses the existing C++ transport, scattering and polarization
solvers to reproduce the `t10_g20_c90_Eu` model in Fernández & Davis (2011),
*ApJ* 730:131, Figure 7 (PDF page 10), Table 1 and Sections 3.4–4.1.

## Physical settings and conventions

| Setting | Value |
| --- | --- |
| Stellar radius | 10 km |
| Compactness | R / r_s = 3; mass = 10 / (3 × 2.9532501) solar masses |
| Polar field | 10^14 G |
| Twist | 1 rad, from `table/bfield_t10.txt` |
| Temperature at infinity | 0.4 keV |
| Emission | Uniform full surface, isotropic intensity, pure E mode |
| Electron momentum distribution | dP/du proportional to u^(-2), u = gamma beta |
| Speed interval | 0.2 < beta < sqrt(3)/2; gamma_max = 2 |
| Electron motion | Along +B, north to south |
| Electron current fraction, mass ratio | 1, 1 |

The initializer samples photon **number** from the Planck distribution and
outgoing flux with p(mu_surface) = 2 mu_surface. The energy is measured at
infinity; it must not be redshifted a second time. The paper's printed equation
(7) has a minus-sign convention inconsistent with its stated alpha = -2. The
decreasing u^(-2) distribution here agrees with the accompanying historical
`magmc_2019/src/physics.F90` routine `disfpw`, which uses u^alpha.

`PhotonEvolution` evaluates the distribution at `charge_sign * beta`.
`FD11Electrons` reflects the distribution argument and bounds while keeping
`charge_sign = -1`, so the physical electron velocity is positive along B.
Neither the magnetic field nor the plotted north/south coordinate is flipped.
This is a local adapter; the shared transport API is unchanged.

Q and U refer to the projected magnetic axis, with the second sky axis chosen
to form a right-handed basis with the photon direction. Individual records
have I = 1. Bin sums, not individual polarization fractions or angles, determine
the displayed observables:

    Pi_L = 100 hypot(sum_Q, sum_U) / N
    chi  = (0.5 atan2(sum_U, sum_Q) in degrees) modulo 180

## Build and run

Run from the repository root. The Makefile uses the installed Homebrew GCC 16
and GSL. It leaves the original `main` target unchanged.

```sh
make fd11_fig7
build/fd11_fig7 --help
```

The following Bash commands emit 1,000,000 photons in four independent batches.
Output files must not already exist; choose a different output path for reruns.
Each process has its own random generator and GSL workspace.

```bash
pids=()
for batch in 1 2 3 4; do
    build/fd11_fig7 --photons 250000 --seed "$((110700 + batch))" \
        --output "data/fd11_fig7_photons_${batch}.txt" &
    pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done
test "$status" -eq 0
```

`--photons` counts emitted photons, including those subsequently absorbed.
`--field PATH` overrides the magnetic table path; its twist must be 1 rad.
`--tight` changes eta_couple and escape_tolerance from 10^-3 to 3×10^-4,
quadrature absolute/relative tolerances from 10^-8 to 10^-9, and n_scan from
8 to 16. Other numerical settings remain those of the existing solver.

Every escaping photon is written, including energies outside the plotted
range. Returning photons are absorbed by the surface, not replaced. Exceptions
fail the run and do not silently discard a photon. A successful output ends
with `# complete = true` and accounting totals. Interrupted files are rejected
by the plotting script. Progress is reported to stderr every 30 seconds.

## Plot and data format

Use a Python environment with NumPy and Matplotlib. On this workstation:

```sh
MPLCONFIGDIR=/tmp/magmc2-mpl XDG_CACHE_HOME=/tmp/magmc2-cache \
  /opt/homebrew/Caskroom/miniforge/base/bin/python py/plot_fd11_fig7.py \
  data/fd11_fig7_photons_1.txt data/fd11_fig7_photons_2.txt \
  data/fd11_fig7_photons_3.txt data/fd11_fig7_photons_4.txt
```

The seven photon columns are:

    photon_id energy_keV mu Q U V n_scatter

Photon IDs are local to each batch and increase strictly; IDs of absorbed
photons do not appear. The file's seed identifies its independent realization.
Metadata includes the physical settings and numerical tolerances. Repeated
paths/seeds, incompatible settings, incomplete files and invalid photon
accounting are rejected. Data are processed in 50,000-row chunks.

The plot uses 120 logarithmic energy bins spanning 0.04–40 keV and 64 uniform
cos(theta) bins spanning -1–1. Total, unscattered and scattered columns share
each row's color scale. Empty bins are black, and zero net linear polarization
has an undefined angle. The count color scale uses actual photon counts; no
normalization to the paper's larger sample is applied. The angle color bar
has extensions to indicate values outside the paper's 90–125 degree display
range. There is no smoothing or low-count threshold.

Outputs, with prefix configurable through `--output`:

- `data/fd11_fig7.png`: nine-panel image.
- `data/fd11_fig7.pdf`: PDF with rasterized heatmaps and vector text.
- `data/fd11_fig7_bins.txt`: every bin for all three groups, including empty
  bins; columns are group, energy edges, mu edges, N, sum_Q, sum_U, sum_V,
  linear_percent and angle_deg. Undefined values are NaN.
- `data/fd11_fig7_validation.txt`: accounting and physical diagnostics from
  the completed production run.

The repository ignores `data/`, `build/` and `table/`; generated results are
local artifacts, not automatically tracked by Git.

## Validation and interpretation

```sh
/opt/homebrew/Caskroom/miniforge/base/bin/python -m unittest discover \
    -s py -p test_fd11_fig7.py -v
```

The deterministic tests cover summing Stokes vectors before deriving the
polarization, atan2 quadrants, cancellation, empty bins, bin endpoints,
Total = Unscattered + Scattered, incomplete files and accounting errors.
The simulator validates the positive electron speed at every scattering and
finite, positive-energy escaping photons with unit Stokes norm.

Before the production run, 1,000 photons with seed 110701 were compared with
the tightened settings: both runs gave 996 escapes, four absorptions and 535
scattered escaping photons, with identical scattering counts per escaping photon.
Mean Q, U and V changed by less than 0.0014, below the sampling uncertainty.
The 95th percentiles of absolute changes in Q, U and V were 0.0028, 0.0034
and 0.0037, respectively. A few individual trajectories changed substantially;
this is an ensemble check, not a claim of uniform per-photon convergence.
A separate repeated 20-photon run matched the corresponding records exactly.

One million initial photons is a preview-quality reproduction. The paper's
120×64 bins and 10^4 seed photons per bin imply 76,800,000 emitted photons;
the two figures therefore have different count scales and Monte Carlo noise.
In particular, sparse tail bins can show spuriously high estimated linear
polarization, and individual polarization angles there are noisy. The full
resolution is retained so this limitation remains visible.

Expected ensemble features are southern enhancement of scattered high-energy
photons, predominantly E-mode unscattered light, scattered polarization of a
few tens of percent, and an angle predominantly dependent on magnetic
colatitude. The implementation retains the current solver's full geodesic
propagation and Magnus Stokes integration, rather than reimplementing the
historical code's approximations.

## Completed production run (2026-10-24)

All four batches completed successfully. The slowest batch took 1,769.7 s
(about 29.5 minutes); no photons were dropped or replaced.

| Quantity | Number |
| --- | ---: |
| Emitted | 1,000,000 |
| Escaped, all energies | 992,706 |
| Absorbed | 7,294 |
| Escaped, 0.04–40 keV | 988,668 |
| Unscattered, 0.04–40 keV | 455,142 |
| Scattered, 0.04–40 keV | 533,526 |

The largest absolute single-photon Stokes norm error was 1.09×10^-12.
The exported 23,040-bin table passed both count partition and Stokes-sum
partition checks, and all finite linear fractions were within 0–100% to
roundoff. The PNG and one-page PDF were generated successfully and the PNG
was visually checked against the paper.

Of 73,542 scattered photons in 4–40 keV, 71.72% escaped southward. In the
2–10 keV band, combining photons over cos(theta) intervals (-0.8,-0.4),
(-0.2,0.2) and (0.4,0.8) gave scattered linear fractions of 29.19%, 24.72%
and 21.22%, respectively; unscattered fractions were 98.74%, 99.99% and
98.76%. Near the equator, the unscattered polarization angle was 118.50
degrees. These are broad-bin diagnostics, not replacements for the figure's
original 120×64 grid. The high-energy tail is visibly noisier than the
paper, as expected for the smaller sample.
