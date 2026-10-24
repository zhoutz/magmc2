"""Compare matched Monte Carlo runs; report numerical changes relative to MC noise."""

import argparse
import json
from pathlib import Path

import numpy as np

from plot_fd11_fig7 import read_photons


def compare(reference, candidate):
    # Reuse accounting, completion and physical-output validation.
    _, am = read_photons([reference])
    _, bm = read_photons([candidate])
    for key in ("requested", "seed", "model", "random_stream", "R_km", "M_solar",
                "B_pole_G", "T_infinity_keV", "field_p", "field_A", "beta_min", "beta_max"):
        if am[0].get(key) != bm[0].get(key):
            raise ValueError(f"Cannot compare different {key}")
    n = int(am[0]["emitted"])
    a = np.loadtxt(reference, ndmin=2)
    b = np.loadtxt(candidate, ndmin=2)
    _, ia, ib = np.intersect1d(a[:, 0], b[:, 0], return_indices=True)
    delta = np.abs(b[ib, 1:6] - a[ia, 1:6])
    report = {
        "emitted": n, "reference_escaped": len(a), "candidate_escaped": len(b),
        "escape_classification_changes": len(a) + len(b) - 2 * len(ia),
        "matched_scattering_count_changes": int(np.count_nonzero(a[ia, 6] != b[ib, 6])),
        "delta_columns": ["energy_keV", "mu", "Q", "U", "V"],
        "absolute_delta_quantiles": {
            str(q): np.quantile(delta, q, axis=0).tolist() for q in (0.5, 0.95, 0.99, 1.0)
        },
        "global": {}, "binned": {},
    }
    # Each emitted photon contributes zero if absorbed. This also accounts for
    # changes in escape classification, rather than conditioning on survivors.
    for name, column in (("I", None), ("energy", 1), ("mu", 2), ("Q", 3), ("U", 4), ("V", 5)):
        wa = np.ones(len(a)) if column is None else a[:, column]
        wb = np.ones(len(b)) if column is None else b[:, column]
        change = (wb.sum() - wa.sum()) / n
        sem = np.sqrt(max(0, (wa * wa).sum() - wa.sum() ** 2 / n)) / n
        report["global"][name] = dict(mean_change=float(change), reference_sem=float(sem),
                                      sigma=float(abs(change) / sem) if sem else None)

    energy = [0, 0.04, 0.4, 1, 2, 4, 10, 40, np.inf]
    mu = np.linspace(-1, 1, 5)
    for group in ("total", "unscattered", "scattered"):
        def select(rows):
            return rows if group == "total" else rows[(rows[:, 6] == 0) if group == "unscattered" else (rows[:, 6] > 0)]
        x, y = select(a), select(b)
        def histogram(rows, weights):
            return np.histogram2d(rows[:, 1], rows[:, 2], bins=(energy, mu), weights=weights)[0]
        counts = histogram(x, np.ones(len(x)))
        for name, column in (("I", None), ("Q", 3), ("U", 4), ("V", 5)):
            wx = np.ones(len(x)) if column is None else x[:, column]
            wy = np.ones(len(y)) if column is None else y[:, column]
            sx, sy = histogram(x, wx), histogram(y, wy)
            variance = histogram(x, wx * wx) - sx * sx / n
            valid = (counts >= 100) & (variance > 1e-20)
            z = np.abs(sy[valid] - sx[valid]) / np.sqrt(variance[valid])
            report["binned"][f"{group}_{name}"] = {
                "bins_with_at_least_100_photons": int(valid.sum()),
                "max_sigma": float(z.max()) if len(z) else None,
            }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = compare(args.reference, args.candidate)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
