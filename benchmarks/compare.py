#!/usr/bin/env python3
"""Collect serial photon batches; compare distributions (NumPy/SciPy)."""
import argparse
import json
import os
from pathlib import Path
import subprocess
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build/perf"
RESULTS = ROOT / "benchmarks/results"
FIELDS = ["seed", "escaped", "n_scattered", "energy_inf", "muk", "Q", "U", "V"]


def collect(stage, name, groups, count, seed, accuracy="fast"):
    blocks = []
    for i in range(groups):
        path = BUILD / f"{name}-{i}.csv"
        with path.open("w") as output:
            command = [str(BUILD / stage / "main"), "--count", str(count), "--seed", str(seed + i * count), "--output", "csv"]
            if accuracy != "fast": command += ["--accuracy", accuracy]
            subprocess.run(command, cwd=ROOT, stdout=output, check=True)
        rows = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
        if rows.shape != (count, 8): raise RuntimeError("Missing photon records")
        if not np.array_equal(rows[:,0], np.arange(seed+i*count,seed+(i+1)*count)):
            raise RuntimeError("Missing or duplicated photon seeds")
        if not np.isfinite(rows).all() or not np.isin(rows[:,1],[0,1]).all():
            raise RuntimeError("Invalid photon record")
        if not ((rows[:,2] >= 0) & (rows[:,2] == np.floor(rows[:,2]))).all():
            raise RuntimeError("Invalid scattering count")
        blocks.append(rows)
        path.unlink()
        print(f"{name}: {i + 1}/{groups}", flush=True)
    np.savez_compressed(BUILD / f"{name}.npz", photons=np.concatenate(blocks), stage=stage, seed=seed, groups=groups, group_size=count, accuracy=accuracy)


def compare(a_name, b_name, name):
    af, bf = np.load(BUILD / f"{a_name}.npz"), np.load(BUILD / f"{b_name}.npz")
    a, b = af["photons"], bf["photons"]
    ae, be = a[a[:, 1] == 1], b[b[:, 1] == 1]
    checks, tests = {}, []
    def test(label, p, **values):
        tests.append({"test": label, "p": float(p), **values})
    def mean(label, x, y, tolerance=None):
        delta = float(y.mean() - x.mean())
        se = float(np.sqrt(x.var(ddof=1)/len(x) + y.var(ddof=1)/len(y)))
        p = 2 * stats.norm.sf(abs(delta)/se) if se else float(delta == 0)
        test(label, p, baseline=float(x.mean()), candidate=float(y.mean()), difference=delta, standard_error=se)
        if tolerance is not None: checks[label] = bool(abs(delta) <= tolerance)
    mean("escape_fraction", a[:, 1], b[:, 1], .01)
    mean("scatter_mean", a[:, 2], b[:, 2])
    for i, field in enumerate(FIELDS[3:], 3):
        mean(field + "_mean", ae[:, i], be[:, i], .01 * abs(ae[:, i].mean()) if i == 3 else .01)
        ks = stats.ks_2samp(ae[:, i], be[:, i], method="asymp")
        test(field + "_KS", ks.pvalue, distance=float(ks.statistic))
    # Merge sparse adjacent cells deterministically using pooled counts.
    def histogram(label, x, y):
        merged, cx, cy = [], 0., 0.
        for u, v in zip(x.ravel(), y.ravel()):
            cx += u; cy += v
            if cx + cy >= 40:
                merged.append([cx, cy]); cx = cy = 0.
        if cx + cy:
            if merged: merged[-1][0] += cx; merged[-1][1] += cy
            else: merged.append([cx, cy])
        if len(merged) > 1:
            result = stats.chi2_contingency(np.array(merged).T, correction=False)
            test(label, result.pvalue, statistic=float(result.statistic), cells=len(merged))
    sedges = np.r_[np.arange(-.5, 10), np.inf]
    histogram("scatter_histogram", np.histogram(a[:, 2], sedges)[0], np.histogram(b[:, 2], sedges)[0])
    edges = np.unique(np.r_[-np.inf, np.quantile(ae[:, 3], np.arange(.1, 1, .1)), np.inf])
    directions = np.linspace(-1, 1, 11)
    histogram("energy_direction_histogram", np.histogram2d(ae[:, 3], ae[:, 4], [edges, directions])[0], np.histogram2d(be[:, 3], be[:, 4], [edges, directions])[0])
    # Conditional polarization: five energy bins by five direction bins.
    ee = np.unique(np.r_[-np.inf, np.quantile(ae[:, 3], [.2,.4,.6,.8]), np.inf])
    mm = np.linspace(-1, 1, 6)
    for j in range(len(ee) - 1):
        for k in range(5):
            aa = ae[(ae[:,3] >= ee[j]) & (ae[:,3] < ee[j+1]) & (ae[:,4] >= mm[k]) & (ae[:,4] <= mm[k+1])]
            bb = be[(be[:,3] >= ee[j]) & (be[:,3] < ee[j+1]) & (be[:,4] >= mm[k]) & (be[:,4] <= mm[k+1])]
            if min(len(aa), len(bb)) >= 40:
                for i in (5, 6, 7): mean(f"conditional_{FIELDS[i]}_{j}_{k}", aa[:,i], bb[:,i])
    # Holm family-wise error control, alpha=.01 over all inferential tests.
    ordered = sorted(tests, key=lambda t: t["p"])
    largest = 0.
    for rank, t in enumerate(ordered):
        largest = max(largest, min(1., t["p"] * (len(tests)-rank)))
        t["holm_p"] = largest
        t["pass"] = largest >= .01
    for label, data, escaped in (("baseline",a,ae), ("candidate",b,be)):
        checks[label + "_finite"] = bool(np.isfinite(data).all())
        checks[label + "_physical"] = bool((escaped[:,3] > 0).all() and (np.abs(escaped[:,4]) <= 1).all() and (np.abs(escaped[:,5:]) <= 1+1e-8).all())
        checks[label + "_stokes_norm"] = bool((np.abs(np.linalg.norm(escaped[:,5:], axis=1)-1) < 1e-7).all())
    metadata = lambda data: {key: data[key].item() for key in ("stage","seed","groups","group_size","accuracy") if key in data.files}
    result = {"baseline": a_name, "candidate": b_name, "datasets": [metadata(af),metadata(bf)], "counts": [len(a),len(b)], "checks": checks, "tests": tests, "passed": all(checks.values()) and all(t["pass"] for t in tests), "energy_quantiles": {"probabilities": [.01,.1,.5,.9,.99,.999], "baseline": np.quantile(ae[:,3],[.01,.1,.5,.9,.99,.999]).tolist(), "candidate": np.quantile(be[:,3],[.01,.1,.5,.9,.99,.999]).tolist()}}
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"{name}.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"passed": result["passed"], "failed_checks": [k for k,v in checks.items() if not v], "failed_tests": [t for t in tests if not t["pass"]]}, indent=2), flush=True)
    return result["passed"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    c = sub.add_parser("collect")
    c.add_argument("stage"); c.add_argument("name")
    c.add_argument("--groups", type=int, default=30)
    c.add_argument("--count", type=int, default=10000)
    c.add_argument("--seed", type=int, default=1000042)
    c.add_argument("--accuracy", choices=["fast","strict"], default="fast")
    c = sub.add_parser("compare")
    c.add_argument("baseline"); c.add_argument("candidate"); c.add_argument("--name", required=True)
    a = p.parse_args()
    if a.command == "collect": collect(a.stage, a.name, a.groups, a.count, a.seed, a.accuracy)
    else: raise SystemExit(0 if compare(a.baseline, a.candidate, a.name) else 1)
