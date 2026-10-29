#!/usr/bin/env python3
"""Serial, reproducible stage builds and end-to-end benchmarks."""
import argparse
import difflib
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = "70b7c2218d9520011f0bf02314d353a97837204e"
BUILD = ROOT / "build/perf"
STAGES = ROOT / "benchmarks/stages"
RESULTS = ROOT / "benchmarks/results"


def original():
    paths = subprocess.check_output(["git", "ls-tree", "-r", "--name-only", BASE, "src"], cwd=ROOT, text=True).splitlines()
    return {p: subprocess.check_output(["git", "show", f"{BASE}:{p}"], cwd=ROOT, text=True) for p in paths}


def save(stage):
    STAGES.mkdir(parents=True, exist_ok=True)
    before = original()
    after = {str(p.relative_to(ROOT)): p.read_text() for p in (ROOT / "src").glob("*") if p.is_file()}
    patch = "".join("".join(difflib.unified_diff(before.get(p, "").splitlines(True), after.get(p, "").splitlines(True), fromfile=f"a/{p}" if p in before else "/dev/null", tofile=f"b/{p}")) for p in sorted(before.keys() | after.keys()))
    (STAGES / f"{stage}.patch").write_text(patch)


def build(stage, diagnostic=False):
    target = BUILD / stage
    if target.exists(): shutil.rmtree(target)
    target.mkdir(parents=True)
    for p, data in original().items():
        out = target / p
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(data)
    if stage != "original":
        with (STAGES / f"{stage}.patch").open() as patch:
            subprocess.run(["patch", "-s", "-p1"], cwd=target, stdin=patch, check=True)
    gsl = subprocess.check_output(["brew", "--prefix", "gsl"], text=True).strip()
    flags = ["-std=c++23", "-O3", "-march=native", "-Wall", "-Wextra", "-Wpedantic"]
    if diagnostic: flags.append("-DMAGMC_DIAGNOSTICS")
    command = [os.environ.get("CXX", "g++-16"), *flags, str(target / "src/main.cpp"), "-o", str(target / "main"), f"-I{gsl}/include", f"-L{gsl}/lib", "-lgsl", "-lgslcblas"]
    subprocess.run(command, cwd=ROOT, check=True)
    (target / "build.json").write_text(json.dumps({"base": BASE, "command": command}, indent=2))


def measure(stages, repeats, count, seed, output, name):
    records = {s: [] for s in stages}
    def run(stage):
        if stage == "original" and (count != 10000 or seed != 42 or output != "text"):
            raise ValueError("Unmodified original only supports the default text workload")
        start = time.perf_counter()
        p = subprocess.run([str(BUILD / stage / "main"), "--count", str(count), "--seed", str(seed), "--output", output, "--summary"], cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, check=True)
        return {"wall_seconds": time.perf_counter() - start, "summary": json.loads(p.stderr) if p.stderr.strip() else None}
    for s in stages: run(s)
    for i in range(repeats):
        order = stages if i % 2 == 0 else stages[::-1]
        for s in order:
            r = run(s)
            records[s].append(r)
            print(s, i + 1, r["wall_seconds"], flush=True)
    medians = {s: statistics.median(r["wall_seconds"] for r in rows) for s, rows in records.items()}
    result = {"base": BASE, "count": count, "seed": seed, "output": output, "repeats": repeats, "records": records, "medians": medians}
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"{name}.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(medians, indent=2), flush=True)


if __name__ == "__main__":
    os.chdir(ROOT)
    p = argparse.ArgumentParser(__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for command in ("save", "build"):
        s = sub.add_parser(command)
        s.add_argument("stage")
        if command == "build": s.add_argument("--diagnostic", action="store_true")
    s = sub.add_parser("measure")
    s.add_argument("stages", nargs="+")
    s.add_argument("--repeats", type=int, default=5)
    s.add_argument("--count", type=int, default=10000)
    s.add_argument("--seed", type=int, default=42)
    s.add_argument("--output", choices=["text", "none"], default="text")
    s.add_argument("--name", required=True)
    a = p.parse_args()
    if a.command == "save": save(a.stage)
    elif a.command == "build": build(a.stage, a.diagnostic)
    else: measure(a.stages, a.repeats, a.count, a.seed, a.output, a.name)
