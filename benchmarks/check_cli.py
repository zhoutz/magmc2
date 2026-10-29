#!/usr/bin/env python3
"""CLI regression checks; optionally verify the untouched reference harness."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(__doc__)
p.add_argument("binary", nargs="?", default="build/main")
p.add_argument("--reference-check", action="store_true")
a = p.parse_args()
binary = str((ROOT / a.binary).resolve())

for args in (["--count","0"], ["--count","-1"], ["--count","1x"],
             ["--count","2147483647","--seed","2"], ["--seed"],
             ["--output","bad"], ["--accuracy","bad"], ["--unknown","1"]):
    result = subprocess.run([binary,*args],cwd=ROOT,capture_output=True,text=True)
    assert result.returncode != 0 and "Simulation failed:" in result.stderr, args
result = subprocess.run([binary,"--count","8","--output","csv","--summary"],cwd=ROOT,capture_output=True,text=True,check=True)
rows = list(csv.DictReader(io.StringIO(result.stdout)))
summary = json.loads(result.stderr)
assert len(rows) == 8 and [int(r["seed"]) for r in rows] == list(range(42,50))
assert summary["count"] == summary["escaped"] + summary["absorbed"] == 8
result = subprocess.run([binary,"--count","4","--seed","99","--accuracy","strict","--output","none","--summary"],cwd=ROOT,capture_output=True,text=True,check=True)
assert result.stdout == "" and json.loads(result.stderr)["count"] == 4
result = subprocess.run([binary,"--count","4"],cwd=ROOT,capture_output=True,text=True,check=True)
assert len(result.stdout.splitlines()) == 4
assert all(line == "Absorbed" or line.startswith("Escaped: energy_inf=") for line in result.stdout.splitlines())
result = subprocess.run([binary,"--help"],cwd=ROOT,capture_output=True,text=True,check=True)
assert "--accuracy" in result.stdout
checks = {"cli": "passed"}
if a.reference_check:
    hashes = {}
    for stage in ("original","stage00"):
        output = subprocess.check_output([str(ROOT / "build/perf" / stage / "main")],cwd=ROOT)
        hashes[stage] = hashlib.sha256(output).hexdigest()
    assert len(set(hashes.values())) == 1, "Reference harness changed baseline output"
    checks["reference_text_sha256"] = hashes
    (ROOT / "benchmarks/results/cli-checks.json").write_text(json.dumps(checks,indent=2)+"\n")
print(json.dumps(checks,indent=2))
