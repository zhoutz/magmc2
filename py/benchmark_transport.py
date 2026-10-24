"""Compare two FD11 executables serially, including exact photon records."""

import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("optimized", type=Path)
    parser.add_argument("--photons", type=int, default=1000)
    parser.add_argument("--seeds", type=int, nargs="+", default=[110701, 110702])
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--tight", action="store_true")
    parser.add_argument("--fast", action="store_true",
                        help="Use fast candidate and matched per-photon streams; reference must support them")
    parser.add_argument("--threads", type=int, default=1, help="Candidate workers in --fast mode")
    args = parser.parse_args()
    if args.photons <= 0 or args.repeat <= 0 or any(s <= 0 for s in args.seeds):
        parser.error("Photons, repeats and seeds must be positive")
    if args.threads <= 0 or (args.threads != 1 and not args.fast):
        parser.error("Positive --threads requires --fast")
    Path("build").mkdir(exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="transport-benchmark-", dir="build"))
    print(f"Results: {directory}", flush=True)
    results = []
    for repeat in range(args.repeat):
        for seed in args.seeds:
            row = dict(repeat=repeat, seed=seed, photons=args.photons, tight=args.tight)
            outputs = []
            for label in ("baseline", "optimized"):
                executable = getattr(args, label).resolve()
                output = directory / f"{repeat}-{seed}-{label}.txt"
                command = [str(executable), "--photons", str(args.photons),
                           "--seed", str(seed), "--output", str(output)]
                if args.tight:
                    command.append("--tight")
                if args.fast:
                    command.append("--photon-streams")
                    if label == "optimized":
                        command.extend(["--fast", "--threads", str(args.threads)])
                start = time.perf_counter()
                run = subprocess.run(command, text=True, capture_output=True)
                row[f"{label}_seconds"] = time.perf_counter() - start
                (directory / f"{repeat}-{seed}-{label}.log").write_text(run.stderr)
                if run.returncode:
                    raise RuntimeError(f"{label} failed: {run.stderr}")
                outputs.append(output)
            row["byte_identical"] = outputs[0].read_bytes() == outputs[1].read_bytes()
            row["speedup"] = row["baseline_seconds"] / row["optimized_seconds"]
            if args.fast:
                from compare_transport import compare
                accuracy = compare(*outputs)
                (directory / f"{repeat}-{seed}-accuracy.json").write_text(
                    json.dumps(accuracy, indent=2) + "\n")
                row["candidate_threads"] = args.threads
            results.append(row)
            (directory / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
            print(json.dumps(row), flush=True)
            if not args.fast and not row["byte_identical"]:
                raise RuntimeError("Output changed; inspect photon records before accepting")


if __name__ == "__main__":
    main()
