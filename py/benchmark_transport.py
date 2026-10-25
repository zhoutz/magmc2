"""Build and compare the single-thread transport against a frozen source tree."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
FLAGS = ['-std=c++23', '-O3', '-march=native', '-Wall', '-Wextra',
         '-Wpedantic', '-I/opt/homebrew/include', '-L/opt/homebrew/lib',
         '-lstdc++exp', '-lgsl', '-lgslcblas']


def build(source, includes, output, defines=()):
    subprocess.run(['g++-16', str(source), '-I' + str(includes), '-o',
                    str(output)] + list(defines) + FLAGS, cwd=ROOT, check=True)


def compare(reference, candidate):
    a, b = reference.splitlines(), candidate.splitlines()
    if not a or len(a) != len(b):
        raise ValueError('Empty output or different photon counts')
    errors = [[] for _ in range(5)]
    seeds = []
    classes, scatters = [], []
    for x, y in zip(a, b):
        x, y = x.split(), y.split()
        if x[0] != y[0]:
            raise ValueError('Different seeds')
        seed = int(x[0])
        if x[2] not in ('E', 'A') or y[2] not in ('E', 'A'):
            raise ValueError('Invalid photon classification')
        if len(x) != (8 if x[2] == 'E' else 3) or len(y) != (8 if y[2] == 'E' else 3):
            raise ValueError('Invalid photon record')
        if x[1] != y[1]:
            scatters.append(seed)
        if x[2] != y[2]:
            classes.append(seed)
        elif x[2] == 'E':
            seeds.append(seed)
            for j, (u, v) in enumerate(zip(x[3:], y[3:])):
                u, v = float(u), float(v)
                if not math.isfinite(u) or not math.isfinite(v) or (j == 0 and u <= 0):
                    err = math.inf
                else:
                    err = abs(u-v) / abs(u) if j == 0 else abs(u-v)
                errors[j].append(err if math.isfinite(err) else math.inf)
    fields = {}
    failed = set(classes)
    for name, values in zip(['energy_relative', 'muk', 'Q', 'U', 'V'], errors):
        if not values:
            continue
        ordered = sorted(values)
        fields[name] = {'max': max(values), 'worst_seed': seeds[values.index(max(values))],
                        'quantiles': {str(q): ordered[min(len(ordered)-1, int(q*(len(ordered)-1)))]
                                      for q in [0.5, 0.95, 0.99]}}
        failed.update(seed for seed, err in zip(seeds, values) if err >= 1e-3)
    return {'photons': len(a), 'classification_changes': classes,
            'scattering_count_changes': scatters, 'fields': fields,
            'failed_seeds': sorted(failed), 'passed': not failed}


def run(binary, first, count, path, production=False):
    start = time.perf_counter()
    with path.open('w') as stream:
        command = [str(binary)] if production else [str(binary), str(first), str(count)]
        subprocess.run(command, stdout=stream,
                       cwd=ROOT, check=True)
    return time.perf_counter() - start


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--baseline', type=Path, default=ROOT/'build/single_thread_optimization/baseline_src')
    parser.add_argument('--output', type=Path, default=ROOT/'build/single_thread_optimization')
    parser.add_argument('--first', type=int, default=42)
    parser.add_argument('--count', type=int, default=2000)
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--no-build', action='store_true')
    parser.add_argument('--skip-reference', action='store_true')
    parser.add_argument('--production', action='store_true', help='Time the unmodified main.cpp entrypoint; accuracy uses the driver separately')
    parser.add_argument('--warmup', action='store_true')
    parser.add_argument('--baseline-revision', help='Create a missing baseline snapshot from this git revision')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.repeats < 1 or args.count < 1 or args.first < 0:
        parser.error('repeats/count must be positive; first must be nonnegative')
    if args.production and (args.first != 42 or args.count != 2000):
        parser.error('Production entrypoint uses the fixed 42/2000 workload')
    if args.baseline_revision:
        if args.baseline.exists():
            parser.error('Refusing to overwrite the existing baseline snapshot')
        args.baseline.mkdir(parents=True)
        paths = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only',
                                         args.baseline_revision, 'src'], cwd=ROOT, text=True).splitlines()
        for path in paths:
            data = subprocess.check_output(['git', 'show', f'{args.baseline_revision}:{path}'], cwd=ROOT)
            if path == 'src/main.cpp':
                data = data.replace(b'int N = 1000;', b'int N = 2000;')
            (args.baseline/Path(path).name).write_bytes(data)
    suffix = '-main' if args.production else ''
    old, new = args.output/('reference'+suffix), args.output/('candidate'+suffix)
    if not args.no_build:
        if not args.skip_reference:
            build(args.baseline/'main.cpp' if args.production else ROOT/'tests/transport_driver.cpp', args.baseline, old)
        build(ROOT/'src/main.cpp' if args.production else ROOT/'tests/transport_driver.cpp', ROOT/'src', new,
              () if args.production else ['-DTRANSPORT_CURRENT'])
    ref, cand = (args.output/f'{args.first}-{args.count}-{kind}{suffix}.txt' for kind in ['reference', 'candidate'])
    if args.warmup:
        for binary in ([new] if args.skip_reference else [old, new]):
            run(binary, args.first, args.count, args.output/'warmup.txt', args.production)
    timings = []
    digests = {}
    for i in range(args.repeats):
        t = {}
        kinds = [('reference', old, ref), ('candidate', new, cand)]
        if i % 2:
            kinds.reverse()
        for kind, binary, path in kinds:
            if kind == 'reference' and args.skip_reference:
                continue
            t[kind] = run(binary, args.first, args.count, path, args.production)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if kind in digests and digests[kind] != digest:
                raise RuntimeError(f'{kind} produced different output on repeated runs')
            digests[kind] = digest
        timings.append(t)
    if args.production:
        a, b = ref.read_bytes(), cand.read_bytes()
        result = {'byte_identical': a == b,
                  'passed': a == b and len(a.splitlines()) == args.count,
                  'rows': len(a.splitlines())}
    else:
        result = compare(ref.read_text(), cand.read_text())
        if result['photons'] != args.count:
            raise ValueError('Output photon count does not match requested count')
    result['timings'] = timings
    result['output_sha256'] = digests
    result['compiler'] = subprocess.check_output(['g++-16', '--version'], text=True).splitlines()[0]
    result['flags'] = FLAGS
    result['source_sha256'] = {
        kind: {str(p.relative_to(tree)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted(tree.glob('*')) if p.is_file()}
        for kind, tree in [('reference', args.baseline), ('candidate', ROOT/'src')]}
    result['table_sha256'] = hashlib.sha256((ROOT/'table/bfield_t10.txt').read_bytes()).hexdigest()
    if not args.skip_reference:
        result['median_seconds'] = {kind: statistics.median(t[kind] for t in timings)
                                    for kind in ['reference', 'candidate']}
        result['speedup'] = result['median_seconds']['reference'] / result['median_seconds']['candidate']
        result['performance_passed'] = result['speedup'] >= 10
        result['target_met'] = result['passed'] and result['performance_passed']
    print(json.dumps(result, indent=2))
    (args.output/f'{args.first}-{args.count}-report{suffix}.json').write_text(json.dumps(result, indent=2))
    if not result['passed'] or (args.production and not result.get('performance_passed', False)):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
