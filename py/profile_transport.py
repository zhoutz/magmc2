"""Count hot-path calls without including instrumentation in speed measurements."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

from benchmark_transport import ROOT, build


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--first', type=int, default=42)
    parser.add_argument('--count', type=int, default=2000)
    args = parser.parse_args()
    directory = ROOT/'build/single_thread_optimization'
    baseline = directory/'baseline_src'
    instrumented = directory/'instrumented_baseline'
    shutil.copytree(baseline, instrumented, dirs_exist_ok=True)
    shutil.copy(ROOT/'src/transport_stats.hpp', instrumented)
    additions = {
        'evolution.hpp': [('  basis(Geodesic::State y) const {', 'basis')],
        'magnus.hpp': [('  do_step(DerivFunc const& derivs) {', 'magnus'),
                       ('  interp(double x) const {', 'polarization')],
        'quad.hpp': [('  qags(F& fn, double a, double b, double epsabs, double epsrel) {', 'quadrature')],
        'resonant.hpp': [('  integration_cuts(double left, double right, BasisAt const& get_basis) const {', 'cuts'),
                         ('    auto integrand = [&](double l) {', 'integrand')],
    }
    for name, changes in additions.items():
        path = instrumented/name
        text = path.read_text()
        text = '#include "transport_stats.hpp"\n' + text
        for marker, counter in changes:
            if text.count(marker) != 1:
                raise ValueError(f'Unrecognized baseline layout in {name}: {marker}')
            text = text.replace(marker, marker+f'\n    TRANSPORT_COUNT({counter});')
        if name == 'resonant.hpp':
            marker = '          cuts.push_back(ridders_root(f, a, b, fa, fb));'
            if text.count(marker) != 1:
                raise ValueError('Unrecognized baseline root finder')
            text = text.replace(marker, '          { TRANSPORT_COUNT(roots); cuts.push_back(ridders_root(f, a, b, fa, fb)); }')
        path.write_text(text)
    counts = {}
    outputs = []
    for name, sources in [('reference', instrumented), ('candidate', ROOT/'src')]:
        binary = directory/f'profile-{name}'
        defines = ['-DTRANSPORT_STATS']
        if name == 'candidate':
            defines.append('-DTRANSPORT_CURRENT')
        build(ROOT/'tests/transport_driver.cpp', sources, binary, defines)
        output, log = directory/f'profile-{name}.txt', directory/f'profile-{name}.log'
        with output.open('w') as stream, log.open('w') as diagnostics:
            subprocess.run([str(binary), str(args.first), str(args.count)],
                           stdout=stream, stderr=diagnostics, cwd=ROOT, check=True)
        counts[name] = {k: int(v) for k, v in
                        (entry.split('=') for entry in log.read_text().splitlines()[-1].split())}
        outputs.append(output.read_bytes())
    report = {'first': args.first, 'count': args.count, 'counts': counts,
              'output_byte_identical': outputs[0] == outputs[1],
              'fallback': 'Active or uncertain intervals use the original QAGS and boundary scan'}
    (directory/'profile-report.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if not report['output_byte_identical']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
