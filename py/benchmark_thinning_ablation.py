#!/usr/bin/env python3
"""Build isolated source copies; serial leave-one-out timing, no src/ edits."""
import argparse, hashlib, json, random, resource, shutil, statistics, subprocess, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'build/thinning/ablation'
FLAGS = ['-std=c++23', '-O3', '-march=native', '-Wall', '-Wextra', '-Wpedantic', '-I/opt/homebrew/include', '-L/opt/homebrew/lib', '-lstdc++exp', '-lgsl', '-lgslcblas']

def replace(src, filename, old, new):
    p=src/filename; text=p.read_text()
    if text.count(old)!=1: raise RuntimeError(f'{filename}: expected one replacement, got {text.count(old)}')
    p.write_text(text.replace(old,new))

def disable(src, name):
    if name=='bounds': replace(src,'evolution.hpp','auto cuts = possible_resonance()','auto cuts = true')
    elif name=='cache': replace(src,'evolution.hpp','double cache_tolerance = 1e-7;','double cache_tolerance = 0;')
    elif name=='dense': replace(src,'evolution.hpp','stepper.do_step(rotation, false)','stepper.do_step(rotation, true)')
    elif name=='initial_step': replace(src,'evolution.hpp','stepper.init(cursor, end, end - cursor, pol.stokes);','stepper.init(cursor, end, 0.1 * (end - cursor), pol.stokes);')
    elif name=='norm':
        original=(ROOT/'tests/reference/magnus.hpp').read_text()
        current=(src/'magnus.hpp').read_text()
        a=current.index('    double angle = a.length();'); b=current.index('    if (!std::isfinite(angle))',a)
        current=current[:a]+'    double angle = std::hypot(a.x, a.y, a.z);'+chr(10)+current[b:]
        a=current.index('    double value = s.length();'); b=current.index(chr(10)+'  }',a)
        current=current[:a]+'    return std::hypot(s.x, s.y, s.z);'+current[b:]
        (src/'magnus.hpp').write_text(current)
    elif name=='quadrature':
        p=src/'resonant.hpp'; p.write_text(p.read_text().replace('quad.smooth(', 'quad.qags('))
    elif name=='cuts': replace(src,'evolution.hpp','resonant.integration_cuts(left, right, get_basis, true)','resonant.integration_cuts(left, right, get_basis, false)')
    elif name=='locator':
        p=src/'resonant.hpp'; text=p.read_text()
        a=text.index('        // Retain cumulative mass'); b=text.index(chr(10)+'      }',a)
        code='''        auto residual = [&](double l) { return target + integral(l); };
        double event = ridders_root(residual, a, b, target, target + mass, 1e-5 * (b-a));
        return std::clamp(event, std::nextafter(a,b), std::nextafter(b,a));'''
        p.write_text(text[:a]+code+text[b:])
    else: raise ValueError(name)

NAMES=['bounds','cache','dense','initial_step','norm','quadrature','cuts','locator']
def build():
    OUT.mkdir(parents=True,exist_ok=True)
    for name in ['baseline','full','core']+NAMES:
        src=OUT/name; shutil.copytree(ROOT/('tests/reference' if name=='baseline' else 'src'),src,dirs_exist_ok=True)
        if name=='baseline':
            replace(src,'main.cpp','main() {','main(int argc, char** argv) {')
            replace(src,'main.cpp','int N = 1000;','int N = std::stoi(argv[1]); int seed = std::stoi(argv[2]); (void)argc;')
            replace(src,'main.cpp','42 + i','seed + i')
        elif name=='core':
            for feature in NAMES: disable(src,feature)
        elif name!='full': disable(src,name)
        subprocess.run(['g++-16',str(src/'main.cpp'),'-o',str(src/'main')]+FLAGS,cwd=ROOT,check=True)
        print('built',name,flush=True)
    (OUT/'manifest.json').write_text(json.dumps({'flags':FLAGS,'compiler':subprocess.check_output(['g++-16','--version'],text=True).splitlines()[0], 'hashes':{str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*/*.hpp')},'input_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in list((ROOT/'src').glob('*'))+[ROOT/'table/bfield_t10.txt'] }},indent=2))

def run(name,n,seed):
    cmd=[str(OUT/name/'main')]+([str(n),str(seed)] if name=='baseline' else ['--transport=thinning',f'--photons={n}',f'--seed={seed}'])
    before=resource.getrusage(resource.RUSAGE_CHILDREN); start=time.perf_counter()
    with (OUT/f'{name}-{n}-{seed}.txt').open('w') as out:
        subprocess.run(cmd,stdout=out,cwd=ROOT,check=True)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(wall=time.perf_counter()-start,user=after.ru_utime-before.ru_utime,system=after.ru_stime-before.ru_stime)

def bench():
    report=[]
    for n,seed in [(1000,42),(3000,1000000)]:
        names=['baseline','full','core']+NAMES
        for name in names: run(name,100,seed)
        samples={name:[] for name in names}; rng=random.Random(4729)
        for rep in range(3):
            order=names.copy();rng.shuffle(order)
            for name in order:
                row=run(name,n,seed);samples[name].append(row);print(n,seed,rep,name,row,flush=True)
        med={name:{k:statistics.median(r[k] for r in rr) for k in ['wall','user','system']} for name,rr in samples.items()}
        report.append(dict(n=n,seed=seed,runs=samples,medians=med, relative_to_full={name:med[name]['user']/med['full']['user'] for name in names}))
        (OUT/'report.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['build','run','all']);args=p.parse_args()
    if args.stage in ('build','all'):build()
    if args.stage in ('run','all'):bench()
