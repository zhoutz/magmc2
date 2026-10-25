#!/usr/bin/env python3
'''Serial benchmark and paired, simultaneous statistical intervals. No dependencies.'''
import argparse, bisect, hashlib, json, math, random, resource, shutil, statistics, subprocess, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'build/thinning/validation'
FLAGS = ['-std=c++23', '-O3', '-march=native', '-Wall', '-Wextra', '-Wpedantic', '-I/opt/homebrew/include', '-L/opt/homebrew/lib', '-lstdc++exp', '-lgsl', '-lgslcblas']
def save(name, data):
    (OUT / (name + '.json')).write_text(json.dumps(data, indent=2))
def build():
    OUT.mkdir(parents=True, exist_ok=True)
    driver = (ROOT / 'tests/reference/main.cpp').read_text()
    driver = driver.replace('main() {', 'main(int argc, char** argv) {').replace('int N = 1000;', 'int N = argc > 1 ? std::stoi(argv[1]) : 1000; int seed = argc > 2 ? std::stoi(argv[2]) : 42;').replace('42 + i', 'seed + i')
    driver = driver.replace('photon.init();', 'photon.init(); int count = 0;').replace('photon.perform_scattering();', 'photon.perform_scattering(); ++count;')
    # Insert diagnostics immediately after the transport loop.
    idx = driver.rfind('  }')
    driver = driver[:idx] + '    if (argc > 3) std::println(stderr, ' + chr(34) + '{} {}' + chr(34) + ', seed + i, count);' + chr(10) + driver[idx:]
    (OUT / 'reference_driver.cpp').write_text(driver)
    robust = OUT / 'robust_reference'
    shutil.copytree(ROOT / 'tests/reference', robust, dirs_exist_ok=True)
    shutil.copyfile(ROOT / 'src/find_root.hpp', robust / 'find_root.hpp')
    for name, source, inc in [('reference', OUT / 'reference_driver.cpp', 'tests/reference'), ('robust', OUT / 'reference_driver.cpp', str(OUT/'robust_reference')), ('main', ROOT / 'src/main.cpp', 'src'), ('test', ROOT / 'tests/test_thinning.cpp', 'src')]:
        subprocess.run(['g++-16', str(source), '-I'+inc, '-o', str(OUT/name)] + FLAGS, cwd=ROOT, check=True)
    subprocess.run([str(OUT/'test')], cwd=ROOT, check=True)
    save('manifest', {'flags': FLAGS, 'hashes': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'src').glob('*'))}, 'table_hash': hashlib.sha256((ROOT/'table/bfield_t10.txt').read_bytes()).hexdigest()})
def run(name, mode, n, seed, extra=(), diagnostics=True):
    cmd = ([str(OUT/('reference' if mode in ('snapshot','robust') else 'robust')),str(n),str(seed)] + (['counts'] if diagnostics else []) if mode in ('snapshot','robust') else [str(OUT/'main'),'--transport='+mode,'--photons='+str(n),'--seed='+str(seed)] + (['--diagnostics'] if diagnostics else []) + list(extra))
    u = resource.getrusage(resource.RUSAGE_CHILDREN); t = time.perf_counter()
    with (OUT/(name+'.txt')).open('w') as out, (OUT/(name+'.log')).open('w') as log:
        subprocess.run(cmd,cwd=ROOT,stdout=out,stderr=log,check=True)
    v = resource.getrusage(resource.RUSAGE_CHILDREN)
    r = dict(wall=time.perf_counter()-t,user=v.ru_utime-u.ru_utime,system=v.ru_stime-u.ru_stime,n=n,seed=seed,mode=mode)
    print(name,r,flush=True); save(name+'-timing',r); return r
def load(name):
    rows=[]; logs=(OUT/(name+'.log')).read_text().splitlines()
    for i,l in enumerate((OUT/(name+'.txt')).read_text().splitlines()):
        count=int(logs[i].split()[1]) if logs else 0
        if l=='Absorbed': rows.append((0,count,0,0,0,0,0)); continue
        v=tuple(float(x.split('=')[1]) for x in l.split()[1:])
        assert len(v)==5 and all(map(math.isfinite,v)) and abs(v[1])<=1 and sum(x*x for x in v[2:])<=1+1e-8
        rows.append((1,count)+v)
    return rows
def quantile(x,p): return x[min(len(x)-1,int(p*len(x)))]
def bins_from(rows):
    e=sorted(r[2] for r in rows if r[0])
    return dict(energy=[quantile(e,j/10) for j in range(1,10)]+[quantile(e,.99)],joint=[quantile(e,j/4) for j in range(1,4)])
def compare(a,b,bins,name):
    assert len(a)==len(b) and len(a)>=1000
    n=len(a); defs=[('escape',lambda r:r[0],None),('scatter_mean',lambda r:r[1],None),('energy_mean',lambda r:r[2],lambda r:r[0])]
    for t in bins['energy']: defs.append((f'energy_cdf_{t:.8g}',lambda r,t=t:r[0]*(r[2]<=t),lambda r:r[0]))
    for t in [-.8,-.6,-.4,-.2,0,.2,.4,.6,.8]: defs.append((f'mu_cdf_{t}',lambda r,t=t:r[0]*(r[3]<=t),lambda r:r[0]))
    for t in range(5): defs.append((f'scatter_cdf_{t}',lambda r,t=t:r[1]<=t,None))
    defs.append(('tail',lambda r:r[0]*(r[2]>bins['energy'][-1]),lambda r:r[0]))
    for ie in range(4):
        for im in range(4):
            def sel(r,ie=ie,im=im): return bool(r[0] and bisect.bisect_right(bins['joint'],r[2])==ie and bisect.bisect_right([-.5,0,.5],r[3])==im)
            defs.append((f'count_{ie}_{im}',sel,None))
            for k,label in enumerate('QUV',4): defs.append((f'{label}_{ie}_{im}',lambda r,sel=sel,k=k:r[k]*sel(r),sel))
    stats=[]; blocks=[]
    for label,value,condition in defs:
        va=[value(r) for r in a]; vb=[value(r) for r in b]
        da=[condition(r) for r in a] if condition else [1]*n
        db=[condition(r) for r in b] if condition else [1]*n
        na,nb=sum(da),sum(db)
        if min(na,nb)<100: raise ValueError('Insufficient counts: '+label)
        ma,mb=sum(va)/na,sum(vb)/nb
        d=[(y-mb*dy)*n/nb-(x-ma*dx)*n/na for x,y,dx,dy in zip(va,vb,da,db)]
        se=math.sqrt(sum(x*x for x in d)/(n*(n-1)))
        blocks.append([sum(d[n*j//100:n*(j+1)//100])/n for j in range(100)])
        stats.append(dict(name=label,reference=ma,candidate=mb,difference=mb-ma,standard_error=se,counts=[na,nb]))
    rng=random.Random(48291); maxima=[]
    for _ in range(1999):
        weights=[0]*100
        for _ in range(100): weights[rng.randrange(100)]+=1
        maxima.append(max(abs(sum((w-1)*v for w,v in zip(weights,block)))/s['standard_error'] if s['standard_error'] else 0 for s,block in zip(stats,blocks)))
    critical=quantile(sorted(maxima),.99)
    for s in stats:
        width=critical*s['standard_error']; delta=s['difference']
        s['simultaneous_99pct_interval']=[delta-width,delta+width]; s['compatible']=abs(delta)<=width+1e-14
    result=dict(n=n,bootstrap_replicates=1999,paired_blocks=100,critical=critical,statistics=stats,passed=all(s['compatible'] for s in stats))
    save(name,result); print(name,'passed',result['passed'],flush=True); return result
def benchmark():
    report=[]
    for n,seed in [(1000,42),(10000,1000000),(10000,2000000),(10000,3000000)]:
        for mode in ['snapshot','thinning']: run(f'warm-{seed}-{mode}',mode,n,seed,diagnostics=False)
        rows={'snapshot':[],'thinning':[]}
        for i in range(5):
            for mode in (['snapshot','thinning'] if i%2==0 else ['thinning','snapshot']): rows[mode].append(run(f'bench-{seed}-{i}-{mode}',mode,n,seed,diagnostics=False))
        med={m:{k:statistics.median(r[k] for r in rr) for k in ['wall','user']} for m,rr in rows.items()}
        speed=med['snapshot']['wall']/med['thinning']['wall']
        report.append(dict(n=n,seed=seed,runs=rows,medians=med,speedup=speed,passed=speed>=10)); save('benchmark',report)
def statistics_suite():
    run('bin-pilot','snapshot',10000,9000000); bins=bins_from(load('bin-pilot')); save('bins',bins)
    for seed in [42,1000000,2000000]:
        for mode in ['robust','thinning']: run(f'stats-{seed}-{mode}',mode,100000,seed)
        compare(load(f'stats-{seed}-robust'),load(f'stats-{seed}-thinning'),bins,f'stats-{seed}-report')
def convergence():
    bins=json.loads((OUT/'bins.json').read_text())
    for mode in ['reference','thinning']:
        for label,extra in [('normal',[]),('tight',['--accuracy=0.1']),('orbit',['--step-scale=0.5']),('nocache',['--no-cache'])]: run(f'convergence-{mode}-{label}',mode,10000,4000000,extra)
        for label in ['tight','orbit','nocache']: compare(load(f'convergence-{mode}-normal'),load(f'convergence-{mode}-{label}'),bins,f'convergence-{mode}-{label}-report')
    compare(load('convergence-reference-tight'),load('convergence-thinning-tight'),bins,'convergence-tight-reference-report')
    for seed in [537,789,696,983,854,1039]:
        for mode in ['reference','thinning']:
            for label,extra in [('normal',[]),('tight',['--accuracy=0.1']),('orbit',['--step-scale=0.5'])]: run(f'regression-{seed}-{mode}-{label}',mode,1,seed,extra)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('stage',choices=['build','benchmark','statistics','convergence','all']); args=p.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    if args.stage in ('build','all'): build()
    if args.stage in ('benchmark','all'): benchmark()
    if args.stage in ('statistics','all'): statistics_suite()
    if args.stage in ('convergence','all'): convergence()

# Each stage records a machine-readable result. A statistics or convergence
# failure must also fail the command, rather than looking like a successful run.
if __name__ == '__main__':
    reports = []
    if args.stage in ('statistics', 'all'):
        reports += [OUT / f'stats-{s}-report.json' for s in [42, 1000000, 2000000]]
    if args.stage in ('convergence', 'all'):
        reports += list(OUT.glob('convergence-*-report.json'))
    passed = all(json.loads(p.read_text())['passed'] for p in reports)
    if args.stage in ('benchmark', 'all'):
        passed &= all(r['passed'] for r in json.loads((OUT / 'benchmark.json').read_text()))
    if not passed:
        raise SystemExit('Validation did not pass; inspect JSON reports.')
