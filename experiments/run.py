"""Execute the declared experiments and record original numerical results."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from recoil_dynamics import PlummerOrbit, core_exit_threshold


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, path):
    fig.savefig(path.with_suffix(".png"), dpi=190)
    fig.savefig(path.with_suffix(".svg"), metadata={"Date": None})
    svg = path.with_suffix(".svg")
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)


def e001(config):
    out = ROOT / "results/e001"
    out.mkdir(parents=True, exist_ok=True)
    rows=[]
    for eta in config['e001']['eta_table']:
        orb = PlummerOrbit(eta, epsabs=config['epsabs'], epsrel=config['epsrel'])
        rows.append({'eta': eta, 'apocenter_over_b': orb.apocenter,
                     'first_return_over_tdyn': orb.period})
    write_csv(out / 'return_times.csv', rows)
    cfg=config['e001']
    eta = np.unique(np.r_[np.linspace(.001,.95,cfg['plot_points']),
                         1-np.geomspace(.05,1-cfg['eta_plot_max'],cfg['plot_points'])])
    periods = [PlummerOrbit(float(e),epsabs=config['epsabs'],epsrel=config['epsrel']).period for e in eta]
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    ax.semilogy(eta,periods,color='#14677d',lw=2,label='Energy quadrature')
    mask=eta>=.8
    ax.semilogy(eta[mask],math.pi/math.sqrt(2)*(1-eta[mask]**2)**(-1.5),
                '--',color='#be6d22',label='Near-escape asymptote')
    ax.axhline(math.pi,color='#888',ls=':',label='Small-kick limit')
    ax.set(xlabel='Kick speed / central escape speed',ylabel='First central return / dynamical time',
           title='E001: ballistic return in a static Plummer potential',xlim=(0,1))
    ax.legend(fontsize=9)
    save_figure(fig,out/'return_time')
    return rows


def e002(config):
    out=ROOT/'results/e002';out.mkdir(parents=True,exist_ok=True)
    cfg=config['e002'];fixed=[];relative=[];cycle=[]
    for eta in cfg['eta']:
        orb=PlummerOrbit(eta,epsabs=config['epsabs'],epsrel=config['epsrel'])
        if orb.period is not None:
            row={'eta':eta,'period_over_tdyn':orb.period,'cycle_density_exposure':orb.cycle_exposure(),
                 'cycle_mean_density':orb.cycle_exposure()/orb.period}
            for c in cfg['core_radii_over_b']:
                row[f'core_fraction_c_{c:g}']=orb.cycle_core_time(c)/orb.period
            cycle.append(row)
        for family,horizons in [('fixed',cfg['fixed_horizons_over_tdyn']),
                                ('relative',cfg['relative_horizons_over_period'])]:
            if family=='relative' and orb.period is None: continue
            for h in horizons:
                H=h if family=='fixed' else h*orb.period
                I=orb.exposure(H)
                # A stationary object has no return gate or lost time.
                wait=max(H-orb.period,0.) if orb.period is not None else H
                row={'eta':eta,'horizon_over_tdyn':H,'horizon_over_period':H/orb.period if orb.period else '',
                     'density_exposure':I,'mean_density':I/H,
                     'central_exposure':H,'wait_then_central_exposure':wait,
                     'wait_error_over_horizon':(wait-I)/H}
                for c in cfg['core_radii_over_b']:
                    row[f'core_fraction_c_{c:g}']=orb.core_time(H,c)/H
                (fixed if family=='fixed' else relative).append(row)
    for name,rows in [('cycles.csv',cycle),('fixed_horizons.csv',fixed),('relative_horizons.csv',relative)]:
        write_csv(out/name,rows)
    ts=np.unique(np.r_[np.geomspace(1e-6,.05,80),
                       np.linspace(.05,cfg['plot_horizon_max_periods'],cfg['plot_points'])])
    fig,axes=plt.subplots(1,3,figsize=(11,3.8),layout='constrained',sharey=True)
    for ax,eta in zip(axes,cfg['plot_eta']):
        orb=PlummerOrbit(eta,epsabs=config['epsabs'],epsrel=config['epsrel'])
        values=np.array([orb.exposure(float(t*orb.period))/(t*orb.period) for t in ts])
        ax.plot(ts,values,color='#14677d',lw=2,label='Continuous orbit')
        ax.plot(ts,np.maximum(ts-1,0)/ts,color='#be6d22',ls='--',label='First return, then central')
        ax.axhline(1,color='#777',ls=':',label='Assumed always central')
        ax.set(title=f'Kick / escape = {eta:g}',xlabel='Horizon / first-return period',ylim=(0,1.04),xlim=(0,10))
    axes[0].set_ylabel('Density exposure / horizon')
    axes[1].legend(loc='upper center',fontsize=7)
    fig.suptitle('E002: continuous exposure versus a first-return approximation',fontsize=13)
    save_figure(fig,out/'exposure_comparison')
    thresholds=np.array([core_exit_threshold(c) for c in cfg['core_radii_over_b']])
    eta_grid=np.unique(np.r_[np.linspace(.001,.95,cfg['plot_points']),
                            1-np.geomspace(.05,.001,cfg['plot_points']),
                            thresholds*(1-1e-5),thresholds,thresholds*(1+1e-5)])
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    for c,col in zip(cfg['core_radii_over_b'],['#14677d','#be6d22','#78539b']):
        values=[]
        for eta in eta_grid:
            orb=PlummerOrbit(float(eta),epsabs=config['epsabs'],epsrel=config['epsrel'])
            values.append(orb.cycle_core_time(c)/orb.period)
        ax.semilogy(eta_grid,values,color=col,lw=2,label=f'Operational radius = {c:g} b')
    ax.set(xlabel='Kick speed / central escape speed',ylabel='Fraction of one period inside the sphere',
           title='E002: central residence of bound radial orbits',xlim=(0,1),ylim=(1e-6,1.3))
    ax.legend(fontsize=9)
    save_figure(fig,out/'core_residence')
    return {'cycles':cycle,'fixed_rows':len(fixed),'relative_rows':len(relative)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment',choices=['all','e001','e002'],default='all')
    args=parser.parse_args()
    config=json.loads((ROOT/'experiments/config.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'recoil-dynamics-v1'})
    summaries={}
    for key,run in [('e001',e001),('e002',e002)]:
        if args.experiment in ['all',key]: summaries[key]=run(config)
    inputs=[ROOT/'experiments/config.json',ROOT/'experiments/run.py',ROOT/'pyproject.toml',ROOT/'requirements.txt',ROOT/'requirements-lock.txt']
    inputs+=sorted((ROOT/'src').rglob('*.py'))
    git=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,text=True,capture_output=True)
    manifest={'researcher':'Venkata Hemanth Athota','research_identity':'Black-hole Recoil Dynamics',
              'run_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),
              'versions':{x:importlib.metadata.version(x) for x in ['numpy','scipy','matplotlib']},
              'source_revision':git.stdout.strip() if git.returncode==0 else None,
              'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
              'output_sha256':{str(p.relative_to(ROOT)):sha(p) for key in summaries
                               for p in sorted((ROOT/'results'/key).iterdir())
                               if p.suffix in ['.csv','.png','.svg']},
              'experiments_run':list(summaries),'summaries':summaries,
              'interpretation':'Controlled collisionless orbital benchmarks; no observed-data inference or merger probabilities.'}
    # A run-specific manifest avoids claiming stale results were regenerated.
    (ROOT/'results'/f'run_{args.experiment}.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'experiments':list(summaries),'e002_counts':{k:v for k,v in summaries.get('e002',{}).items() if k.endswith('rows')}},indent=2))


if __name__=='__main__': main()
