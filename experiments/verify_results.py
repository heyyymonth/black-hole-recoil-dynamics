"""Recompute recorded table values and verify their source/output provenance."""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recoil_dynamics import PlummerOrbit


def read_csv(relative):
    with (ROOT/relative).open() as f: return list(csv.DictReader(f))


def equal(actual,expected,label):
    if not math.isclose(float(actual),expected,rel_tol=2e-9,abs_tol=2e-10):
        raise AssertionError(f'{label}: recorded {actual}, recomputed {expected}')


def main():
    cfg=json.loads((ROOT/'experiments/config.json').read_text())
    meta=json.loads((ROOT/'results/run_all.json').read_text())
    if meta['experiments_run']!=['e001','e002']: raise AssertionError('Expected complete recorded run')
    for group in ['source_sha256','output_sha256']:
        for relative,expected in meta[group].items():
            actual=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
            if actual!=expected: raise AssertionError(f'{group} mismatch: {relative}')
    rows=read_csv('results/e001/return_times.csv')
    assert len(rows)==len(cfg['e001']['eta_table'])
    checked=len(rows)
    for row in rows:
        orbit=PlummerOrbit(float(row['eta']),epsabs=cfg['epsabs'],epsrel=cfg['epsrel'])
        equal(row['first_return_over_tdyn'],orbit.period,'first return')
        equal(row['apocenter_over_b'],orbit.apocenter,'apocenter')
    rows=read_csv('results/e002/cycles.csv')
    assert len(rows)==sum(e>0 for e in cfg['e002']['eta'])
    checked+=len(rows)
    for row in rows:
        orbit=PlummerOrbit(float(row['eta']),epsabs=cfg['epsabs'],epsrel=cfg['epsrel'])
        equal(row['period_over_tdyn'],orbit.period,'cycle period')
        equal(row['cycle_density_exposure'],orbit.cycle_exposure(),'cycle exposure')
        equal(row['cycle_mean_density'],orbit.cycle_exposure()/orbit.period,'cycle mean')
        for c in cfg['e002']['core_radii_over_b']:
            equal(row[f'core_fraction_c_{c:g}'],orbit.cycle_core_time(c)/orbit.period,'cycle core')
    for name,expected_rows in [('fixed_horizons',len(cfg['e002']['eta'])*len(cfg['e002']['fixed_horizons_over_tdyn'])),
                               ('relative_horizons',sum(e>0 for e in cfg['e002']['eta'])*len(cfg['e002']['relative_horizons_over_period']))]:
        rows=read_csv(f'results/e002/{name}.csv')
        assert len(rows)==expected_rows
        checked+=len(rows)
        for row in rows:
            orbit=PlummerOrbit(float(row['eta']),epsabs=cfg['epsabs'],epsrel=cfg['epsrel'])
            h=float(row['horizon_over_tdyn']);exposure=orbit.exposure(h)
            equal(row['density_exposure'],exposure,'finite-window exposure')
            equal(row['mean_density'],exposure/h,'finite-window mean')
            equal(row['central_exposure'],h,'central comparator')
            wait=max(h-orbit.period,0) if orbit.period is not None else h
            equal(row['wait_then_central_exposure'],wait,'wait comparator')
            equal(row['wait_error_over_horizon'],(wait-exposure)/h,'wait error')
            for c in cfg['e002']['core_radii_over_b']:
                equal(row[f'core_fraction_c_{c:g}'],orbit.core_time(h,c)/h,'finite-window core')
    print(f'Verified {checked} result rows and all recorded source/output hashes.')


if __name__=='__main__': main()
