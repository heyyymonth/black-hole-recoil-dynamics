"""Recompute every E003 table field and check source/output hashes."""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recoil_dynamics.encounters import experiment_rows


def close(actual, expected, label):
    if not math.isclose(float(actual), float(expected), rel_tol=2e-9, abs_tol=1e-16):
        raise AssertionError(f'{label}: recorded {actual}, recomputed {expected}')


def main():
    config = json.loads((ROOT/'experiments/e003.json').read_text())
    meta = json.loads((ROOT/'results/run_e003.json').read_text())
    assert meta['experiment'] == 'e003'
    for group in ['source_sha256', 'output_sha256']:
        for relative, expected in meta[group].items():
            actual = hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
            if actual != expected:
                raise AssertionError(f'{group} mismatch: {relative}')
    scales, expected_rows = experiment_rows(config)
    with (ROOT/'results/e003/encounters.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(expected_rows) == meta['rows']
    for index, (row, expected) in enumerate(zip(rows, expected_rows)):
        assert row.keys() == expected.keys()
        for key, value in expected.items():
            close(row[key], value, f'row {index}, {key}')
    summary = json.loads((ROOT/'results/e003/summary.json').read_text())
    assert summary['scales'].keys() == scales.keys()
    for key, value in scales.items():
        close(summary['scales'][key], value, key)
    close(summary['proxy_relative_error_tolerance'], config['proxy_relative_error_tolerance'], 'tolerance')
    close(summary['maximum_absolute_proxy_relative_error'],
          max(abs(r['proxy_relative_error']) for r in expected_rows), 'maximum error')
    assert summary['nonzero_kicks_exceeding_tolerance'] == [
        r['eta'] for r in expected_rows if r['eta'] and r['exceeds_tolerance']]
    print(f'Verified {len(rows)} E003 rows, physical scales, decision summary, and all recorded hashes.')


if __name__ == '__main__':
    main()
