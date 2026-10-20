"""Recompute the external-snapshot audit and verify all recorded provenance."""
from pathlib import Path
import csv
import hashlib
import json
import math

from run_e004 import calculate

ROOT = Path(__file__).resolve().parents[1]


def check(actual, expected, label):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), label
        for key, value in expected.items():
            check(actual[key], value, label+'/'+key)
    elif isinstance(expected, list):
        assert len(actual) == len(expected), label
        for index, (a, e) in enumerate(zip(actual, expected)):
            check(a, e, label+f'/{index}')
    elif expected is None:
        assert actual in (None, ''), label
    elif isinstance(expected, bool):
        assert actual is expected, label
    elif isinstance(expected, (int, float)):
        assert math.isclose(float(actual), expected, rel_tol=2e-8, abs_tol=2e-12), (label, actual, expected)
    else:
        assert actual == expected, label


def main():
    meta = json.loads((ROOT/'results/run_e004.json').read_text())
    assert meta['experiment'] == 'e004'
    for group in ['source_sha256', 'output_sha256']:
        for name, expected in meta[group].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, (group, name)
    config = json.loads((ROOT/'experiments/e004.json').read_text())
    assert meta['external_revision'] == config['source_revision']
    assert meta['external_sha256'] == {k: v['sha256'] for k, v in config['files'].items()}
    tables, summary = calculate(config)
    for name, expected in tables.items():
        with (ROOT/'results/e004'/f'{name}.csv').open() as stream:
            actual = list(csv.DictReader(stream))
        check(actual, expected, name)
        assert len(actual) == meta['rows'][name]
    check(json.loads((ROOT/'results/e004/summary.json').read_text()), summary, 'summary')
    print('Verified E004: external input hashes, 5 population rows, 12 trajectory rows, '
          '160 profile points, summary and source/output hashes.')


if __name__ == '__main__':
    main()
