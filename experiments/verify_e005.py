"""Reproduce E005 and independently check its primary scalar census."""
from pathlib import Path
import csv
import hashlib
import json
import math
import statistics

from run_e005 import calculate

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
        assert actual is expected or actual == str(expected), label
    elif isinstance(expected, (int, float)):
        assert math.isclose(float(actual), expected, rel_tol=1e-12, abs_tol=1e-12), (label, actual, expected)
    else:
        assert actual == expected, label


def independent_census(paths, config, tables, summary):
    """No production parser/reducer: extract scalars, select by a separate scan."""
    histories = {}
    for name, path in paths.items():
        fields = dict(line.split() for line in path.read_text().splitlines() if line.strip())
        histories.setdefault(name.split('/')[1], []).append(
            (float(fields['age']), name, int(fields['Nbin_14_14'])))
    selected = {}
    for model, states in histories.items():
        best = None
        for age, name, count in states:
            distance = abs(age-config['comparison_age_gyr'])
            if distance > config['age_tolerance_gyr']:
                continue
            candidate = (distance, age, name, count)
            if best is None or candidate[:3] < best[:3]:
                best = candidate
        if best is not None:
            selected[model] = best
    assert {r['model']: (abs(r['age_offset_gyr']), r['age_gyr'], r['path'], r['bbh_count'])
            for r in tables['selected']} == selected
    counts = [s[3] for s in selected.values()]
    assert len(selected) == summary['selected_models']
    assert statistics.median(counts) == summary['selected_bbh_counts']['median']
    assert sum(c == 0 for c in counts) == summary['selected_zero_bbh_models']
    assert sum(c >= config['minimum_bbh_count_screen'] for c in counts) == summary['selected_meeting_count_screen_models']
    assert sorted(set(histories)-set(selected)) == summary['models_without_age_coverage']


def main():
    meta = json.loads((ROOT/'results/run_e005.json').read_text())
    assert meta['experiment'] == 'e005'
    for group in ['source_sha256', 'output_sha256']:
        for name, expected in meta[group].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, (group, name)
    config = json.loads((ROOT/'experiments/e005.json').read_text())
    assert meta['external_revision'] == config['source_revision']
    tables, summary, paths = calculate(config)
    assert meta['external_sha256'].keys() == paths.keys()
    for name, path in paths.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == meta['external_sha256'][name], name
    for name, expected in tables.items():
        with (ROOT/'results/e005'/f'{name}.csv').open() as stream:
            actual = list(csv.DictReader(stream))
        check(actual, expected, name)
        assert len(actual) == meta['rows'][name]
    check(json.loads((ROOT/'results/e005/summary.json').read_text()), summary, 'summary')
    independent_census(paths, config, tables, summary)
    print(f"Verified E005: {len(paths)} exact external inputs, {len(tables['models'])} histories, "
          f"{len(tables['selected'])} selected states; all tables, source/output hashes, E004 cross-checks "
          "and independent scalar census.")


if __name__ == '__main__':
    main()
