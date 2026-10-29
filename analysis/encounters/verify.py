"""Check recorded finite-angle calculations without running physical integrals."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/finite-angle-encounters'
PROTOCOL = ROOT/'analysis/encounters/protocol.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(record=False):
    p = json.loads(PROTOCOL.read_text())
    summary = json.loads((OUT/'summary.json').read_text())
    all_rows = []
    input_files = [PROTOCOL, *(ROOT/'analysis/encounters').glob('*.py')]
    output_files = [OUT/'summary.json', OUT/'states.csv']
    maxima = []
    for cluster in p['clusters']:
        host = ROOT/f'results/cluster-comparison/{cluster}/reference/host.npz'
        file = OUT/f'{cluster}.json'
        r = json.loads(file.read_text())
        input_files.append(host); output_files.append(file)
        assert r['host_sha256'] == digest(host)
        assert r['protocol_sha256'] == digest(PROTOCOL)
        for name, expected in r['implementation_sha256'].items():
            assert digest(ROOT/'analysis/encounters'/name) == expected, name
        assert len(r['states']) == len(r['tensors']) == 27
        assert r['all_numerical_checks_passed']
        for row, tensor in zip(r['states'], r['tensors'], strict=True):
            assert row['accepted'] and row['state'] == tensor['state']
            A, D = np.asarray(tensor['A']), np.asarray(tensor['D'])
            np.testing.assert_allclose(D, D.T, atol=1e-14)
            assert np.linalg.eigvalsh(D).min() >= 0
            np.testing.assert_allclose(np.trace(D), row['trace_D_kms2_per_Myr'])
            np.testing.assert_allclose(row['allowed_rate_per_Myr']+row['excluded_rate_per_Myr'], row['parent_rate_per_Myr'],rtol=1e-12)
            np.testing.assert_allclose(row['bmax_pc']/row['local_scale_pc'], row['bmax_over_local_scale'])
            assert 0 <= row['b90_above_L_trace_fraction'] <= 1+1e-14
            assert 0 <= row['excluded_flux_fraction'] <= 1
            np.testing.assert_allclose(row['rms_increment_over_sigma'],
                np.sqrt(np.trace(D)/row['allowed_rate_per_Myr'])/row['sigma_BH_1d_kms'])
            a = p['acceptance']
            assert max(row[k] for k in ('A_refinement','D_refinement','rate_refinement','excluded_rate_refinement')) < a['moment_relative_refinement']
            assert row['diagnostic_fraction_refinement'] < a['diagnostic_fraction_absolute_refinement']
            assert row['number_density_relative_error'] < a['number_density_relative_error']
            if row['speed_over_sigma'] == 0:
                assert np.linalg.norm(A) < a['zero_drift_relative_scale']*np.trace(D)/row['sigma_BH_1d_kms']
        s = r['sampling']
        expected, sampled, se = [np.asarray(s[k]) for k in ('expected_conditional_moments','sampled_conditional_moments','standard_errors')]
        assert np.all(se>0)
        np.testing.assert_allclose(abs(sampled-expected)/se,s['standardized_errors'])
        assert max(s['standardized_errors']) < p['sampling']['standard_error_limit']
        maxima.append(max(s['standardized_errors']))
        all_rows.extend(r['states'])
    assert summary['states'] == summary['accepted_states'] == len(all_rows) == 162
    assert summary['clusters'] == len(p['clusters']) == 6
    assert summary['sampled_elastic_events'] == 6*p['sampling']['events_per_host']
    assert summary['all_numerical_checks_passed']
    assert summary['transport_simulations_performed'] is False and summary['residence_response_evaluated'] is False
    np.testing.assert_allclose(summary['max_sample_standardized_error'],max(maxima))
    with (OUT/'states.csv').open() as f:
        csv_rows=list(csv.DictReader(f))
    assert len(csv_rows)==len(all_rows)
    for csv_row,row in zip(csv_rows,all_rows,strict=True):
        assert csv_row=={k:str(v) for k,v in row.items()}
    provenance = {'revision':p['revision'],
                  'scope':p['interpretation'],
                  'input_sha256':{str(f.relative_to(ROOT)):digest(f) for f in sorted(set(input_files))},
                  'output_sha256':{str(f.relative_to(ROOT)):digest(f) for f in sorted(output_files)}}
    path=OUT/'provenance.json'
    if record:
        path.write_text(json.dumps(provenance,indent=2)+'\n')
    else:
        assert json.loads(path.read_text())==provenance, 'Stale encounter provenance'
    print('Verified 162 conditional finite-angle states and 600000 sampled increments; no residence results.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true')
    verify(parser.parse_args().record)
