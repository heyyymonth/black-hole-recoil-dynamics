"""Post-pilot descriptive concentration of the estimated sampling variance.

This diagnostic was specified after the completed efficiency pilot. It does
not alter scores, sample allocation, the six-ratio family or the summary.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.precision.statistics import estimate


def calculate():
    protocol_path = ROOT/'analysis/precision/protocol.json'
    p = json.loads(protocol_path.read_text())
    files = [protocol_path, Path(__file__), ROOT/'analysis/precision/statistics.py']
    rows = []
    for launch in p['launches']:
        for method in p['methods']:
            paths = [ROOT/f"results/precision/{launch}-{method['name']}-block{block}.npz"
                     for block in range(p['blocks'])]
            scores = []; labels = []
            for path in paths:
                with np.load(path, allow_pickle=False) as data:
                    scores.append(data['response_sample'])
                    labels.append(data['insertion_stratum'])
            z = np.concatenate(scores); strata = np.concatenate(labels)
            stratified = method['selection']=='stratified'
            contribution = np.zeros_like(z)
            for label in (range(4) if stratified else [-1]):
                mask = strata==label; values = z[mask]; n = len(values)
                weight = .25 if stratified else 1.
                contribution[mask] = weight**2*(values-values.mean())**2/(n*(n-1))
            total = float(contribution.sum())
            expected = estimate(z, strata if stratified else None)['variance_of_mean']
            np.testing.assert_allclose(total, expected, rtol=1e-12, atol=1e-20)
            ordered = np.sort(contribution)[::-1]; count = int(np.ceil(.01*len(z)))
            rows.append({'launch':launch, 'method':method['name'], 'n':len(z),
                'variance_of_mean':total, 'largest_one_percent_count':count,
                'largest_single_share':float(ordered[0]/total) if total else None,
                'largest_one_percent_share':float(ordered[:count].sum()/total) if total else None})
            files.extend(paths)
    return {'scope':'Post-pilot descriptive quadratic contributions to estimated mean variance; no population-tail or confidence claim and no trimming',
        'formula':'For stratum s, c_i = w_s^2 (Z_i - mean_s)^2 / [n_s (n_s-1)]; sum_i c_i estimates Var(mean). Uniform samples have one stratum with weight 1.',
        'rows':rows,
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in files}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    path=ROOT/'results/precision/tail-diagnostic.json';result=calculate()
    if args.verify:
        stored=json.loads(path.read_text())
        from analysis.response.verify import compare_tree
        compare_tree(stored,result)
        print('Verified post-pilot variance-concentration diagnostic; no trajectories rerun.')
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print('Recorded post-pilot variance concentration without changing the pilot summary.')


if __name__=='__main__':main()
