"""Summarize published derived draws without inventing raw-chain row alignment."""
from pathlib import Path
import hashlib,json,csv,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
CLUSTERS=['NGC0104','NGC1851','NGC3201','NGC5139','NGC6121','NGC6397']

def quantiles(a):
    a=np.asarray(a,dtype=float)
    if not np.isfinite(a).all():raise ValueError('Nonfinite derived quantity')
    return dict(zip(['p16','median','p84'],np.quantile(a,[.16,.5,.84]).tolist())) if len(a) else None

def summarize(number,mass,radius):
    number,mass,radius=map(np.asarray,[number,mass,radius])
    if number.shape!=mass.shape or radius.shape!=mass.shape or mass.ndim!=1:raise ValueError('Incompatible derived rows')
    if not np.isfinite(number).all() or not np.isfinite(mass).all() or np.any(number<0) or np.any(mass<0):raise ValueError('Invalid BH quantity')
    zero=mass==0
    if not np.array_equal(zero,number==0):raise ValueError('Inconsistent absent BH population')
    if np.any(radius[~zero]<=0) or not np.isfinite(radius[~zero]).all():raise ValueError('Invalid occupied BH radius')
    return {'derived_rows':len(mass),'zero_BH_rows':int(zero.sum()),'zero_BH_fraction':float(zero.mean()),
            'BH_expected_number':quantiles(number),'BH_mass_Msun':quantiles(mass),
            'BH_half_mass_radius_pc_given_present':quantiles(radius[~zero]),
            'probe_to_BH_mass_ratio_given_present':quantiles(60/mass[~zero]),
            'rows_BH_mass_at_most_probe':int(((mass>0)&(mass<=60)).sum())}

def main():
    import h5py
    pins={x['file']:x for x in json.loads((ROOT/'data/observed-clusters/sources.json').read_text())['files']}
    rows=[]
    for c in CLUSTERS:
        path=ROOT/f'data/observed-clusters/raw/gcfit_{c}_sampler.hdf'
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if digest!=pins[path.name]['sha256']:raise ValueError('Source hash mismatch')
        with h5py.File(path) as f:
            group=f['model_output/quantities'];q=summarize(group['BH_num'][:],group['BH_mass'][:],group['BH_rh'][:])
            logw=f['nested/logwt'][:];w=np.exp(logw-np.max(logw));w/=w.sum()
            samples=f['nested/samples'][:]
            if samples.shape!=(len(w),13) or not np.isfinite(samples).all():raise ValueError('Invalid joint chain')
            units={k:str(group[k].attrs['unit']) for k in ['BH_num','BH_mass','BH_rh']}
            if units!={'BH_num':'','BH_mass':'solMass','BH_rh':'pc'}:raise ValueError('Unknown derived units')
        rows.append({'cluster':c,**q,'joint_rows':len(w),'weight_ESS':float(1/(w@w)),
                     'source_sha256':digest,'source_url':pins[path.name]['url'],'source_units':units,
                     'preferred_model_chain_available':c!='NGC6397'})
    out=ROOT/'results/cluster-comparison';out.mkdir(exist_ok=True,parents=True)
    report={'revision':'six-cluster-published-derived-1','interpretation':'Empirical quantiles of released 4000 derived model realizations; not six new fits or direct BH detections. Historical compatible GCfit constructs these from equal-weight chains. Original draw mapping and generating random seed are not released; no alignment to raw nested rows is assumed.',
            'zero_policy':'Zero-BH rows retained; BH subsystem radius and probe/reservoir mass ratio undefined there and excluded only from explicitly conditional summaries.',
            'clusters':rows}
    (out/'posterior-inventory.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    flat=[]
    for r in rows:
        a={k:r[k] for k in ['cluster','derived_rows','zero_BH_rows','zero_BH_fraction','joint_rows','weight_ESS']}
        for k in ['BH_expected_number','BH_mass_Msun','BH_half_mass_radius_pc_given_present','probe_to_BH_mass_ratio_given_present']:
            for stat,v in r[k].items():a[f'{k}_{stat}']=v
        flat.append(a)
    with (out/'posterior-inventory.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=flat[0],lineterminator='\n');w.writeheader();w.writerows(flat)
    print(json.dumps([{k:r[k] for k in ['cluster','BH_expected_number','zero_BH_fraction']} for r in rows],indent=2))
if __name__=='__main__':main()
