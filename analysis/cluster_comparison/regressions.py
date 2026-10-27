"""Reproduce the detected outer-interpolation failure on two fixed pilot meshes.

This is a numerical regression diagnostic. The unbounded interpolant is never
used for accepted host or transport results. Regenerate pilot meshes first.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'analysis/host_reconstruction'))
from host import Host

def main():
    rows=[]
    for label in ['pilot-02','pilot-05']:
        path=ROOT/'tmp/cluster-comparison/NGC3201'/label/'host.npz'
        arrays=dict(np.load(path,allow_pickle=False));fixed=Host(arrays);original=Host(arrays)
        for profile in [original.rho_spline,original.pr_spline,original.pt_spline]:profile.shape_fallback[:]=False
        star=int(np.flatnonzero(fixed.types=='MS')[-1]);grid=np.r_[0,np.geomspace(1e-10,fixed.rt,100000)]
        old=original.rho_spline(grid)[star];new=fixed.rho_spline(grid)[star]
        q_old=[original.project(0,n)[:,star] for n in [128,256,512]]
        q_new=[fixed.project(0,n)[:,star] for n in [128,256,512]]
        rows.append({'cluster':'NGC3201','label':label,'host_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                     'original_log_cubic_max_density_over_nodal_max':float(old.max()/fixed.rhoj[star].max()),
                     'shape_preserving_max_density_over_nodal_max':float(new.max()/fixed.rhoj[star].max()),
                     'original_maximum_radius_pc':float(grid[old.argmax()]),
                     'original_central_projection_128_256_relative_change':float(np.max(np.abs(q_old[1]/q_old[0]-1))),
                     'original_central_projection_256_512_relative_change':float(np.max(np.abs(q_old[2]/q_old[1]-1))),
                     'corrected_central_projection_128_256_relative_change':float(np.max(np.abs(q_new[1]/q_new[0]-1))),
                     'corrected_central_projection_256_512_relative_change':float(np.max(np.abs(q_new[2]/q_new[1]-1)))})
    record={'interpretation':'Deliberate reproduction of the original unbounded log-cubic behavior versus the shape-preserving interval replacement, on unchanged reconstructed meshes. These are regression checks, not additional cluster experiments.','rows':rows}
    p=ROOT/'results/cluster-comparison/numerical-regressions.json';p.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
