"""Reconstruct and independently validate one published joint host sample."""
from pathlib import Path
import csv
import hashlib
import json
import os
import sys
import warnings

ROOT=Path(__file__).resolve().parents[2]
os.environ['GCFIT_DIR']=str(ROOT/'tmp/host-reconstruction/gcfit-cache')
import gcfit
import h5py
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simpson, cumulative_simpson
from scipy.interpolate import CubicSpline
from scipy.special import gammainc

OUT=ROOT/'results/host-reconstruction'

def values(x):
    return np.asarray(getattr(x,'value',x))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

from host import Host


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)

def run():
    protocol=json.loads((Path(__file__).with_name('protocol.json')).read_text())
    for path,digest in protocol['sources'].items():
        if sha(ROOT/path)!=digest:raise ValueError('Source hash mismatch: '+path)
    obs_path=ROOT/'data/observed-clusters/raw/gcfit_NGC3201_observations.hdf'
    observations=gcfit.Observations(str(obs_path),standardize_name=False,restrict_to='local')
    sampler=ROOT/'data/observed-clusters/raw/gcfit_NGC3201_sampler.hdf'
    with h5py.File(sampler) as f:
        logw=f['nested/logwt'][:];i=int(np.argmax(logw));row=f['nested/samples'][i]
        weight=np.exp(logw-logw.max());weight/=weight.sum()
        rpub=f['model_output/metadata/r'][:]
        published={key:f['model_output/profiles/'+key][0] for key in ['LOS','pm_R','pm_T','rho_tot']}
        mjpub=f['model_output/metadata/mj'][:]
    theta=dict(zip(protocol['parameter_order'],row))
    if '--reuse-solution' in sys.argv:
        arrays=dict(np.load(OUT/'host.npz',allow_pickle=False))
        previous=json.loads((OUT/'summary.json').read_text())
        if previous['theta']!=theta:raise ValueError('Cached solution has different parameters')
        warning_text=previous['warnings']
        star_bin=int(np.where(arrays['types']=='MS')[0][-1])
    else:
        print('Constructing fixed joint row',i,flush=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            model=gcfit.FittableModel(theta,observations,ode_rtol=protocol['host']['ode_rtol'],ode_maxstep=protocol['host']['ode_maxstep'])
        arrays={key:values(getattr(model,key)) for key in ['r','mj','Mj','rhoj','mcj','s2j','raj','A','v2rj','v2tj','G','M','rh','rt','g','v2pj','v2Rj','v2Tj','Sigmaj']}
        arrays['psi']=values(model.phi);arrays['types']=np.asarray(model.star_types)
        warning_text=sorted(set(str(w.message) for w in caught))
        star_bin=model.nms-1
    for key,array in arrays.items():
        if key != "types" and not np.isfinite(array).all():
            raise ValueError("Nonfinite host array: "+key)
    if np.any(arrays["mj"]<=0) or np.any(arrays["Mj"]<0):
        raise ValueError("Invalid species masses")
    host=Host(arrays)
    grid=np.r_[0,np.geomspace(host.rt*1e-9,host.rt,protocol['validation']['radial_grid_points']-1)]
    rho=host.rho_spline(grid);density=rho.sum(axis=0)
    cumulative=cumulative_simpson(4*np.pi*grid**2*density,x=grid,initial=0)
    mass=float(cumulative[-1]);masserr=abs(mass/host.M-1)
    rh_check=float(np.interp(mass/2,cumulative,grid))
    inner=(grid>.001*host.rh)&(grid<min(10*host.rh,.8*host.rt))
    force_expected=host.G*cumulative[inner]/grid[inner]**2
    force_relative=np.abs(-host.psi_spline(grid[inner],1)/force_expected-1)
    grad=np.zeros_like(grid);grad[1:]=host.G*cumulative[1:]/grid[1:]**2
    integ=cumulative_simpson(grad,x=grid,initial=0)
    psi_integral=integ[-1]-integ
    potential_error=float(np.max(np.abs(psi_integral-host.psi(grid)))/host.psi(0))
    kinetic=2*np.pi*simpson(grid**2*(host.pr_spline(grid)+host.pt_spline(grid)).sum(axis=0),x=grid)
    potential=-4*np.pi*host.G*simpson(grid*density*cumulative,x=grid)
    virial_error=float(abs(2*kinetic/potential+1))
    dfrows=[]
    for fraction in protocol['validation']['DF_radius_fractions_rh']:
        radius=fraction*host.rh;cutoff=np.sqrt(2*host.psi(radius))
        expected_rho=host.rho_spline(radius)
        expected_r=host.pr_spline(radius)/expected_rho
        expected_t=host.pt_spline(radius)/expected_rho
        previous=None
        for order in protocol['validation']['velocity_quadrature_orders']:
            x,w=leggauss(order);mu,wm=leggauss(order)
            speed=(x+1)*cutoff/2;wv=w*cutoff/2
            v,cos=np.meshgrid(speed,mu,indexing='ij')
            u=np.column_stack([(v*np.sqrt(1-cos*cos)).ravel(),np.zeros(v.size),(v*cos).ravel()])
            f=host.number_df(radius,u)*host.mj[:,None]
            weights=(2*np.pi*(wv*speed**2)[:,None]*wm[None,:]).ravel()
            moments=np.array([f@weights,(f*(u[:,2]**2))@weights,(f*(u[:,0]**2))@weights])
            error=np.max(np.abs(moments/np.array([expected_rho,expected_rho*expected_r,expected_rho*expected_t])-1))
            refinement=None if previous is None else float(np.max(np.abs(moments/previous-1)))
            dfrows.append({'r_over_rh':fraction,'order':order,'max_density_and_raw_moment_relative_error':float(error),'refinement_relative_change':refinement})
            previous=moments
    projected=[];projection_err=[];projection_ref=[]
    for R in host.rh*np.array([0.,.1,.3,1.,3.]):
        p128=host.project(R,128);p256=host.project(R,256)
        native=np.array([CubicSpline(host.r,arrays[key],axis=1)(R) for key in ['Sigmaj','v2pj','v2Rj','v2Tj']])
        err=float(np.max(np.abs(p256[:,star_bin]/native[:,star_bin]-1)))
        refinement=float(np.max(np.abs(p256[:,star_bin]/p128[:,star_bin]-1)))
        projection_err.append(err);projection_ref.append(refinement)
        projected.append({'R_pc':float(R),'sigma_LOS_kms':float(np.sqrt(p256[1,star_bin])),'sigma_PMR_kms':float(np.sqrt(p256[2,star_bin])),'sigma_PMT_kms':float(np.sqrt(p256[3,star_bin])),'observed_tracer_native_relative_difference':err,'all_bins_native_relative_difference':float(np.max(np.abs(p256/native-1))),'quadrature_relative_change':refinement})
    comparison=[]
    for k,pub in published.items():
        selection=(rpub>=.01*host.rh)&(rpub<=3*host.rh)
        if k=='rho_tot':prediction=host.rho_spline(rpub[selection]).sum(axis=0)
        else:
            native={'LOS':'v2pj','pm_R':'v2Rj','pm_T':'v2Tj'}[k]
            prediction=np.sqrt(CubicSpline(host.r,arrays[native][star_bin])(rpub[selection]))
        low,median,high=pub[4,selection],pub[2,selection],pub[0,selection]
        for R,p,l,med,h in zip(rpub[selection],prediction,low,median,high):
            comparison.append({'profile':k,'R_pc':float(R),'reconstructed':float(p),'published_lower_95_45':float(l),'published_median':float(med),'published_upper_95_45':float(h),'within_envelope':bool(l<=p<=h)})
    species=[]
    for j in range(len(host.mj)):
        species.append({'bin':j,'type':str(host.types[j]),'m_Msun':float(host.mj[j]),'M_Msun':float(host.Mj[j]),'expected_number':float(host.Mj[j]/host.mj[j]),'central_density_Msun_pc3':float(host.rhoj[j,0]),'central_sigma_r_kms':float(np.sqrt(host.v2rj[j,0])),'ra_pc':float(host.raj[j]),'included_in_scattering':bool(host.physical[j])})
    checks={'integrated_mass_relative_error':masserr,'half_mass_radius_relative_error':abs(rh_check/host.rh-1),'force_max_relative_error':float(force_relative.max()),'potential_max_error_over_central_relative_potential':potential_error,'virial_relative_error':virial_error,'DF_max_relative_error':max(d['max_density_and_raw_moment_relative_error'] for d in dfrows if d['order']==64),'projection_max_relative_difference':max(projection_err),'projection_max_refinement_change':max(projection_ref)}
    # A failed check remains in the result and blocks downstream adoption.
    tests={key:bool(value<(.002 if key.startswith('projection') else .001)) for key,value in checks.items()}
    result={'revision':'host-reconstruction-1','cluster':'NGC3201','protocol_sha256':sha(Path(__file__).with_name('protocol.json')),'joint_row_index':i,'normalized_posterior_weight':float(weight[i]),'theta':theta,'observations_metadata':{k:float(observations.mdata[k]) for k in ['age','FeH','vesc','Ndot']},'summary':{'M_Msun':host.M,'rh_pc':host.rh,'rt_pc':host.rt,'distance_kpc':theta['d'],'G_pc_kms2_per_Msun':host.G,'central_density_Msun_pc3':float(host.rhoj[:,0].sum()),'central_cutoff_kms':float(np.sqrt(2*host.psi(0))),'central_escape_to_infinity_kms':float(np.sqrt(2*(host.psi(0)+host.G*host.M/host.rt))),'physical_species_bins':int(sum(host.physical)),'tracer_mass_fraction':float(host.Mj[~host.physical].sum()/host.M),'total_BH_expected_number':float((host.Mj/host.mj)[host.types=='BH'].sum()),'total_BH_mass_Msun':float(host.Mj[host.types=='BH'].sum()),'published_tracer_masses_Msun':mjpub.tolist(),'selected_tracer_mass_Msun':float(host.mj[star_bin])},'checks':checks,'checks_passed':tests,'all_numerical_checks_passed':all(tests.values()),'warnings':warning_text,'interpretation':'One compatible reconstruction, not an exact historical refit or uncertainty-marginalized host. Published profile envelopes are descriptive; large truncation radius is not an empirical escape boundary.'}
    OUT.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(OUT/'host.npz',**arrays)
    (OUT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    write_csv(OUT/'species.csv',species);write_csv(OUT/'df-checks.csv',dfrows);write_csv(OUT/'projection-checks.csv',projected);write_csv(OUT/'published-profile-comparison.csv',comparison)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':run()
