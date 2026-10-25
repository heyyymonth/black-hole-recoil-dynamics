"""Inspect originals and harmonize variables; no fits, potential evaluation or orbits."""
from collections import Counter
from io import StringIO
from pathlib import Path
import hashlib
import json
import re
import zipfile

import h5py
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/observed-clusters"
RAW, OUT = DATA / "raw", DATA / "derived"
IDS = ["NGC0104", "NGC1851", "NGC3201", "NGC5139", "NGC6121", "NGC6397"]
NAMES = dict(zip(IDS, ["47 Tuc", "NGC 1851", "NGC 3201", "Omega Centauri", "M4", "NGC 6397"]))
PARAMS = ["W0", "M", "rh", "ra", "g", "delta", "s2", "F", "a1", "a2", "a3", "BHret", "d"]
G = 4.30091727003628e-3  # pc (km/s)^2 / Msun
PC_PER_KMS_MYR = 0.9777922216807892
PM_FACTOR = 4.740470463533349  # km/s per (kpc mas/yr); AU and Julian year
ARCSEC_PER_RAD = 206264.80624709636


def plain(value):
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        return [plain(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {k: plain(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    return value


def write_json(name, obj):
    (OUT/name).write_text(json.dumps(plain(obj), indent=2, allow_nan=False)+"\n")


def catalog(filename):
    soup = BeautifulSoup((RAW/filename).read_text(), "html.parser")
    rows=[]
    for tr in soup.find_all("tr"):
        cells=tr.find_all("td")
        if len(cells)!=27:
            continue
        text=[c.get_text(" ",strip=True).replace("\xa0"," ") for c in cells]
        match=re.match(r"NGC\s+(\d+)",text[0])
        if not match:
            continue
        cid=f"NGC{int(match[1]):04d}"
        if cid not in IDS:
            continue
        def pair(cell):
            numbers=re.findall(r"[-+]?\d+(?:\.\d+)?",cell)
            return [float(numbers[0]),float(numbers[1]) if "±" in cell else None]
        d,de=pair(text[3])
        mass_text=cells[7]
        exponent=int(mass_text.find("sup").get_text()) if mass_text.find("sup") else 0
        m,me=pair(text[7])
        mass=m*10**exponent
        rh=float(text[12])
        row=dict(cluster=cid,name=NAMES[cid],source_file=filename,
                 distance_kpc=d,distance_error_kpc=de,mass_msun=mass,
                 mass_error_msun=me*10**exponent if me is not None else None,
                 core_spitzer_pc=float(text[10]),half_light_projected_pc=float(text[11]),
                 half_mass_3d_pc=rh,tidal_catalog_pc=float(text[13]),
                 log10_core_density_msun_pc3=float(text[14]),
                 log10_relaxation_years=float(text[16]),
                 sigma0_1d_kms=float(text[21]),escape0_catalog_kms=float(text[22]),
                 rotation_amplitude_kms=pair(text[25])[0],
                 rotation_error_kms=pair(text[25])[1],
                 mass_inferred_from_luminosity=bool(cells[7].find("i")),
                 provenance_class="catalog_model_estimate",
                 radii_and_escape_uncertainty_available=False)
        row["scale_v_kms"]=(G*mass/rh)**0.5
        row["scale_t_myr"]=rh/row["scale_v_kms"]*PC_PER_KMS_MYR
        rows.append(row)
    assert set(r["cluster"] for r in rows)==set(IDS)
    assert len(rows)==6
    return pd.DataFrame(rows).set_index("cluster").loc[IDS].reset_index()


def numeric_lines(path, ncols):
    """Author table data lines must have exactly the documented numeric width."""
    rows=[]
    for line in path.read_text().splitlines():
        fields=line.split()
        if not fields:
            continue
        # Byte-range header tokens such as "1-" are descriptions, not data.
        try:
            float(fields[0])
        except ValueError:
            continue
        try:
            numeric = [float(v) for v in fields]
        except ValueError as exc:
            raise ValueError(f"Malformed numeric row in {path.name}: {line}") from exc
        assert len(numeric)==ncols, (path.name, 'numeric row width', len(numeric), ncols)
        rows.append(numeric)
    arr=np.asarray(rows)
    assert arr.ndim==2 and arr.shape[1]==ncols and np.isfinite(arr).all(),path
    return arr


def hst_profiles():
    rows=[]
    expected_bins=dict(zip(IDS,[28,20,6,40,6,7]))
    for cid in IDS:
        if cid=="NGC5139":
            filename="h25_proper_motion_dispersion_log_bins.dat"
            array=numeric_lines(RAW/filename,15)
            assert len(array)==expected_bins[cid]
            for b,a in enumerate(array,1):
                for component,start in [("PM_COMBINED",4),("PM_R",7),("PM_T",10)]:
                    rows.append(dict(cluster=cid,source="H25",source_file=filename,
                        bin_id=b,n_stars=int(a[3]),radius_arcsec=a[1],radius_lower_arcsec=a[0],
                        radius_upper_arcsec=a[2],component=component,sigma_masyr=a[start],
                        error_upper_masyr=a[start+1],error_lower_masyr=a[start+2],
                        tracer="full magnitude range; not matched to HACKS giant selection",
                        uncertainty="posterior 16/50/84 asymmetric offsets",
                        covariance_group=f"{cid}:H25:{b}",independent_third_component=False,
                        likelihood_ready=False))
        else:
            filename=f"hacks_{cid}.txt"
            array=numeric_lines(RAW/filename,9)
            assert len(array)==expected_bins[cid]
            assert np.array_equal(array[:,0],np.arange(1,len(array)+1))
            for a in array:
                for component,start in [("PM_COMBINED",3),("PM_R",5),("PM_T",7)]:
                    rows.append(dict(cluster=cid,source="HACKS",source_file=filename,
                        bin_id=int(a[0]),n_stars=int(a[1]),radius_arcsec=a[2],
                        radius_lower_arcsec=None,radius_upper_arcsec=None,
                        component=component,sigma_masyr=a[start],
                        error_upper_masyr=a[start+1],error_lower_masyr=a[start+1],
                        tracer="SGB/RGB brighter than main-sequence turnoff",
                        uncertainty="average of 16/84 offsets; bias-corrected dispersion",
                        covariance_group=f"{cid}:HACKS:{int(a[0])}",
                        independent_third_component=False,likelihood_ready=False))
    frame=pd.DataFrame(rows)
    assert not frame.duplicated(["cluster","source","bin_id","component"]).any()
    assert (frame.sigma_masyr>0).all()
    assert (frame[["error_upper_masyr","error_lower_masyr"]]>=0).all().all()
    return frame


def compiled_profiles():
    rows,surface,flags=[],[],[]
    for cid in IDS:
        filename=f"uq_{cid}_kin.txt"
        for line_no,line in enumerate((RAW/filename).read_text().splitlines(),1):
            if line.startswith("#") or not line.strip():continue
            bits=line.split()
            assert len(bits)==7
            _,r,s,up,lo,typ,src=bits
            r,s,up,lo=map(float,(r,s,up,lo))
            key=f"{cid}:{line_no}"
            reason=[]
            if r<=0:reason.append("nonpositive_radius_needs_original_mapping")
            if s<=0 or up<=0 or lo<=0:reason.append("nonpositive_dispersion_or_uncertainty")
            if src=="L24" and s==up==lo:reason.append("identical_dispersion_and_errors_unreconciled")
            row=dict(cluster=cid,source_file=filename,source_line=line_no,
                     radius_arcsec=r,sigma_kms=s,error_upper_kms=up,error_lower_kms=lo,
                     component="LOS" if typ=="1" else "PM_UNSPECIFIED_COMBINED",
                     literature_tag=src,quality_flags=";".join(reason),
                     quarantine=bool(reason),role="audit_only_do_not_stack_with_originals")
            rows.append(row)
            if reason:flags.append(dict(row_id=key,reasons=reason))
        a=[]
        for line_no,line in enumerate((RAW/f"uq_{cid}_sb.txt").read_text().splitlines(),1):
            b=line.split()
            if len(b)!=4 or not b[0].lower().startswith("ngc"):continue
            surface.append(dict(cluster=cid,source_file=f"uq_{cid}_sb.txt",source_line=line_no,
                log10_radius_arcsec=float(b[1]),radius_arcsec=10**float(b[1]),
                native_log_surface=float(b[2]),error_log_surface=float(b[3]),
                normalization="compiler stitched counts/brightness; absolute mass normalization unresolved",
                mass_density_ready=False))
    df=pd.DataFrame(rows)
    keys=["cluster","radius_arcsec","component","literature_tag"]
    duplicated=df[df.duplicated(keys,keep=False)]
    return df,pd.DataFrame(surface),flags,duplicated


def gaia_curves():
    members,rows=[],[]
    with zipfile.ZipFile(RAW/"gaia_edr3_v2.zip") as z:
        readme=z.read("!readme.txt")
        (OUT/"gaia_profile_README.txt").write_bytes(readme)
        for cid in IDS:
            names=[n for n in z.namelist() if re.match(rf"profiles/NGC_{int(cid[3:])}(?:_|\.txt)",n)]
            assert len(names)==1,(cid,names)
            name=names[0]; content=z.read(name)
            members.append(dict(cluster=cid,member=name,bytes=len(content),
                sha256=hashlib.sha256(content).hexdigest()))
            a=np.loadtxt(StringIO(content.decode()))
            assert a.shape==(101,12)
            for i,row in enumerate(a):
                rows.append(dict(cluster=cid,member=name,grid_row=i,
                    radius_degree=row[0],radius_arcsec=row[0]*3600,
                    sigma_p023_masyr=row[1],sigma_p159_masyr=row[2],sigma_median_masyr=row[3],
                    sigma_p841_masyr=row[4],sigma_p977_masyr=row[5],
                    rotation_p023_masyr=row[6],rotation_p159_masyr=row[7],
                    rotation_median_masyr=row[8],rotation_p841_masyr=row[9],
                    rotation_p977_masyr=row[10],undocumented_column_12=row[11],
                    record_type="correlated_fitted_curve_not_independent_bin",
                    covariance_available=False,likelihood_ready=False))
    return pd.DataFrame(rows),members



def gaia_stars():
    """Descriptive source audit only; does not estimate intrinsic dispersions."""
    rows=[]
    columns=["source_id","ra","dec","x","y","plx","pmra","pmdec","plxe",
             "pmrae","pmdece","pmcorr","g_mag","bp_rp","Sigma","qflag","memberprob"]
    with zipfile.ZipFile(RAW/"gaia_edr3_v2.zip") as z:
        for cid in IDS:
            names=[n for n in z.namelist() if re.match(rf"catalogues/NGC_{int(cid[3:])}(?:_|\.txt)",n)]
            assert len(names)==1
            content=z.read(names[0])
            df=pd.read_csv(StringIO(content.decode()),sep=r"\s+",skiprows=1,names=columns,
                           dtype={"source_id":str})
            assert df.shape[1]==17
            assert isinstance(df.index,pd.RangeIndex), "Unexpected extra catalog fields"
            assert df.source_id.str.fullmatch(r"\d+").all(), "Invalid Gaia source identifier"
            assert df.memberprob.between(0,1).all()
            assert df.pmcorr.between(-1,1).all()
            assert df.qflag.isin([0,1,2,3]).all()
            assert (df[["pmrae","pmdece","plxe"]]>0).all().all()
            quality=(df.qflag.astype(int)&2)!=0
            selected=quality & (df.memberprob>=0.9)
            r_arcmin=np.hypot(df.x,df.y)*60
            rows.append(dict(cluster=cid,archive_member=names[0],
                member_sha256=hashlib.sha256(content).hexdigest(),rows=len(df),
                duplicate_source_ids=int(df.source_id.duplicated().sum()),
                astrometric_quality_pass=int(quality.sum()),
                high_probability_quality_diagnostic=int(selected.sum()),
                diagnostic_selection="(qflag & 2) != 0 and memberprob >=0.9; not adopted fit selection",
                missing_bp_rp=int(df.bp_rp.isna().sum()),
                missing_g_mag=int(df.g_mag.isna().sum()),
                min_projected_radius_arcmin=float(r_arcmin[selected].min()),
                max_projected_radius_arcmin=float(r_arcmin[selected].max()),
                median_pmra_measurement_error_masyr=float(df.loc[selected,"pmrae"].median()),
                likelihood_ready=False))
    return pd.DataFrame(rows)


def weighted_quantile(x,w,qs=(0.16,0.5,0.84)):
    order=np.argsort(x)
    return np.interp(qs,np.cumsum(w[order]),x[order])


def posterior_products():
    stats,diagnostics,inventories=[],[],[]
    for cid in IDS:
        summary=(RAW/f"gcfit_{cid}_summary.txt").read_text()
        reported={}
        for line in summary.splitlines():
            m=re.match(r"\s*(\w+)\s*=\s*([-.\d]+)\s*\(±([-.\d]+)\)\s*\|\s*([-.\d]+)\s*\(±([-.\d]+)\)",line)
            if m:reported[m[1]]=tuple(float(m[k]) for k in (2,4,3,5))
        assert list(reported)==PARAMS,"Parameter order changed; verify authoritative model source."
        with h5py.File(RAW/f"gcfit_{cid}_sampler.hdf") as f:
            assert f["metadata"].attrs["cluster"]==cid
            samples=f["nested/samples"][:];logw=f["nested/logwt"][:]
            assert samples.shape==(len(logw),13) and np.isfinite(samples).all()
            assert np.isfinite(logw).all()
            w=np.exp(logw-np.max(logw));w/=w.sum()
            mean=np.sum(samples*w[:,None],axis=0)
            centered=samples-mean
            cov=(centered*w[:,None]).T@centered
            sd=np.sqrt(np.diag(cov));corr=cov/np.outer(sd,sd)
            ess=1/np.sum(w*w)
            assert abs(ess-float(re.search(r"ESS\s*=\s*([\d.]+)",summary).group(1)))<=0.0051
            mismatches=[]
            for j,param in enumerate(PARAMS):
                lo,med,hi=weighted_quantile(samples[:,j],w)
                pubmean,pubsd,pubmeanerr,pubsderr=reported[param]
                # Summaries have three decimals; no extra tolerance hides mismatches.
                if abs(mean[j]-pubmean)>0.00051 or abs(sd[j]-pubsd)>0.00051:
                    mismatches.append(param)
                factor=1e6 if param=="M" else 1.0
                units={"M":"Msun","rh":"pc","d":"kpc","ra":"log10 dimensionless anisotropy radius",
                       "BHret":"percent model retention parameter; not BH mass fraction"}
                stats.append(dict(cluster=cid,parameter=param,unit=units.get(param,"model convention"),
                    mean=mean[j]*factor,standard_deviation=sd[j]*factor,
                    p16=lo*factor,median=med*factor,p84=hi*factor,
                    provenance="published standard-fit posterior; not an adopted host",
                    reported_mean=pubmean*factor,reported_standard_deviation=pubsd*factor,
                    reported_mean_numerical_uncertainty=pubmeanerr*factor,
                    reported_standard_deviation_numerical_uncertainty=pubsderr*factor))
            meta=dict(f["metadata"].attrs)
            meta.pop("GCFIT_DIR",None)
            priors=dict(f["metadata/specified_priors"].attrs)
            di=dict(cluster=cid,rows=len(samples),parameters=PARAMS,posterior_ess=ess,
                reported_ess=float(re.search(r"ESS\s*=\s*([\d.]+)",summary).group(1)),
                sum_posterior_weights=float(w.sum()),max_weight=float(w.max()),
                summary_rounding_differences=mismatches,metadata=meta,
                specified_priors=priors,
                fixed_parameters=dict(f["metadata/fixed_params"].attrs),
                excluded_likelihoods=dict(f["metadata/excluded_likelihoods"].attrs),
                input_hash_in_chain=False,software_revision_in_chain=False,
                exact_model_data_pairing="unverified",
                moment_estimator="direct_weighted_empirical_posterior",
                summary_estimator="historical GCfit average of resampled-run moments; three-decimal output",
                posterior_covariance_native=cov.tolist(),posterior_correlation=corr.tolist(),
                covariance_parameter_units=["dimensionless","1e6 Msun","pc","log10 dimensionless ra",
                    "dimensionless","dimensionless","number-density likelihood units",
                    "dimensionless","IMF slope","IMF slope","IMF slope","percent","kpc"],
                covariance_scale_to_parameter_csv=[1,1e6,1,1,1,1,1,1,1,1,1,1,1],
                intended_sampling="preserve joint weighted rows; no independent Gaussians",
                standard_fit_is_final_selected_model=False if cid=="NGC6397" else None)
            diagnostics.append(di)
            inventory=[]
            def visit(n,obj):
                if isinstance(obj,h5py.Dataset):
                    inventory.append(dict(path=n,shape=obj.shape,attributes=dict(obj.attrs)))
            f["model_output"].visititems(visit)
            inventories.append(dict(cluster=cid,model_products=inventory))
    return pd.DataFrame(stats),diagnostics,inventories


def observations_inventory():
    out=[]
    for cid in IDS:
        with h5py.File(RAW/f"gcfit_{cid}_observations.hdf") as f:
            groups=[]
            def visit(n,obj):
                if not isinstance(obj,h5py.Group):return
                datasets=[dict(name=k,shape=v.shape,unit=v.attrs.get("unit")) for k,v in obj.items()
                          if isinstance(v,h5py.Dataset)]
                if datasets:groups.append(dict(path=n,datasets=datasets,source=obj.attrs.get("source"),
                                               notes=obj.attrs.get("notes")))
            f.visititems(visit)
            out.append(dict(cluster=cid,groups=groups,compatibility="candidate only; private run inputs not identified"))
    return out


def analyze():
    OUT.mkdir(parents=True,exist_ok=True)
    sources=json.loads((DATA/"sources.json").read_text())
    for source in sources["files"]:
        assert hashlib.sha256((RAW/source["file"]).read_bytes()).hexdigest()==source["sha256"]
    cats=catalog("uq_parameters.html");old=catalog("uq_root_parameters.html")
    hst=hst_profiles()
    values=cats.set_index("cluster")
    for col in ["distance_kpc","distance_error_kpc","half_mass_3d_pc","scale_v_kms"]:
        hst[col]=hst.cluster.map(values[col])
    hst["radius_pc_at_catalog_distance"]=hst.radius_arcsec/ARCSEC_PER_RAD*1000*hst.distance_kpc
    hst["sigma_kms_at_catalog_distance"]=PM_FACTOR*hst.distance_kpc*hst.sigma_masyr
    for side in ["upper","lower"]:
        hst[f"measurement_error_{side}_kms"]=PM_FACTOR*hst.distance_kpc*hst[f"error_{side}_masyr"]
    hst["d_sigma_kms_d_distance_kpc"]=PM_FACTOR*hst.sigma_masyr
    hst["d_radius_pc_d_distance_kpc"]=hst.radius_arcsec/ARCSEC_PER_RAD*1000
    hst["projected_radius_over_rh"]=hst.radius_pc_at_catalog_distance/hst.half_mass_3d_pc
    hst["sigma_over_scale_v"]=hst.sigma_kms_at_catalog_distance/hst.scale_v_kms
    hst["distance_basis"]="UQ newdata point estimate; illustrative conversion, shared uncertainty retained"
    assert np.allclose(hst.sigma_kms_at_catalog_distance/PM_FACTOR/hst.distance_kpc,hst.sigma_masyr)
    assert np.allclose(hst.radius_pc_at_catalog_distance/1000/hst.distance_kpc*ARCSEC_PER_RAD,hst.radius_arcsec)
    compiled,surface,flags,dupes=compiled_profiles()
    gaia,members=gaia_curves()
    stars=gaia_stars()
    stars.to_csv(OUT/"gaia_star_catalog_quality.csv",index=False,float_format="%.12g")
    stats,posterior,modelinv=posterior_products()
    for name,frame in [("catalog_scalars.csv",cats),("root_catalog_comparison.csv",old),
        ("hst_native_and_common_units.csv",hst),("compiled_kinematics_audit.csv",compiled),
        ("surface_profiles.csv",surface),("gaia_correlated_curves.csv",gaia),
        ("posterior_parameters.csv",stats)]:
        frame.to_csv(OUT/name,index=False,float_format="%.12g")
    write_json("posterior_diagnostics.json",posterior)
    write_json("model_product_inventory.json",modelinv)
    write_json("observation_bundle_inventory.json",observations_inventory())
    write_json("gaia_archive_members.json",members)
    inventory=[]
    for cid in IDS:
        h=hst[(hst.cluster==cid)&(hst.component=="PM_COMBINED")]
        c=compiled[compiled.cluster==cid]
        inventory.append(dict(cluster=cid,hst_radial_bins=len(h),
            hst_scalar_components=len(h)*3,compiled_bins=len(c),
            compiled_quarantined_rows=int(c.quarantine.sum()),
            surface_bins=int((surface.cluster==cid).sum()),
            gaia_correlated_grid_points=int((gaia.cluster==cid).sum()),
            hst_R_over_rh_min=h.projected_radius_over_rh.min(),
            hst_R_over_rh_max=h.projected_radius_over_rh.max(),
            hst_native_radius_min_arcsec=h.radius_arcsec.min(),
            hst_native_radius_max_arcsec=h.radius_arcsec.max()))
    inv=pd.DataFrame(inventory);inv.to_csv(OUT/"coverage.csv",index=False,float_format="%.12g")
    quality=dict(source_files=len(sources["files"]),source_bytes=sum(x["bytes"] for x in sources["files"]),
        clusters=6,gaia_stars_inspected=int(stars.rows.sum()),hst_unique_bins=len(hst)//3,hst_component_records=len(hst),
        compiled_kinematic_rows=len(compiled),surface_rows=len(surface),
        gaia_curve_samples=len(gaia),compiled_flagged_rows=len(flags),
        compiled_duplicate_candidate_rows=len(dupes),
        duplicate_groups=dupes[["cluster","radius_arcsec","component","literature_tag"]].drop_duplicates().to_dict("records"),
        flags=flags,intersection_hst_median_radius_spans_R_over_rh=[inv.hst_R_over_rh_min.max(),inv.hst_R_over_rh_max.min()],
        posterior_summary_rounding_differences={d["cluster"]:d["summary_rounding_differences"] for d in posterior},
        common_units_ready=True,common_tracer_selection_ready=False,
        full_local_species_distribution_ready=False,dynamics_enabled=False,
        checks=["source SHA256","Gaia provider MD5 during retrieval","exact source table widths",
                "six unique canonical IDs","nonnegative HST dispersions/errors","HST component key uniqueness",
                "native/physical unit round trips","nested sample/weight alignment",
                "posterior normalized weights"],
        comparisons_performed=["all 13 means/SDs vs published summary rounding; resampled-summary estimator differs; rounding differences recorded"])
    write_json("quality.json",quality)
    print(json.dumps({k:v for k,v in quality.items() if k not in ("flags","checks","duplicate_groups")},indent=2))
    return quality


if __name__=="__main__":
    analyze()
