"""Data-audit figures; plots existing observations and published posteriors."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"data/observed-clusters/derived"
IDS=["NGC0104","NGC1851","NGC3201","NGC5139","NGC6121","NGC6397"]
NAMES=["47 Tuc","NGC 1851","NGC 3201","Omega Centauri","M4","NGC 6397"]


def save_figure(fig, stem):
    fig.savefig(OUT/f"{stem}.png")
    # Exclude wall-clock metadata from the reproducible vector artifact.
    fig.savefig(OUT/f"{stem}.pdf", metadata={"CreationDate": None, "ModDate": None})


def create():
    plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,
                         "axes.titleweight":"bold","savefig.dpi":160})
    hst=pd.read_csv(OUT/"hst_native_and_common_units.csv")
    cats=pd.read_csv(OUT/"catalog_scalars.csv").set_index("cluster").loc[IDS]
    stats=pd.read_csv(OUT/"posterior_parameters.csv")
    fig,axes=plt.subplots(2,3,figsize=(12,8.2),sharex=True,sharey=True)
    for ax,cid,name in zip(axes.flat,IDS,NAMES):
        sub=hst[hst.cluster==cid]
        for component,color,marker,label in [
            ("PM_R","#245b82","o","Radial PM"),("PM_T","#b17618","s","Tangential PM")]:
            d=sub[sub.component==component]
            errors=np.stack([d.measurement_error_lower_kms,d.measurement_error_upper_kms])/d.scale_v_kms.to_numpy()
            ax.errorbar(d.projected_radius_over_rh,d.sigma_over_scale_v,yerr=errors,
                        color=color,fmt=marker,ms=3,capsize=1.5,elinewidth=.7,label=label,alpha=.9)
        ax.set_title(f"{name} · {sub.source.iloc[0]}",fontsize=11)
        ax.set_xscale("log");ax.set_xlim(.003,2);ax.set_ylim(0,.95)
        ax.grid(alpha=.14,which="major")
    axes[0,0].legend(frameon=False,fontsize=9,loc="lower left")
    for ax in axes[-1]:ax.set_xlabel(r"Projected $R/r_h$")
    for ax in axes[:,0]:ax.set_ylabel(r"$\sigma_{\rm PM}/\sqrt{GM/r_h}$")
    fig.suptitle("HST kinematics in common comparison units",y=.985,fontsize=16)
    fig.text(.07,.03,"Catalog distance, mass and radius are point scales, not fitted here. Error bars show conditional measurement errors only.\nHACKS uses bright SGB/RGB tracers; H25 uses the full magnitude range. No common tracer distribution is assumed.",fontsize=10)
    fig.tight_layout(rect=[0,.10,1,.95])
    save_figure(fig, "hst_common_scales")
    plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(12,5.4),sharey=True)
    y=np.arange(6)
    specs=[("M","mass_msun","mass_error_msun","Mass (solar masses)"),
           ("rh","half_mass_3d_pc",None,"3D half-mass radius (pc)"),
           ("d","distance_kpc","distance_error_kpc","Distance (kpc)")]
    for ax,(param,col,errcol,label) in zip(axes,specs):
        d=stats[stats.parameter==param].set_index("cluster").loc[IDS]
        ax.errorbar(d["median"],y-.12,xerr=np.stack([d["median"]-d.p16,d.p84-d["median"]]),
                    fmt="o",color="#245b82",ms=5,capsize=3,label="GCfit marginal: 16–84%")
        ax.errorbar(cats[col],y+.12,xerr=cats[errcol] if errcol else None,fmt="s",
                    mfc="white",mec="#444444",ecolor="#777777",ms=5,capsize=3,label="UQ catalog estimate")
        ax.set_xlabel(label);ax.grid(axis="x",alpha=.16)
        if param=="M":ax.set_xscale("log")
    axes[0].set_yticks(y,NAMES);axes[0].invert_yaxis()
    fig.suptitle("Published host estimates depend on the model and data pairing",fontsize=15,y=.99)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.56,.93),ncol=2,frameon=False,fontsize=10)
    fig.text(.08,.025,"GCfit standard-fit chains are not adopted hosts; the NGC 6397 no-BH alternative lacks a released chain.\nCatalog radius errors are unavailable. Differences are model/data comparisons, not temporal evolution or independent detections.",fontsize=10)
    fig.tight_layout(rect=[0,.12,1,.86])
    save_figure(fig, "host_estimate_comparison")
    plt.close(fig)

    cov=pd.read_csv(OUT/"coverage.csv").set_index("cluster").loc[IDS]
    fig,ax=plt.subplots(figsize=(9,4.7))
    for i,cid in enumerate(IDS):
        lo,hi=cov.loc[cid,["hst_R_over_rh_min","hst_R_over_rh_max"]]
        ax.plot([lo,hi],[i,i],"-o",color="#245b82",lw=2,ms=5)
        ax.annotate(f"{int(cov.loc[cid,'hst_radial_bins'])} bins",(hi,i),xytext=(7,0),
                    textcoords="offset points",va="center",fontsize=10)
    lo=cov.hst_R_over_rh_min.max();hi=cov.hst_R_over_rh_max.min()
    ax.axvline(lo,color="#555555",ls="--",lw=1)
    ax.axvline(hi,color="#555555",ls="--",lw=1)
    ax.set_yticks(y,NAMES);ax.invert_yaxis();ax.set_xscale("log");ax.set_xlim(.003,3)
    ax.set_xlabel(r"Projected median-bin radius $R/r_h$");ax.grid(axis="x",alpha=.15)
    ax.set_title("Radial sampling differs across the six clusters",pad=14)
    fig.text(.11,.02,f"Dashed bounds: intersection of median-radius spans, {lo:.3f}–{hi:.3f}.\nThese are not annular edges; HACKS bin edges are not supplied. No interpolation or shared likelihood is constructed.",fontsize=9)
    fig.tight_layout(rect=[0,.12,1,1])
    save_figure(fig, "radial_sampling")
    plt.close(fig)
    return ["hst_common_scales","host_estimate_comparison","radial_sampling"]


if __name__=="__main__":
    print("Rendered data-audit figures:",", ".join(create()))
