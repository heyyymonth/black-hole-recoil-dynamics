"""Publication figures from recorded cluster comparisons; no host evaluation."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

IDS = ["NGC0104", "NGC1851", "NGC3201", "NGC5139", "NGC6121", "NGC6397"]
NAMES = ["47 Tuc", "NGC 1851", "NGC 3201", r"$\omega$ Cen", "M4", "NGC 6397"]


def render(data: Path, output: Path):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "axes.labelsize": 8, "xtick.labelsize": 7.5,
                         "ytick.labelsize": 7.5, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    hst = pd.read_csv(data / "hst_native_and_common_units.csv")
    cats = pd.read_csv(data / "catalog_scalars.csv").set_index("cluster").loc[IDS]
    stats = pd.read_csv(data / "posterior_parameters.csv")
    coverage = pd.read_csv(data / "coverage.csv").set_index("cluster").loc[IDS]
    output.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        fig.savefig(output / name, metadata={"CreationDate": None, "ModDate": None})
        plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.5), sharex=True, sharey=True)
    for ax, cid, name in zip(axes.flat, IDS, NAMES):
        sub = hst[hst.cluster == cid]
        for component, color, marker, label in [
                ("PM_R", "#245b82", "o", "Radial PM"),
                ("PM_T", "#b17618", "s", "Tangential PM")]:
            d = sub[sub.component == component]
            err = np.stack([d.measurement_error_lower_kms,
                            d.measurement_error_upper_kms]) / d.scale_v_kms.to_numpy()
            ax.errorbar(d.projected_radius_over_rh, d.sigma_over_scale_v, yerr=err,
                        color=color, fmt=marker, ms=2.5, capsize=1.2, elinewidth=.6, label=label)
        ax.set_title(name, fontsize=9)
        ax.set_xscale("log")
        ax.set_xlim(.003, 2)
        ax.set_ylim(0, .95)
        ax.grid(alpha=.15)
    axes[0, 0].legend(frameon=False, fontsize=7, loc="lower left")
    for ax in axes[-1]:
        ax.set_xlabel(r"Projected $R/r_h$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\sigma_{\rm PM}/\sqrt{GM/r_h}$")
    fig.tight_layout(pad=.7, h_pad=1.0, w_pad=.8)
    save(fig, "cluster-hst.pdf")

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.2), sharey=True)
    y = np.arange(6)
    for ax, (param, column, errcolumn, label) in zip(axes, [
            ("M", "mass_msun", "mass_error_msun", r"Mass ($M_\odot$)"),
            ("rh", "half_mass_3d_pc", None, r"$r_h$ (pc)"),
            ("d", "distance_kpc", "distance_error_kpc", "Distance (kpc)")]):
        d = stats[stats.parameter == param].set_index("cluster").loc[IDS]
        ax.errorbar(d["median"], y-.12, xerr=np.stack([d["median"]-d.p16, d.p84-d["median"]]),
                    fmt="o", color="#245b82", ms=3.5, capsize=2, elinewidth=.8,
                    label="Multimass fit: 16–84% marginals")
        ax.errorbar(cats[column], y+.12, xerr=cats[errcolumn] if errcolumn else None,
                    fmt="s", mfc="white", mec="#444444", ecolor="#777777", ms=3.5,
                    capsize=2, label="Structural catalog")
        ax.set_xlabel(label)
        ax.grid(axis="x", alpha=.15)
        if param == "M":
            ax.set_xscale("log")
    axes[0].set_yticks(y, NAMES)
    axes[0].invert_yaxis()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=[0, 0, 1, .91], pad=.7, w_pad=1.0)
    save(fig, "cluster-host-comparison.pdf")

    fig, ax = plt.subplots(figsize=(3.4, 2.45))
    for i, cid in enumerate(IDS):
        lo, hi = coverage.loc[cid, ["hst_R_over_rh_min", "hst_R_over_rh_max"]]
        ax.plot([lo, hi], [i, i], "-o", color="#245b82", lw=1.4, ms=3)
        ax.annotate(str(int(coverage.loc[cid, "hst_radial_bins"])), (hi, i), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=7)
    for bound in [coverage.hst_R_over_rh_min.max(), coverage.hst_R_over_rh_max.min()]:
        ax.axvline(bound, color="#555555", ls="--", lw=.7)
    ax.set_yticks(y, NAMES)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(.003, 4)
    ax.set_xlabel(r"Projected median-bin radius $R/r_h$")
    ax.grid(axis="x", alpha=.15)
    fig.tight_layout(pad=.6)
    save(fig, "cluster-coverage.pdf")
