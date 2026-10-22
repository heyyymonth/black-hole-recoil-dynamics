"""Render vector manuscript figures from the declared orbital experiments.

Run from any directory with the repository environment::

    .venv/bin/python paper/figures.py

This changes presentation only: the existing model, quadrature tolerances,
parameter grids, and comparator definitions are taken from the experiment
configuration. The original result tables and figures are left untouched.
"""
from pathlib import Path
import csv
import json
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from recoil_dynamics import PlummerOrbit, core_exit_threshold


OUTPUT = ROOT / "paper" / "figures"


def orbit(eta, config):
    return PlummerOrbit(float(eta), epsabs=config["epsabs"],
                        epsrel=config["epsrel"])


def save(fig, filename):
    """Retain vector text/curves, fixed page dimensions, and stable metadata."""
    fig.savefig(OUTPUT / filename, metadata={
        "Creator": "Black-hole Recoil Dynamics manuscript figures",
        "CreationDate": None,
        "ModDate": None,
    })
    plt.close(fig)


def return_time(config):
    cfg = config["e001"]
    eta = np.unique(np.r_[
        np.linspace(.001, .95, cfg["plot_points"]),
        1 - np.geomspace(.05, 1 - cfg["eta_plot_max"], cfg["plot_points"]),
    ])
    periods = [orbit(e, config).period for e in eta]
    fig, ax = plt.subplots(figsize=(3.4, 2.75), layout="constrained")
    ax.semilogy(eta, periods, color="black", label="Orbit quadrature")
    mask = eta >= .8
    ax.semilogy(eta[mask], math.pi / math.sqrt(2) *
                (1 - eta[mask] ** 2) ** (-1.5), color="0.35", ls="--",
                label="Near-escape limit")
    ax.axhline(math.pi, color="0.5", ls=":", label="Small-kick limit")
    ax.set(xlabel=r"Kick ratio $\eta = v_k/v_{\rm esc}(0)$",
           ylabel=r"First return (units of $t_{\rm d}$)", xlim=(0, 1))
    ax.legend(loc="upper left", frameon=False)
    save(fig, "return-time.pdf")


def core_residence(config):
    cfg = config["e002"]
    thresholds = np.array([core_exit_threshold(c)
                           for c in cfg["core_radii_over_b"]])
    eta = np.unique(np.r_[
        np.linspace(.001, .95, cfg["plot_points"]),
        1 - np.geomspace(.05, .001, cfg["plot_points"]),
        thresholds * (1 - 1e-5), thresholds, thresholds * (1 + 1e-5),
    ])
    fig, ax = plt.subplots(figsize=(3.4, 2.75), layout="constrained")
    for c, color, style in zip(cfg["core_radii_over_b"],
                               ["black", "0.3", "0.5"], ["-", "--", ":"]):
        values = []
        for e in eta:
            o = orbit(e, config)
            values.append(o.cycle_core_time(c) / o.period)
        ax.semilogy(eta, values, color=color, ls=style, label=rf"$c={c:g}$")
    ax.set(xlabel=r"Kick ratio $\eta = v_k/v_{\rm esc}(0)$",
           ylabel=r"Core residence $C(P,c)/P$", xlim=(0, 1),
           ylim=(1e-6, 1.6))
    ax.legend(loc="lower left", frameon=False)
    save(fig, "core-residence.pdf")


def exposure_comparison(config):
    cfg = config["e002"]
    horizons = np.unique(np.r_[
        np.geomspace(1e-6, .05, 80),
        np.linspace(.05, cfg["plot_horizon_max_periods"], cfg["plot_points"]),
    ])
    fig, axes = plt.subplots(1, 3, figsize=(7, 2.8), sharey=True)
    fig.subplots_adjust(left=.075, right=.99, bottom=.20, top=.80, wspace=.13)
    for ax, eta, panel in zip(axes, cfg["plot_eta"], "abc"):
        o = orbit(eta, config)
        means = np.array([o.exposure(float(t * o.period)) / (t * o.period)
                          for t in horizons])
        ax.plot(horizons, means, color="black", label="Orbit")
        ax.plot(horizons, np.maximum(horizons - 1, 0) / horizons,
                color="0.3", ls="--", label="Wait, then central")
        ax.axhline(1, color="0.5", ls=":", label="Always central")
        ax.text(.04, .86, rf"({panel}) $\eta={eta:g}$", transform=ax.transAxes)
        ax.set(xlabel=r"Horizon $H/P$", ylim=(-.035, 1.06),
               xlim=(0, cfg["plot_horizon_max_periods"]))
        ax.set_xticks([0, 5, 10])
    axes[0].set_ylabel(r"Mean exposure $I(H)/H$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(.53, .99), columnspacing=1.8, handlelength=2.6)
    save(fig, "exposure-comparison.pdf")


def encounter_proxy():
    """Present E003's committed outputs; do not rerun its rate calculation."""
    rows = np.genfromtxt(ROOT / "results/e003/encounters.csv", delimiter=",",
                         names=True)
    summary = json.loads((ROOT / "results/e003/summary.json").read_text())
    positions = np.arange(len(rows))
    labels = [f"{eta:g}" for eta in rows["eta"]]
    labels[0] = "0\n(control)"
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.15))
    fig.subplots_adjust(left=.095, right=.985, bottom=.22, top=.84,
                        wspace=.31)
    left, right = axes
    left.plot(positions - .065, rows["density_ratio"], "o", color="black",
              ms=5, label=r"Density proxy $D=I(H)/H$")
    left.plot(positions + .065, rows["encounter_ratio"], "s", color="0.35",
              markerfacecolor="white", ms=5,
              label=r"Passage count $Q=\tau(\eta,H)/\tau(0,H)$")
    left.set(yscale="log", ylim=(.02, 2.5),
             ylabel="Ratio to stationary control")
    left.set_yticks([.03, .1, .3, 1], ["0.03", "0.1", "0.3", "1"])
    left.legend(loc="upper right", fontsize=8.5, frameon=False,
                handlelength=1, handletextpad=.5)

    errors = 100 * rows["proxy_relative_error"]
    right.bar(positions, errors, width=.52, facecolor="0.82",
              edgecolor="black", linewidth=.8)
    tolerance = 100 * summary["proxy_relative_error_tolerance"]
    right.axhline(tolerance, color="0.3", ls="--", linewidth=1,
                  label=f"Declared {tolerance:g}% screen")
    for x, error in zip(positions, errors):
        right.text(x, error + 2.5, f"{error:.1f}%", ha="center", va="bottom",
                   fontsize=9)
    right.set(ylabel=r"Density-proxy error (%)",
              ylim=(-2, 119), yticks=[0, 25, 50, 75, 100])
    right.legend(loc="upper left", fontsize=8.5, frameon=False)
    for ax, panel in zip(axes, "ab"):
        ax.set_xticks(positions, labels)
        ax.set(xlabel=r"Kick ratio $\eta$", xlim=(-.48, 3.48))
        ax.text(.02, 1.055, f"({panel})", transform=ax.transAxes,
                 ha="left", va="bottom")
    fig.text(.52, .985, r"Common window $T=10\,t_{\rm d}$; imposed Plummer tracers",
             ha="center", va="top", fontsize=10)
    save(fig, "encounter-proxy.pdf")


def snapshot_comparison():
    """Present the saved E004 population, potential, and ballistic results."""
    profiles = np.genfromtxt(ROOT / "results/e004/profiles.csv", delimiter=",",
                             names=True)
    trajectories = np.genfromtxt(ROOT / "results/e004/trajectories.csv",
                                 delimiter=",", names=True)
    summary = json.loads((ROOT / "results/e004/summary.json").read_text())
    finest = trajectories[trajectories["shell_count"] ==
                          trajectories["shell_count"].max()]
    radius = profiles["radius_over_rh"]
    fig, axes = plt.subplots(1, 3, figsize=(7, 3.15))
    fig.subplots_adjust(left=.067, right=.99, bottom=.23, top=.85,
                        wspace=.47)
    a, b, c = axes
    a.plot(radius, profiles["host_mass_cdf"], color="black",
            label="Host mass")
    a.plot(radius, profiles["plummer_mass_cdf"], color="0.4", ls="--",
            label="Matched Plummer")
    # This is the discrete CDF sampled on the committed plotting grid. A
    # step preserves the small-sample structure; no continuous fit is made.
    a.step(radius, profiles["bbh_number_cdf"], where="post", color="0.3",
            ls=":", linewidth=1.7,
            label=rf"BH--BH number ($N={summary['binary_black_hole_count']}$)")
    a.set(xscale="log", xlim=(radius[0], radius[-1]), ylim=(-.025, 1.045),
           xlabel=r"Radius $r/r_{\rm h}$", ylabel="Cumulative fraction")
    a.legend(loc="lower right", fontsize=8, frameon=True,
              facecolor="white", edgecolor="none", framealpha=1,
              handlelength=1.6, handletextpad=.4, borderaxespad=.2)

    b.plot(radius, profiles["snapshot_potential_over_GM_b"], color="black",
            label="Frozen snapshot")
    b.plot(radius, profiles["plummer_potential_over_GM_b"], color="0.4",
            ls="--", label="Matched Plummer")
    b.set(xscale="log", xlim=(radius[0], radius[-1]), ylim=(-1.43, .015),
           xlabel=r"Radius $r/r_{\rm h}$", ylabel=r"Potential $\Phi/(GM/b)$")
    b.legend(loc="upper left", fontsize=8, frameon=False,
              handlelength=1.6, handletextpad=.4, borderaxespad=.2)

    kicks = finest["kick_kms"]
    # Discrete points at the same physical kicks, not separately normalized
    # escape fractions. No line interpolation or error bars are implied.
    c.plot(kicks, finest["snapshot_return_myr"], "o", color="black", ms=4.5,
            label="Frozen snapshot")
    c.plot(kicks, finest["plummer_return_myr"], "s", color="0.4", ms=4.5,
            markerfacecolor="white", label="Matched Plummer")
    c.set(xlabel=r"Same physical kick $v_k$ (km s$^{-1}$)",
           ylabel="First return (Myr)", ylim=(0, 7.1),
           xlim=(kicks.min()-1.6, kicks.max()+1.6))
    c.set_xticks(kicks, [f"{v:.2f}" for v in kicks], rotation=45, ha="right")
    c.legend(loc="upper left", fontsize=8, frameon=False,
              handlelength=1, handletextpad=.4, borderaxespad=.2)
    for ax, panel in zip(axes, "abc"):
        ax.tick_params(axis="both", labelsize=8.5)
        ax.xaxis.label.set_size(9)
        ax.yaxis.label.set_size(9)
        ax.text(.02, 1.045, f"({panel})", transform=ax.transAxes,
                 ha="left", va="bottom")
    for ax in (a, b):
        ax.set_xticks([.001, .01, .1, 1, 10])
        ax.tick_params(axis="x", which="minor", bottom=False)
    fig.text(.52, .985, "One evolved CMC snapshot: population and frozen potential",
             ha="center", va="top", fontsize=10)
    save(fig, "snapshot-comparison.pdf")


def catalog_counts():
    """Display the committed common-age census without smoothing or fits."""
    with (ROOT / 'results/e005/selected.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    summary = json.loads((ROOT / 'results/e005/summary.json').read_text())
    bbh = np.array([int(row['bbh_count']) for row in rows])
    bh = np.array([int(row['bh_count']) for row in rows])
    if len(rows) != summary['selected_models']:
        raise ValueError('Selected catalog rows disagree with the saved census')
    if len({row['model'] for row in rows}) != len(rows):
        raise ValueError('The common-age comparison must contain one state per model')
    bins = np.arange(bbh.max() + 1)
    frequency = np.bincount(bbh, minlength=len(bins))
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.0))
    fig.subplots_adjust(left=.075, right=.97, bottom=.20, top=.84, wspace=.29)
    left, right = axes
    left.bar(bins, frequency, width=.78, facecolor='0.82', edgecolor='black',
             linewidth=.8)
    for count, number in zip(bins, frequency):
        if number:
            left.text(count, number + .7, str(number), ha='center', va='bottom',
                      fontsize=8.5)
    left.set(xlabel=r'Instantaneous BBH count $N_{\rm BBH}$',
             ylabel='Number of selected models', xlim=(-.6, bins[-1]+.6),
             ylim=(0, frequency.max()*1.16), xticks=bins)
    # Exact integer coordinates are retained. Coincident states overlap; the
    # left panel accounts for every selected state without arbitrary jitter.
    right.scatter(bh, bbh, s=18, facecolor='none', edgecolor='black',
                  linewidth=.7, alpha=.65)
    right.set(xlabel=r'Total BH-object count $N_{\rm BH}$',
              ylabel=r'Instantaneous BBH count $N_{\rm BBH}$',
              xlim=(-50, 2050), ylim=(-.45, bbh.max()+.65),
              xticks=[0, 500, 1000, 1500, 2000], yticks=np.arange(0, bbh.max()+1, 2))
    for ax, panel in zip(axes, 'ab'):
        ax.tick_params(axis='both', labelsize=9)
        ax.text(.02, 1.045, f'({panel})', transform=ax.transAxes,
                ha='left', va='bottom')
    fig.text(.52, .985,
             f"One state per available model near 12 Gyr ($N={len(rows)}$)",
             ha='center', va='top', fontsize=10)
    save(fig, 'catalog-counts.pdf')


def main():
    config = json.loads((ROOT / "experiments" / "config.json").read_text())
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "STIXGeneral",
        "mathtext.fontset": "stix",
        "font.size": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "lines.linewidth": 1.2,
        "axes.linewidth": .6,
        "xtick.major.width": .6,
        "ytick.major.width": .6,
        "xtick.minor.width": .5,
        "ytick.minor.width": .5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    })
    return_time(config)
    core_residence(config)
    exposure_comparison(config)
    encounter_proxy()
    snapshot_comparison()
    catalog_counts()
    print("Rendered six vector manuscript figures in paper/figures/.")


if __name__ == "__main__":
    main()
