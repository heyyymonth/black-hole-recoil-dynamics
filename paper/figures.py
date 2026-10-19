"""Render vector manuscript figures from the declared orbital experiments.

Run from any directory with the repository environment::

    .venv/bin/python paper/figures.py

This changes presentation only: the existing model, quadrature tolerances,
parameter grids, and comparator definitions are taken from the experiment
configuration. The original result tables and figures are left untouched.
"""
from pathlib import Path
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
    print("Rendered three vector manuscript figures in paper/figures/.")


if __name__ == "__main__":
    main()
