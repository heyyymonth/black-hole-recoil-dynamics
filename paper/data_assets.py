"""Prepare or verify paper assets from the recorded cluster-data revision."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/observed-clusters/derived"
PAPER = ROOT / "paper"
MANIFEST = PAPER / "data-assets.json"
IDS = ["NGC0104", "NGC1851", "NGC3201", "NGC5139", "NGC6121", "NGC6397"]
NAMES = dict(zip(IDS, ["47 Tuc", "NGC 1851", "NGC 3201", r"$\omega$ Cen", "M4", "NGC 6397"]))
SOURCES = ["coverage.csv", "catalog_scalars.csv", "posterior_parameters.csv",
           "posterior_diagnostics.json", "hst_native_and_common_units.csv",
           "gaia_star_catalog_quality.csv"]
INPUTS = ["paper/data_assets.py", "paper/data_figures.py",
          "analysis/cluster_data/requirements-lock.txt", "data/observed-clusters/audit.json",
          *["data/observed-clusters/derived/"+name for name in SOURCES]]
OUTPUTS = ["tables/cluster-coverage.tex", "tables/cluster-scales.tex", "tables/cluster-posterior.tex",
           "figures/cluster-hst.pdf", "figures/cluster-host-comparison.pdf", "figures/cluster-coverage.pdf"]


def rows(name):
    with (DATA / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def tables():
    coverage = {r["cluster"]: r for r in rows("coverage.csv")}
    catalog = {r["cluster"]: r for r in rows("catalog_scalars.csv")}
    pars = {(r["cluster"], r["parameter"]): r for r in rows("posterior_parameters.csv")}
    diagnostics = {d["cluster"]: d for d in json.loads((DATA / "posterior_diagnostics.json").read_text())}
    stars = {r["cluster"]: r for r in rows("gaia_star_catalog_quality.csv")}

    def environment(body, caption, label, columns, wide=False):
        env = "table*" if wide else "table"
        return (f"\\begin{{{env}}}[t]\n\\caption{{{caption}}}\\label{{{label}}}\n"
                f"\\centering\\begin{{tabular}}{{{columns}}}\n\\toprule\n"+body+
                f"\n\\bottomrule\n\\end{{tabular}}\n\\end{{{env}}}\n")

    body = r"Cluster & Gaia field rows & HST bins & Flagged rows \\" + "\n\\midrule\n"
    for cid in IDS:
        body += f"{NAMES[cid]} & {int(stars[cid]['rows']):,} & {coverage[cid]['hst_radial_bins']} & {coverage[cid]['compiled_quarantined_rows']} \\\\\n"
    totals = [sum(int(stars[c]["rows"]) for c in IDS),
              sum(int(coverage[c]["hst_radial_bins"]) for c in IDS),
              sum(int(coverage[c]["compiled_quarantined_rows"]) for c in IDS)]
    body += "\\midrule\nTotal & "+f"{totals[0]:,} & {totals[1]} & {totals[2]}"+r" \\"
    compiled_total = sum(int(coverage[c]["compiled_bins"]) for c in IDS)
    counts = environment(body,
        "Source coverage. Gaia field rows include nonmembers. HST bins refer to the original proper-motion profiles; the combined, radial and tangential quantities share stars. Flagged rows belong to the separate "
        +str(compiled_total)+"-row compilation and are not corrected or used as likelihood inputs.",
        "tab:cluster-coverage", "lrrr")
    body = (r"Cluster & $D$ (kpc) & $M/(10^5 M_\odot)$ & $r_h$ (pc) & $v_{\rm esc,0}$ & $v_0$ & $t_0$ (Myr) \\"+"\n"
            +r" & & & & (km\,s$^{-1}$) & (km\,s$^{-1}$) & \\"+"\n\\midrule\n")
    for cid in IDS:
        c = {k: float(v) for k,v in catalog[cid].items() if k in [
            "distance_kpc","distance_error_kpc","mass_msun","mass_error_msun","half_mass_3d_pc",
            "escape0_catalog_kms","scale_v_kms","scale_t_myr"]}
        body += (NAMES[cid]+" & $"+f"{c['distance_kpc']:.2f}\\pm{c['distance_error_kpc']:.2f}$"
                 +" & $"+f"{c['mass_msun']/1e5:.3f}\\pm{c['mass_error_msun']/1e5:.3f}$"
                 +f" & {c['half_mass_3d_pc']:.2f} & {c['escape0_catalog_kms']:.1f}"
                 +f" & {c['scale_v_kms']:.2f} & {c['scale_t_myr']:.3f} \\\\\n")
    scales = environment(body.rstrip(),
        r"Structural-catalog estimates and derived display scales. The selected catalog path supplies the quoted mass and distance errors, but no errors for $r_h$ or $v_{\rm esc,0}$; missing uncertainty is not zero. Model systematics and joint scale uncertainty are not included. The gravitational scales $v_0$ and $t_0$ are defined in Eq.~\eqref{eq:comparison-scales}; $t_0$ is not a friction time.",
        "tab:cluster-scales", "lrrrrrr", wide=True)
    body = (r"Cluster & $M/(10^5 M_\odot)$ & $r_h$ (pc) & $D$ (kpc) & $\mathrm{Corr}(M,r_h)$ \\"
            +"\n\\midrule\n")
    for cid in IDS:
        cells = []
        for key, factor in [("M", 1e5), ("rh", 1), ("d", 1)]:
            p = pars[cid,key]
            cells.append("$"+f"{float(p['median'])/factor:.3f}\\;[{float(p['p16'])/factor:.3f},\\,{float(p['p84'])/factor:.3f}]$")
        corr = diagnostics[cid]["posterior_correlation"][1][2]
        body += NAMES[cid]+" & "+" & ".join(cells)+f" & {corr:.3f} \\\\\n"
    posterior = environment(body.rstrip(),
        r"Direct weighted marginal medians and 16th--84th percentiles from the published standard-fit joint samples of Dickson et al.~\cite{Dickson2023,Dickson2024}. These are conditional model estimates, not new fits. The NGC 6397 row is the standard fit; its distinct no-black-hole alternative lacks a released chain in the inspected product. Correlations refer to the complete weighted samples.",
        "tab:cluster-posterior", "lrrrr", wide=True)
    return dict(zip(OUTPUTS[:3], [counts, scales, posterior]))


def checksums(paths, base):
    return {path: hashlib.sha256((base/path).read_bytes()).hexdigest() for path in paths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="read-only check; no plotting dependencies")
    args = parser.parse_args()
    subprocess.run([sys.executable, str(ROOT/"analysis/cluster_data/verify.py")], check=True)
    rendered = tables()
    if args.verify:
        manifest = json.loads(MANIFEST.read_text())
        if manifest["input_sha256"] != checksums(INPUTS, ROOT):
            raise SystemExit("Paper data inputs changed; explicitly regenerate and review the assets.")
        if manifest["output_sha256"] != checksums(OUTPUTS, PAPER):
            raise SystemExit("Paper data assets differ from their manifest.")
        if any((PAPER/name).read_text() != value for name,value in rendered.items()):
            raise SystemExit("Paper tables differ from recorded evidence.")
        print("Verified three evidence tables and three vector figures; no new dynamics.")
        return
    for name,value in rendered.items():
        (PAPER/name).parent.mkdir(parents=True, exist_ok=True)
        (PAPER/name).write_text(value)
    from data_figures import render
    render(DATA, PAPER/"figures")
    manifest = {"data_revision":"cluster-data-v0.1.0",
                "scope":"presentation of recorded observations and published posterior summaries",
                "transport_simulations_performed":False,
                "input_sha256":checksums(INPUTS, ROOT),
                "output_sha256":checksums(OUTPUTS, PAPER)}
    MANIFEST.write_text(json.dumps(manifest, indent=2)+"\n")
    print("Prepared three evidence tables and three vector figures from recorded data.")


if __name__ == "__main__":
    main()
