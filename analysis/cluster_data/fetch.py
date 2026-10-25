"""Download pinned source products for a six-cluster data audit, not dynamics.

First capture writes a manifest. Later runs require matching bytes and never
silently update it. Third-party originals are cached locally, outside Git.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/observed-clusters"
RAW = DATA / "raw"
CLUSTERS = ["NGC0104", "NGC1851", "NGC3201", "NGC5139", "NGC6121", "NGC6397"]
RESULTS_REV = "c45a03b1681570b6ee96be76f362894c6628e1f1"
OBS_REV = "ca63e5b465d01bb5b006dce7d40e28d0a95643ed"


def specifications():
    rows = []
    def add(name, url, category, role, license_note="Provider rights; not relicensed as MIT"):
        rows.append(dict(file=name, url=url, category=category, role=role, license=license_note))
    uq = "https://people.smp.uq.edu.au/HolgerBaumgardt/globular/"
    for name, path in [("uq_parameters.html", "newdata/parameter.html"),
                       ("uq_root_parameters.html", "parameter.html"),
                       ("uq_profiles.html", "newdata/indivdata/velsurf.html")]:
        add(name, uq+path, "documentation", "source definitions and catalog version audit")
    for cluster in CLUSTERS:
        cid = "ngc"+str(int(cluster[3:]))
        add(f"uq_{cluster}_kin.txt", uq+f"newdata/indivdata/{cid}_kin.txt",
            "compiled_profile", "audit only; source overlap and suspect rows retained")
        add(f"uq_{cluster}_sb.txt", uq+f"newdata/indivdata/{cluster.lower()}_sb.txt",
            "projected_surface_profile", "native log profile; normalization unresolved")
        if cluster != "NGC5139":
            add(f"hacks_{cluster}.txt",
                "https://archive.stsci.edu/hlsps/hacks/velocity_dispersions/"
                f"hlsp_hacks_hst_acs-wfc3_{cluster.lower()}_multi_v1.0_vdisp.txt",
                "original_hst_binned_profile", "native mas/yr; combined and components share stars", "CC BY 4.0; cite Libralato et al. 2022 and MAST HACKS")
        base = f"https://raw.githubusercontent.com/nmdickson/GCfit-results/{RESULTS_REV}/"
        for suffix, remote, category in [
            ("sampler.hdf", f"clusters/{cluster}/{cluster}_sampler.hdf", "published_joint_posterior"),
            ("summary.txt", f"clusters/{cluster}/summary.txt", "published_fit_summary"),
            ("gaia.hdf", f"Gaia_profiles/{cluster}_profiles.hdf", "gcfit_gaia_reanalysis")]:
            add(f"gcfit_{cluster}_{suffix}", base+remote, category,
                "pinned published model product; exact run/data compatibility must be audited")
        add(f"gcfit_{cluster}_observations.hdf",
            f"https://raw.githubusercontent.com/nmdickson/GCfit/{OBS_REV}/gcfit/resources/{cluster}.hdf",
            "gcfit_observation_bundle", "historical v1.0 candidate pairing; not a proven chain reproduction")
    add("gcfit_results_README.md", f"https://raw.githubusercontent.com/nmdickson/GCfit-results/{RESULTS_REV}/README.md",
        "documentation", "posterior release instructions")
    add("gcfit_NGC6397_no_BHs_summary.txt",
        f"https://raw.githubusercontent.com/nmdickson/GCfit-results/{RESULTS_REV}/clusters/NGC6397/no_BHs/summary.txt",
        "published_fit_summary", "alternative model; corresponding chain absent from inspected release")
    add("gaia_edr3_v2.zip", "https://zenodo.org/records/4891252/files/clusters.zip?download=1",
        "original_gaia_catalog", "full original archive, membership and correlated fitted profile curves",
        "CC BY 4.0; Vasiliev and Baumgardt 2021; DOI 10.5281/zenodo.4891252")
    for filename in ["proper_motion_dispersion_log_bins.dat", "los_profile.dat",
                     "energy_equipartition_profiles.dat", "energy_equipartition_massbins.dat"]:
        add("h25_"+filename, f"https://zenodo.org/records/14978551/files/{filename}",
            "original_hst_profile" if filename.startswith("proper") else "original_h25_support",
            "Omega Centauri original release; no pooling alternative binnings",
            "CC BY 4.0; Haberle et al. 2025; DOI 10.5281/zenodo.14978551")
    return rows


def fetch_one(spec):
    path = RAW / spec["file"]
    if not path.exists():
        part = path.with_suffix(path.suffix + ".part")
        subprocess.run(["curl", "--fail", "--silent", "--show-error", "--location",
                        "--retry", "2", "--connect-timeout", "20", "--max-time", "240",
                        spec["url"], "--output", str(part)], check=True)
        part.replace(path)
    body = path.read_bytes()
    if path.suffix in (".hdf", ".zip"):
        signature = b"\x89HDF" if path.suffix == ".hdf" else b"PK"
        if not body.startswith(signature):
            raise ValueError(f"Unexpected binary signature: {path.name}")
    if spec["file"] == "gaia_edr3_v2.zip":
        assert hashlib.md5(body).hexdigest() == "3a07093be9a1fcbfb0e3c545ca175c99"
    return {**spec, "bytes":len(body), "sha256":hashlib.sha256(body).hexdigest()}


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    path = DATA / "sources.json"
    prior = json.loads(path.read_text()) if path.exists() else None
    specs = specifications()
    with ThreadPoolExecutor(max_workers=5) as pool:
        rows = list(pool.map(fetch_one, specs))
    if prior:
        assert rows == prior["files"], "Source bytes/definitions changed; do not overwrite this release."
    else:
        manifest = {"capture_date_local":"2026-10-24", "timezone":"America/Los_Angeles",
                    "retrieved_utc":datetime.now(timezone.utc).isoformat(),
                    "scope":"original source retrieval and data analysis; no host fit or dynamics",
                    "clusters":CLUSTERS, "gcfit_results_commit":RESULTS_REV,
                    "gcfit_observation_candidate_commit":OBS_REV, "files":rows}
        path.write_text(json.dumps(manifest, indent=2)+"\n")
    print(f"Verified {len(rows)} originals, {sum(r['bytes'] for r in rows):,} bytes. No physical experiment run.")


if __name__ == "__main__":
    main()
