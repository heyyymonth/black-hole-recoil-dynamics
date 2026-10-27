"""Compile the formulation and data-analysis manuscript from verified assets.

The source archive contains the article, tables, vector figures, supporting
data evidence, bibliography, license, and pinned REVTeX distribution. The build
checks existing records and typesets them; it does not rerun data analysis,
generate figures, evaluate a host, or integrate a remnant trajectory.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import tempfile
import zipfile

from verify import (
    ARCHIVE_SOURCE_FILES, DATA_REVISION, FIGURE_ASSETS, INPUT_FILES, ROOT,
    TABLE_ASSETS, VERSION, inspect_pdf, source_metadata, verify_evidence,
)


PAPER = ROOT / "paper"
OUT = ROOT / "output/pdf"
TMP = ROOT / "tmp/latex"
EPOCH = datetime(2026, 10, 26, tzinfo=timezone.utc)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    verify_evidence(ROOT)
    source = source_metadata(ROOT)
    toolchain = json.loads((PAPER / "toolchain.json").read_text())
    for relative in INPUT_FILES:
        if not (ROOT / relative).is_file():
            raise SystemExit(f"Missing manuscript input: {relative}")
    engine = shutil.which("tectonic")
    if not engine:
        raise SystemExit("Install Tectonic 0.16.9; see paper/README.md.")
    engine_version = subprocess.check_output([engine, "--version"], text=True).strip()
    if engine_version != f"Tectonic {toolchain['tectonic_version']}":
        raise SystemExit(f"Unexpected LaTeX engine: {engine_version}")

    TMP.mkdir(parents=True, exist_ok=True)
    revtex_archive = ROOT / "tmp/revtex-4.2f.zip"
    if not revtex_archive.exists():
        source_archive = OUT / "manuscript-source.zip"
        if source_archive.exists():
            with zipfile.ZipFile(source_archive) as bundle:
                revtex_archive.write_bytes(bundle.read("vendor/revtex-4.2f.zip"))
        else:
            subprocess.run(
                ["curl", "--fail", "--silent", "--show-error", "--location",
                 "--connect-timeout", "15", "--max-time", "120", "--retry", "2",
                 toolchain["revtex_archive_url"], "--output", str(revtex_archive)],
                check=True,
            )
    if sha256(revtex_archive) != toolchain["revtex_archive_sha256"]:
        raise SystemExit("REVTeX dependency checksum mismatch")

    env = os.environ.copy()
    env["SOURCE_DATE_EPOCH"] = str(int(EPOCH.timestamp()))
    with tempfile.TemporaryDirectory(prefix="manuscript-source-", dir=TMP) as stage_name:
        stage = Path(stage_name)
        with zipfile.ZipFile(revtex_archive) as distribution:
            for item in distribution.infolist():
                if item.filename.startswith(("tex/latex/revtex/", "bibtex/bst/revtex/")) and not item.is_dir():
                    (stage / Path(item.filename).name).write_bytes(distribution.read(item))
        for name in ("manuscript.tex", "references.bib", *TABLE_ASSETS, *FIGURE_ASSETS):
            (stage / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(PAPER / name, stage / name)
        command = [engine, "--untrusted", "--keep-logs", "--keep-intermediates",
                   "--outdir", str(stage), "manuscript.tex"]
        completed = subprocess.run(command, cwd=stage, env=env, text=True, capture_output=True)
        (TMP / "build-console.txt").write_text(completed.stdout + completed.stderr)
        for name in ("manuscript.log", "manuscript.blg"):
            if (stage / name).exists():
                shutil.copyfile(stage / name, TMP / name)
        if completed.returncode:
            raise SystemExit(completed.stdout + completed.stderr)

        log = (stage / "manuscript.log").read_text(errors="replace")
        blg = (stage / "manuscript.blg").read_text(errors="replace")
        for problem in (r"Overfull \\[hv]box", r"Citation .+ undefined", r"Reference .+ undefined",
                        r"There were undefined", r"Missing character:", r"! LaTeX Error"):
            if re.search(problem, log):
                raise SystemExit(f"LaTeX validation failed: {problem}; inspect tmp/latex/manuscript.log")
        style_diagnostics = ["Warning--jnrlst (dependency: not reversed) set 1"]
        entry_warnings = [line for line in blg.splitlines()
                          if line.startswith("Warning--") and line not in style_diagnostics]
        if entry_warnings or "I couldn't open" in blg:
            raise SystemExit("BibTeX warnings found; inspect tmp/latex/manuscript.blg")
        pdf_info = inspect_pdf(stage / "manuscript.pdf", source["title"])

        archive_path = stage / "manuscript-source.zip"
        with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
            def add(name, content):
                info = zipfile.ZipInfo(name, date_time=(2026, 10, 24, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                bundle.writestr(info, content)

            for member, relative in ARCHIVE_SOURCE_FILES.items():
                add(member, (ROOT / relative).read_bytes())
            add("manuscript.bbl", (stage / "manuscript.bbl").read_bytes())
            add("vendor/revtex-4.2f.zip", revtex_archive.read_bytes())

        OUT.mkdir(parents=True, exist_ok=True)
        pdf = OUT / "central-residence.pdf"
        archive = OUT / "manuscript-source.zip"
        shutil.copyfile(stage / "manuscript.pdf", pdf)
        shutil.copyfile(archive_path, archive)

    manifest = {
        "paper": "recoil-transport-and-environmental-exposure",
        "version": VERSION,
        "title": source["title"],
        "pages": pdf_info["pages"],
        "scope": "Two-step transport formulation, six-cluster source audit, six conditional reference hosts, joint-row sensitivity, local kernels and environmental diagnostics",
        "transport_simulations_performed": False,
        "local_calculations": {"host": "host-reconstruction-1", "coefficients": "local-transport-kernels-1", "boundary": "outer-boundary-1", "galactic_orbit": "galactic-orbit-1", "cross_cluster_comparison": "six-cluster-comparison-1", "recomputed_during_build": False},
        "executed_during_build": [],
        "data_analysis": {"revision": DATA_REVISION, "mode": "reuse-verified-recorded-analysis",
                          "recomputed_during_build": False},
        "figure_generation": {"mode": "reuse-verified-assets", "generated_during_build": False},
        "format": "REVTeX 4.2 APS/PRD reprint; US Letter; two columns; txfonts",
        "engine": engine_version,
        "source_date_epoch": env["SOURCE_DATE_EPOCH"],
        "revtex_log": next(line.strip() for line in log.splitlines()
                           if line.startswith("Document Class: revtex")),
        "versions": {"pypdf": importlib.metadata.version("pypdf")},
        "checks": {
            "undefined_references": 0, "overfull_boxes": 0,
            "bibliography_entry_warnings": 0, "visible_experiment_identifiers": 0,
            "embedded_vector_fonts": pdf_info["embedded_vector_fonts"],
            "raster_objects": 0,
            "allowed_bst_diagnostics": [line for line in style_diagnostics if line in blg],
        },
        "input_sha256": {name: sha256(ROOT / name) for name in INPUT_FILES},
        "output_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (pdf, archive)},
    }
    (OUT / "build.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Built {pdf_info['pages']}-page formulation and data manuscript {VERSION} "
          "and portable source archive. Reused verified data assets; no transport simulations performed.")


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__).parse_args()
    build()
