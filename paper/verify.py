"""Verify the formulation and observed-cluster manuscript without rerunning research."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import re
import subprocess
import sys
import unicodedata
import zipfile

from pypdf import PdfReader


ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.8.0-manuscript.2"
TAG = "manuscript-v0.8.0.2"
DATA_REVISION = "cluster-data-v0.1.0"
CITATION_URL = f"https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/{TAG}"
TABLE_ASSETS = (
    "tables/cluster-coverage.tex", "tables/cluster-scales.tex",
    "tables/cluster-posterior.tex", "tables/host-validation.tex", "tables/local-kernels.tex", "tables/outer-boundary.tex", "tables/galactic-orbit.tex",
)
FIGURE_ASSETS = (
    "figures/cluster-hst.pdf", "figures/cluster-host-comparison.pdf",
    "figures/cluster-coverage.pdf", "figures/host-projections.pdf", "figures/scattering-composition.pdf", "figures/outer-boundary.pdf", "figures/galactic-orbit.pdf",
)
DATA_EVIDENCE = (
    "README.md", "sources.json", "audit.json", "candidate-hosts.json",
    "derived/catalog_scalars.csv", "derived/root_catalog_comparison.csv",
    "derived/coverage.csv", "derived/hst_native_and_common_units.csv",
    "derived/compiled_kinematics_audit.csv", "derived/surface_profiles.csv",
    "derived/gaia_correlated_curves.csv", "derived/gaia_star_catalog_quality.csv",
    "derived/posterior_parameters.csv", "derived/posterior_diagnostics.json",
    "derived/quality.json", "derived/model_product_inventory.json",
    "derived/observation_bundle_inventory.json",
)
CALCULATION_RECORD = json.loads((ROOT / "results/calculations.json").read_text())
CALCULATION_FILES = tuple(sorted({"results/calculations.json", *CALCULATION_RECORD["input_sha256"], *CALCULATION_RECORD["output_sha256"]}))
BOUNDARY_RECORD = json.loads((ROOT / "results/outer-boundary/provenance.json").read_text())
BOUNDARY_FILES = tuple(sorted({"results/outer-boundary/provenance.json", *BOUNDARY_RECORD["input_sha256"], *BOUNDARY_RECORD["output_sha256"]}))
ORBIT_RECORD = json.loads((ROOT / "results/orbit/provenance.json").read_text())
ORBIT_FILES = tuple(sorted({"results/orbit/provenance.json", *ORBIT_RECORD["input_sha256"], *ORBIT_RECORD["output_sha256"]}))
ARCHIVE_SOURCE_FILES = {
    "manuscript.tex": "paper/manuscript.tex",
    "references.bib": "paper/references.bib",
    "CITATION.bib": "paper/CITATION.bib",
    "SOURCE_README.txt": "paper/SOURCE_README.txt",
    "LICENSE": "LICENSE",
    **{name: f"paper/{name}" for name in (*TABLE_ASSETS, *FIGURE_ASSETS)},
    "evidence/paper-data-assets.json": "paper/data-assets.json",
    "evidence/paper-calculation-assets.json": "paper/calculation-assets.json",
    "evidence/paper-boundary-assets.json": "paper/boundary-assets.json",
    "evidence/paper-orbit-assets.json": "paper/orbit-assets.json",
    **{f"evidence/orbit/{name}": name for name in ORBIT_FILES if name not in {*CALCULATION_FILES, *BOUNDARY_FILES}},
    **{f"evidence/boundary/{name}": name for name in BOUNDARY_FILES
       if name not in {*CALCULATION_FILES, *(f"data/observed-clusters/{item}" for item in DATA_EVIDENCE)}},
    **{f"evidence/calculations/{name}": name for name in CALCULATION_FILES},
    **{f"evidence/observed-clusters/{name}": f"data/observed-clusters/{name}"
       for name in DATA_EVIDENCE},
}
DATA_AUDIT_FILES = tuple(item["path"] for item in json.loads(
    (ROOT / "data/observed-clusters/audit.json").read_text())["files"])
INPUT_FILES = tuple(sorted({
    *ARCHIVE_SOURCE_FILES.values(),
    *DATA_AUDIT_FILES,
    "paper/data_assets.py", "paper/data_figures.py", "paper/calculation_assets.py", "paper/calculation_figures.py", "docs/local-transport-calculations.md", "docs/host-reconstruction-protocol.md",
    "paper/boundary_assets.py", "docs/outer-boundary.md", "analysis/boundary/README.md",
    "paper/orbit_assets.py", "docs/galactic-orbit.md", "analysis/orbit/README.md", "data/orbit/README.md",
    "paper/build.py", "paper/verify.py", "paper/requirements.txt",
    "paper/toolchain.json", "paper/working.json", "paper/CLAIMS.md",
    "docs/problem.md", "docs/step1-transport.md", "docs/step2-exposure.md",
    "docs/recoil-transport-theory.md", "docs/calculation-handoff.md",
    "formulation/inputs.json",
}))
OUTPUT_FILES = {"output/pdf/central-residence.pdf", "output/pdf/manuscript-source.zip"}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def compact(text):
    return re.sub(r"\s+", "", unicodedata.normalize("NFKC", text)).lower()


def check_citation(citation, readme):
    require(re.findall(r"revision\s+([^},\s]+)", citation) == [VERSION],
            "Citation revision mismatch")
    require(re.findall(r"\burl\s*=\s*\{([^{}]+)\}", citation) == [CITATION_URL],
            "Citation source tag mismatch")
    require(f"revision {VERSION}" in readme and CITATION_URL in readme,
            "Source README revision or source tag mismatch")


def source_metadata(root):
    paper = root / "paper"
    state = json.loads((paper / "working.json").read_text())
    require(state.get("compile_enabled") is True,
            "Working manuscript source is not ready for compilation (paper/working.json)")
    require(state.get("status") == "theory-data-and-local-calculations-draft", "Expected a formulation and data-analysis draft")
    require(state.get("data_analysis_performed") is True, "Completed data analysis must be recorded")
    require(state.get("transport_simulations_performed") is False, "No transport simulation is reported")
    require(state.get("host_reconstruction_performed") is True, "Missing reconstructed-host status")
    require(state.get("local_collision_coefficients_evaluated") is True, "Missing collision-evaluation status")
    require(state.get("outer_boundary_diagnostics_evaluated") is True, "Missing boundary-evaluation status")
    require(state.get("galactic_orbits_integrated") is True and state.get("host_re_equilibrated") is False, "Incorrect orbit or host-evolution status")
    manuscript = (paper / "manuscript.tex").read_text()
    title_match = re.search(r"\\title\{([^{}]+)\}", manuscript)
    require(title_match is not None, "Expected a plain-text manuscript title")
    title = re.sub(r"\s+", " ", title_match.group(1)).strip()
    require(title == state["title"], "Working metadata and LaTeX titles differ")
    document_class = re.search(r"\\documentclass\[([^]]+)\]\{revtex4-2\}", manuscript)
    require(document_class is not None and "reprint" in document_class.group(1).split(","),
            "Use the two-column REVTeX reprint format")
    require(re.search(r"\\usepackage(?:\[[^]]*\])?\{txfonts\}", manuscript),
            "Expected the declared vector font package")
    asset_commands = re.findall(r"\\(input|include|includegraphics)\s*(?:\[[^]]*\])?\s*\{([^{}]+)\}", manuscript)
    require(len(asset_commands) == len(re.findall(r"\\(?:input|include|includegraphics)\b", manuscript)),
            "Unrecognized external asset command")
    require(sorted(path for kind, path in asset_commands if kind == "input") == sorted(TABLE_ASSETS),
            "The article must use exactly the declared cluster-data tables")
    require(sorted(path for kind, path in asset_commands if kind == "includegraphics") == sorted(FIGURE_ASSETS),
            "The article must use exactly the declared cluster-data figures")
    require(not any(kind == "include" for kind, _ in asset_commands), "Undeclared section include")
    require(not re.search(r"E00\d+", manuscript, flags=re.IGNORECASE),
            "Internal experiment identifiers do not belong in the article")
    bibliography = (paper / "references.bib").read_text()
    require("github.com/heyyymonth" not in bibliography.lower(), "The article must not cite its own repository")
    require(not re.search(r"\\section\*?\{[^}]*data and code availability[^}]*\}", manuscript, flags=re.I),
            "No data and code availability section is requested")
    check_citation((paper / "CITATION.bib").read_text(), (paper / "SOURCE_README.txt").read_text())
    return {"title": title, "manuscript": manuscript}


def inspect_pdf(path, title):
    reader = PdfReader(path, strict=True)
    require(len(reader.pages) > 0, "Empty manuscript PDF")
    text = unicodedata.normalize("NFKC", "\n".join(page.extract_text() for page in reader.pages))
    normalized = compact(text)
    for required in (title, "Venkata Hemanth Athota", "bharadwaj.hemanth.002@gmail.com",
                     "Step 1", "Step 2", "NGC 3201", "NGC 6397"):
        require(compact(required) in normalized, f"Missing manuscript text: {required}")
    require(not re.search(r"E00\d+", re.sub(r"\s+", "", text), flags=re.IGNORECASE),
            "Visible internal experiment identifier")
    for forbidden in ("BHRD-TN-001", "Not peer reviewed", "Published as a conference paper", "DRAFT"):
        require(forbidden not in text, f"Unexpected publication label: {forbidden}")

    fonts = set()
    visited = set()

    def inspect_resources(reference):
        resources = reference.get_object()
        marker = id(resources)
        if marker in visited:
            return
        visited.add(marker)
        font_dictionary = resources.get("/Font")
        if font_dictionary is not None:
            for reference in font_dictionary.get_object().values():
                font = reference.get_object()
                require(font.get("/Subtype") != "/Type3", "Bitmap/Type3 font found")
                for face in font.get("/DescendantFonts", [font]):
                    face = face.get_object()
                    descriptor = face.get("/FontDescriptor")
                    require(descriptor is not None, f"Missing font descriptor: {face.get('/BaseFont')}")
                    descriptor = descriptor.get_object()
                    require(any(name in descriptor for name in ("/FontFile", "/FontFile2", "/FontFile3")),
                            "Font is not embedded")
                    fonts.add(str(face.get("/BaseFont")))
        objects = resources.get("/XObject")
        if objects is not None:
            for reference in objects.get_object().values():
                obj = reference.get_object()
                require(obj.get("/Subtype") != "/Image", "Raster figure or equation found")
                if "/Resources" in obj:
                    inspect_resources(obj["/Resources"])

    for page in reader.pages:
        require(abs(float(page.mediabox.width) - 612) < 0.1
                and abs(float(page.mediabox.height) - 792) < 0.1, "Expected US Letter pages")
        inspect_resources(page["/Resources"])
    require(fonts, "No embedded fonts found")
    return {"pages": len(reader.pages), "embedded_vector_fonts": len(fonts)}


def verify_evidence(root=ROOT):
    """Read-only checks; neither verifier downloads data or evaluates a host."""
    for command in ([sys.executable, str(root / "analysis/cluster_data/verify.py")],
                    [sys.executable, str(root / "paper/data_assets.py"), "--verify"],
                    [sys.executable, str(root / "analysis/host_reconstruction/verify.py")],
                    [sys.executable, str(root / "paper/calculation_assets.py"), "--verify"],
                    [sys.executable, str(root / "analysis/boundary/verify.py")],
                    [sys.executable, str(root / "paper/boundary_assets.py"), "--verify"],
                    [sys.executable, str(root / "analysis/orbit/verify.py")],
                    [sys.executable, str(root / "paper/orbit_assets.py"), "--verify"]):
        subprocess.run(command, cwd=root, check=True)


def verify(root=ROOT):
    verify_evidence(root)
    meta = json.loads((root / "output/pdf/build.json").read_text())
    require(meta["version"] == VERSION, "Unexpected manuscript version")
    require(set(meta["input_sha256"]) == set(INPUT_FILES), "Unexpected manuscript input inventory")
    require(set(meta["output_sha256"]) == OUTPUT_FILES, "Unexpected manuscript output inventory")
    require(meta.get("transport_simulations_performed") is False, "No transport simulation is reported")
    require(meta.get("local_calculations") == {"host": "host-reconstruction-1", "coefficients": "local-transport-kernels-1", "boundary": "outer-boundary-1", "galactic_orbit": "galactic-orbit-1", "recomputed_during_build": False}, "Incorrect calculation provenance")
    require(meta.get("executed_during_build") == [], "Compilation must not execute research calculations")
    require(meta.get("data_analysis") == {
        "revision": DATA_REVISION, "mode": "reuse-verified-recorded-analysis",
        "recomputed_during_build": False}, "Incorrect data-analysis provenance")
    require(meta.get("figure_generation") == {
        "mode": "reuse-verified-assets", "generated_during_build": False},
        "Compilation must reuse the verified scientific assets")
    for group in ("input_sha256", "output_sha256"):
        for relative, expected in meta[group].items():
            path = PurePosixPath(relative)
            require(not path.is_absolute() and ".." not in path.parts, f"Unsafe manifest path: {relative}")
            require(re.fullmatch(r"[0-9a-f]{64}", expected), f"Invalid checksum: {relative}")
            actual = hashlib.sha256((root / relative).read_bytes()).hexdigest()
            require(actual == expected, f"{group} mismatch: {relative}; rebuild the paper")

    source = source_metadata(root)
    require(meta["title"] == source["title"], "Build metadata and source titles differ")
    toolchain = json.loads((root / "paper/toolchain.json").read_text())
    require(meta["engine"] == f"Tectonic {toolchain['tectonic_version']}", "Incorrect build engine")
    require(len(set(ARCHIVE_SOURCE_FILES.values())) == len(ARCHIVE_SOURCE_FILES),
            "A source file is duplicated in the portable archive")
    expected_members = set(ARCHIVE_SOURCE_FILES) | {"manuscript.bbl", "vendor/revtex-4.2f.zip"}
    with zipfile.ZipFile(root / "output/pdf/manuscript-source.zip") as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)) and set(names) == expected_members,
                "Unexpected or duplicate portable source archive member")
        require(archive.testzip() is None, "Source archive integrity failure")
        for member, relative in ARCHIVE_SOURCE_FILES.items():
            require(archive.read(member) == (root / relative).read_bytes(), f"Archive mismatch: {member}")
        require(archive.read("manuscript.bbl").strip(), "Empty compiled bibliography")
        require(hashlib.sha256(archive.read("vendor/revtex-4.2f.zip")).hexdigest()
                == toolchain["revtex_archive_sha256"], "Bundled REVTeX checksum mismatch")
        check_citation(archive.read("CITATION.bib").decode(), archive.read("SOURCE_README.txt").decode())

    pdf_info = inspect_pdf(root / "output/pdf/central-residence.pdf", source["title"])
    require(pdf_info["pages"] == meta["pages"], "Recorded PDF page count mismatch")
    require(pdf_info["embedded_vector_fonts"] == meta["checks"]["embedded_vector_fonts"],
            "Recorded embedded font count mismatch")
    print(f"Verified formulation and data manuscript {VERSION}: {pdf_info['pages']} US Letter pages, "
          f"{pdf_info['embedded_vector_fonts']} embedded vector fonts, current source and archive hashes, "
          "and matching citation metadata. No transport simulations performed.")


if __name__ == "__main__":
    verify()
