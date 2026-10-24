"""Verify the current theoretical manuscript and its portable source archive."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import re
import unicodedata
import zipfile

from pypdf import PdfReader


ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.4.0-manuscript.1"
TAG = "manuscript-v0.4.0.1"
CITATION_URL = f"https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/{TAG}"
ARCHIVE_SOURCE_FILES = {
    "manuscript.tex": "paper/manuscript.tex",
    "references.bib": "paper/references.bib",
    "CITATION.bib": "paper/CITATION.bib",
    "SOURCE_README.txt": "paper/SOURCE_README.txt",
    "LICENSE": "LICENSE",
}
INPUT_FILES = tuple(sorted({
    *ARCHIVE_SOURCE_FILES.values(),
    "paper/build.py", "paper/verify.py", "paper/requirements.txt",
    "paper/toolchain.json", "paper/working.json", "paper/CLAIMS.md",
    "docs/problem.md", "docs/e001-ballistic-return.md", "docs/e002-core-exposure.md",
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
    require(state.get("status") == "theory-draft", "Expected a complete theoretical draft")
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
    require(not re.search(r"\\(?:input|include|includegraphics)\b", manuscript),
            "The theory manuscript must be self-contained and use no external result assets")
    require(not re.search(r"E00\d+", manuscript, flags=re.IGNORECASE),
            "Internal experiment identifiers do not belong in the article")
    check_citation((paper / "CITATION.bib").read_text(), (paper / "SOURCE_README.txt").read_text())
    return {"title": title, "manuscript": manuscript}


def inspect_pdf(path, title):
    reader = PdfReader(path, strict=True)
    require(len(reader.pages) > 0, "Empty manuscript PDF")
    text = unicodedata.normalize("NFKC", "\n".join(page.extract_text() for page in reader.pages))
    normalized = compact(text)
    for required in (title, "Venkata Hemanth Athota", "bharadwaj.hemanth.002@gmail.com", "Step 1", "Step 2"):
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


def verify(root=ROOT):
    meta = json.loads((root / "output/pdf/build.json").read_text())
    require(meta["version"] == VERSION, "Unexpected manuscript version")
    require(set(meta["input_sha256"]) == set(INPUT_FILES), "Unexpected theory input inventory")
    require(set(meta["output_sha256"]) == OUTPUT_FILES, "Unexpected manuscript output inventory")
    require(meta.get("physical_calculations_performed") is False, "Theory compilation cannot claim physical calculations")
    require(meta.get("executed_experiments") == [], "This manuscript reports no executed experiments")
    require(meta.get("figure_generation") == {"mode": "none", "physical_calculations_performed": False},
            "Theory compilation must not generate scientific figures")
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
    print(f"Verified theory manuscript {VERSION}: {pdf_info['pages']} US Letter pages, "
          f"{pdf_info['embedded_vector_fonts']} embedded vector fonts, current source and archive hashes, "
          "and matching citation metadata. No physical calculations performed.")


if __name__ == "__main__":
    verify()
