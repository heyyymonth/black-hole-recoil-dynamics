"""Validate the two-step specification without evaluating a physical model."""
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]
ACTIVE_DOCS = sorted({
    str(p.relative_to(ROOT))
    for pattern in ('*.md', 'docs/*.md', 'paper/*.md',
                    'analysis/*/README.md', 'data/*/README.md')
    for p in ROOT.glob(pattern)
})
INPUT_GROUPS = ("host", "launch", "collisions", "termination", "exposure", "numerics")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def verify():
    working = json.loads((ROOT / "paper/working.json").read_text())
    require(working["status"] == "theory-data-and-controlled-dynamics-draft", "Unexpected manuscript scope")
    require(working["compile_enabled"] is True, "Theory typesetting must be enabled")
    require(working["transport_simulations_performed"] is True
            and working["collisionless_trajectories_performed"] is True
            and working["collisional_trajectories_performed"] is False,
            "Only collisionless reference trajectories are supplied")
    require(working["data_analysis_performed"] is True, "Completed cluster-data analysis is included")
    require(working["active_steps"] == ["transport", "exposure"], "Unexpected paper scope")
    source = (ROOT / "paper/manuscript.tex").read_text()
    title = re.search(r"\\title\{([^{}]*)\}", source).group(1)
    abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", source, re.S).group(1).strip()
    require(title == working["title"], "Working title and LaTeX differ")
    normalize = lambda s: re.sub(r"\s+", " ", s).strip()
    require(normalize(abstract) == normalize(working["abstract"]), "Working abstract and LaTeX differ")
    readme = (ROOT / "README.md").read_text()
    require(readme.startswith("# " + title + "\n"), "Repository title differs")
    require(working["abstract"] in readme, "Working abstract missing from README")
    for name in ("docs/problem.md", "CITATION.cff", "CITATION.bib", "paper/CITATION.bib"):
        require(title in (ROOT / name).read_text(), f"Title mismatch: {name}")
    require(working["manuscript_version"] in (ROOT / "paper/CITATION.bib").read_text(),
            "Manuscript citation revision differs")
    require(working["manuscript_tag"] in (ROOT / "paper/CITATION.bib").read_text(),
            "Manuscript citation tag differs")
    for name in ("CITATION.cff", "CITATION.bib"):
        require(working["formulation_tag"] in (ROOT / name).read_text(), "Formulation citation lacks tag")

    inputs = json.loads((ROOT / "formulation/inputs.json").read_text())
    require(inputs["status"] == "specification-only", "Unexpected calculation status")
    require(inputs["execution_enabled"] is False, "Physics execution remains disabled")
    require(inputs["active_steps"] == ["transport", "exposure"], "Scope must remain two steps")
    for group in INPUT_GROUPS:
        require(inputs[group] and all(value is None for value in inputs[group].values()),
                f"Input assigned without a calculation specification: {group}")
    for stage in inputs["completed_calculations"]:
        require((ROOT / stage["protocol"]).is_file() and (ROOT / stage["result"]).is_file(), "Missing defined stage")
    require(working.get("host_reconstruction_performed") is True and working.get("local_collision_coefficients_evaluated") is True, "Completed calculation status missing")
    require((ROOT / inputs["handoff"]).is_file(), "Missing calculation handoff")

    links = 0
    displays = 0
    for name in ACTIVE_DOCS:
        path = ROOT / name
        text = path.read_text()
        require(len(re.findall(r"^\x60\x60\x60", text, re.M)) % 2 == 0,
                f"Unclosed fence: {name}")
        for expression in re.findall(r"\x60\x60\x60math\n(.*?)\n\x60\x60\x60", text, re.S):
            level = 0
            for brace in re.findall(r"(?<!\\)[{}]", expression):
                level += 1 if brace == "{" else -1
                require(level >= 0, f"Unbalanced display equation: {name}")
            require(level == 0, f"Unbalanced display equation: {name}")
            displays += 1
        for target in re.findall(r"\]\(([^)\s]+)\)", text):
            if "://" in target or target.startswith(("#", "mailto:")):
                continue
            require((path.parent / target.split("#", 1)[0]).exists(),
                    f"Broken link in {name}: {target}")
            links += 1
    require("Step 1:" in (ROOT / "docs/step1-transport.md").read_text(), "Missing Step 1")
    require("Step 2:" in (ROOT / "docs/step2-exposure.md").read_text(), "Missing Step 2")
    require(not re.search(r"E00[1-5]", source), "Internal experiment identifier in article")
    expected_tables = {"tables/collisionless-residence.tex", "tables/finite-angle-encounters.tex", "tables/cluster-coverage.tex", "tables/cluster-scales.tex",
                       "tables/cluster-posterior.tex", "tables/cluster-transport-comparison.tex", "tables/host-validation.tex", "tables/local-kernels.tex", "tables/outer-boundary.tex", "tables/galactic-orbit.tex"}
    require(set(re.findall(r"\\input\{([^}]+)\}", source)) == expected_tables,
            "Expected the ten data and calculation evidence tables")
    require(not any(key in source for key in ("tables/return.tex", "figures/return-time.pdf")),
            "Historical reference dynamics do not validate the active transport problem")
    print(f"Verified two-step foundation: {len(ACTIVE_DOCS)} documents, {links} local links, "
          f"{displays} display equations; inputs unassigned and title/abstract/citations consistent. "
          "Recorded data, coefficients and collisionless controls included; collisional residence remains unevaluated.")


if __name__ == "__main__":
    verify()
