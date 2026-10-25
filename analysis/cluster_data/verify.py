"""Check released data-audit artifacts without downloads or scientific execution.

Normal invocation is read-only.  --write-lock explicitly records a reviewed
artifact revision, after the same structural and consistency checks pass.
"""
import argparse
import base64
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "analysis/cluster_data"
DATA = ROOT / "data/observed-clusters"
DERIVED = DATA / "derived"
LOCK = DATA / "audit.json"
IDS = {"NGC0104", "NGC1851", "NGC3201", "NGC5139", "NGC6121", "NGC6397"}
DOCS = [ROOT / "docs/cluster-data-assessment.md", ROOT / "docs/cluster-host-schema.md",
        DATA / "README.md"]
PARAMS = ["W0", "M", "rh", "ra", "g", "delta", "s2", "F", "a1", "a2", "a3", "BHret", "d"]
UNASSIGNED = ["mass_resolved_distribution_functions", "consistent_potential",
              "external_field_and_evolution", "merger_launch_law", "remnant_mass_law",
              "collision_cutoffs_and_closure", "region_and_horizon"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reject_constant(value):
    raise ValueError(f"Nonstandard JSON constant: {value}")


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path):
    result = json.loads(path.read_text(), parse_constant=reject_constant,
                        object_pairs_hook=unique_object)
    def finite(value):
        if isinstance(value, float):
            require(math.isfinite(value), f"Nonfinite JSON value: {path.name}")
        elif isinstance(value, dict):
            for item in value.values():
                finite(item)
        elif isinstance(value, list):
            for item in value:
                finite(item)
    finite(result)
    return result


def read_csv(name):
    path = DERIVED / name
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)),
                f"Missing or duplicate CSV headers: {name}")
        rows = list(reader)
    require(rows and all(None not in row and None not in row.values() for row in rows),
            f"Empty or inconsistent CSV rows: {name}")
    return rows


def number(value):
    result = float(value)
    require(math.isfinite(result), f"Nonfinite numeric value: {value}")
    return result


def integer(value):
    result = number(value)
    require(result.is_integer(), f"Expected integer: {value}")
    return int(result)


def same_number(left, right, label):
    require(math.isclose(number(left), number(right), rel_tol=1e-10, abs_tol=1e-11),
            f"Inconsistent {label}: {left} versus {right}")


def keyed(rows, field="cluster"):
    require(len(rows) == 6 and {row[field] for row in rows} == IDS,
            f"Expected exactly six canonical IDs in {field}")
    return {row[field]: row for row in rows}


def scoped_files():
    paths = set(DOCS)
    for base in (ANALYSIS, DATA):
        for path in base.rglob("*"):
            relative = path.relative_to(base)
            if ((base == ANALYSIS and "__pycache__" in relative.parts)
                    or (base == DATA and relative.parts[0] == "raw")):
                continue
            if path == LOCK or not path.is_file():
                continue
            require(not path.is_symlink(), f"Audited artifact is a symlink: {path}")
            paths.add(path)
    require(all(path.is_file() for path in paths), "An audited document is missing")
    return sorted(paths)


def validate_sources(quality):
    sources = read_json(DATA / "sources.json")
    require(len(sources["clusters"]) == 6 and set(sources["clusters"]) == IDS,
            "Source manifest cluster list differs")
    files = sources["files"]
    require(len(files) == 51, "Expected 51 captured source records")
    require(len({f["file"] for f in files}) == 51, "Duplicate source file")
    require(len({f["url"] for f in files}) == 51, "Duplicate source URL")
    for item in files:
        require(Path(item["file"]).name == item["file"], "Source name must be a basename")
        uri = urlsplit(item["url"])
        require(uri.scheme in {"https", "http"} and bool(uri.netloc), "Invalid source URL")
        require(re.fullmatch(r"[0-9a-f]{64}", item["sha256"]), "Invalid source SHA256")
        require(type(item["bytes"]) is int and item["bytes"] > 0, "Invalid source size")
    require(quality["source_files"] == len(files), "Source count differs from quality report")
    require(quality["source_bytes"] == sum(f["bytes"] for f in files),
            "Source bytes differ from quality report")
    return {f["file"] for f in files}


def validate_tables(quality):
    coverage = keyed(read_csv("coverage.csv"))
    catalog = keyed(read_csv("catalog_scalars.csv"))
    keyed(read_csv("root_catalog_comparison.csv"))
    stars = keyed(read_csv("gaia_star_catalog_quality.csv"))
    specifications = [("hst_native_and_common_units.csv", "hst_component_records"),
                      ("compiled_kinematics_audit.csv", "compiled_kinematic_rows"),
                      ("surface_profiles.csv", "surface_rows"),
                      ("gaia_correlated_curves.csv", "gaia_curve_samples")]
    tables = {}
    for filename, metric in specifications:
        rows = read_csv(filename)
        require({r["cluster"] for r in rows} == IDS, f"Cluster coverage differs: {filename}")
        require(len(rows) == quality[metric], f"Count differs: {filename}")
        tables[filename] = rows
    hst = tables["hst_native_and_common_units.csv"]
    groups = Counter((r["cluster"], r["source"], r["bin_id"]) for r in hst)
    require(set(groups.values()) == {3}, "HST bins must each retain three component records")
    require(len({(r["cluster"], r["source"], r["bin_id"], r["component"]) for r in hst}) == len(hst),
            "Duplicate HST component record")
    for group in groups:
        require({r["component"] for r in hst if (r["cluster"], r["source"], r["bin_id"]) == group}
                == {"PM_COMBINED", "PM_R", "PM_T"}, "Unexpected HST component set")
    require(len(groups) == quality["hst_unique_bins"], "HST unique-bin count differs")
    compiled = tables["compiled_kinematics_audit.csv"]
    require(all(r["quarantine"] in {"True", "False"} for r in compiled), "Invalid quarantine flag")
    flagged = [r for r in compiled if r["quarantine"] == "True"]
    require(len(flagged) == quality["compiled_flagged_rows"] == len(quality["flags"]),
            "Quarantine counts differ")
    require({f'{r["cluster"]}:{r["source_line"]}' for r in flagged}
            == {r["row_id"] for r in quality["flags"]}, "Quarantine row IDs differ")
    for cid in IDS:
        row = coverage[cid]
        own = [r for r in hst if r["cluster"] == cid]
        bins = len(own) // 3
        require(integer(row["hst_radial_bins"]) == bins, f"HST coverage differs: {cid}")
        require(integer(row["hst_scalar_components"]) == len(own), f"HST records differ: {cid}")
        for filename, key in [("compiled_kinematics_audit.csv", "compiled_bins"),
                              ("surface_profiles.csv", "surface_bins"),
                              ("gaia_correlated_curves.csv", "gaia_correlated_grid_points")]:
            require(integer(row[key]) == sum(r["cluster"] == cid for r in tables[filename]),
                    f"Coverage {key} differs: {cid}")
        require(integer(row["compiled_quarantined_rows"]) == sum(r["cluster"] == cid for r in flagged),
                f"Quarantine coverage differs: {cid}")
        require(integer(stars[cid]["duplicate_source_ids"]) == 0, f"Duplicate Gaia IDs: {cid}")
    require(sum(integer(r["rows"]) for r in stars.values()) == quality["gaia_stars_inspected"],
            "Gaia stellar-row total differs")
    parameters = read_csv("posterior_parameters.csv")
    require(len(parameters) == 78 and len({(r["cluster"], r["parameter"]) for r in parameters}) == 78,
            "Expected 78 unique posterior parameter records")
    for cid in IDS:
        require({r["parameter"] for r in parameters if r["cluster"] == cid} == set(PARAMS),
                f"Posterior parameter names differ: {cid}")
    return coverage, catalog


def validate_candidates(coverage, catalog, source_names):
    candidates = read_json(DATA / "candidate-hosts.json")
    hosts = keyed(candidates["hosts"], "cluster_id")
    posterior = keyed(read_json(DERIVED / "posterior_diagnostics.json"))
    for cid, host in hosts.items():
        require(host["execution_enabled"] is False, f"Execution enabled: {cid}")
        require(all(key in host and host[key] is None for key in UNASSIGNED),
                f"Physical input assigned or absent: {cid}")
        require(host["native_hst_radial_bins"] == integer(coverage[cid]["hst_radial_bins"]),
                f"Candidate HST count differs: {cid}")
        require(set(host["catalog_display_scales"]) == {
            "distance_kpc", "distance_error_kpc", "mass_msun", "mass_error_msun",
            "half_mass_3d_pc", "scale_v_kms", "scale_t_myr"},
            f"Candidate display-scale fields differ: {cid}")
        for key, value in host["catalog_display_scales"].items():
            same_number(value, catalog[cid][key], f"candidate {cid}/{key}")
        proposal = host["candidate_joint_host_posterior"]
        diag = posterior[cid]
        require(proposal["parameters"] == diag["parameters"] == PARAMS,
                f"Posterior ordering differs: {cid}")
        require(proposal["joint_rows"] == diag["rows"], f"Posterior row count differs: {cid}")
        same_number(proposal["effective_sample_size"], diag["posterior_ess"], f"ESS {cid}")
        require(0 < number(diag["posterior_ess"]) <= integer(diag["rows"]), f"Invalid ESS: {cid}")
        for relative in [host["catalog_display_scale_source"], proposal["source_file"]]:
            require(relative.startswith("raw/") and relative[4:] in source_names,
                    f"Candidate source absent from manifest: {relative}")


def validate_notebook():
    notebook = read_json(ANALYSIS / "cluster_data_audit.ipynb")
    cells = [c for c in notebook["cells"] if c["cell_type"] == "code"]
    require(cells, "Notebook has no code cells")
    images = 0
    for cell in cells:
        require(type(cell.get("execution_count")) is int and cell["execution_count"] > 0,
                "Notebook contains an unexecuted code cell")
        for output in cell.get("outputs", []):
            require(output.get("output_type") != "error", "Notebook contains an execution error")
            png = output.get("data", {}).get("image/png")
            if png is not None:
                encoded = "".join(png) if isinstance(png, list) else png
                payload = base64.b64decode(re.sub(r"\s+", "", encoded), validate=True)
                require(payload.startswith(b"\x89PNG\r\n\x1a\n"), "Invalid embedded PNG")
                images += 1
    require(images >= 3, "Notebook must embed at least three PNG outputs")
    return len(cells), images


def validate_links():
    count = 0
    for document in DOCS:
        for target in re.findall(r"!?\[[^\]]*\]\(([^\s)]+)(?:\s+\"[^\"]*\")?\)", document.read_text()):
            uri = urlsplit(target.strip("<>"))
            if uri.scheme or uri.netloc or not uri.path:
                continue
            path = (document.parent / unquote(uri.path)).resolve()
            require(path.is_relative_to(ROOT) and path.exists(),
                    f"Broken or external local link in {document.relative_to(ROOT)}: {target}")
            count += 1
    return count


def validate_semantics(paths):
    for path in paths:
        if path.suffix in {".json", ".ipynb"}:
            read_json(path)
    quality = read_json(DERIVED / "quality.json")
    require(quality["clusters"] == 6 and quality["dynamics_enabled"] is False,
            "Unexpected scope or enabled dynamics")
    require(quality["full_local_species_distribution_ready"] is False
            and quality["common_tracer_selection_ready"] is False,
            "Audit claims an approved distribution or matched tracers")
    sources = validate_sources(quality)
    coverage, catalog = validate_tables(quality)
    validate_candidates(coverage, catalog, sources)
    code_cells, images = validate_notebook()
    links = validate_links()
    return code_cells, images, links


def artifact_records(paths):
    return [{"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size,
             "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in paths]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-lock", action="store_true",
                        help="explicitly record current hashes after consistency validation")
    args = parser.parse_args()
    paths = scoped_files()
    code_cells, images, links = validate_semantics(paths)
    records = artifact_records(paths)
    if args.write_lock:
        lock = {"schema_version": 1,
                "data_revision": read_json(DATA / "candidate-hosts.json")["data_revision"],
                "scope": "artifact integrity and consistency; no scientific reproduction",
                "files": records}
        LOCK.write_text(json.dumps(lock, indent=2, allow_nan=False) + "\n")
        action = "Wrote explicitly requested lock for"
    else:
        require(LOCK.is_file(), "Audit lock missing; normal verification never creates or refreshes it")
        lock = read_json(LOCK)
        require(lock.get("schema_version") == 1, "Unsupported audit lock schema")
        require(lock.get("data_revision") == read_json(DATA / "candidate-hosts.json")["data_revision"],
                "Lock data revision differs from candidate records")
        require(lock.get("files") == records,
                "Artifact hashes, sizes, or inventory differ from audit.json; review changes before an explicit --write-lock")
        action = "Verified"
    print(f"{action} {len(records)} artifacts: six clusters, 51 source records, "
          f"{code_cells} executed code cells, {images} embedded PNGs, {links} local links. "
          "No downloads or physical calculations performed.")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError, csv.Error) as error:
        print(f"Data audit verification failed: {error}", file=sys.stderr)
        sys.exit(1)
