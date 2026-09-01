#!/usr/bin/env python3
"""Map NHE3 into the existing pooled whole-worker transcriptomes.

This reuses the repository's exact 31-mer probe method.  The comparison has
one untreated pooled whole-worker library per species, so the output is a
descriptive check for a gross constitutive difference, not differential
expression.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import sys
from pathlib import Path


def load_expression_module(project: Path):
    path = project / "scripts/analyze_constitutive_transcript_expression.py"
    spec = importlib.util.spec_from_file_location("repo_expression", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    module = load_expression_module(args.project)
    module.CATEGORIES = {"NHE3": r"sodium/hydrogen exchanger 3(?:,|\b)"}
    module.PATTERNS = {"NHE3": re.compile(module.CATEGORIES["NHE3"], re.I)}

    rows = []
    for species, paths in module.SPECIES.items():
        rows.extend(module.analyze_species(species, paths))
    if len(rows) != 2:
        raise RuntimeError(f"Expected one NHE3 annotation per species, found {len(rows)} rows")

    detail_path = args.output / "nhe3_constitutive_expression.tsv"
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    by_species = {row["species"]: row for row in rows}
    lab = by_species["Apis_laboriosa"]
    dor = by_species["Apis_dorsata"]
    if not lab["mapped"] or not dor["mapped"]:
        raise RuntimeError("NHE3 did not pass the exact-probe mapping threshold in both species")
    lab_fpkm, dor_fpkm = float(lab["fpkm"]), float(dor["fpkm"])
    lab_cpm, dor_cpm = float(lab["cpm"]), float(dor["cpm"])
    summary = [{
        "laboriosa_fpkm": f"{lab_fpkm:.6f}",
        "dorsata_fpkm": f"{dor_fpkm:.6f}",
        "laboriosa_to_dorsata_fpkm_ratio": f"{lab_fpkm / dor_fpkm:.6f}" if dor_fpkm else "NA",
        "laboriosa_cpm": f"{lab_cpm:.6f}",
        "dorsata_cpm": f"{dor_cpm:.6f}",
        "laboriosa_to_dorsata_cpm_ratio": f"{lab_cpm / dor_cpm:.6f}" if dor_cpm else "NA",
        "design_warning": "one pooled untreated whole-worker library per species; descriptive only",
    }]
    with (args.output / "nhe3_constitutive_expression_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary)
    print(summary[0])


if __name__ == "__main__":
    main()
