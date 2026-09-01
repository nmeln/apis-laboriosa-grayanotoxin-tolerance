#!/usr/bin/env python3
"""Map known grayanotoxin-relevant Nav sites onto the current Bombus reference."""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path


ADDENDUM_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ADDENDUM_ROOT.parent
RESULTS = ADDENDUM_ROOT / "results"
CURRENT_PROTEOME = (
    ADDENDUM_ROOT
    / "inputs/bombus_terrestris_current/GCF_910591885.1_iyBomTerr1.2_protein.faa.gz"
)


def load_nav_module():
    path = PROJECT_ROOT / "scripts/analyze_nav.py"
    spec = importlib.util.spec_from_file_location("repo_analyze_nav", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    nav = load_nav_module()
    current, records = nav.longest_para(CURRENT_PROTEOME)
    old, _ = nav.longest_para(nav.PROTEOME_SOURCES["Bombus_terrestris_old_refseq"])
    rat = nav.read_fasta(PROJECT_ROOT / "references/rat_Nav1.4_NP_037310.2.faa")[0]
    aligner = nav.make_aligner()
    rat_sequence = str(rat.seq)
    current_sequence = str(current.seq).replace("*", "")
    old_sequence = str(old.seq).replace("*", "")
    rat_alignment = aligner.align(rat_sequence, current_sequence)[0]
    old_current_alignment = aligner.align(old_sequence, current_sequence)[0]
    aligned, matches, identity = nav.alignment_stats(old_current_alignment, old_sequence, current_sequence)

    rows = []
    for target in nav.TARGETS:
        current_position = nav.map_position(rat_alignment, target.position)
        current_residue = current_sequence[current_position - 1] if current_position else "-"
        old_row = next(
            row
            for row in csv.DictReader(
                (PROJECT_ROOT / "results/gtx_target_residues.tsv").open(), delimiter="\t"
            )
            if int(row["rat_position"]) == target.position
        )
        rows.append(
            {
                "rat_position": target.position,
                "rat_residue": target.expected,
                "region": target.region,
                "current_bombus_position": current_position or "gap",
                "current_bombus_residue": current_residue,
                "old_bombus_residue": old_row["Bombus_terrestris_old_refseq_residue"],
                "current_matches_old": current_residue == old_row["Bombus_terrestris_old_refseq_residue"],
                "current_matches_apis_laboriosa": current_residue == old_row["Apis_laboriosa_residue"],
                "current_matches_apis_mellifera": current_residue == old_row["Apis_mellifera_refseq_residue"],
            }
        )

    RESULTS.mkdir(parents=True, exist_ok=True)
    with (RESULTS / "current_bombus_para_gtx_sites.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    summary = [
        {
            "current_accession": current.id,
            "current_length_aa": len(current_sequence),
            "matching_para_models_in_current_proteome": len(records),
            "old_accession": old.id,
            "old_length_aa": len(old_sequence),
            "old_current_aligned_aa": aligned,
            "old_current_identity_over_aligned": f"{identity:.8f}",
            "known_gtx_sites_checked": len(rows),
            "known_gtx_sites_matching_old": sum(bool(row["current_matches_old"]) for row in rows),
            "known_gtx_sites_matching_laboriosa": sum(bool(row["current_matches_apis_laboriosa"]) for row in rows),
            "known_gtx_sites_matching_mellifera": sum(bool(row["current_matches_apis_mellifera"]) for row in rows),
        }
    ]
    with (RESULTS / "current_bombus_para_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary)
    print(summary[0])


if __name__ == "__main__":
    main()
