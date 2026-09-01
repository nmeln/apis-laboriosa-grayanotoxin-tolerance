#!/usr/bin/env python3
"""Assert the quantitative claims reported by the comparative addendum."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def rows(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def one(name: str, **criteria: str) -> dict[str, str]:
    matches = [
        row for row in rows(name)
        if all(row[key] == value for key, value in criteria.items())
    ]
    if len(matches) != 1:
        raise AssertionError(f"Expected one {name} row matching {criteria}, found {len(matches)}")
    return matches[0]


def close(actual: str | float, expected: float, tolerance: float = 1e-9) -> None:
    if not math.isclose(float(actual), expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"Expected {expected}, observed {actual}")


def main() -> None:
    summary = json.loads((RESULTS / "analysis_summary.json").read_text())
    expected = {
        "orthogroups_total": 10003,
        "complete_single_copy_orthogroups": 8273,
        "callable_amino_acid_sites": 4902144,
        "laboriosa_bombus_strict_sites": 2024,
        "laboriosa_bombus_strict_site_orthogroups": 1311,
        "para_complete_single_copy_orthogroups": 1,
        "para_laboriosa_bombus_strict_sites": 0,
        "strict_copy_number_hits_laboriosa_bombus": 1,
        "alignment_source": "OrthoFinder final MultipleSequenceAlignments",
    }
    for key, value in expected.items():
        if summary[key] != value:
            raise AssertionError(f"{key}: expected {value!r}, observed {summary[key]!r}")

    primary = rows("primary_proteome_summary_current.tsv")
    if sum(int(row["primary_gene_representatives"]) for row in primary) != 60216:
        raise AssertionError("Unexpected representative-protein total")
    overall: dict[str, str] = {}
    for line in (RESULTS / "orthofinder_key/Statistics_Overall.tsv").read_text().splitlines():
        fields = line.split("\t")
        if len(fields) == 2:
            overall[fields[0]] = fields[1]
    if overall.get("Number of genes") != "60216":
        raise AssertionError("Unexpected OrthoFinder gene total")
    if overall.get("Number of genes in orthogroups") != "59381":
        raise AssertionError("Unexpected assigned-gene total")
    if overall.get("Number of orthogroups with all species present") != "8683":
        raise AssertionError("Unexpected all-species orthogroup total")

    dorsata_control = one("strict_site_aggregate.tsv", focal_apis="Apis_dorsata")
    if dorsata_control["strict_sites"] != "2036":
        raise AssertionError("Unexpected A. dorsata internal-control site count")

    para = json.loads((RESULTS / "para_strict_sharing_summary.json").read_text())
    if para["callable_amino_acid_sites"] != 2050:
        raise AssertionError("Unexpected callable Para residue count")
    if para["laboriosa_bombus_strict_shared_sites"] != 0:
        raise AssertionError("Unexpected strict Para state sharing")

    barrier = one(
        "candidate_site_enrichment.tsv",
        focal_apis="Apis_laboriosa",
        candidate_set="barrier_or_detox",
    )
    core = one(
        "candidate_site_enrichment.tsv",
        focal_apis="Apis_laboriosa",
        candidate_set="core_toxicokinetic",
    )
    close(barrier["matched_density_difference_per_1000_sites"], -0.061762679504032775)
    if barrier["candidate_strict_sites"] != "58":
        raise AssertionError("Unexpected barrier-or-detox strict-site count")
    close(barrier["matched_permutation_p_two_sided"], 0.4458277086145693)
    close(barrier["matched_permutation_bh_q"], 0.5390286041253494)
    close(core["matched_density_difference_per_1000_sites"], -0.06902361786650257)
    if core["candidate_strict_sites"] != "32":
        raise AssertionError("Unexpected core-toxicokinetic strict-site count")
    close(core["matched_permutation_p_two_sided"], 0.48512574371281436)
    close(core["matched_permutation_bh_q"], 0.5390286041253494)

    copy_hits = [
        row for row in rows("strict_copy_number_hits.tsv")
        if row["focal_apis"] == "Apis_laboriosa"
    ]
    if len(copy_hits) != 1 or copy_hits[0]["orthogroup"] != "OG0000196":
        raise AssertionError("Unexpected laboriosa plus Bombus copy-number result")
    if copy_hits[0]["candidate_any"] != "False":
        raise AssertionError("Copy-number hit unexpectedly carries a candidate annotation")

    exchanger = one(
        "candidate_category_fisher.tsv",
        focal_apis="Apis_laboriosa",
        category="sodium_hydrogen_exchange",
    )
    if exchanger["category_strict_sites"] != "4":
        raise AssertionError("Unexpected exchanger strict-site count")
    close(exchanger["site_fisher_p_two_sided"], 0.07586731170946386)
    close(exchanger["site_fisher_bh_q"], 0.3652870563789001)

    panel = {int(row["laboriosa_position_1based"]): row for row in rows("nhe3_bombus_panel_summary.tsv")}
    if set(panel) != {44, 159, 232, 353}:
        raise AssertionError("Unexpected exchanger positions")
    if panel[159]["additional_bombus_sequences_matching_focal_state"] != "7":
        raise AssertionError("Unexpected Bombus support at S159")
    if panel[232]["additional_bombus_sequences_matching_focal_state"] != "8":
        raise AssertionError("Unexpected Bombus support at I232")

    external = {int(row["laboriosa_position_1based"]): row for row in rows("nhe3_external_bee_summary.tsv")}
    for position in (44, 159, 232):
        if external[position]["external_bees_matching_focal_state"] != "0":
            raise AssertionError(f"Unexpected external-bee focal state at {position}")
    if external[353]["external_bees_matching_focal_state"] != "5":
        raise AssertionError("Unexpected external-bee state at 353")

    expression = rows("nhe3_constitutive_expression_summary.tsv")
    if len(expression) != 1:
        raise AssertionError("Unexpected exchanger expression summary row count")
    close(expression[0]["laboriosa_fpkm"], 20.54)
    close(expression[0]["dorsata_fpkm"], 24.97)
    close(expression[0]["laboriosa_to_dorsata_fpkm_ratio"], 0.822587, 1e-6)

    transcript = rows("nhe3_transcript_residue_validation.tsv")
    if len(transcript) != 8 or any(row["state_matches"] != "True" for row in transcript):
        raise AssertionError("The eight transcript-validated exchanger states do not all match")

    para_sites = rows("current_bombus_para_gtx_sites.tsv")
    if len(para_sites) != 16 or any(row["current_matches_apis_laboriosa"] != "True" for row in para_sites):
        raise AssertionError("Current Bombus does not match laboriosa at all 16 mapped Para sites")

    print("comparative addendum claims validated")


if __name__ == "__main__":
    main()
