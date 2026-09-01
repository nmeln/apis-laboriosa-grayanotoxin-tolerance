#!/usr/bin/env python3
"""Assert the README's headline findings from generated result tables."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def rows(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def one(items: list[dict[str, str]], key: str, value: str) -> dict[str, str]:
    matches = [item for item in items if item[key] == value]
    assert len(matches) == 1, f"Expected one {key}={value}, found {len(matches)}"
    return matches[0]


def close(observed: str, expected: float, tolerance: float = 1e-6) -> None:
    assert math.isclose(float(observed), expected, rel_tol=tolerance, abs_tol=tolerance), (
        observed,
        expected,
    )


def validate_para_sites() -> None:
    table = rows("gtx_target_residues.tsv")
    assert len(table) == 16
    core = [
        "Apis_laboriosa_residue",
        "Apis_dorsata_residue",
        "Apis_mellifera_refseq_residue",
        "Apis_cerana_residue",
        "Apis_florea_residue",
        "Bombus_terrestris_old_refseq_residue",
    ]
    for row in table:
        assert len({row[column] for column in core}) == 1, row
    print("PASS Para: 16 of 16 mapped sites invariant across six reference bees")


def validate_para_sequence() -> None:
    identity = one(
        rows("laboriosa_pairwise_identity.tsv"),
        "comparison",
        "Apis_laboriosa_vs_Apis_dorsata",
    )
    close(identity["identity_over_aligned"], 0.999025)
    unique = rows("laboriosa_unique_genomewide.tsv")
    observed = {
        (int(row["laboriosa_position"]), row["laboriosa_residue"], row["topology"])
        for row in unique
    }
    assert observed == {(95, "N", "non-transmembrane"), (465, "V", "unmapped")}
    print("PASS Para: 99.9025% identity to A. dorsata; unique residues N95 and V465")


def validate_isoforms() -> None:
    summary = one(rows("para_isoform_species_summary.tsv"), "species", "Apis_laboriosa")
    assert int(summary["full_length_para_isoforms"]) == 17
    assert int(summary["variable_known_gtx_positions"]) == 0
    assert int(summary["diii_s3_s4_haplotypes"]) == 3
    overlap = rows("para_diii_s3_s4_haplotype_overlap.tsv")
    for species in ["Apis_dorsata", "Apis_mellifera", "Apis_cerana", "Apis_florea"]:
        row = one(overlap, "species", species)
        assert row["all_laboriosa_haplotypes_present"] == "True"
    print("PASS Para isoforms: 17 laboriosa isoforms; all three haplotypes shared")


def validate_transcript() -> None:
    row = one(
        rows("para_transcriptome_alignment_summary.tsv"),
        "species",
        "Apis_laboriosa",
    )
    assert int(row["reference_cds_length_nt"]) == 6129
    close(row["reference_coverage"], 1.0)
    assert int(row["single_base_mismatches"]) == 4
    assert int(row["nonsynonymous_mismatches"]) == 0
    print("PASS Para transcript: complete CDS; four mismatches; zero nonsynonymous")


def validate_population_mapping() -> None:
    para = rows("para_eastern_scaffold_probe_hits.tsv")
    assert len(para) == 1
    assert para[0]["eastern_accession"] == "GWHAOTM00000105"
    assert para[0]["header_metadata"].startswith("OriSeqID=scaffold_105")
    assert int(para[0]["probe_hits"]) == int(para[0]["total_probes"]) == 29
    candidates = rows("population_sodium_candidate_crossmap.tsv")
    assert {row["likely_refseq_gene"] for row in candidates} == {
        "LOC122714529",
        "LOC122714475",
        "LOC122718769",
    }
    assert {row["eastern_original_scaffold"] for row in candidates} == {
        "scaffold_8",
        "scaffold_25",
    }
    print("PASS population mapping: Para on scaffold 105; candidates on 8 and 25")


def validate_dsc1() -> None:
    table = rows("selected_60e_gtx_site_residues.tsv")
    assert len(table) == 16
    residue_columns = [column for column in table[0] if column.endswith("_60e_residue")]
    assert len(residue_columns) == 6
    for row in table:
        assert len({row[column] for column in residue_columns}) == 1, row
    print("PASS DSC1/60E: 16 of 16 mapped homologous sites invariant")


def validate_chitin_qc() -> None:
    table = rows("chitin_synthase_locus_qc.tsv")
    assert len(table) == 3
    assert sum(row["full_length_chitin_synthase_like"] == "True" for row in table) == 2
    fragment = one(table, "gene_name", "LOC122719633")
    assert int(fragment["max_protein_length_aa"]) == 114
    assert "nested in" in fragment["overlap_relation"]
    print("PASS chitin QC: two full-length loci; one nested 114-aa fragment")


def validate_expression() -> None:
    table = rows("constitutive_expression_contrasts.tsv")
    expected = {
        "peritrophic_matrix_chitin_binding": (492.587966, 383.627193),
        "chitin_synthase": (19.428144, 46.469687),
        "ABCC_MRP": (3.186328, 3.024383),
        "organic_anion_transport": (3.880838, 5.087768),
        "ABCB_P_glycoprotein": (0.742297, 0.604117),
        "ABCG_half_transporter": (0.921179, 0.806465),
        "Para": (0.875786, 1.079238),
    }
    for category, (fpkm, cpm) in expected.items():
        row = one(table, "category", category)
        close(row["laboriosa_to_dorsata_fpkm_ratio"], fpkm)
        close(row["laboriosa_to_dorsata_cpm_ratio"], cpm)
        assert "descriptive only" in row["design_warning"]
    print("PASS expression: README ratios match descriptive pooled-worker outputs")


def main() -> None:
    validate_para_sites()
    validate_para_sequence()
    validate_isoforms()
    validate_transcript()
    validate_population_mapping()
    validate_dsc1()
    validate_chitin_qc()
    validate_expression()
    print("All headline claims validated from generated tables.")


if __name__ == "__main__":
    main()
