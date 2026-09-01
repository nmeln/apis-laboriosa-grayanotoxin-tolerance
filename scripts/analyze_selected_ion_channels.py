#!/usr/bin/env python3
"""Stress-test sodium-channel candidates from the population-selection signal.

Table S6 of Cao et al. (2023) places sodium-transport candidates on eastern
assembly scaffolds 8 and 25.  Cross-assembly mapping identifies likely RefSeq
orthologs as DSC1/NaCP60E, a DEG/ENaC Nach-like channel, and COMMD3-like.

Only DSC1 is homologous enough to the four-domain Nav architecture to map the
experimentally implicated grayanotoxin positions.  Nach is a distinct DEG/ENaC
family channel, so for it this script reports ortholog-level conservation only.
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / ".deps"))

from Bio import SeqIO  # type: ignore  # noqa: E402
from Bio.Align import PairwiseAligner  # type: ignore  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
PROTEOMES = {
    "Apis_laboriosa": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz",
    "Apis_dorsata": ROOT / "genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz",
    "Apis_mellifera": ROOT / "genomes/apis_mellifera/GCF_003254395.2_protein.faa.gz",
    "Apis_cerana": ROOT / "genomes/apis_cerana/GCF_029169275.2_protein.faa.gz",
    "Apis_florea": ROOT / "genomes/apis_florea/GCF_048593485.1_protein.faa.gz",
    "Bombus_terrestris": ROOT / "genomes/bombus_terrestris/GCF_000214255.1_protein.faa.gz",
}
GTX_SITES = [237, 243, 246, 248, 249, 250, 251, 433, 434, 437, 784, 1276, 1463, 1575, 1579, 1586]
LAB_NACH_ACCESSIONS = {"XP_043800019.1", "XP_043800020.1"}
NACH_PATTERN = re.compile(r"sodium channel protein Nach|pickpocket protein 28", re.I)


def records(path: Path):
    with gzip.open(path, "rt") as handle:
        return list(SeqIO.parse(handle, "fasta"))


def aligner() -> PairwiseAligner:
    result = PairwiseAligner()
    result.mode = "global"
    result.match_score = 2.0
    result.mismatch_score = -1.0
    result.open_gap_score = -8.0
    result.extend_gap_score = -0.5
    return result


def map_position(alignment, reference_position_1based: int):
    index = reference_position_1based - 1
    reference_blocks, query_blocks = alignment.aligned
    for (r0, r1), (q0, q1) in zip(reference_blocks, query_blocks):
        if r0 <= index < r1:
            return int(q0 + index - r0) + 1
    return None


def alignment_stats(alignment, reference: str, query: str):
    aligned = matches = 0
    reference_blocks, query_blocks = alignment.aligned
    for (r0, r1), (q0, q1) in zip(reference_blocks, query_blocks):
        left, right = reference[r0:r1], query[q0:q1]
        aligned += len(left)
        matches += sum(a == b for a, b in zip(left, right))
    return aligned, matches, matches / aligned if aligned else 0.0


def main() -> None:
    proteomes = {species: records(path) for species, path in PROTEOMES.items()}
    tool = aligner()

    # DSC1 / NaCP60E comparison.
    dsc1 = {}
    for species, species_records in proteomes.items():
        candidates = [record for record in species_records if "sodium channel protein 60E" in record.description]
        if not candidates:
            raise RuntimeError(f"No DSC1/60E record for {species}")
        dsc1[species] = max(candidates, key=lambda record: len(record.seq))

    lab = dsc1["Apis_laboriosa"]
    lab_seq = str(lab.seq).replace("*", "")
    identity_rows = []
    lab_alignments = {}
    for species, record in dsc1.items():
        sequence = str(record.seq).replace("*", "")
        alignment = tool.align(lab_seq, sequence)[0]
        lab_alignments[species] = alignment
        aligned, matches, identity = alignment_stats(alignment, lab_seq, sequence)
        identity_rows.append(
            {
                "species": species,
                "accession": record.id,
                "length_aa": len(sequence),
                "aligned_aa_to_laboriosa": aligned,
                "matches_to_laboriosa": matches,
                "identity_over_aligned": f"{identity:.6f}",
            }
        )
    with (OUT / "selected_60e_pairwise_identity.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=identity_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(identity_rows)

    with (ROOT / "references/rat_Nav1.4_NP_037310.2.faa").open() as handle:
        rat = next(SeqIO.parse(handle, "fasta"))
    rat_seq = str(rat.seq)
    rat_alignments = {
        species: tool.align(rat_seq, str(record.seq).replace("*", ""))[0]
        for species, record in dsc1.items()
    }
    site_rows = []
    for site in GTX_SITES:
        row = {"rat_nav1.4_position": site, "rat_residue": rat_seq[site - 1]}
        for species, record in dsc1.items():
            sequence = str(record.seq).replace("*", "")
            mapped = map_position(rat_alignments[species], site)
            row[f"{species}_60e_position"] = mapped if mapped is not None else "gap"
            row[f"{species}_60e_residue"] = sequence[mapped - 1] if mapped is not None else "-"
        site_rows.append(row)
    with (OUT / "selected_60e_gtx_site_residues.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=site_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(site_rows)

    # Species-specific DSC1 positions relative to the other Apis consensus.
    apis_comparators = ["Apis_dorsata", "Apis_mellifera", "Apis_cerana", "Apis_florea"]
    lab_site_positions = [map_position(rat_alignments["Apis_laboriosa"], site) for site in GTX_SITES]
    unique_rows = []
    for lab_position, lab_residue in enumerate(lab_seq, start=1):
        other_residues = []
        mapped_positions = {}
        for species in apis_comparators:
            mapped = map_position(lab_alignments[species], lab_position)
            mapped_positions[species] = mapped
            sequence = str(dsc1[species].seq).replace("*", "")
            other_residues.append(sequence[mapped - 1] if mapped is not None else "-")
        if "-" in other_residues or len(set(other_residues)) != 1 or lab_residue == other_residues[0]:
            continue
        nearest_index = min(
            range(len(GTX_SITES)),
            key=lambda index: abs((lab_site_positions[index] or 10**9) - lab_position),
        )
        unique_rows.append(
            {
                "laboriosa_60e_position": lab_position,
                "laboriosa_residue": lab_residue,
                "other_apis_consensus": other_residues[0],
                "nearest_rat_gtx_site": GTX_SITES[nearest_index],
                "distance_to_mapped_gtx_site_aa": lab_position - int(lab_site_positions[nearest_index] or 0),
                **{f"{species}_position": mapped_positions[species] for species in apis_comparators},
            }
        )
    with (OUT / "selected_60e_laboriosa_unique_positions.tsv").open("w", newline="") as handle:
        fieldnames = list(unique_rows[0].keys()) if unique_rows else [
            "laboriosa_60e_position",
            "laboriosa_residue",
            "other_apis_consensus",
            "nearest_rat_gtx_site",
            "distance_to_mapped_gtx_site_aa",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(unique_rows)

    # Nach is a DEG/ENaC, not a four-domain Nav.  Select the closest annotated
    # family member in each proteome to the mapped A. laboriosa locus.
    lab_nach_candidates = [record for record in proteomes["Apis_laboriosa"] if record.id in LAB_NACH_ACCESSIONS]
    lab_nach = max(lab_nach_candidates, key=lambda record: len(record.seq))
    lab_nach_seq = str(lab_nach.seq).replace("*", "")
    nach_rows = []
    for species, species_records in proteomes.items():
        candidates = [record for record in species_records if NACH_PATTERN.search(record.description)]
        best = None
        for record in candidates:
            sequence = str(record.seq).replace("*", "")
            alignment = tool.align(lab_nach_seq, sequence)[0]
            aligned, matches, identity = alignment_stats(alignment, lab_nach_seq, sequence)
            candidate = (tool.score(lab_nach_seq, sequence), identity, aligned, record, matches)
            if best is None or candidate[:3] > best[:3]:
                best = candidate
        if best is None:
            raise RuntimeError(f"No Nach/PPK candidate for {species}")
        _score, identity, aligned, record, matches = best
        nach_rows.append(
            {
                "species": species,
                "best_accession": record.id,
                "length_aa": len(record.seq),
                "aligned_aa": aligned,
                "matches": matches,
                "identity_to_laboriosa_candidate": f"{identity:.6f}",
                "description": record.description,
            }
        )
    with (OUT / "selected_nach_best_family_match.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=nach_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(nach_rows)

    invariant_sites = sum(
        len({row[f"{species}_60e_residue"] for species in dsc1}) == 1 for row in site_rows
    )
    print(f"DSC1/60E GTX-site homologues invariant across species: {invariant_sites}/{len(GTX_SITES)}")
    print(f"A. laboriosa DSC1/60E-specific consensus substitutions: {len(unique_rows)}")
    if unique_rows:
        nearest = min(abs(int(row["distance_to_mapped_gtx_site_aa"])) for row in unique_rows)
        print(f"Nearest unique DSC1/60E substitution to a mapped GTX site: {nearest} aa")
    print(f"A. laboriosa mapped Nach candidate: {lab_nach.id} ({len(lab_nach.seq)} aa; DEG/ENaC family)")


if __name__ == "__main__":
    main()
