#!/usr/bin/env python3
"""Compare annotated Para splice diversity across Apis species.

The central question is whether the A. laboriosa domain-III S3/S4 alternative
segment is species-specific.  All isoforms are independently aligned to rat Nav1.4
so homologous positions can be compared without assuming equal protein lengths.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".deps"))

from Bio import SeqIO  # type: ignore  # noqa: E402
from Bio.Align import PairwiseAligner  # type: ignore  # noqa: E402


SOURCES = {
    "Apis_laboriosa": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz",
    "Apis_dorsata": ROOT / "genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz",
    "Apis_mellifera": ROOT / "genomes/apis_mellifera/GCF_003254395.2_protein.faa.gz",
    "Apis_cerana": ROOT / "genomes/apis_cerana/GCF_029169275.2_protein.faa.gz",
    "Apis_florea": ROOT / "genomes/apis_florea/GCF_048593485.1_protein.faa.gz",
    "Bombus_terrestris": ROOT / "genomes/bombus_terrestris/GCF_000214255.1_protein.faa.gz",
}
PARA_RE = re.compile(r"sodium (?:voltage-gated channel paralytic|channel protein (?:para|paralytic))", re.I)
GTX_SITES = {237, 243, 246, 248, 249, 250, 251, 433, 434, 437, 784, 1276, 1463, 1575, 1579, 1586}
S3S4_START = 1100
S3S4_END = 1140


def read_records(path: Path):
    with gzip.open(path, "rt") as handle:
        return [r for r in SeqIO.parse(handle, "fasta") if PARA_RE.search(r.description)]


def make_aligner():
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -8.0
    a.extend_gap_score = -0.5
    return a


def ref_to_query_map(alignment):
    result = {}
    ref_blocks, query_blocks = alignment.aligned
    for (r0, r1), (q0, q1) in zip(ref_blocks, query_blocks):
        for offset in range(r1 - r0):
            result[int(r0 + offset)] = int(q0 + offset)
    return result


def main():
    with (ROOT / "references/rat_Nav1.4_NP_037310.2.faa").open() as handle:
        rat_record = next(SeqIO.parse(handle, "fasta"))
    rat = str(rat_record.seq).replace("*", "")
    a = make_aligner()

    uniprot = json.loads((ROOT / "references/uniprot_P15390.json").read_text())
    tm = []
    for feature in uniprot["features"]:
        if feature["type"] == "Transmembrane":
            tm.append(
                (
                    int(feature["location"]["start"]["value"]),
                    int(feature["location"]["end"]["value"]),
                    feature.get("description", "transmembrane"),
                )
            )

    def topology(rat_pos):
        for start, end, label in tm:
            if start <= rat_pos <= end:
                return label
        return "non-transmembrane"

    summary_rows = []
    haplotype_rows = []
    tm_variable_rows = []
    all_species_haplotypes = {}

    for species, path in SOURCES.items():
        records = read_records(path)
        if not records:
            raise SystemExit(f"No Para records in {path}")
        records = [r for r in records if len(r.seq) >= 1500]
        reference = max(records, key=lambda r: len(r.seq))
        ref = str(reference.seq).replace("*", "")

        rat_maps = {}
        haplotypes = defaultdict(list)
        for record in records:
            seq = str(record.seq).replace("*", "")
            mapping = ref_to_query_map(a.align(rat, seq)[0])
            rat_maps[record.id] = mapping
            hap = "".join(seq[mapping[pos - 1]] if pos - 1 in mapping else "-" for pos in range(S3S4_START, S3S4_END + 1))
            haplotypes[hap].append(record.id)

        all_species_haplotypes[species] = set(haplotypes)
        for hap, accessions in sorted(haplotypes.items(), key=lambda item: (-len(item[1]), item[0])):
            haplotype_rows.append(
                {
                    "species": species,
                    "rat_region": f"{S3S4_START}-{S3S4_END}",
                    "haplotype": hap,
                    "isoform_count": len(accessions),
                    "accessions": ",".join(accessions),
                }
            )

        variable_rat_positions = []
        for rat_pos in range(1, len(rat) + 1):
            states = []
            for record in records:
                seq = str(record.seq).replace("*", "")
                q = rat_maps[record.id].get(rat_pos - 1)
                states.append(seq[q] if q is not None else "-")
            if len(set(states)) > 1:
                variable_rat_positions.append(rat_pos)
                top = topology(rat_pos)
                if top != "non-transmembrane":
                    tm_variable_rows.append(
                        {
                            "species": species,
                            "rat_position": rat_pos,
                            "rat_residue": rat[rat_pos - 1],
                            "topology": top,
                            "known_gtx_site": rat_pos in GTX_SITES,
                            "state_counts": ";".join(f"{k}:{v}" for k, v in sorted(Counter(states).items())),
                        }
                    )

        tm_positions = [p for p in variable_rat_positions if topology(p) != "non-transmembrane"]
        summary_rows.append(
            {
                "species": species,
                "full_length_para_isoforms": len(records),
                "longest_accession": reference.id,
                "longest_length_aa": len(ref),
                "variable_rat_mapped_positions": len(variable_rat_positions),
                "variable_transmembrane_positions": len(tm_positions),
                "variable_known_gtx_positions": sum(p in GTX_SITES for p in variable_rat_positions),
                "diii_s3_s4_haplotypes": len(haplotypes),
                "diii_s3_s4_isoforms_with_deletion": sum("-" in h for h in haplotypes for _ in haplotypes[h]),
            }
        )

    lab_haps = all_species_haplotypes["Apis_laboriosa"]
    overlap_rows = []
    for species, haps in all_species_haplotypes.items():
        overlap = lab_haps & haps
        overlap_rows.append(
            {
                "species": species,
                "laboriosa_haplotypes": len(lab_haps),
                "species_haplotypes": len(haps),
                "shared_exact_haplotypes": len(overlap),
                "all_laboriosa_haplotypes_present": lab_haps <= haps,
                "shared_haplotypes": " | ".join(sorted(overlap)),
            }
        )

    out = ROOT / "results"
    for filename, rows in [
        ("para_isoform_species_summary.tsv", summary_rows),
        ("para_diii_s3_s4_haplotypes.tsv", haplotype_rows),
        ("para_isoform_transmembrane_variation.tsv", tm_variable_rows),
        ("para_diii_s3_s4_haplotype_overlap.tsv", overlap_rows),
    ]:
        with (out / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    for row in summary_rows:
        print(row)
    print("\nExact DIII S3/S4 haplotype overlap with A. laboriosa:")
    for row in overlap_rows:
        print(row["species"], row["shared_exact_haplotypes"], row["all_laboriosa_haplotypes_present"])


if __name__ == "__main__":
    main()
