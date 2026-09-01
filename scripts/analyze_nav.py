#!/usr/bin/env python3
"""Compare bee Para sodium channels at experimentally mapped GTX-contact sites.

The reference coordinates are from rat skeletal-muscle Nav1.4 (mu1/SCN4A,
NP_037310.2).  Pairwise affine-gap alignments map those coordinates to each
bee sequence.  The output is deliberately residue-level: whole-protein
identity is not treated as evidence for or against toxin binding.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / ".deps"))

from Bio import SeqIO  # type: ignore  # noqa: E402
from Bio.Align import PairwiseAligner  # type: ignore  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class Target:
    position: int
    expected: str
    region: str
    evidence: str


TARGETS = [
    Target(237, "K", "DI S4-S5 linker", "alanine scan: moderate GTX effect"),
    Target(243, "L", "DI S4-S5 linker", "synergistic determinant"),
    Target(246, "S", "DI S4-S5 linker", "alanine scan: moderate GTX effect"),
    Target(248, "K", "DI S4-S5 linker", "alanine scan: moderate GTX effect"),
    Target(249, "K", "DI S4-S5 linker", "alanine scan: moderate GTX effect"),
    Target(250, "L", "DI S4-S5 linker", "alanine scan: moderate GTX effect"),
    Target(251, "S", "DI S4-S5 linker", "isoform sensitivity determinant"),
    Target(433, "I", "DI S6", "Lys substitution abolished GTX response"),
    Target(434, "N", "DI S6", "Lys substitution abolished GTX response"),
    Target(437, "L", "DI S6", "Lys substitution abolished GTX response"),
    Target(784, "N", "DII S6", "binding/unbinding determinant"),
    Target(1276, "S", "DIII S6", "binding/unbinding determinant"),
    Target(1463, "T", "DIV S4-S5 linker", "alanine scan: moderate GTX effect"),
    Target(1575, "I", "DIV S6", "Ala substitution abolished GTX effect"),
    Target(1579, "F", "DIV S6", "access-gate kinetics determinant"),
    Target(1586, "Y", "DIV S6", "GTX-affinity determinant; Lys abolished effect"),
]


PROTEOME_SOURCES = {
    "Apis_laboriosa": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz",
    "Apis_dorsata": ROOT / "genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz",
    "Apis_mellifera_refseq": ROOT / "genomes/apis_mellifera/GCF_003254395.2_protein.faa.gz",
    "Apis_cerana": ROOT / "genomes/apis_cerana/GCF_029169275.2_protein.faa.gz",
    "Apis_florea": ROOT / "genomes/apis_florea/GCF_048593485.1_protein.faa.gz",
    "Bombus_terrestris_old_refseq": ROOT / "genomes/bombus_terrestris/GCF_000214255.1_protein.faa.gz",
}

PARA_PATTERN = re.compile(
    r"sodium (?:voltage-gated channel paralytic|channel protein (?:para|paralytic))",
    re.IGNORECASE,
)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_fasta(path: Path):
    with open_text(path) as handle:
        return list(SeqIO.parse(handle, "fasta"))


def longest_para(path: Path):
    records = [r for r in read_fasta(path) if PARA_PATTERN.search(r.description)]
    if not records:
        raise RuntimeError(f"No Para/Nav record found in {path}")
    return max(records, key=lambda r: len(r.seq)), records


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -8.0
    aligner.extend_gap_score = -0.5
    return aligner


def map_position(alignment, ref_pos_1based: int):
    """Map a 1-based reference residue through a pairwise alignment."""
    ref_index = ref_pos_1based - 1
    ref_blocks, query_blocks = alignment.aligned
    for (r0, r1), (q0, q1) in zip(ref_blocks, query_blocks):
        if r0 <= ref_index < r1:
            if (r1 - r0) != (q1 - q0):
                raise AssertionError("Aligned blocks must be ungapped diagonals")
            return int(q0 + (ref_index - r0)) + 1
    return None


def map_query_position(alignment, query_pos_1based: int):
    """Map a 1-based query residue back to the reference."""
    query_index = query_pos_1based - 1
    ref_blocks, query_blocks = alignment.aligned
    for (r0, r1), (q0, q1) in zip(ref_blocks, query_blocks):
        if q0 <= query_index < q1:
            if (r1 - r0) != (q1 - q0):
                raise AssertionError("Aligned blocks must be ungapped diagonals")
            return int(r0 + (query_index - q0)) + 1
    return None


def alignment_stats(alignment, ref: str, query: str):
    ref_blocks, query_blocks = alignment.aligned
    aligned = 0
    matches = 0
    for (r0, r1), (q0, q1) in zip(ref_blocks, query_blocks):
        rs = ref[r0:r1]
        qs = query[q0:q1]
        aligned += len(rs)
        matches += sum(a == b for a, b in zip(rs, qs))
    return aligned, matches, matches / aligned if aligned else 0.0


def seq_window(seq: str, pos_1based: int | None, radius: int = 8):
    if pos_1based is None:
        return "-"
    i = pos_1based - 1
    lo = max(0, i - radius)
    hi = min(len(seq), i + radius + 1)
    return f"{lo + 1}-{hi}:{seq[lo:hi]}"


def main() -> None:
    rat_record = read_fasta(ROOT / "references/rat_Nav1.4_NP_037310.2.faa")[0]
    rat = str(rat_record.seq)
    for target in TARGETS:
        observed = rat[target.position - 1]
        if observed != target.expected:
            raise RuntimeError(
                f"Reference mismatch at {target.position}: expected {target.expected}, got {observed}"
            )

    selected = {}
    all_counts = {}
    all_para_records = {}
    for species, path in PROTEOME_SOURCES.items():
        record, records = longest_para(path)
        selected[species] = record
        all_counts[species] = len(records)
        all_para_records[species] = records

    functional = {r.id: r for r in read_fasta(ROOT / "references/bee_nav_functional_panel.faa")}
    selected.update(
        {
            "Apis_mellifera_functional": functional["AMB38675.1"],
            "Bombus_terrestris_functional": functional["XP_012167116.1"],
            "Apis_dorsata_functional": functional["XP_006613070.1"],
            "Apis_florea_functional": functional["XP_012347667.1"],
        }
    )

    aligner = make_aligner()
    rat_alignments = {}
    summary_rows = []
    for species, record in selected.items():
        seq = str(record.seq).replace("*", "")
        aln = aligner.align(rat, seq)[0]
        rat_alignments[species] = aln
        aligned, matches, identity = alignment_stats(aln, rat, seq)
        summary_rows.append(
            {
                "taxon_or_sequence": species,
                "accession": record.id,
                "length_aa": len(seq),
                "matching_para_models_in_proteome": all_counts.get(species, "curated panel"),
                "aligned_to_rat_aa": aligned,
                "identity_to_rat_over_aligned": f"{identity:.5f}",
                "x_residues": seq.count("X"),
                "description": record.description,
            }
        )

    with (OUT / "sequence_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    target_rows = []
    for target in TARGETS:
        row = {
            "rat_position": target.position,
            "rat_residue": target.expected,
            "region": target.region,
            "functional_evidence": target.evidence,
        }
        for species, record in selected.items():
            seq = str(record.seq).replace("*", "")
            mapped = map_position(rat_alignments[species], target.position)
            row[f"{species}_position"] = mapped if mapped is not None else "gap"
            row[f"{species}_residue"] = seq[mapped - 1] if mapped is not None else "-"
        target_rows.append(row)

    with (OUT / "gtx_target_residues.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=target_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(target_rows)

    isoform_rows = []
    full_length_records = {
        species: [r for r in records if len(r.seq) >= 1500]
        for species, records in all_para_records.items()
    }
    for target in TARGETS:
        row = {
            "rat_position": target.position,
            "rat_residue": target.expected,
            "region": target.region,
        }
        for species, records in full_length_records.items():
            residues = set()
            mapped_count = 0
            for record in records:
                seq = str(record.seq).replace("*", "")
                aln = aligner.align(rat, seq)[0]
                mapped = map_position(aln, target.position)
                if mapped is not None:
                    mapped_count += 1
                    residues.add(seq[mapped - 1])
            row[f"{species}_full_length_isoforms"] = len(records)
            row[f"{species}_mapped_isoforms"] = mapped_count
            row[f"{species}_residue_set"] = ",".join(sorted(residues)) or "-"
        isoform_rows.append(row)
    with (OUT / "gtx_target_isoform_residue_sets.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=isoform_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(isoform_rows)

    lab_seq = str(selected["Apis_laboriosa"].seq).replace("*", "")
    lab_target_positions = {
        target.position: map_position(rat_alignments["Apis_laboriosa"], target.position)
        for target in TARGETS
    }

    windows_lines = []
    for target in TARGETS:
        windows_lines.append(
            f"Rat {target.expected}{target.position} | {target.region} | {target.evidence}"
        )
        for species, record in selected.items():
            seq = str(record.seq).replace("*", "")
            mapped = map_position(rat_alignments[species], target.position)
            residue = seq[mapped - 1] if mapped is not None else "-"
            windows_lines.append(
                f"  {species}\t{record.id}\t{residue}{mapped or '-'}\t{seq_window(seq, mapped)}"
            )
    (OUT / "gtx_target_windows.txt").write_text("\n".join(windows_lines) + "\n")

    # Find A. laboriosa-specific substitutions within +/- 12 aa of mapped GTX sites.
    # A site is called unique only when every non-laboriosa Apis sequence maps there,
    # all non-laboriosa Apis agree, and A. laboriosa differs from that consensus.
    apis_comparators = [
        "Apis_dorsata",
        "Apis_mellifera_functional",
        "Apis_cerana",
        "Apis_florea",
    ]
    lab_alignments = {}
    for species in apis_comparators + ["Bombus_terrestris_functional"]:
        query = str(selected[species].seq).replace("*", "")
        lab_alignments[species] = aligner.align(lab_seq, query)[0]

    lab_summary_rows = []
    for species, aln in lab_alignments.items():
        query = str(selected[species].seq).replace("*", "")
        aligned, matches, identity = alignment_stats(aln, lab_seq, query)
        lab_summary_rows.append(
            {
                "comparison": f"Apis_laboriosa_vs_{species}",
                "laboriosa_accession": selected["Apis_laboriosa"].id,
                "query_accession": selected[species].id,
                "aligned_aa": aligned,
                "matches": matches,
                "identity_over_aligned": f"{identity:.6f}",
                "laboriosa_coverage": f"{aligned / len(lab_seq):.6f}",
                "query_coverage": f"{aligned / len(query):.6f}",
            }
        )
    with (OUT / "laboriosa_pairwise_identity.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=lab_summary_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(lab_summary_rows)

    near_positions = set()
    for mapped in lab_target_positions.values():
        if mapped is not None:
            near_positions.update(range(max(1, mapped - 12), min(len(lab_seq), mapped + 12) + 1))

    unique_rows = []
    for lab_pos in sorted(near_positions):
        lab_aa = lab_seq[lab_pos - 1]
        states = {}
        mapped_positions = {}
        for species in apis_comparators + ["Bombus_terrestris_functional"]:
            qpos = map_position(lab_alignments[species], lab_pos)
            mapped_positions[species] = qpos
            qseq = str(selected[species].seq).replace("*", "")
            states[species] = qseq[qpos - 1] if qpos is not None else "-"
        apis_states = [states[s] for s in apis_comparators]
        if "-" not in apis_states and len(set(apis_states)) == 1 and lab_aa != apis_states[0]:
            nearest = min(TARGETS, key=lambda t: abs((lab_target_positions[t.position] or 10**9) - lab_pos))
            unique_rows.append(
                {
                    "laboriosa_position": lab_pos,
                    "laboriosa_residue": lab_aa,
                    "other_apis_consensus": apis_states[0],
                    "bombus_residue": states["Bombus_terrestris_functional"],
                    "nearest_rat_gtx_site": f"{nearest.expected}{nearest.position}",
                    "nearest_region": nearest.region,
                    "distance_from_mapped_site_aa": lab_pos - (lab_target_positions[nearest.position] or 0),
                    **{f"{s}_position": mapped_positions[s] for s in mapped_positions},
                }
            )

    fields = [
        "laboriosa_position",
        "laboriosa_residue",
        "other_apis_consensus",
        "bombus_residue",
        "nearest_rat_gtx_site",
        "nearest_region",
        "distance_from_mapped_site_aa",
        *[f"{s}_position" for s in apis_comparators + ["Bombus_terrestris_functional"]],
    ]
    with (OUT / "laboriosa_unique_near_gtx_sites.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(unique_rows)

    uniprot = json.loads((ROOT / "references/uniprot_P15390.json").read_text())
    transmembranes = []
    for feature in uniprot["features"]:
        if feature["type"] != "Transmembrane":
            continue
        transmembranes.append(
            (
                int(feature["location"]["start"]["value"]),
                int(feature["location"]["end"]["value"]),
                feature.get("description", "transmembrane"),
            )
        )

    def topology_at(rat_pos):
        if rat_pos is None:
            return "unmapped"
        for start, end, name in transmembranes:
            if start <= rat_pos <= end:
                return name
        return "non-transmembrane"

    genomewide_rows = []
    for lab_pos, lab_aa in enumerate(lab_seq, start=1):
        states = {}
        positions = {}
        for species in apis_comparators + ["Bombus_terrestris_functional"]:
            qpos = map_position(lab_alignments[species], lab_pos)
            positions[species] = qpos
            qseq = str(selected[species].seq).replace("*", "")
            states[species] = qseq[qpos - 1] if qpos is not None else "-"
        apis_states = [states[s] for s in apis_comparators]
        if "-" in apis_states or len(set(apis_states)) != 1 or lab_aa == apis_states[0]:
            continue
        rat_pos = map_query_position(rat_alignments["Apis_laboriosa"], lab_pos)
        nearest = None
        if rat_pos is not None:
            nearest = min(TARGETS, key=lambda t: abs(t.position - rat_pos))
        genomewide_rows.append(
            {
                "laboriosa_position": lab_pos,
                "laboriosa_residue": lab_aa,
                "other_apis_consensus": apis_states[0],
                "bombus_residue": states["Bombus_terrestris_functional"],
                "laboriosa_matches_bombus": lab_aa == states["Bombus_terrestris_functional"],
                "rat_position": rat_pos if rat_pos is not None else "gap",
                "rat_residue": rat[rat_pos - 1] if rat_pos is not None else "-",
                "topology": topology_at(rat_pos),
                "nearest_rat_gtx_site": f"{nearest.expected}{nearest.position}" if nearest else "unmapped",
                "distance_to_nearest_gtx_site": rat_pos - nearest.position if rat_pos is not None and nearest else "unmapped",
                **{f"{s}_residue": states[s] for s in apis_comparators},
            }
        )

    genomewide_fields = [
        "laboriosa_position",
        "laboriosa_residue",
        "other_apis_consensus",
        "bombus_residue",
        "laboriosa_matches_bombus",
        "rat_position",
        "rat_residue",
        "topology",
        "nearest_rat_gtx_site",
        "distance_to_nearest_gtx_site",
        *[f"{s}_residue" for s in apis_comparators],
    ]
    with (OUT / "laboriosa_unique_genomewide.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=genomewide_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(genomewide_rows)

    print(f"Rat reference: {rat_record.id}, {len(rat)} aa")
    for row in summary_rows:
        print(
            f"{row['taxon_or_sequence']}: {row['accession']}, {row['length_aa']} aa, "
            f"identity-to-rat={row['identity_to_rat_over_aligned']}, X={row['x_residues']}"
        )
    print(f"Mapped {len(TARGETS)} experimentally supported GTX positions")
    print(f"Laboriosa-specific substitutions within GTX-site windows: {len(unique_rows)}")
    print(
        "Laboriosa-specific substitutions genome-wide: "
        f"{len(genomewide_rows)}; in mapped transmembrane helices: "
        f"{sum(row['topology'] not in ('non-transmembrane', 'unmapped') for row in genomewide_rows)}"
    )


if __name__ == "__main__":
    main()
