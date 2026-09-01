#!/usr/bin/env python3
"""Locate every protein-level difference among annotated A. laboriosa Para isoforms."""

from __future__ import annotations

import csv
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".deps"))

from Bio import SeqIO  # type: ignore  # noqa: E402
from Bio.Align import PairwiseAligner  # type: ignore  # noqa: E402


GTX_RAT_SITES = {237, 243, 246, 248, 249, 250, 251, 433, 434, 437, 784, 1276, 1463, 1575, 1579, 1586}
PARA_RE = re.compile(r"sodium channel protein para", re.I)


def read_fasta(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return list(SeqIO.parse(handle, "fasta"))


def aligner():
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -8.0
    a.extend_gap_score = -0.5
    return a


def map_reference_positions(alignment):
    mapping = {}
    ref_blocks, query_blocks = alignment.aligned
    for (r0, r1), (q0, q1) in zip(ref_blocks, query_blocks):
        for offset in range(r1 - r0):
            mapping[int(r0 + offset)] = int(q0 + offset)
    return mapping


def map_one(alignment, ref_pos_1based: int):
    return map_reference_positions(alignment).get(ref_pos_1based - 1)


def main():
    proteome = ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz"
    records = [r for r in read_fasta(proteome) if PARA_RE.search(r.description)]
    records.sort(key=lambda r: r.id)
    reference = next(r for r in records if r.id == "XP_043795192.1")
    ref = str(reference.seq).replace("*", "")

    rat_record = read_fasta(ROOT / "references/rat_Nav1.4_NP_037310.2.faa")[0]
    rat = str(rat_record.seq).replace("*", "")
    a = aligner()
    rat_to_ref = a.align(rat, ref)[0]
    ref_to_rat = {}
    rblocks, qblocks = rat_to_ref.aligned
    for (rat0, rat1), (ref0, ref1) in zip(rblocks, qblocks):
        for offset in range(rat1 - rat0):
            ref_to_rat[int(ref0 + offset)] = int(rat0 + offset)

    uniprot = json.loads((ROOT / "references/uniprot_P15390.json").read_text())
    transmembranes = []
    for feature in uniprot["features"]:
        if feature["type"] == "Transmembrane":
            transmembranes.append(
                (
                    int(feature["location"]["start"]["value"]),
                    int(feature["location"]["end"]["value"]),
                    feature.get("description", "transmembrane"),
                )
            )

    def topology(rat_pos):
        if rat_pos is None:
            return "unmapped"
        for start, end, label in transmembranes:
            if start <= rat_pos <= end:
                return label
        return "non-transmembrane"

    maps = {}
    seqs = {}
    for record in records:
        seq = str(record.seq).replace("*", "")
        seqs[record.id] = seq
        maps[record.id] = map_reference_positions(a.align(ref, seq)[0])

    rows = []
    variable_positions = []
    for ref_i, ref_aa in enumerate(ref):
        states = []
        state_by_isoform = {}
        for record in records:
            q_i = maps[record.id].get(ref_i)
            aa = seqs[record.id][q_i] if q_i is not None else "-"
            states.append(aa)
            state_by_isoform[record.id] = aa
        if len(set(states)) == 1:
            continue
        variable_positions.append(ref_i + 1)
        rat_i = ref_to_rat.get(ref_i)
        rat_pos = rat_i + 1 if rat_i is not None else None
        counts = Counter(states)
        rows.append(
            {
                "reference_position": ref_i + 1,
                "reference_residue": ref_aa,
                "rat_position": rat_pos if rat_pos is not None else "gap",
                "rat_residue": rat[rat_i] if rat_i is not None else "-",
                "topology": topology(rat_pos),
                "known_gtx_site": rat_pos in GTX_RAT_SITES if rat_pos is not None else False,
                "state_counts": ";".join(f"{k}:{v}" for k, v in sorted(counts.items())),
                **state_by_isoform,
            }
        )

    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    with (out / "para_isoform_variable_positions.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    blocks = []
    if variable_positions:
        start = previous = variable_positions[0]
        for pos in variable_positions[1:] + [None]:
            if pos is not None and pos == previous + 1:
                previous = pos
                continue
            block_rows = [r for r in rows if start <= r["reference_position"] <= previous]
            mapped_rat = [r["rat_position"] for r in block_rows if isinstance(r["rat_position"], int)]
            topologies = sorted(set(r["topology"] for r in block_rows))
            blocks.append(
                {
                    "reference_start": start,
                    "reference_end": previous,
                    "length_aa": previous - start + 1,
                    "rat_start": min(mapped_rat) if mapped_rat else "gap",
                    "rat_end": max(mapped_rat) if mapped_rat else "gap",
                    "topologies": ";".join(topologies),
                    "contains_transmembrane_residue": any(t not in {"non-transmembrane", "unmapped"} for t in topologies),
                    "contains_known_gtx_site": any(r["known_gtx_site"] for r in block_rows),
                    "reference_sequence": ref[start - 1 : previous],
                }
            )
            if pos is None:
                break
            start = previous = pos

    with (out / "para_isoform_variable_blocks.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=blocks[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(blocks)

    summary = {
        "isoforms": len(records),
        "reference": reference.id,
        "reference_length_aa": len(ref),
        "isoform_lengths_aa": {r.id: len(seqs[r.id]) for r in records},
        "variable_reference_positions": len(rows),
        "variable_blocks": len(blocks),
        "variable_positions_in_transmembrane_helices": sum(
            r["topology"] not in {"non-transmembrane", "unmapped"} for r in rows
        ),
        "variable_positions_at_known_gtx_sites": sum(bool(r["known_gtx_site"]) for r in rows),
        "blocks_touching_transmembrane_helices": sum(bool(b["contains_transmembrane_residue"]) for b in blocks),
    }
    (out / "para_isoform_variation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
