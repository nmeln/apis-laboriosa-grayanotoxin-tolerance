#!/usr/bin/env python3
"""Map the RefSeq Apis laboriosa Para locus to the eastern-Yunnan assembly.

This uses exact genomic probes sampled across the Para CDS. The
assemblies are 99.48% ANI, and dozens of probes are sampled from independent Para CDS
exons.  Concordant hits are enough to identify the eastern-Yunnan scaffold and test
whether the population-study candidate scaffolds (8 and 25) can contain Para.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def fasta_records(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        name = None
        chunks: list[str] = []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks).upper()


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def para_cds(gff_path: Path, protein_id: str):
    rows = []
    with gzip.open(gff_path, "rt") as handle:
        for line in handle:
            if line.startswith("#") or f"protein_id={protein_id}" not in line:
                continue
            parts = line.rstrip().split("\t")
            if len(parts) == 9 and parts[2] == "CDS":
                rows.append((parts[0], int(parts[3]), int(parts[4]), parts[6]))
    if not rows:
        raise SystemExit(f"No CDS rows found for {protein_id}")
    scaffolds = {r[0] for r in rows}
    if len(scaffolds) != 1:
        raise SystemExit(f"Unexpected CDS scaffolds: {sorted(scaffolds)}")
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseq-genome", type=Path, required=True)
    parser.add_argument("--refseq-gff", type=Path, required=True)
    parser.add_argument("--eastern-genome", type=Path, required=True)
    parser.add_argument("--protein", default="XP_043795192.1")
    parser.add_argument("--probe-length", type=int, default=51)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cds = para_cds(args.refseq_gff, args.protein)
    ref_scaffold = cds[0][0]
    ref_sequence = None
    for header, sequence in fasta_records(args.refseq_genome):
        if header.split()[0] == ref_scaffold:
            ref_sequence = sequence
            break
    if ref_sequence is None:
        raise SystemExit(f"RefSeq scaffold {ref_scaffold} not found")

    probes: list[str] = []
    for _, start, end, _ in cds:
        segment = ref_sequence[start - 1 : end]
        if len(segment) < args.probe_length or "N" in segment:
            continue
        offset = (len(segment) - args.probe_length) // 2
        probe = segment[offset : offset + args.probe_length]
        if probe not in probes:
            probes.append(probe)

    hits = []
    for header, sequence in fasta_records(args.eastern_genome):
        matched = []
        for i, probe in enumerate(probes, start=1):
            if probe in sequence or revcomp(probe) in sequence:
                matched.append(i)
        if matched:
            hits.append((header, len(matched), matched))

    hits.sort(key=lambda row: row[1], reverse=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as out:
        out.write("eastern_accession\theader_metadata\tprobe_hits\ttotal_probes\tprobe_indexes\n")
        for header, count, matched in hits:
            header_fields = header.split()
            accession = header_fields[0]
            metadata = " ".join(header_fields[1:])
            out.write(
                f"{accession}\t{metadata}\t{count}\t{len(probes)}\t{','.join(map(str, matched))}\n"
            )

    print(f"RefSeq Para scaffold: {ref_scaffold}")
    print(f"Independent CDS probes: {len(probes)} x {args.probe_length} bp")
    for header, count, _ in hits[:10]:
        print(f"{count}/{len(probes)} probes -> {header}")


if __name__ == "__main__":
    main()
