#!/usr/bin/env python3
"""Run every analysis in dependency order."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def command(script: str, *arguments: str) -> list[str]:
    return [sys.executable, str(ROOT / "scripts" / script), *arguments]


PIPELINE = [
    command("analyze_nav.py"),
    command("analyze_para_isoforms.py"),
    command("compare_para_isoforms.py"),
    command("analyze_para_transcriptomes.py"),
    command("extract_population_table_s6.py"),
    command(
        "map_para_to_eastern.py",
        "--refseq-genome",
        "genomes/apis_laboriosa/GCF_014066325.1_genomic.fna.gz",
        "--refseq-gff",
        "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz",
        "--eastern-genome",
        "genomes/apis_laboriosa/eastern_yunnan/GWHAOTM00000000.genome.fasta.gz",
        "--output",
        "results/para_eastern_scaffold_probe_hits.tsv",
    ),
    command("map_refseq_genes_to_selected_scaffolds.py"),
    command("analyze_selected_ion_channels.py"),
    command("analyze_detox_families.py"),
    command("analyze_barrier_clearance.py"),
    command("qc_chitin_synthase_loci.py"),
    command("analyze_constitutive_transcript_expression.py"),
]


def run(step: list[str]) -> None:
    display = " ".join(Path(part).name if part.startswith(str(ROOT)) else part for part in step)
    print(f"\n$ {display}", flush=True)
    subprocess.run(step, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-input-verification", action="store_true",
        help="run without checking data_sources.tsv first",
    )
    parser.add_argument(
        "--skip-output-verification", action="store_true",
        help="do not compare generated files with results.sha256",
    )
    args = parser.parse_args()

    if not args.skip_input_verification:
        run(command("verify_project.py", "--inputs"))
    for step in PIPELINE:
        run(step)
    run(command("validate_claims.py"))
    if not args.skip_output_verification:
        run(command("verify_project.py", "--results"))


if __name__ == "__main__":
    main()
