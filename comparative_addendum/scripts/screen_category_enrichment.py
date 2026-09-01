#!/usr/bin/env python3
"""Exploratory category-level calibration of strict amino-acid sharing.

The tests are unstratified Fisher screens over annotation-derived categories.
They are reported only to calibrate individual leads and should not replace the
length/divergence-matched aggregate test in screen_strict_convergence.py.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from scipy.stats import fisher_exact


APIS = [
    "Apis_laboriosa",
    "Apis_dorsata",
    "Apis_mellifera",
    "Apis_cerana",
    "Apis_florea",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bh_adjust(pvalues: list[float]) -> list[float]:
    order = sorted(range(len(pvalues)), key=pvalues.__getitem__)
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    for reverse_rank, index in enumerate(reversed(order), start=1):
        rank = len(pvalues) - reverse_rank + 1
        running = min(running, pvalues[index] * len(pvalues) / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = read_tsv(args.summary)
    categories = sorted({
        category
        for row in rows
        for category in row["candidate_categories"].split(";")
        if category
    })

    results = []
    for focal in APIS:
        hit_key = f"hits_{focal}"
        for category in categories:
            selected = [
                row for row in rows
                if category in row["candidate_categories"].split(";")
            ]
            if len(selected) < 2:
                continue
            background = [
                row for row in rows
                if category not in row["candidate_categories"].split(";")
            ]
            selected_hits = sum(int(row[hit_key]) for row in selected)
            background_hits = sum(int(row[hit_key]) for row in background)
            selected_sites = sum(int(row["callable_sites"]) for row in selected)
            background_sites = sum(int(row["callable_sites"]) for row in background)
            site_or, site_p = fisher_exact(
                [
                    [selected_hits, selected_sites - selected_hits],
                    [background_hits, background_sites - background_hits],
                ],
                alternative="two-sided",
            )
            selected_genes_hit = sum(int(row[hit_key]) > 0 for row in selected)
            background_genes_hit = sum(int(row[hit_key]) > 0 for row in background)
            gene_or, gene_p = fisher_exact(
                [
                    [selected_genes_hit, len(selected) - selected_genes_hit],
                    [background_genes_hit, len(background) - background_genes_hit],
                ],
                alternative="two-sided",
            )
            results.append({
                "focal_apis": focal,
                "category": category,
                "category_orthogroups": len(selected),
                "category_callable_sites": selected_sites,
                "category_strict_sites": selected_hits,
                "category_orthogroups_with_hit": selected_genes_hit,
                "site_fisher_odds_ratio": site_or,
                "site_fisher_p_two_sided": site_p,
                "gene_fisher_odds_ratio": gene_or,
                "gene_fisher_p_two_sided": gene_p,
                "warning": "unstratified exploratory test; annotation categories include false positives",
            })
    for row, qvalue in zip(
        results, bh_adjust([float(row["site_fisher_p_two_sided"]) for row in results])
    ):
        row["site_fisher_bh_q"] = qvalue
    for row, qvalue in zip(
        results, bh_adjust([float(row["gene_fisher_p_two_sided"]) for row in results])
    ):
        row["gene_fisher_bh_q"] = qvalue

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(results[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(results)
    print(f"tests={len(results)}")


if __name__ == "__main__":
    main()
