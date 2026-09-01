#!/usr/bin/env python3
"""Fetch and verify every public input used by the analysis.

The analyzed byte snapshot is defined by data_sources.tsv. Files are downloaded
to temporary paths and published only after their size and SHA-256 hash match.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import os
import re
import shutil
import sys
import tarfile
import tempfile
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data_sources.tsv"
USER_AGENT = "apis-laboriosa-grayanotoxin-tolerance/1.0"
SNAPSHOT_NAME = "apis-laboriosa-inputs-v1.tar"
SNAPSHOT_CHECKSUM_FILE = ROOT / "input_snapshot.sha256"
SNAPSHOT_URL = (
    "https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/"
    "releases/download/input-snapshot-v1/apis-laboriosa-inputs-v1.tar"
)

ASSEMBLY_FILES = {
    "GCF_014066325.1": {
        "genomic.fna.gz": "genomes/apis_laboriosa/GCF_014066325.1_genomic.fna.gz",
        "genomic.gff.gz": "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz",
        "protein.faa.gz": "genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz",
        "rna.fna.gz": "genomes/apis_laboriosa/refseq/GCF_014066325.1_rna.fna.gz",
    },
    "GCF_000469605.1": {
        "genomic.gff.gz": "genomes/apis_dorsata/GCF_000469605.1_genomic.gff.gz",
        "protein.faa.gz": "genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz",
        "rna.fna.gz": "genomes/apis_dorsata/GCF_000469605.1_rna.fna.gz",
    },
    "GCF_003254395.2": {
        "genomic.gff.gz": "genomes/apis_mellifera/GCF_003254395.2_genomic.gff.gz",
        "protein.faa.gz": "genomes/apis_mellifera/GCF_003254395.2_protein.faa.gz",
    },
    "GCF_029169275.2": {
        "genomic.gff.gz": "genomes/apis_cerana/GCF_029169275.2_genomic.gff.gz",
        "protein.faa.gz": "genomes/apis_cerana/GCF_029169275.2_protein.faa.gz",
    },
    "GCF_048593485.1": {
        "genomic.gff.gz": "genomes/apis_florea/GCF_048593485.1_genomic.gff.gz",
        "protein.faa.gz": "genomes/apis_florea/GCF_048593485.1_protein.faa.gz",
    },
    "GCF_000214255.1": {
        "genomic.gff.gz": "genomes/bombus_terrestris/GCF_000214255.1_genomic.gff.gz",
        "protein.faa.gz": "genomes/bombus_terrestris/GCF_000214255.1_protein.faa.gz",
    },
}

ASSEMBLY_BASES = {
    "GCF_014066325.1": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/066/325/"
        "GCF_014066325.1_ASM1406632v1"
    ),
    "GCF_000469605.1": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/469/605/"
        "GCF_000469605.1_Apis_dorsata_1.3"
    ),
    "GCF_003254395.2": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/"
        "GCF_003254395.2_Amel_HAv3.1"
    ),
    "GCF_029169275.2": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/169/275/"
        "GCF_029169275.2_AcerK_1.0"
    ),
    "GCF_048593485.1": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/593/485/"
        "GCF_048593485.1_ASM4859348v1"
    ),
    "GCF_000214255.1": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/"
        "GCF_000214255.1_Bter_1.0"
    ),
}

DIRECT_FILES = {
    "references/GCF_014066325.1_gene_ontology.gaf.gz": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/066/325/"
        "GCF_014066325.1_ASM1406632v1/"
        "GCF_014066325.1_ASM1406632v1_gene_ontology.gaf.gz"
    ),
    "genomes/apis_laboriosa/eastern_yunnan/GWHAOTM00000000.genome.fasta.gz": (
        "https://download.cncb.ac.cn/gwh/Animals/"
        "Apis_laboriosa_A.laboriosa_scaffold_GWHAOTM00000000/"
        "GWHAOTM00000000.genome.fasta.gz"
    ),
    "references/transcriptome_2019/GSM3757258_AL.unigene.fasta.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3757nnn/GSM3757258/suppl/"
        "GSM3757258_AL.unigene.fasta.gz"
    ),
    "references/transcriptome_2019/GSM3757258_AL.Readcount_FPKM.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3757nnn/GSM3757258/suppl/"
        "GSM3757258_AL.Readcount_FPKM.txt.gz"
    ),
    "references/transcriptome_2019/GSM3757259_AD.unigene.fasta.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3757nnn/GSM3757259/suppl/"
        "GSM3757259_AD.unigene.fasta.gz"
    ),
    "references/transcriptome_2019/GSM3757259_AD.Readcount_FPKM.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3757nnn/GSM3757259/suppl/"
        "GSM3757259_AD.Readcount_FPKM.txt.gz"
    ),
    "references/uniprot_P15390.json": (
        "https://rest.uniprot.org/uniprotkb/P15390.json"
    ),
}

SUPPLEMENT_PAGE_URL = (
    "https://oup.silverchair-cdn.com/article-minimal/7044694"
)
SUPPLEMENT_TARGET = (
    "references/population_2023/extracted/Supplemental_Tables.xlsx"
)

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
RAT_ACCESSION = "NP_037310.2"
FUNCTIONAL_PANEL = [
    "AMB38675.1",
    "XP_012167116.1",
    "XP_006613070.1",
    "XP_012347667.1",
]


class SnapshotMismatch(RuntimeError):
    """A completed download did not match the pinned bytes."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest() -> dict[str, dict[str, str]]:
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["relative_path"]: row for row in rows}


def load_snapshot_hash() -> str:
    line = SNAPSHOT_CHECKSUM_FILE.read_text().strip()
    digest, filename = line.split("  ", 1)
    if filename != SNAPSHOT_NAME or len(digest) != 64:
        raise RuntimeError(f"Invalid {SNAPSHOT_CHECKSUM_FILE.name}")
    return digest


def exact_match(path: Path, row: dict[str, str]) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == int(row["size_bytes"])
        and sha256(path) == row["sha256"]
    )


def all_inputs_match(manifest: dict[str, dict[str, str]]) -> bool:
    return all(
        exact_match(ROOT / relative_path, row)
        for relative_path, row in manifest.items()
    )


def urlopen(url: str):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=120)


def download_to_temp(url: str, directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(1, 4):
        temp = tempfile.NamedTemporaryFile(
            prefix="download-", suffix=".part", dir=directory, delete=False
        )
        temp_path = Path(temp.name)
        try:
            with temp, urlopen(url) as response:
                shutil.copyfileobj(response, temp, length=1024 * 1024)
            return temp_path
        except Exception as error:  # network errors vary by platform
            last_error = error
            temp_path.unlink(missing_ok=True)
            if attempt < 3:
                time.sleep(attempt * 2)
    raise RuntimeError(f"Download failed after three attempts: {url}") from last_error


def install_verified(temp_path: Path, relative_path: str, manifest) -> None:
    target = ROOT / relative_path
    row = manifest[relative_path]
    observed_size = temp_path.stat().st_size
    observed_hash = sha256(temp_path)
    expected_size = int(row["size_bytes"])
    expected_hash = row["sha256"]
    if observed_size != expected_size or observed_hash != expected_hash:
        temp_path.unlink(missing_ok=True)
        raise SnapshotMismatch(
            f"Snapshot mismatch for {relative_path}: "
            f"size {observed_size} vs {expected_size}; "
            f"sha256 {observed_hash} vs {expected_hash}"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    os.replace(temp_path, target)
    print(f"verified  {relative_path}")


def fetch_verified(relative_path: str, url: str, manifest, force: bool) -> None:
    target = ROOT / relative_path
    if not force and exact_match(target, manifest[relative_path]):
        print(f"cached    {relative_path}")
        return
    last_error: SnapshotMismatch | None = None
    for attempt in range(1, 4):
        temp_path = download_to_temp(url, target.parent)
        try:
            install_verified(temp_path, relative_path, manifest)
            return
        except SnapshotMismatch as error:
            last_error = error
            if attempt < 3:
                print(
                    f"retrying  {relative_path} after checksum mismatch "
                    f"({attempt}/3)",
                    file=sys.stderr,
                )
    raise last_error


def fetch_assemblies(manifest, force: bool) -> None:
    for accession, files in ASSEMBLY_FILES.items():
        base = ASSEMBLY_BASES[accession]
        basename = base.rsplit("/", 1)[-1]
        for suffix, relative_path in files.items():
            url = f"{base}/{basename}_{suffix}"
            fetch_verified(relative_path, url, manifest, force)


def efetch_url(accessions: list[str]) -> str:
    query = urllib.parse.urlencode(
        {
            "db": "protein",
            "id": ",".join(accessions),
            "rettype": "fasta",
            "retmode": "text",
        }
    )
    return f"{EFETCH_URL}?{query}"


def fetch_supplement(manifest, force: bool) -> None:
    target = ROOT / SUPPLEMENT_TARGET
    if not force and exact_match(target, manifest[SUPPLEMENT_TARGET]):
        print(f"cached    {SUPPLEMENT_TARGET}")
        return
    article_page = download_to_temp(SUPPLEMENT_PAGE_URL, ROOT / ".cache")
    try:
        article_html = article_page.read_text(encoding="utf-8")
    finally:
        article_page.unlink(missing_ok=True)
    match = re.search(
        r'href="([^"]*evad025_supplementary_data\.zip[^"]*)"',
        article_html,
        flags=re.IGNORECASE,
    )
    if not match:
        raise RuntimeError("Cao et al. supplementary ZIP link missing from OUP page")
    supplement_url = html.unescape(match.group(1))
    archive = download_to_temp(supplement_url, ROOT / ".cache")
    try:
        with zipfile.ZipFile(archive) as bundle:
            members = [
                name
                for name in bundle.namelist()
                if Path(name).name.lower() == "supplemental_tables.xlsx"
            ]
            if len(members) != 1:
                raise RuntimeError(
                    "Expected one Supplemental_Tables.xlsx in the Cao et al. supplement; "
                    f"found {members}"
                )
            target.parent.mkdir(parents=True, exist_ok=True)
            with bundle.open(members[0]) as source:
                temp = tempfile.NamedTemporaryFile(
                    prefix="supplement-", suffix=".xlsx.part",
                    dir=target.parent, delete=False
                )
                with temp:
                    shutil.copyfileobj(source, temp)
                temp_path = Path(temp.name)
        install_verified(temp_path, SUPPLEMENT_TARGET, manifest)
    finally:
        archive.unlink(missing_ok=True)


def fetch_snapshot(url: str, manifest: dict[str, dict[str, str]]) -> None:
    archive_path = download_to_temp(url, ROOT / ".cache")
    try:
        observed_hash = sha256(archive_path)
        expected_hash = load_snapshot_hash()
        if observed_hash != expected_hash:
            raise SnapshotMismatch(
                f"Input archive SHA-256 mismatch: {observed_hash} vs {expected_hash}"
            )

        with tarfile.open(archive_path, mode="r:") as archive:
            members = archive.getmembers()
            expected_names = set(manifest)
            names = [member.name for member in members]
            if len(names) != len(set(names)):
                raise RuntimeError("Input archive contains duplicate paths")
            if set(names) != expected_names:
                missing = sorted(expected_names - set(names))
                extra = sorted(set(names) - expected_names)
                raise RuntimeError(
                    f"Input archive path mismatch; missing={missing}; extra={extra}"
                )
            if any(not member.isfile() for member in members):
                raise RuntimeError("Input archive contains a non-file member")

            members_by_name = {member.name: member for member in members}
            for relative_path in manifest:
                member = members_by_name[relative_path]
                source = archive.extractfile(member)
                if source is None:
                    raise RuntimeError(f"Cannot read archive member: {relative_path}")
                target = ROOT / relative_path
                target.parent.mkdir(parents=True, exist_ok=True)
                temporary = tempfile.NamedTemporaryFile(
                    prefix="snapshot-",
                    suffix=".part",
                    dir=target.parent,
                    delete=False,
                )
                temporary_path = Path(temporary.name)
                try:
                    with source, temporary:
                        shutil.copyfileobj(source, temporary, length=1024 * 1024)
                    install_verified(temporary_path, relative_path, manifest)
                except Exception:
                    temporary_path.unlink(missing_ok=True)
                    raise
    finally:
        archive_path.unlink(missing_ok=True)


def verify_complete(manifest) -> None:
    errors = []
    for relative_path, row in manifest.items():
        path = ROOT / relative_path
        if not exact_match(path, row):
            errors.append(relative_path)
    if errors:
        raise RuntimeError(f"Input verification failed: {errors}")
    print(f"All {len(manifest)} input files match the analyzed snapshot.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--force", action="store_true",
        help="redownload files that already pass verification",
    )
    source_group = parser.add_mutually_exclusive_group()
    source_group.add_argument(
        "--official-only",
        action="store_true",
        help="skip the GitHub input snapshot and use official repositories",
    )
    source_group.add_argument(
        "--snapshot-only",
        action="store_true",
        help="require the pinned GitHub input snapshot",
    )
    parser.add_argument(
        "--snapshot-url",
        default=SNAPSHOT_URL,
        help="override the pinned input snapshot URL",
    )
    args = parser.parse_args()
    manifest = load_manifest()

    if not args.force and all_inputs_match(manifest):
        print(f"All {len(manifest)} input files already match the analyzed snapshot.")
        return

    if not args.official_only:
        if not args.snapshot_url:
            if args.snapshot_only:
                raise RuntimeError("No input snapshot URL is configured")
        else:
            try:
                fetch_snapshot(args.snapshot_url, manifest)
                verify_complete(manifest)
                return
            except Exception as error:
                if args.snapshot_only:
                    raise
                print(
                    f"WARNING: input snapshot unavailable ({error}); "
                    "trying official repositories",
                    file=sys.stderr,
                )

    fetch_assemblies(manifest, args.force)
    for relative_path, url in DIRECT_FILES.items():
        fetch_verified(relative_path, url, manifest, args.force)
    fetch_supplement(manifest, args.force)
    fetch_verified(
        "references/rat_Nav1.4_NP_037310.2.faa",
        efetch_url([RAT_ACCESSION]),
        manifest,
        args.force,
    )
    fetch_verified(
        "references/bee_nav_functional_panel.faa",
        efetch_url(FUNCTIONAL_PANEL),
        manifest,
        args.force,
    )
    verify_complete(manifest)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
