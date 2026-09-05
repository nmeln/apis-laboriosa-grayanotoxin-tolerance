"""Shared paths and deterministic tabular I/O. No network access."""
from pathlib import Path
import csv
import subprocess

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
INPUT = HERE / 'inputs'
WORK = HERE / 'work'
OUT = HERE / 'results'
MM = ROOT / '.tools/minimap2-2.31_x64-linux/minimap2'
MP = ROOT / '.tools/miniprot-0.18_x64-linux/miniprot'
GENOMES = {
    'AD_Malaysia': INPUT / 'GCF_000469605.1_genomic.fna.gz',
    'AD_Thailand': INPUT / 'GCA_009792835.1_genomic.fna.gz',
    'AL_Shangrila': ROOT / 'genomes/apis_laboriosa/GCF_014066325.1_genomic.fna.gz',
    'AL_Pingbian': ROOT / 'genomes/apis_laboriosa/eastern_yunnan/GWHAOTM00000000.genome.fasta.gz',
}
CANDIDATES = {'AL': 'XP_043798878.1', 'AD': 'XP_006621542.1'}


def table(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(name, rows, fields=None):
    rows = list(rows)
    if fields is None:
        if not rows:
            raise ValueError('Explicit columns required for empty table: ' + name)
        fields = list(rows[0])
    OUT.mkdir(exist_ok=True)
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def align(tool, args, name):
    WORK.mkdir(exist_ok=True)
    with (WORK / name).open('w') as out, (WORK / (name + '.log')).open('w') as err:
        subprocess.run([str(tool), *map(str, args)], stdout=out, stderr=err, check=True)
