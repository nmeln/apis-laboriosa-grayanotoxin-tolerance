"""Paths, deterministic tables and fixed biological coordinates."""
from pathlib import Path
import csv
import hashlib
import json

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
INPUT = HERE / 'inputs'
WORK = HERE / 'work'
OUT = HERE / 'results'
MM = ROOT / '.tools/minimap2-2.31_x64-linux/minimap2'
GENOMES = {
    'AL_Shangrila': ROOT / 'genomes/apis_laboriosa/GCF_014066325.1_genomic.fna.gz',
    'AD_Malaysia': ROOT / 'candidate_followup/inputs/GCF_000469605.1_genomic.fna.gz',
    'AL_Pingbian': ROOT / 'genomes/apis_laboriosa/eastern_yunnan/GWHAOTM00000000.genome.fasta.gz',
    'AD_Thailand': ROOT / 'candidate_followup/inputs/GCA_009792835.1_genomic.fna.gz',
}
SITES = {254: 'primary', 549: 'primary', 1134: 'primary',
         452: 'within_species_comparison', 822: 'within_species_comparison',
         474: 'secondary', 717: 'secondary', 727: 'secondary', 766: 'secondary'}


def table(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(name, rows, fields=None):
    rows = list(rows)
    if fields is None:
        if not rows:
            raise ValueError('Need explicit fields for an empty table: ' + name)
        fields = list(rows[0])
    OUT.mkdir(exist_ok=True)
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        while data := f.read(1024 * 1024):
            h.update(data)
    return h.hexdigest()


def json_write(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
