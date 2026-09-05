"""Map preserved read pairs against both complete species references."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import pysam
from common import INPUT, WORK, OUT, MM, GENOMES, digest, json_write, table

REFERENCES = ('AL_Shangrila', 'AD_Malaysia')


def merge_reads(runs, destination):
    paths = [destination / f'merged_{m}.fastq' for m in (1, 2)]
    with paths[0].open('w') as o1, paths[1].open('w') as o2:
        for run in runs:
            source = INPUT / 'selected_reads' / run
            evidence = json.loads((source / 'collection.json').read_text())
            assert evidence['complete_source_files_verified'] and evidence['selected_records_copied_from_original_fastqs']
            for item in evidence['selected_files']:
                assert digest(source/item['name']) == item['sha256']
            count = 0
            with (source / f'{run}_1.fastq').open() as f1, (source / f'{run}_2.fastq').open() as f2:
                while h1 := f1.readline():
                    h2 = f2.readline()
                    a, b = [f1.readline() for _ in range(3)], [f2.readline() for _ in range(3)]
                    n1, n2 = h1.split()[0][1:], h2.split()[0][1:]
                    if n1.endswith('/1'): n1 = n1[:-2]
                    if n2.endswith('/2'): n2 = n2[:-2]
                    assert n1 == n2, (run, n1, n2)
                    for out, record in ((o1, a), (o2, b)):
                        out.write('@' + run + ':' + n1 + '\n' + ''.join(record))
                    count += 1
                assert not f2.readline() and count == evidence['selected_pairs']
    return paths


def map_reference(reference, runs, fastqs, destination):
    cache = Path(tempfile.gettempdir()) / 'bee_population_references'
    cache.mkdir(exist_ok=True)
    fasta = cache / (reference + '.fna')
    source_hash = digest(GENOMES[reference])
    stamp = cache / (reference + '.source.sha256')
    if not fasta.exists() or not stamp.exists() or stamp.read_text().strip() != source_hash:
        with gzip.open(GENOMES[reference], 'rb') as src, fasta.open('wb') as out:
            shutil.copyfileobj(src, out)
        stamp.write_text(source_hash+'\n')
    sam = destination / f'{reference}.sam'
    log = destination / f'{reference}.mapping.log'
    command = [str(MM), '-ax', 'sr', '--secondary=yes', '--MD', '-t', '2', str(fasta), *map(str, fastqs)]
    with sam.open('w') as out, log.open('w') as err:
        subprocess.run(command, stdout=out, stderr=err, check=True)
    tagged = destination / f'{reference}.tagged.bam'
    with pysam.AlignmentFile(str(sam), 'r') as src:
        header = src.header.to_dict()
        header['RG'] = [dict(ID=r, SM=r, LB=r, PL='ILLUMINA') for r in runs]
        with pysam.AlignmentFile(str(tagged), 'wb', header=header) as out:
            for read in src:
                run = read.query_name.split(':', 1)[0]
                assert run in runs
                read.set_tag('RG', run)
                out.write(read)
    named, fixed, sorted_bam = [destination / f'{reference}.{suffix}.bam' for suffix in ('name', 'fixmate', 'sorted')]
    bam = destination / f'{reference}.bam'
    pysam.sort('-n', '-@', '2', '-o', str(named), str(tagged))
    pysam.fixmate('-m', str(named), str(fixed))
    pysam.sort('-@', '2', '-o', str(sorted_bam), str(fixed))
    pysam.markdup('--use-read-groups', '-s', '-f', str(destination / f'{reference}.duplicates.txt'), str(sorted_bam), str(bam))
    pysam.index(str(bam))
    summary = json.loads(pysam.flagstat('-O', 'json', str(bam)))
    json_write(destination/f'{reference}.mapping_summary.json', summary)
    for p in (sam, tagged, named, fixed, sorted_bam): p.unlink()
    return bam


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--runs', default='all')
    p.add_argument('--label', default='population')
    args = p.parse_args()
    available = [r['run_accession'] for r in table(OUT/'sample_inventory.tsv')]
    runs = sorted(available if args.runs == 'all' else args.runs.split(','))
    assert set(runs) <= set(available)
    destination = WORK / args.label
    destination.mkdir(parents=True, exist_ok=True)
    fastqs = merge_reads(runs, destination)
    for reference in REFERENCES:
        map_reference(reference, runs, fastqs, destination)
        print('Mapped', len(runs), 'individuals to', reference, flush=True)
    json_write(destination/'mapping_inputs.json', dict(runs=runs, minimap2_version=subprocess.check_output([str(MM), '--version'],text=True).strip(), pysam_version=pysam.__version__, samtools_version=pysam.__samtools_version__, references={k:digest(GENOMES[k]) for k in REFERENCES}))

if __name__ == '__main__': main()
