"""Reference swaps use all primary host CDSs, retaining competing genes."""
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from common import ROOT, WORK, MM, table, align


def main():
    WORK.mkdir(exist_ok=True)
    source = ROOT / 'transcriptomic_addendum/work/pooled_primary_cds.fna'
    metadata = table(ROOT / 'transcriptomic_addendum/results/raw_reference_metadata.tsv')
    records = list(SeqIO.parse(source, 'fasta'))
    assert len(records) == len(metadata) == 19632
    for ref in ['AL', 'AD']:
        SeqIO.write([r for r in records if r.id.startswith(ref + '|')], WORK / (ref + '_primary_cds.fna'), 'fasta')
    def run(pair):
        sample, ref = pair
        accession = {'AL': 'SRR9034695', 'AD': 'SRR9034696'}[sample]
        reads = ROOT / 'transcriptomic_addendum/inputs' / (accession + '.first_2000000.R1.fastq.gz')
        name = sample + '_vs_' + ref + '.paf'
        align(MM, ['-x', 'sr', '-c', '--secondary=yes', '-N', '20', '-p', '0.8', '-t', '2', WORK / (ref + '_primary_cds.fna'), reads], name)
        print(name, flush=True)
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(run, [(s, r) for s in ['AL', 'AD'] for r in ['AL', 'AD']]))


if __name__ == '__main__':
    main()
