"""Recover the candidate independently from four complete assembly files.

Miniprot is run with both focal-species proteins against each assembly. CDS
coordinates must agree between queries. Translation is derived from genomic
bases, including split codons across exon boundaries. No predicted query
residue is substituted for a genomic base.
"""
import gzip
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq
from common import ROOT, WORK, OUT, GENOMES, MP, MM, CANDIDATES, align, write, table


def attrs(s):
    return dict(x.split('=', 1) for x in s.split(';') if '=' in x)


def prepare():
    metadata = table(ROOT / 'transcriptomic_addendum/results/raw_reference_metadata.tsv')
    paths = {'AL': 'genomes/apis_laboriosa/refseq/GCF_014066325.1',
             'AD': 'genomes/apis_dorsata/GCF_000469605.1'}
    proteins, rnas = {}, {}
    for sp, prefix in paths.items():
        tid = next(r['transcript'] for r in metadata if r['protein'] == CANDIDATES[sp])
        for suffix, target, output in [('protein.faa.gz', CANDIDATES[sp], proteins), ('rna.fna.gz', tid, rnas)]:
            with gzip.open(ROOT / (prefix + '_' + suffix), 'rt') as f:
                record = next(r for r in SeqIO.parse(f, 'fasta') if r.id == target)
            output[sp + '|' + target] = str(record.seq)
    for name, records in [('candidate_queries.faa', proteins), ('candidate_queries.fna', rnas)]:
        with (WORK / name).open('w') as f:
            for rid, s in records.items():
                f.write('>' + rid + '\n' + s + '\n')
    return proteins


def main():
    WORK.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
    queries = prepare()
    def run(item):
        label, genome = item
        align(MP, ['--gff', '-t', '2', genome, WORK / 'candidate_queries.faa'], label + '.miniprot.gff')
        align(MM, ['-x', 'splice:hq', '-uf', '-c', '--cs=long', '-t', '2', genome, WORK / 'candidate_queries.fna'], label + '.rna.paf')
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(run, GENOMES.items()))
    summaries, exons, recovered = [], [], {}
    rna_rows = []
    for label, genome in GENOMES.items():
        models, features = {}, defaultdict(list)
        with (WORK / (label + '.miniprot.gff')).open() as f:
            for line in f:
                if line.startswith('#'): continue
                p = line.rstrip().split('\t'); a = attrs(p[8])
                if p[2] == 'mRNA': models[a['ID']] = (p, a)
                elif p[2] == 'CDS': features[a['Parent']].append(p)
        assert len(models) == 2, (label, 'Expected one locus per query; inspect ambiguous hits')
        seqids = {p[0] for p, a in models.values()}
        assert len(seqids) == 1, (label, seqids)
        with gzip.open(genome, 'rt') as f:
            contigs = {r.id: r.seq for r in SeqIO.parse(f, 'fasta') if r.id in seqids}
        reconstructed = []
        coordinate_sets = []
        for model, (p, a) in models.items():
            cds = sorted(features[model], key=lambda p: int(p[3]), reverse=p[6] == '-')
            coordinate_sets.append([(c[0], c[3], c[4], c[6], c[7]) for c in cds])
            sequence = ''
            for idx, c in enumerate(cds, 1):
                fragment = contigs[c[0]][int(c[3])-1:int(c[4])]
                if c[6] == '-': fragment = fragment.reverse_complement()
                assert int(c[7]) == (3 - len(sequence) % 3) % 3, (label, idx, 'CDS phase')
                sequence += str(fragment)
                exons.append(dict(assembly=label, query=a['Target'].split()[0], exon=idx,
                                  scaffold=c[0], start_1based=int(c[3]), end_1based=int(c[4]), strand=c[6], phase=int(c[7])))
            assert len(sequence) % 3 == 0 and set(sequence.upper()) <= set('ACGT'), label
            protein = str(Seq(sequence).translate())
            assert protein.startswith('M') and protein.endswith('*') and '*' not in protein[:-1], label
            reconstructed.append((sequence[:-3], protein[:-1]))
            qid, qstart, qend = a['Target'].split()
            assert int(qstart) == 1 and int(qend) == len(queries[qid]), label
            summaries.append(dict(assembly=label, query=qid, scaffold=p[0],
                                  start_1based=int(p[3]), end_1based=int(p[4]), strand=p[6],
                                  coding_exons=len(cds), recovered_protein_aa=len(protein)-1,
                                  miniprot_identity=a['Identity'], query_start_1based=qstart,
                                  query_end_1based=qend, internal_stops=0, ambiguous_coding_bases=0))
        assert coordinate_sets[0] == coordinate_sets[1], (label, 'Query-dependent exon structure')
        assert reconstructed[0] == reconstructed[1]
        recovered[label] = reconstructed[0]
        if label in ['AL_Shangrila', 'AD_Malaysia']:
            query = next(s for q, s in queries.items() if q.startswith(label[:2] + '|'))
            assert recovered[label][1] == query, (label, 'Known annotation reconstruction failed')
        with (WORK / (label + '.rna.paf')).open() as f:
            for line in f:
                p = line.rstrip().split('\t')
                rna_rows.append(dict(assembly=label, query=p[0], query_nt=int(p[1]),
                                     aligned_query_nt=int(p[3])-int(p[2]), scaffold=p[5],
                                     target_start_0based=p[7], target_end_exclusive=p[8],
                                     matched_bases=int(p[9]), alignment_block_nt=int(p[10]),
                                     identity=float(p[9])/int(p[10]), mapq=p[11]))
    write('genome_candidate_models.tsv', summaries)
    write('genome_candidate_exons.tsv', exons)
    write('genome_rna_alignment_summary.tsv', rna_rows)
    for suffix, index in [('fna', 0), ('faa', 1)]:
        with (OUT / ('genome_recovered_cds.' + suffix)).open('w') as f:
            for label, sequences in recovered.items():
                f.write('>' + label + '\n' + sequences[index] + '\n')
    print('Recovered matching exon models from both queries in all four assemblies.')


if __name__ == '__main__': main()
