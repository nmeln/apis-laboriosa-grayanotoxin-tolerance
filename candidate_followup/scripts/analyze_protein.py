"""Compare the full candidate, independent genomes and ABC family references.

The eight-bee alignment uses the comparative addendum's deterministic FAMSA
settings. Human functional features are projected through a separate global
protein alignment and remain homology-based annotations. Similarity scores
are not a gene tree and do not establish one-to-one fly/human orthology.
"""
import gzip
import json
import re
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from pyfamsa import Aligner, Sequence
from common import ROOT, INPUT, OUT, WORK, table, write

SOURCES = {
    'AL_Shangrila': ('genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz', 'XP_043798878.1'),
    'AD_Malaysia': ('genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz', 'XP_006621542.1'),
    'AM': ('genomes/apis_mellifera/GCF_003254395.2_protein.faa.gz', 'XP_393750.5'),
    'AC': ('genomes/apis_cerana/GCF_029169275.2_protein.faa.gz', 'XP_061938405.1'),
    'AF': ('genomes/apis_florea/GCF_048593485.1_protein.faa.gz', 'XP_086761534.1'),
    'BT': ('candidate_followup/inputs/GCF_910591885.1_protein.faa.gz', 'XP_048266277.1'),
}


def pair_aligner(mode='global'):
    return Align.PairwiseAligner(mode=mode, substitution_matrix=substitution_matrices.load('BLOSUM62'),
                                open_gap_score=-10, extend_gap_score=-.5)


def projection(target, query):
    a = pair_aligner().align(target, query)[0]
    mapping = {}; t = q = 0
    for x, y in zip(str(a[0]), str(a[1])):
        if x != '-': t += 1
        if y != '-': q += 1
        if x != '-' and y != '-': mapping[t] = q
    return a, mapping


def main():
    group = next(r for r in table(ROOT / 'comparative_addendum/results/orthofinder_key/Orthogroups.tsv') if r['Orthogroup'] == 'OG0000499')
    accessions = set(group.values()) - {'OG0000499'}
    seqs = {}; labels = []
    for label, (path, accession) in SOURCES.items():
        assert accession in accessions
        with gzip.open(ROOT / path, 'rt') as f:
            r = next(r for r in SeqIO.parse(f, 'fasta') if r.id == accession)
        seqs[label] = str(r.seq)
        labels.append(dict(label=label, accession=accession, length_aa=len(r), source='RefSeq primary protein in OG0000499'))
    for r in SeqIO.parse(OUT / 'genome_recovered_cds.faa', 'fasta'):
        if r.id in seqs: assert seqs[r.id] == str(r.seq)
        else:
            seqs[r.id] = str(r.seq)
            labels.append(dict(label=r.id, accession=r.id, length_aa=len(r), source='Translation of genomic CDS reconstructed with miniprot'))
    ordered = list(seqs)
    msa = Aligner(threads=1, refine=False).align([Sequence(k.encode(), seqs[k].encode()) for k in ordered])
    aligned = {s.id.decode(): s.sequence.decode() for s in msa}
    assert set(aligned) == set(seqs)
    with (OUT / 'candidate_eight_bee_alignment.faa').open('w') as f:
        for k in ordered: f.write('>' + k + '\n' + aligned[k] + '\n')
    write('candidate_protein_members.tsv', labels)
    al = seqs['AL_Shangrila']
    human = json.loads((INPUT / 'Human_ABCC_reviewed.json').read_text())['results']
    fly = json.loads((INPUT / 'Dmel_ABCC_domain.json').read_text())['results']
    rank_rows = []
    for species, panel in [('Homo_sapiens', human), ('Drosophila_melanogaster', fly)]:
        for r in panel:
            protein = r['sequence']['value']
            if len(protein) < 1000: continue
            a, mapping = projection(al, protein)
            left, right = str(a[0]), str(a[1]); paired = sum(x != '-' and y != '-' for x,y in zip(left,right))
            matches = sum(x == y and x != '-' for x,y in zip(left,right))
            rank_rows.append(dict(species=species, accession=r['primaryAccession'], gene=r.get('genes',[{}])[0].get('geneName',{}).get('value',''),
                                  reference_aa=len(protein), global_blosum62_score=a.score, paired_columns=paired,
                                  identity_in_paired_columns=matches/paired, identical_aa=matches,
                                  candidate_paired_fraction=paired/len(al)))
    write('abc_reference_similarity.tsv', sorted(rank_rows, key=lambda r: (r['species'], -r['global_blosum62_score'], r['accession'])))
    mrp4 = next(r for r in human if r['primaryAccession'] == 'O15439')
    hseq = mrp4['sequence']['value']; pair, mapping = projection(al, hseq)
    with (OUT / 'candidate_human_ABCC4_alignment.faa').open('w') as f:
        f.write('>AL_XP_043798878.1\n' + str(pair[0]) + '\n>Human_O15439\n' + str(pair[1]) + '\n')
    feature_rows = []
    for feature in mrp4['features']:
        if feature['type'] not in ['Domain','Transmembrane','Binding site','Motif']: continue
        start, end = feature['location']['start']['value'], feature['location']['end']['value']
        positions = [p for p, h in mapping.items() if start <= h <= end]
        feature_rows.append(dict(reference='O15439', type=feature['type'], description=feature.get('description',''),
                                 ligand=feature.get('ligand',{}).get('name',''), human_start=start, human_end=end,
                                 AL_mapped_positions=';'.join(map(str, positions)),
                                 AL_mapped_min=min(positions) if positions else None, AL_mapped_max=max(positions) if positions else None,
                                 annotation_basis='Human UniProt feature projected by protein alignment; bee function unmeasured'))
    write('human_ABCC4_feature_projection.tsv', feature_rows)
    changes = []; pairwise = []
    positions = {k: 0 for k in ordered}
    for c in range(len(aligned['AL_Shangrila'])):
        residues = {k: aligned[k][c] for k in ordered}
        for k in ordered:
            if residues[k] != '-': positions[k] += 1
        if residues['AL_Shangrila'] == residues['AD_Malaysia']: continue
        pos = positions['AL_Shangrila']; human_pos = mapping.get(pos)
        features = []
        for r in feature_rows:
            if human_pos is not None and r['human_start'] <= human_pos <= r['human_end']:
                features.append(r['type'] + ':' + (r['description'] or r['ligand']))
        changes.append(dict(alignment_column_1based=c+1, AL_position_1based=pos,
                            kind='substitution' if residues['AL_Shangrila'] != '-' and residues['AD_Malaysia'] != '-' else 'gap_column',
                            **residues,
                            consistent_across_two_assemblies_each=residues['AL_Shangrila'] == residues['AL_Pingbian'] and residues['AD_Malaysia'] == residues['AD_Thailand'],
                            AL_state_absent_from_other_four_bees=all(residues['AL_Shangrila'] != residues[k] for k in ['AM','AC','AF','BT']),
                            human_ABCC4_position=human_pos, projected_human_features=';'.join(features)))
    write('candidate_coding_differences.tsv', changes)
    for k in ordered:
        x, y = aligned['AL_Shangrila'], aligned[k]
        paired = sum(a != '-' and b != '-' for a,b in zip(x,y)); matches = sum(a == b and a != '-' for a,b in zip(x,y))
        pairwise.append(dict(reference='AL_Shangrila', comparator=k, paired_columns=paired, identical_aa=matches,
                             identity_in_paired_columns=matches/paired, gap_columns=sum((a=='-') != (b=='-') for a,b in zip(x,y))))
    write('candidate_bee_pairwise.tsv', pairwise)
    # Search the complete submitted focal transcript in all six reading frames.
    with gzip.open(ROOT / 'references/transcriptome_2019/GSM3757258_AL.unigene.fasta.gz', 'rt') as f:
        rna = next(r for r in SeqIO.parse(f, 'fasta') if r.id == 'AL|c41541_g1')
    exact = []
    for strand, s in [('+',rna.seq),('-',rna.seq.reverse_complement())]:
        for frame in range(3):
            protein = str(s[frame:len(s)-(len(s)-frame)%3].translate())
            offset = protein.find(al)
            if offset >= 0: exact.append((strand, frame+3*offset, str(s[frame+3*offset:frame+3*(offset+len(al))])))
    assert len(exact) == 1, 'Expected one complete exact AL transcript ORF'
    strand, offset, cds = exact[0]
    write('candidate_transcript_support.tsv', [dict(contig=rna.id, strand=strand, oriented_cds_start_0based=offset,
                                                   coding_nt=len(cds), translated_aa=len(al), exact_reference_protein_match=True)])
    with (OUT / 'candidate_AL_transcript_cds.fna').open('w') as f: f.write('>AL_c41541_g1_CDS\n'+cds+'\n')
    print('Stable two-assembly differences:', [(r['AL_position_1based'],r['kind'],r['AL_Shangrila'],r['AD_Malaysia']) for r in changes if r['consistent_across_two_assemblies_each']])


if __name__ == '__main__': main()
