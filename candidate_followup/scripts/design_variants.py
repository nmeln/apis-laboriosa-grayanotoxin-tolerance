"""Define an untested substitution panel, using natural dorsata codons.

These sequences have not been synthesized, expressed or assayed. They are
provided to make the coding hypothesis precise. CDS records exclude the
terminal stop codon, as do the recovered wild-type CDS records.
"""
from Bio import SeqIO
from Bio.Seq import Seq
from common import OUT, table, write


def main():
    changes = [r for r in table(OUT / 'candidate_coding_differences.tsv')
               if r['kind'] == 'substitution' and r['consistent_across_two_assemblies_each'] == 'True'
               and r['AL_state_absent_from_other_four_bees'] == 'True']
    assert [int(r['AL_position_1based']) for r in changes] == [254,549,1134]
    cds = {r.id: str(r.seq) for r in SeqIO.parse(OUT / 'genome_recovered_cds.fna', 'fasta')}
    proteins = {r.id: str(r.seq) for r in SeqIO.parse(OUT / 'genome_recovered_cds.faa', 'fasta')}
    msa = {r.id: str(r.seq) for r in SeqIO.parse(OUT / 'candidate_eight_bee_alignment.faa', 'fasta')}
    variants = {'AL_reference': cds['AL_Shangrila'], 'AD_reference': cds['AD_Malaysia']}
    instructions = []; triple = list(cds['AL_Shangrila'])
    for r in changes:
        pos = int(r['AL_position_1based']); column = int(r['alignment_column_1based']) - 1
        ad_pos = sum(x != '-' for x in msa['AD_Malaysia'][:column+1])
        al_codon = cds['AL_Shangrila'][3*(pos-1):3*pos]
        ad_codon = cds['AD_Malaysia'][3*(ad_pos-1):3*ad_pos]
        assert str(Seq(al_codon).translate()) == r['AL_Shangrila']
        assert str(Seq(ad_codon).translate()) == r['AD_Malaysia']
        name = 'AL_' + r['AL_Shangrila'] + str(pos) + r['AD_Malaysia']
        variant = list(cds['AL_Shangrila']); variant[3*(pos-1):3*pos] = ad_codon
        variants[name] = ''.join(variant)
        triple[3*(pos-1):3*pos] = ad_codon
        instructions.append(dict(design=name, AL_position_1based=pos, AD_position_1based=ad_pos,
                                 original_codon=al_codon, replacement_codon=ad_codon,
                                 original_aa=r['AL_Shangrila'], replacement_aa=r['AD_Malaysia'],
                                 status='Untested sequence design; no function or transport result'))
    variants['AL_three_substitutions'] = ''.join(triple)
    write('untested_variant_definitions.tsv', instructions)
    with (OUT / 'untested_variant_panel.fna').open('w') as nt, (OUT / 'untested_variant_panel.faa').open('w') as aa:
        for name, sequence in variants.items():
            protein = str(Seq(sequence).translate())
            assert '*' not in protein
            if name.startswith('AL_'):
                differences = sum(a != b for a,b in zip(proteins['AL_Shangrila'],protein))
                expected = 0 if name == 'AL_reference' else 3 if name == 'AL_three_substitutions' else 1
                assert len(protein) == 1336 and differences == expected
            nt.write('>' + name + ' CDS_without_terminal_stop;untested\n' + sequence + '\n')
            aa.write('>' + name + ' untested\n' + protein + '\n')


if __name__ == '__main__': main()
