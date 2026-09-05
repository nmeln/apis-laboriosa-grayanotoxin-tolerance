"""Map codons through coding exons; make allele-masked, low-complexity-filtered baits."""
from collections import Counter
import gzip
import itertools
import math
from Bio import SeqIO
from Bio.Seq import Seq
from common import ROOT, OUT, GENOMES, SITES, table, write, json_write


def canonical(s):
    return min(s, str(Seq(s).reverse_complement()))


def useful(s):
    counts = Counter(s)
    if set(s) - set('ACGT') or len(counts) < 3:
        return False
    entropy = -sum((n / len(s)) * math.log2(n / len(s)) for n in counts.values())
    longest = max(sum(1 for _ in g) for _, g in itertools.groupby(s))
    return entropy >= 1.5 and longest <= 8


def main():
    alignment = {r.id: str(r.seq) for r in SeqIO.parse(
        ROOT / 'candidate_followup/results/candidate_eight_bee_alignment.faa', 'fasta')}
    exons = table(ROOT / 'candidate_followup/results/genome_candidate_exons.tsv')
    models = table(ROOT / 'candidate_followup/results/genome_candidate_models.tsv')
    columns = {}
    position = 0
    for column, aa in enumerate(alignment['AL_Shangrila']):
        if aa != '-':
            position += 1
            if position in SITES:
                columns[position] = column
    assert len(columns) == len(SITES)
    rows, regions, masked, capture_checks = [], [], {}, []
    baits = set()
    full_regions = {}
    for label, genome in GENOMES.items():
        model = next(r for r in models if r['assembly'] == label)
        selected = [r for r in exons if r['assembly'] == label and r['query'] == model['query']]
        selected.sort(key=lambda r: int(r['exon']))
        assert len(selected) == 14
        strand, scaffold = model['strand'], model['scaffold']
        coords = []
        for exon in selected:
            start, end = int(exon['start_1based']) - 1, int(exon['end_1based'])
            coords += list(range(start, end)) if strand == '+' else list(range(end - 1, start - 1, -1))
        with gzip.open(genome, 'rt') as f:
            record = next(r for r in SeqIO.parse(f, 'fasta') if r.id == scaffold)
        sequence = str(record.seq).upper()
        start, end = max(0, min(coords) - 1000), min(len(sequence), max(coords) + 1001)
        full_regions[label] = sequence[start:end]
        masked_sequence = list(sequence[start:end])
        for focal in sorted(SITES):
            column = columns[focal]
            residue = alignment[label][column]
            assert residue != '-'
            protein_pos = len(alignment[label][:column + 1].replace('-', ''))
            positions = coords[3 * (protein_pos - 1):3 * protein_pos]
            codon = ''.join(sequence[i] for i in positions)
            if strand == '-':
                codon = str(Seq(codon).complement())
            assert str(Seq(codon).translate()) == residue, (label, focal, codon, residue)
            rows.append(dict(assembly=label, focal_AL_position=focal, role=SITES[focal],
                             scaffold=scaffold, strand=strand, protein_position=protein_pos,
                             genomic_codon_positions_1based=','.join(str(i + 1) for i in positions),
                             reference_codon=codon, reference_amino_acid=residue))
            for i in positions:
                masked_sequence[i - start] = 'N'
        masked[label] = ''.join(masked_sequence)
        regions.append(dict(assembly=label, scaffold=scaffold, start_1based=start + 1,
                            end_1based=end, flank_nt=1000, masked_codon_count=len(SITES)))
        for i in range(len(masked[label]) - 20):
            s = masked[label][i:i + 21]
            if useful(s):
                baits.add(canonical(s))
    OUT.mkdir(exist_ok=True)
    with (OUT / 'capture_baits_21.fna').open('w') as f:
        for i, s in enumerate(sorted(baits), 1):
            f.write('>bait_' + str(i) + '\n' + s + '\n')
    with (OUT / 'masked_capture_regions.fna').open('w') as f:
        for label, s in masked.items():
            f.write('>' + label + '\n' + s + '\n')
    # Exhaustively check every fully codon-spanning, error-free 150-nt read,
    # in each assembly, after replacing that codon by all 64 possible codons.
    # This measures local capture sensitivity independently of the real reads.
    for row in rows:
        label = row['assembly']
        region = next(r for r in regions if r['assembly'] == label)
        offset = region['start_1based'] - 1
        p = [int(v) - 1 - offset for v in row['genomic_codon_positions_1based'].split(',')]
        for coding_codon in map(''.join, itertools.product('ACGT', repeat=3)):
            seq = list(full_regions[label])
            bases = coding_codon if row['strand'] == '+' else str(Seq(coding_codon).complement())
            for pos, base in zip(p, bases):
                seq[pos] = base
            seq = ''.join(seq)
            tested = captured = 0
            for start in range(max(p) - 149, min(p) + 1):
                read = seq[start:start + 150]
                assert len(read) == 150
                tested += 1
                captured += any(canonical(read[i:i + 21]) in baits for i in range(130))
            capture_checks.append(dict(assembly=label, focal_AL_position=row['focal_AL_position'],
                                       coding_codon=coding_codon, tested_read_starts=tested,
                                       captured_read_starts=captured))
    write('target_coordinates.tsv', rows)
    write('capture_regions.tsv', regions)
    write('capture_synthetic_sensitivity.tsv', capture_checks)
    json_write(OUT / 'capture_design_summary.json', dict(k=21, distinct_canonical_baits=len(baits),
               masked_codons_per_assembly=len(SITES), entropy_minimum=1.5, maximum_homopolymer=8,
               synthetic_read_starts=sum(r['tested_read_starts'] for r in capture_checks),
               synthetic_captured=sum(r['captured_read_starts'] for r in capture_checks),
               limitation='Error-free 150-nt codon-spanning reads in four assembly backgrounds; does not guarantee capture of unknown highly divergent alleles.'))
    assert all(r['tested_read_starts'] == r['captured_read_starts'] for r in capture_checks), 'Investigate capture gaps'
    print('Mapped', len(rows), 'codons; constructed', len(baits), 'masked baits; all synthetic codon states captured.')


if __name__ == '__main__':
    main()
