"""Exercise the actual BBDuk invocation on known positive and negative pairs."""
import itertools
import random
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from common import WORK, OUT, table, json_write
from setup_tools import classpath
from collect_reads import restore_original_records


def main():
    WORK.mkdir(exist_ok=True)
    regions = {r['assembly']: r for r in table(OUT / 'capture_regions.tsv')}
    targets = table(OUT / 'target_coordinates.tsv')
    masked = {r.id: str(r.seq) for r in SeqIO.parse(OUT / 'masked_capture_regions.fna', 'fasta')}
    baits = {str(r.seq) for r in SeqIO.parse(OUT / 'capture_baits_21.fna', 'fasta')}
    def hit(s):
        return any(min(s[i:i+21], str(Seq(s[i:i+21]).reverse_complement())) in baits for i in range(130))
    rng = random.Random(4701)
    def background():
        while True:
            s = ''.join(rng.choices('ACGT', k=150))
            if not hit(s):
                return s
    inputs = [WORK / 'capture_test_1.fastq', WORK / 'capture_test_2.fastq']
    outputs = [WORK / 'capture_test_selected_1.fastq', WORK / 'capture_test_selected_2.fastq']
    original, expected = {}, set()
    counter = 0
    with inputs[0].open('w') as f1, inputs[1].open('w') as f2:
        for row in targets:
            label = row['assembly']
            offset = int(regions[label]['start_1based']) - 1
            positions = [int(v)-1-offset for v in row['genomic_codon_positions_1based'].split(',')]
            for codon in map(''.join, itertools.product('ACGT', repeat=3)):
                bases = codon if row['strand'] == '+' else str(Seq(codon).complement())
                sequence = list(masked[label])
                for pos, base in zip(positions, bases):
                    sequence[pos] = base
                sequence = ''.join(sequence)
                for displacement in [0, 70, 147]:
                    start = min(positions) - displacement
                    read = sequence[start:start+150]
                    assert len(read) == 150 and hit(read)
                    counter += 1
                    name = 'positive_' + str(counter)
                    mate = background()
                    if counter % 2:
                        read = str(Seq(read).reverse_complement())
                    for f, end, s in [(f1, 1, read), (f2, 2, mate)]:
                        record = f'@{name}/{end}\n{s}\n+\n' + 'I'*150 + '\n'
                        f.write(record)
                        original[(name, end)] = record
                    expected.add(name)
        for n in range(100):
            for f, end in [(f1, 1), (f2, 2)]:
                f.write(f'@negative_{n}/{end}\n{background()}\n+\n' + 'I'*150 + '\n')
    args = [f'in={inputs[0]}', f'in2={inputs[1]}', f'outm={outputs[0]}', f'outm2={outputs[1]}',
            f'ref={OUT / "capture_baits_21.fna"}', 'k=21', 'hdist=0', 'maskmiddle=f',
            'rcomp=t', 'minkmerhits=1', 'ordered=t', 'threads=2', 'overwrite=t',
            'qtrim=f', 'minlength=1', 'maxns=-1', 'usejni=f', 'pigz=f', 'unpigz=f',
            'qin=33', 'qout=33', 'changequality=f']
    with (WORK / 'capture_tool_test.log').open('w') as f:
        subprocess.run(['java', '-Xmx1g', '-cp', str(classpath()), 'jgi.BBDuk', *args],
                       stdout=f, stderr=subprocess.STDOUT, check=True)
    restored = [WORK / 'capture_test_restored_1.fastq', WORK / 'capture_test_restored_2.fastq']
    for source, captured, output in zip(inputs, outputs, restored):
        restore_original_records(source, captured, output)
    observed = []
    for end, path in enumerate(restored, 1):
        names = set()
        with path.open() as f:
            while first := f.readline():
                rest = [f.readline() for _ in range(3)]
                name = first[1:].strip().rsplit('/', 1)[0]
                assert first + ''.join(rest) == original[(name, end)], (name, 'Sequence or quality changed')
                names.add(name)
        assert names == expected, (end, len(names), len(expected))
        observed.append(len(names))
    json_write(OUT / 'capture_tool_validation.json', dict(bbduk_version='39.91',
               positive_pairs=counter, negative_pairs=100, retained_R1=observed[0], retained_R2=observed[1],
               all_positive_pairs_retained=True, all_negative_pairs_excluded=True,
               unmatched_mates_preserved=True, read_sequences_and_qualities_preserved=True,
               limitation='Synthetic method validation; not biological replication or a guarantee for unknown divergent sequence.'))
    print('Actual BBDuk control passed:', counter, 'positive pairs; 100 negative pairs; both mates unchanged.')


if __name__ == '__main__':
    main()
