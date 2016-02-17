
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import six

hg19_fasta = '/projects/ps-yeolab/genomes/hg19/gencode/v19/GRCh37.p13.genome.fa'

def overlap(x, y):
    return not ((x.start > y.stop) or (x.stop < y.start))

seqrecords = []


splice_type_isoform1_exons = {'SE': ('exon1', 'exon3'), 'MXE': ('exon1', 'exon3', 'exon4')}
splice_type_isoform2_exons = {'SE': ('exon1', 'exon2','exon3'), 'MXE': ('exon1', 'exon2', 'exon4')}

for i, (event_id, row) in enumerate(splicing_feature_data.iterrows()):
    if (i+1) % 1000 == 0:
        six.print_('\t Features completed', i+1)

    isoform1_exons = splice_type_isoform1_exons[row.splice_type]
    isoform2_exons = splice_type_isoform2_exons[row.splice_type]

    isoform_to_exons = {'isoform1': map(lambda x: v19db[row[x]], isoform1_exons),
                        'isoform2': map(lambda x: v19db[row[x]], isoform2_exons)}

    isoform_to_seq = {'isoform1': [], 'isoform2': []}

    reverse = isoform_to_exons['isoform1'][0].strand == '-'

    for isoform, exons in isoform_to_exons.items():
        if reverse:
            sequence = Seq(''.join(exon.sequence(hg19_fasta)[::-1] for exon in exons), generic_dna).complement()
        else:
            sequence = Seq(''.join(exon.sequence(hg19_fasta) for exon in exons), generic_dna)
        transcribed = sequence.transcribe()
        seqrecord = SeqRecord(transcribed, id='{0}|{1}'.format(event_id, isoform))
        seqrecords.append(seqrecord)


with open(transcribed_fasta, 'w') as f:
    SeqIO.write(seqrecords, f, 'fasta')
