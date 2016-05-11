
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna

import pandas as pd

import numpy as np

hg19_fasta = '/projects/ps-yeolab/genomes/hg19/gencode/v19/GRCh37.p13.genome.fa'

splice_type_isoform_exons = {'SE': {'isoform1': ('exon1', 'exon3'),
                                    'isoform2': ('exon1', 'exon2', 'exon3')},
                             'MXE': {'isoform1': ('exon1', 'exon3', 'exon4'),
                                     'isoform2': ('exon1', 'exon2', 'exon4')}
                            }

best_tags = 'appris_principal', 'appris_candidate', 'CCDS', 'exp_conf', 'basic'

class SplicingTranslator(object):
    pass


def overlap(x, y):
    return not ((x.start > y.stop) or (x.stop < y.start))


def filter_on_tags(features, tags):
    if len(features) == 0:
        return None

    if len(features) == 1:
        return features.pop()

    feature_tags = map(lambda x: (x, tuple(x.attributes['tag']))
                    if 'tag' in x.attributes else (x, ()), features)
    for tag in tags:
        for feature, feature_tag in feature_tags:
            if tag in feature_tag:
                return feature

    return np.random.choice(list(features))

seqrecords = []

for i, (event_id, row) in enumerate(splicing_feature_data.iterrows()):
#     if i > 10:
#         break

    if (i+1) % 1000 == 0:
        print i+1
#     exon1 = v19db[row.exon1]
#     exon2 = v19db[row.exon2]
#     exon3 = v19db[row.exon3]
#     print event_id

#     exon_trio = exon1, exon2, exon3

    isoform1_exons = splice_type_isoform_exons[row.splice_type]['isoform1']
    isoform2_exons = splice_type_isoform_exons[row.splice_type]['isoform2']

    isoform_to_exons = {'isoform1': map(lambda x: v19db[row[x]], isoform1_exons),
                        'isoform2': map(lambda x: v19db[row[x]], isoform2_exons)}

    isoform_to_transcripts = dict((k, set.intersection(*map(lambda x: set(v19db.parents(x, featuretype='transcript')), v)))
                                  for k, v in isoform_to_exons.items())

#     isoform_to_exons = {'isoform1': (exon1, exon3), 'isoform2': exon_trio}


    # Make sure isoform1 transcripts are distinct from isoform2
    isoform1_transcripts = isoform_to_transcripts['isoform1'] - isoform_to_transcripts['isoform2']
#     print '\tisoform1_transcripts', ','.join(map(lambda x: x.id, isoform1_transcripts))
#     print '\t', '\t'.join(map(lambda x: str(list(x.attributes['tag'])) if 'tag' in x.attributes else '-', isoform1_transcripts))

    # Make sure isoform2 transcripts are distinct from isoform1
    isoform2_transcripts = isoform_to_transcripts['isoform2']
    if row.splice_type == 'MXE':
        isoform2_transcripts = isoform2_transcripts - isoform1_transcripts
#     print '\tisoform2_transcripts', ','.join(map(lambda x: x.id, isoform2_transcripts))
#     print '\t', '\t'.join(map(lambda x: str(list(x.attributes['tag'])) if 'tag' in x.attributes else '-', isoform2_transcripts))

    isoform1_transcript = filter_on_tags(isoform1_transcripts, best_tags)
    isoform2_transcript = filter_on_tags(isoform2_transcripts, best_tags)

#     if isoform1_transcript is not None:
#         print '\t--isoform1', isoform1_transcript.id,
#         if 'tag' in isoform1_transcript.attributes:
#             print isoform1_transcript.attributes['tag']
#         else:
#             print
#     if isoform2_transcript is not None:
#         print '\t--isoform2', isoform2_transcript.id,
#         if 'tag' in isoform2_transcript.attributes:
#             print isoform2_transcript.attributes['tag']
#         else:
#             print

#     isoforms = {'isoform1': isoform1_transcripts,
#                 'isoform2': isoform2_transcripts}
    isoforms = {'isoform1': isoform1_transcript,
                'isoform2': isoform2_transcript}


#     isoform_to_cds = {'isoform1': [], 'isoform2': []}
    for isoform, transcript in isoforms.items():
        if transcript is None:
#             print '\t\t', isoform, 'is None'
            continue
        exons = isoform_to_exons[isoform]
#         print '\t', '\t'.join(map(lambda x: str(list(x.attributes['tag'])) if 'tag' in x.attributes else '-', exons))

#         for transcript in transcripts:


        reverse = transcript.strand == '-'
        cdss = db.children(transcript, featuretype='CDS', order_by='start',
                              reverse=reverse)

        cdss = filter(lambda cds: any(map(lambda exon: overlap(cds, exon), exons)), cdss)
#         print cdss
        if len(cdss) == len(exons):
            cds_str = '@'.join(map(lambda x: x.id, cdss))

            if reverse:
                coding_sequence = Seq(
                    ''.join(cds.sequence(hg19_fasta)[::-1]
                            for cds in cdss), generic_dna).complement()
            else:
                coding_sequence = Seq(''.join(cds.sequence(hg19_fasta)
                                              for cds in cdss), generic_dna)
            coding_sequence = coding_sequence[int(cdss[0].frame):]
            translated = coding_sequence.translate()
            seqrecord = SeqRecord(translated, id='{0}|{1}|{2}'.format(
                event_id, cds_str, isoform))
            seqrecords.append(seqrecord)


with open(translated_fasta, 'w') as f:
    SeqIO.write(seqrecords, f, 'fasta')


translated_df = pd.DataFrame(columns=['isoform1', 'isoform2'])


with open(translated_fasta) as f:
    for record in SeqIO.parse(f, 'fasta'):
        event_id, cds_str, isoform = record.id.split('|')
        translated_df.loc[event_id, isoform] = str(record.seq)

translated_df.columns = [x + '_translation' for x in translated_df]

translated_df.to_csv(
    '{}/protein_translations.csv'.format(alternative_feature_folder))
