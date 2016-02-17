from outrigger.region import Region

def get_gene_attributes(df, exon_cols, db):
    df['gencode_id'] = df[exon_cols].apply(
        lambda x: ','.join(set(itertools.chain(*[db[i].attributes['gene_id'] for i in x]))), axis=1)
    df['transcript_id'] = df[exon_cols].apply(
        lambda x: ','.join(set(itertools.chain(*[db[i].attributes['transcript_id'] for i in x]))), axis=1)
    df['gene_name'] = df[exon_cols].apply(
        lambda x: ','.join(set(itertools.chain(*[db[i].attributes['gene_name'] for i in x]))), axis=1)
    df['ensembl_id'] = df[exon_cols].apply(
        lambda x: ','.join(set(itertools.chain(*[map(lambda y: y.split('.')[0],
                                                     db[i].attributes['gene_id']) for i in x]))),
        axis=1)

    for exon_col in exon_cols:
        df['{}_region'.format(exon_col)] = df[exon_col].map(Region)
        df['{}_length'.format(exon_col)] = df['{}_region'.format(exon_col)].map(len)
    df['strand'] = df[exon_cols[0]].str[-1]
    return df
