
import pandas as pd


def get_first_last_internal(exon_id, db):
    exon_i = None

    exon = db[exon_id]
    reverse = exon.strand=='-'

    first_internals = []
    last_internals = []
    for transcript in db.parents(exon, featuretype='transcript'):
        children = db.children(transcript, featuretype='exon',
                                  order_by='start', reverse=reverse)
        for i, child in enumerate(children):
            if child == exon:
                exon_i = i
        first_internal = (exon_i == 1)
        last_internal = (exon_i == (i - 1))
        first_internals.append(first_internal)
        last_internals.append(last_internal)
    return pd.Series(dict(last_internal=any(last_internals),
                          first_internal=any(first_internals)))
