def miso_exon_to_gencode_exon(exon):
    """Convert a single miso exon to one or more gffutils database exon id

    >>> # Skipped exon (SE) or Mutually exclusive exon (MXE) ID
    >>> miso_exon_to_gencode_exon('chr2:9624561:9624679:+')
    'exon:chr2:9624561-9624679:+'
    >>> # Alt 5' splice site (pick first of alternative)
    >>> miso_exon_to_gencode_exon('chr15:42565276:42565087|42565161:-')
    'exon:chr15:42565276-42565087:-'
    >>> # Alt 3' splice site (pick first of alternative)
    >>> miso_exon_to_gencode_exon('chr2:130914199|130914248:130914158:-')
    'exon:chr2:130914199-130914158:-'
    >>> # Retained intron: start, stop separated by '-' instead of ':'
    >>> miso_exon_to_gencode_exon('chr1:906259-906386:+')
    'exon:chr1:906259-906386:+'
    """
    return 'exon:{}:{}-{}:{}'.format(*miso_exon_to_coords(exon))


def miso_id_to_exon_ids(miso_id):
    """Convert a MISO-style alternative event ID to a gffutils exon id of
    all exons in all possible transcripts

    Split on the pipe ("|") to account for Alt 5'/3' splice site events

    # A skipped exon (SE) ID
    >>> miso_id_to_exon_ids('chr2:9624561:9624679:+@chr2:9627585:9627676:+@chr2:9628276:9628591:+')
    ['exon:chr2:9624561-9624679:+', 'exon:chr2:9627585-9627676:+', 'exon:chr2:9628276-9628591:+']
    >>> # A mutually exclusive (MXE) ID
    >>> miso_id_to_exon_ids('chr16:89288500:89288591:+@chr16:89289565:89289691:+@chr16:89291127:89291210:+@chr16:89291963:89292039:+')
    ['exon:chr16:89288500-89288591:+', 'exon:chr16:89289565-89289691:+', 'exon:chr16:89291127-89291210:+', 'exon:chr16:89291963-89292039:+']
    >>> # An Alt 5' splice site (A5SS) ID
    >>> miso_id_to_exon_ids("chr15:42565276:42565087|42565161:-@chr15:42564261:42564321:-")
    ['exon:chr15:42565276-42565161:-', 'exon:chr15:42564261-42564321:-']
    >>> # An Alt 3' splice site (A3SS) ID
    >>> miso_id_to_exon_ids('chr2:130914824:130914969:-@chr2:130914199|130914248:130914158:-')
    ['exon:chr2:130914824-130914969:-', 'exon:chr2:130914199-130914158:-']
    >>> # A retained intron (RI) ID
    >>> miso_id_to_exon_ids('chr1:906066-906138:+@chr1:906259-906386:+')
    ['exon:chr1:906066-906138:+', 'exon:chr1:906259-906386:+', 'exon:chr1:906066-906386:+']
    """
    return map(miso_exon_to_gencode_exon, miso_id.split('@'))


def miso_exon_to_coords(exon):
    """Convert a miso exon to gffutils coordinates

    >>> miso_exon_to_coords('chr2:130914824:130914969:-')
    ('chr2', '130914824', '130914969', '-')
    >>> # Alt 5' SS - pick the first of the alternative ends
    >>> miso_exon_to_coords('chr15:42565276:42565087|42565161:-')
    ('chr15', '42565276', '42565087', '-')
    >>> # Alt 3' SS - pick the first of the alternative starts
    >>> miso_exon_to_coords('chr2:130914199|130914248:130914158:-')
    ('chr2', '130914199', '130914158', '-')
    >>> # Retained intron
    >>> miso_exon_to_coords('chr1:906066-906138:+')
    ('chr1', '906066', '906138', '+')
    """
    strand = exon[-1]
    coords = map(lambda x: x.split('|')[0],
                 exon.split(':'))
    if '-' in coords[1]:
        start, stop = coords[1].split('-')
        coords = coords[0], start, stop, strand
    return coords[0], coords[1], coords[2], strand
