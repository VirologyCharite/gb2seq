from Bio.Seq import Seq

from dark.aligners import mafft
from dark.reads import AARead, Reads

DEBUG = False
SLICE = slice(300)


def translate(seq, name):
    """
    Translate a sequence.

    @param seq: A C{str} nucelotide sequence.
    @param name: A C{str} feature name.
    @return: A translated C{str} amino acid sequence.
    """
    if name == 'ORF1ab polyprotein':
        # See Fields Virology (figure 10.6a on page 421, 7th edition or
        # figure 28.7a on page 836, 6th edition) plus
        # https://www.ncbi.nlm.nih.gov/nuccore/NC_045512 for details of
        # what happens below. Note that the nucelotide sequence we are
        # passed is the one that's made from the alignment with the
        # reference ORF1ab nucleotide sequence (in sequence.py) and so is
        # just that ORF and does not include the leading ~265 nucleotides
        # of the 5' UTR. As a result, the offset used to begin the search
        # for the slippery sequence is 13000, which is chosen to be a bit
        # before 13468 - 265.  There are various occurrences of the
        # slippery sequence in the reference genome (and hence probably in
        # other CoV genomes), but only one in this region and with a stop
        # codon shortly (<20 nt) downstream.
        slipperySeq = 'TTTAAAC'
        slipperyLen = len(slipperySeq)
        offset = seq.find(slipperySeq, 13000)
        stop = seq.find('TAA', offset + slipperyLen)
        if DEBUG:
            print(f'LEN: {len(seq)}, OFFSET: {offset}, STOP: {stop}')
        assert offset > -1 and stop > -1 and stop - offset < 20
        seq = seq[:offset + slipperyLen] + seq[offset:]

    # Pad with 'N' to avoid a 'BiopythonWarning: Partial codon' warning.
    remainder = len(seq) % 3
    seq += 'N' * (3 - remainder if remainder else 0)

    return Seq(seq).translate()


class Alignment:
    """
    Hold information about a matched feature.

    @param sequence: A C{dark.reads.DNARead} instance holding the sequence
        of the feature of interest. This has been aligned with C{reference}
        and may contain '-' gap characters.
    @param reference: A C{dark.reads.DNARead} instance holding the
        corresponding sequence extracted from a reference genome. This has been
        aligned with C{reference} and may contain '-' gap characters.
    @param features: A features C{dict}.
    @param name: The C{str} name of the feature.
    @param sequenceOffset: The C{int} offset of where the feature was found in
        C{sequence}.
    """
    def __init__(self, sequence, reference, features, name, sequenceOffset):
        self.sequence = sequence
        self.reference = reference
        self.features = features
        self.name = name
        self.sequenceOffset = sequenceOffset

    def ntSequences(self):
        """
        Get the aligned nucelotide sequences.

        @return: A 2-C{tuple} of C{dark.reads.DNARead} instances, holding
            the nucleotides of the sequence that is being examined and those
            of the reference sequence.
        """
        return self.sequence, self.reference

    def aaSequences(self):
        """
        Match the sequence and the reference at the amino acid level.

        @return: A 2-C{tuple} of C{dark.reads.AARead} instances, holding
            the amino acids of the sequence that is being examined and those
            of the reference sequence.
        """
        name = self.name

        sequenceRead = AARead(
            self.sequence.id,
            translate(self.sequence.sequence.replace('-', ''), name))

        features = self.features[name]
        referenceTranslation = features.get(
            'translation', translate(features['sequence'], name))
        referenceRead = AARead(self.reference.id, referenceTranslation)

        if DEBUG:
            print(f'AA MATCH {name}:')
            print('seq  nt:', self.sequence.sequence[SLICE])
            print('ref  nt:', self.reference.sequence[SLICE])
            print('seq  aa:', sequenceRead.sequence[SLICE])
            print('ref  aa:', referenceRead.sequence[SLICE])

        return align(sequenceRead, referenceRead, nt=False)


def align(sequenceRead, referenceRead, nt=True):
    """
    Align two sequences to find where the reference read fits in the sequence
    read.

    @param sequenceRead: A C{dark.reads.Read} instance.
    @param referenceRead: A C{dark.reads.Read} instance.
    @param nt: If C{True} the sequences are nucleotide. Else protein.
    @return: A 2-C{tuple} of aligned C{dark.reads.Read} instances, with the
        sequence result and the reference result. Both may contain gaps ('-').
    """
    alignment = mafft(
        Reads([sequenceRead, referenceRead]),
        options='--anysymbol --preservecase' + (
            ' --nuc' if nt else ' --amino'))

    sequenceResult, referenceResult = list(alignment)
    if DEBUG:
        print('process alignment')

    if (sequenceResult.sequence.startswith('-') or
            sequenceResult.sequence.endswith('-')):
        if DEBUG:
            print('Sequence result has leading/trailing gaps!')

    # One but not both of the sequence and reference can start with a
    # gap.  Same goes for ends. The reason is that an aligner has no
    # reason to put gaps at the extremes of both sequences at once.
    assert not (sequenceResult.sequence.startswith('-') and
                referenceResult.sequence.startswith('-'))
    assert not (sequenceResult.sequence.endswith('-') and
                referenceResult.sequence.endswith('-'))

    # In the alignment, the reference will have a bunch of leading and
    # trailing '-' chars. The length of these gives us the offsets in
    # the original full-length sequence where the match is.

    # offset = max(
    #     (len(referenceResult) - len(referenceResult.sequence.lstrip('-'))),
    #     (len(sequenceResult) - len(sequenceResult.sequence.lstrip('-'))))
    # length = len(referenceResult.sequence.strip('-'))
    # assert length == len(referenceRead)

    # sequenceTrimmed = sequenceResult[offset:offset + length]

    if DEBUG:
        print('seq   al', sequenceResult.sequence[SLICE])
        print('ref   al', referenceResult.sequence[SLICE])
        with open('/tmp/align-%s.fasta' % ('nt' if nt else 'aa'), 'w') as fp:
            print('seq', sequenceResult.sequence, file=fp)
            print('res', referenceResult.sequence, file=fp)

    if referenceResult.sequence.startswith('-'):
        offset = (len(referenceResult) -
                  len(referenceResult.sequence.lstrip('-')))
        sequenceResult = sequenceResult[offset:]
        referenceResult = referenceResult[offset:]

        if DEBUG:
            print(f'CLIPPED {offset} on left')
            print('seq   al', sequenceResult.sequence[SLICE])
            print('ref   al', referenceResult.sequence[SLICE])

    if referenceResult.sequence.endswith('-'):
        offset = (len(referenceResult) -
                  len(referenceResult.sequence.rstrip('-')))

        if DEBUG:
            print(f'CLIPPING {offset} on right')
            print('seq   al', sequenceResult.sequence)
            print('ref   al', referenceResult.sequence)

        sequenceResult = sequenceResult[:len(referenceResult) - offset]
        referenceResult = referenceResult[:len(referenceResult) - offset]

        if DEBUG:
            print(f'CLIPPED {offset} on right')
            print('seq   al', sequenceResult.sequence[SLICE])
            print('ref   al', referenceResult.sequence[SLICE])

    if DEBUG:
        print('seq  end', sequenceResult.sequence[SLICE])
        print('ref  end', referenceResult.sequence[SLICE])
        print('END process alignment')

    return sequenceResult, referenceResult
