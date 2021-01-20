from Bio.Seq import Seq


class TranslationError(Exception):
    'No slippery sequence could be found in a genome.'


class NoSlipperySequenceError(TranslationError):
    'No slippery sequence could be found in a genome.'


class NoStopCodonError(TranslationError):
    'No stop codon was found downstream from the slippery sequence.'


class StopCodonTooDistantError(TranslationError):
    'The stop codon following the slippery sequence was too far away.'


# The maximum difference (number of nucleotides) to allow between the
# offset of the start of the slippery sequence and the downstream stop
# codon.
_MAX_DISTANCE_TO_STOP = 20


def translate(seq, name=None):
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
        # codon shortly (up to _MAX_DISTANCE_TO_STOP nucleotides) downstream.
        slipperySeq = 'TTTAAAC'
        slipperyLen = len(slipperySeq)
        offset = seq.find(slipperySeq, 13000)
        stop = seq.find('TAA', offset + slipperyLen)
        if offset == -1:
            raise NoSlipperySequenceError('No slippery sequence found.')
        if stop == -1:
            raise NoStopCodonError(
                f'Could not find a stop codon downstream from the start of '
                f'the slippery sequence at location {offset + 1}.')
        if stop - offset > _MAX_DISTANCE_TO_STOP:
            raise StopCodonTooDistantError(
                f'The stop codon was too far ({stop - offset} nucleotides) '
                f'downstream (max allowed distance is {_MAX_DISTANCE_TO_STOP} '
                f'from the start of the slippery sequence at location '
                f'{offset + 1}.')
        seq = seq[:offset + slipperyLen] + seq[offset:]

    # Pad with 'N' to avoid a 'BiopythonWarning: Partial codon' warning.
    remainder = len(seq) % 3
    seq += 'N' * (3 - remainder if remainder else 0)

    return Seq(seq).translate()
