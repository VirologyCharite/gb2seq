from Bio.Seq import Seq
from itertools import groupby

from dark.aa import CODONS, STOP_CODONS


class TranslationError(Exception):
    'Error when using custom translation of sequences.'


class NoSlipperySequenceError(TranslationError):
    'No slippery sequence could be found in a genome.'


class NoStopCodonError(TranslationError):
    'No stop codon was found downstream from the slippery sequence.'


class StopCodonTooDistantError(TranslationError):
    'The stop codon following the slippery sequence was too far away.'


class TranslatedSequenceLengthError(TranslationError):
    'The translated reference and genome have different lengths.'


class TranslatedReferenceAndGenomeLengthError(TranslationError):
    'The translated reference and genome have different lengths.'


class TranslatedGapLengthError(TranslationError):
    'The sequence to translate has stretches of gaps that are not '
    'multiples of three.'


codons = dict([(codon, aa) for aa, cods in CODONS.items() for codon in cods] +
              [(codon, '*') for codon in STOP_CODONS] + [('---', '-')])


# The maximum difference (number of nucleotides) to allow between the
# offset of the start of the slippery sequence and the downstream stop
# codon.
_MAX_DISTANCE_TO_STOP = 20

SLIPPERY_SEQUENCE = 'TTTAAAC'

_SLIPPERY_LEN = len(SLIPPERY_SEQUENCE)


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
        # before 13462 - 265.  There are various occurrences of the
        # slippery sequence in the reference genome (and hence probably in
        # other CoV genomes), but only one in this region and with a stop
        # codon shortly (up to _MAX_DISTANCE_TO_STOP nucleotides) downstream.
        offset = seq.find(SLIPPERY_SEQUENCE, 13000)
        stop = seq.find('TAA', offset + _SLIPPERY_LEN)
        if offset == -1:
            raise NoSlipperySequenceError('No slippery sequence found.')
        if stop == -1:
            raise NoStopCodonError(
                f'Could not find a stop codon downstream from the start of '
                f'the slippery sequence at location {offset + 1}.')
        if stop - offset > _MAX_DISTANCE_TO_STOP:
            raise StopCodonTooDistantError(
                f'The stop codon was too far ({stop - offset} nucleotides) '
                f'downstream (max allowed distance is '
                f'{_MAX_DISTANCE_TO_STOP}) from the start of the slippery '
                f'sequence at location {offset + 1}.')

        seq = seq[:offset + _SLIPPERY_LEN] + seq[offset + _SLIPPERY_LEN - 1:]

    # Pad with 'N' to avoid a 'BiopythonWarning: Partial codon' warning.
    remainder = len(seq) % 3
    seq += 'N' * (3 - remainder if remainder else 0)

    return Seq(seq).translate()


def translateSpike(seq):
    """
    Translate a Spike sequence, taking into account gaps introduced in the
    nucleotide alignment. This means that the amino acid sequences do not
    have to be re-aligned after translating, which avoids the introduction of
    gaps at the amino acid level that may be different from gaps at the
    nucleotide level.

    @param seq: A C{str} nucelotide sequence.
    @return: A translated C{str} amino acid sequence that retains the gaps.
    """
    sequence = ''
    current = 0
    seqLen = len(seq)

    if not seqLen % 3 == 0:
        raise TranslatedSequenceLengthError(
            f'The length of a sequence to be translated must '
            f'be a multiple of 3 but is {seqLen!r}.')

    groups = groupby(seq)
    result = [(label, len(list(group))) for label, group in groups if
              label == '-']

    if any(length % 3 != 0 for g, length in result):
        raise TranslatedGapLengthError(
            'Length of stretch of gaps not divisible by 3.')

    while current + 3 <= seqLen:
        codon = seq[current:current + 3]
        codonGaps = codon.count('-')
        if codon in codons:
            # This is a codon that corresponds to a normal amino acid.
            sequence += codons[codon]
            current += 3
        elif codon not in codons and codonGaps == 0:
            # This is a codon that contains ambiguities.
            sequence += 'X'
            current += 3
        elif codonGaps > 0:
            # Count how many gaps there are after the current codon.
            c = 3
            while seq[current + c] == '-':
                c += 1
            subsequentGaps = c - 3

            # Find the next nucleotide after the gap.
            index = 3 + codonGaps
            nextNt = seq[current + subsequentGaps + 3:
                         current + subsequentGaps + index]

            newCodon = codon.strip('-') + nextNt

            sequence += codons.get(newCodon, 'X')

            # Add the correct number of gaps after the amino acid.
            totalGaps = subsequentGaps + codonGaps
            gapAA = totalGaps // 3
            sequence += (gapAA * '-')

            # Set current so it starts at the right place after the gap.
            current += (subsequentGaps + index)

    return sequence


def checkSpikeInsertions(accession, seq):
    """
    Check and parse out known insertions in the Spike protein.

    @param accession: The C{str} accession number of the sequence.
    @param seq: A C{str} amino acid sequence.

    @return: The C{str} corrected translated sequence.
    """
    seqLen = len(seq)
    if seqLen == 1274:
        # There are no insertions
        return seq
    if seqLen > 1274:
        if seq.find('QTKGIALSPR', 650, 700) > -1:
            return seq[:679] + seq[683:]
        elif seq.find('NLVRTDRDLPQ', 200, 230) > -1:
            return seq[:214] + seq[217:]
        elif seq.find('LVRAAGYLPQ', 200, 230) > -1:
            return seq[:214] + seq[217:]
        else:
            raise TranslatedSequenceLengthError(
                f'Sequence with accession {accession} is too long '
                f'(length: {seqLen}) and does not have a known insertion.',
            )
    raise TranslatedSequenceLengthError(
        f'Sequence with accession {accession} is too short '
        f'(length: {seqLen}.'
    )
