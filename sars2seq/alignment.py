from Bio.Seq import Seq

from dark.aligners import mafft
from dark.reads import AARead, Reads

DEBUG = False
PREFIX = 300


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
    """
    def __init__(self, sequence, reference, features, name):
        self.sequence = sequence
        self.reference = reference
        self.features = features
        self.name = name

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
        if DEBUG:
            print('AA MATCH:')
            print('seq  nt:', self.sequence.sequence[:PREFIX])
            print('ref  nt:', self.reference.sequence[:PREFIX])

        features = self.features[self.name]
        try:
            referenceTranslation = features['translation']
        except KeyError:
            referenceTranslation = Seq(features['sequence']).translate()

        # featureLen = len(features['sequence'])
        # assert featureLen % 3 == 0  # Not true for orf1ab

        sequenceRead = AARead(
            self.sequence.id,
            Seq(self.sequence.sequence.replace('-', '')).translate())

        referenceRead = AARead(self.reference.id, referenceTranslation)

        if DEBUG:
            print('seq  aa:', sequenceRead.sequence[:PREFIX])
            print('ref  aa:', referenceRead.sequence[:PREFIX])

        alignment = mafft(Reads([sequenceRead, referenceRead]),
                          options='--preservecase --amino')

        sequenceResult, referenceResult = list(alignment)

        if (sequenceResult.sequence.startswith('-') or
                sequenceResult.sequence.endswith('-')):
            if DEBUG:
                print('Sequence result has leading/trailing gaps!')

        # In the alignment, the reference will have a bunch of leading and
        # trailing '-' chars. The length of these gives us the offsets in
        # the original full-length sequence where the match is.
        offset = min(
            (len(referenceResult) - len(referenceResult.sequence.lstrip('-'))),
            (len(sequenceResult) - len(sequenceResult.sequence.lstrip('-'))))
        length = len(referenceResult.sequence.strip('-'))

        sequenceTrimmed = sequenceResult[offset:offset + length]

        if DEBUG:
            print('seq res:', sequenceTrimmed.sequence[:PREFIX])
            print('ref res:', referenceRead.sequence[:PREFIX])
            print(f'offset, length = {offset}, {length}')
            print('AA MATCH END')

        assert len(sequenceTrimmed) == len(referenceRead), (
            f'{len(sequenceTrimmed)} != {len(referenceRead)}')

        return sequenceTrimmed, referenceRead
