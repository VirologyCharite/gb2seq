from dark.reads import DNARead

from sars2seq.alignment import Alignment, align

DEBUG = False
SLICE = slice(300)

# Alias names should have a lower case key.
ALIASES = {
    '2': "2'-O-ribose methyltransferase",
    '3clpro': '3C-like proteinase',
    '3utr': "3'UTR",
    '5utr': "5'UTR",
    'e': 'envelope protein',
    'endornase': 'endoRNAse',
    'envelope': 'envelope protein',
    'exon': "3'-to-5' exonuclease",
    'exonuclease': "3'-to-5' exonuclease",
    'leader': 'leader protein',
    'm': 'membrane glycoprotein',
    'membrane': 'membrane glycoprotein',
    'mpro': '3C-like proteinase',
    'n': 'nucleocapsid phosphoprotein',
    'nsp2': 'nsp2',
    'nsp3': 'nsp3',
    'nsp4': 'nsp4',
    'nsp6': 'nsp6',
    'nsp7': 'nsp7',
    'nsp8': 'nsp8',
    'nsp9': 'nsp9',
    'nsp10': 'nsp10',
    'nsp11': 'nsp11',
    'orf1a': 'ORF1a polyprotein',
    'orf1ab': 'ORF1ab polyprotein',
    'orf3a': 'ORF3a protein',
    'orf6': 'ORF6 protein',
    'orf7a': 'ORF7a protein',
    'orf7b': 'ORF7b',
    'orf8': 'ORF8 protein',
    'orf10': 'ORF10 protein',
    'rdrp': 'RNA-dependent RNA polymerase',
    's': 'surface glycoprotein',
    'spike': 'surface glycoprotein',
    'stem loop 1': 'stem loop 1',
    'stem loop 2': 'stem loop 2',
    'stem loop 3': 'stem loop 3',
    'stem loop 4': 'stem loop 4',
    'stem loop 5': 'stem loop 5',
    'surface glycoprotein': 'surface glycoprotein',
}


class SARS2Sequence:
    """
    Methods for working with SARS-CoV-2 sequences.

    @param sequence: A C{dark.reads.Read} instance.
    @param features: An C{Features} instance.
    """
    def __init__(self, sequence, features, window=500):
        self.sequence = sequence
        self.features = features
        self.window = window
        self._features = features.featuresDict()
        self._results = {}

    def featureNames(self):
        """
        Get the list of feature names.

        @return: A C{list} of C{str} feature names.
        """
        return list(self._features)

    def feature(self, name):
        """
        Get a feature (e.g., nsp2).

        @param name: A C{str} feature name.
        @raise KeyError: If the feature is unknown.
        @return: A C{dict} with information about the feature.
        """
        if name not in self._features:
            alt = ALIASES.get(name.lower())
            if alt is None:
                raise KeyError(name)
            else:
                name = alt

        offsets = self._features[name]['start'], self._features[name]['stop']

        try:
            result = self._results[offsets]
        except KeyError:
            result = self._results[offsets] = self._findFeature(name)

        return result

    def _findFeature(self, name):
        """
        Find a reference feature in the sequence.

        @param name: A C{str} feature name.
        """
        sequenceSubsequence = self.sequence.sequence[
            max(0, self._features[name]['start'] - self.window):
            self._features[name]['stop'] + self.window]

        sequenceRead = DNARead(self.sequence.id + f' ({name})',
                               sequenceSubsequence)
        referenceRead = DNARead(self.features.id + f' ({name})',
                                self._features[name]['sequence'])
        assert len(sequenceRead) >= len(referenceRead)

        if DEBUG:
            print('NT MATCH:')
            print('seq  nt:', sequenceRead.sequence[SLICE])
            print('ref  nt:', referenceRead.sequence[SLICE])

        sequenceResult, referenceResult = align(sequenceRead, referenceRead)
        sequenceNoGaps = sequenceResult.sequence.replace('-', '')
        sequenceOffset = self.sequence.sequence.find(sequenceNoGaps)
        assert sequenceOffset > -1

        return Alignment(sequenceResult, referenceResult, self._features, name,
                         sequenceOffset)
