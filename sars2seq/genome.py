from sars2seq.alignment import Alignment

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
    'orf10': 'ORF10 protein',
    'orf1a': 'ORF1a polyprotein',
    'orf1ab': 'ORF1ab polyprotein',
    'orf3a': 'ORF3a protein',
    'orf6': 'ORF6 protein',
    'orf7a': 'ORF7a protein',
    'orf7b': 'ORF7b',
    'orf8': 'ORF8 protein',
    'rdrp': 'RNA-dependent RNA polymerase',
    's': 'surface glycoprotein',
    'sl1': 'stem loop 1',
    'sl2': 'stem loop 2',
    'sl3': 'stem loop 3',
    'sl4': 'stem loop 4',
    'sl5': 'stem loop 5',
    'spike': 'surface glycoprotein',
    'surface glycoprotein': 'surface glycoprotein',
}


class SARS2Genome:
    """
    Methods for working with SARS-CoV-2 genomes.

    @param genome: A C{dark.reads.Read} instance.
    @param features: An C{Features} instance.
    """
    def __init__(self, genome, features):
        self._genome = genome
        self.features = features
        self._featuresDict = features.featuresDict()

    def feature(self, name):
        """
        Get a feature (e.g., nsp2).

        @param name: A C{str} feature name.
        @raise KeyError: If the feature is unknown.
        @return: An C{Alignment} instance, providing access to the found
            feature in the genome and its corresponding sequence in the
            reference.
        """
        if name not in self._featuresDict:
            nameLower = name.lower()
            for featureName in self._featuresDict:
                if nameLower == featureName.lower():
                    name = featureName
                    break
            else:
                alt = ALIASES.get(nameLower)
                if alt is None:
                    raise KeyError(name)
                else:
                    name = alt

        return Alignment(self._genome, self.features.id,
                         self._featuresDict[name])
