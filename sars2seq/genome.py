from dark.reads import DNARead

from sars2seq.alignment import Alignment
from sars2seq.variants import VARIANTS


class SARS2Genome:
    """
    Methods for working with SARS-CoV-2 genomes.

    @param genome: A C{dark.reads.Read} instance.
    @param features: An C{Features} instance.
    """
    def __init__(self, genome, features):
        # Geneious uses ? to indicate unknown nucleotides, at least when it
        # exports a consensus / alignment. Replace with N.
        self.genome = DNARead(genome.id, genome.sequence.replace('?', 'N'))
        self.features = features
        self._featuresDict = features.featuresDict()

    def feature(self, name, window=500):
        """
        Get a feature (e.g., nsp2).

        @param name: A C{str} feature name.
        @param window: The C{int} number of nucleotides to include in the
            region of the genome that is searched for the reference
            sequence. This is in addition to the range of the feature in
            the reference. So e.g., if the feature is in positions
            1000-1300 in the reference and the window is 500, then the
            alignment call that tries to find the reference in the genome
            will look in the region 500-1800 of the genome. This is to
            allow for indels in the genome that cause its feature to be
            located at quite different offsets than the location of the
            feature in the reference.  It greatly speeds up the alignment
            and also prevents cases where a short prefix of the reference
            matches very early in the genome, followed by a huge gap and
            then the rest of the reference matching the genome. The default
            value below was arrived at by experimentation.
        @raise KeyError: If the feature is unknown.
        @return: An C{Alignment} instance, providing access to the found
            feature in the genome and its corresponding sequence in the
            reference.

        """
        return Alignment(self.genome, self.features.id,
                         self.features.getFeature(name), window=window)

    def checkVariant(self, variant, window=500):
        """
        Check that a set of changes in different features all happened as
        expected.

        @param variant: The C{str} name of a key in the VARIANTS C{dict}.
        @param window: The C{int} number of nucleotides to include in the
            region of the genome that is searched for the reference
            sequence. This is in addition to the range of the feature in
            the reference. So e.g., if the feature is in positions
            1000-1300 in the reference and the window is 500, then the
            alignment call that tries to find the reference in the genome
            will look in the region 500-1800 of the genome. This is to
            allow for indels in the genome that cause its feature to be
            located at quite different offsets than the location of the
            feature in the reference.  It greatly speeds up the alignment
            and also prevents cases where a short prefix of the reference
            matches very early in the genome, followed by a huge gap and
            then the rest of the reference matching the genome. The default
            value below was arrived at by experimentation.
        @return: A 3-C{tuple} with the number of checks done, the number of
            errors, and a C{dict} keyed by changes in C{changes}, with values
            a 2-C{tuple} of Booleans to indicate success or failure of the
            check for the reference and the genome respectively.
        """
        result = {}
        testCountTotal = errorCountTotal = 0
        changeDict = VARIANTS[variant]['changes']

        for featureName in changeDict:
            alignment = self.feature(featureName, window=window)
            result[featureName] = {}
            for what in 'aa', 'nt':
                try:
                    changes = changeDict[featureName][what]
                except KeyError:
                    pass
                else:
                    result[featureName][what] = {}
                    testCount, errorCount, result_ = alignment.check(
                        changes, what == 'nt')
                    testCountTotal += testCount
                    errorCountTotal += errorCount
                    for change, values in result_.items():
                        result[featureName][what][change] = values

        return testCountTotal, errorCountTotal, result
