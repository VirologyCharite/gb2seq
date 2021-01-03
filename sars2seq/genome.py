from sars2seq.alignment import Alignment


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
        return Alignment(self._genome, self.features.id,
                         self.features.getFeature(name))

    def checkVariant(self, variant):
        """
        Check that a set of changes in different features all happened as
        expected.

        @param variant: A C{dict}, as in variants.py
        @return: A 3-C{tuple} with the number of checks done, the number of
            errors, and a C{dict} keyed by changes in C{changes}, with values
            a 2-C{tuple} of Booleans to indicate success or failure of the
            check for the reference and the genome respectively.
        """
        result = {}
        testCountTotal = errorCountTotal = 0

        for featureName in variant:
            alignment = self.feature(featureName)
            result[featureName] = {}
            for what in 'aa', 'nt':
                try:
                    changes = variant[featureName][what]
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
