from unittest import TestCase

from sars2seq.features import Features
from sars2seq.variants import spikeDeletion, VOC_20201201_UK, N501Y

VARIANTS = spikeDeletion, VOC_20201201_UK, N501Y


class TestVariants(TestCase):
    """
    Test the variants.
    """
    def testKnownFeatures(self):
        """
        Only known feature names are allowed.
        """
        features = Features()
        for variant in VARIANTS:
            for featureName in variant:
                self.assertIsInstance(features.getFeature(featureName), dict)

    def testOnlyNtOrAa(self):
        """
        Only 'aa' or 'nt' are permitted as sub-keys.
        """
        for variant in VARIANTS:
            for featureName, specification in variant.items():
                self.assertEqual(set(), set(specification) - {'aa', 'nt'})
