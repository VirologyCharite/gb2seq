from unittest import TestCase

from gb2seq.features import Features
from gb2seq.variants import VARIANTS


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
            for featureName in VARIANTS[variant]["changes"]:
                self.assertIsInstance(features[featureName], dict)

    def testKeys(self):
        """
        All variants must have 'changes' and 'description' keys.
        """
        for info in VARIANTS.values():
            self.assertIn("changes", info)
            self.assertIn("description", info)

    def testOnlyNtOrAa(self):
        """
        Only 'aa' or 'nt' are permitted as sub-keys.
        """
        for variant in VARIANTS:
            for featureName in VARIANTS[variant]["changes"]:
                specification = VARIANTS[variant]["changes"][featureName]
                self.assertEqual(set(), set(specification) - {"aa", "nt"})
