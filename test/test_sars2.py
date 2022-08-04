from unittest import TestCase

from gb2seq.sars2 import SARS_COV_2_ALIASES


class TestSars2Aliases(TestCase):
    """
    Test the SARS_COV_2_ALIASES dict.
    """

    def testAliasKeysLowerCase(self):
        """
        Alphanumeric alias keys must be lower case in order to be found.
        """
        self.assertTrue(
            all(key.islower() for key in SARS_COV_2_ALIASES if key.isalpha())
        )
