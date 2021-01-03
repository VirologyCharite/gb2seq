from unittest import TestCase

from sars2seq.features import Features

_FEATURES = Features()


class TestFeatures(TestCase):
    """
    Test the Features class.
    """
    def testGetFeatures(self):
        """
        The getFeature method must return a dict.
        """
        self.assertIsInstance(_FEATURES.getFeature('spike'), dict)

    def testUnknownFeature(self):
        """
        If an unknown feature is asked for, a KeyError must be raised.
        """
        self.assertRaisesRegex(KeyError, "^'xx'$", _FEATURES.getFeature, 'xx')

    def testAliasLookup(self):
        """
        A getFeature request that results in a alias lookup must return
        the expected dict.
        """
        featuresDict = _FEATURES.featuresDict()
        self.assertNotIn('membrane', featuresDict)
        self.assertEqual(_FEATURES.getFeature('membrane'),
                         _FEATURES.getFeature('membrane glycoprotein'))
