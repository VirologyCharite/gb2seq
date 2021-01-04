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

    def testMembrane(self):
        """
        Test the return value when requesting the membrane protein.
        """
        self.assertEqual(
            _FEATURES.getFeature('membrane'),
            {
                'name': 'membrane glycoprotein',
                'start': 26522,
                'stop': 27191,
                'note': 'ORF5; structural protein',
                'product': 'membrane glycoprotein',
                'sequence': (
                    'ATGGCAGATTCCAACGGTACTATTACCGTTGAAGAGCTTAAAAAGCTCCTTGAAC'
                    'AATGGAACCTAGTAATAGGTTTCCTATTCCTTACATGGATTTGTCTTCTACAATT'
                    'TGCCTATGCCAACAGGAATAGGTTTTTGTATATAATTAAGTTAATTTTCCTCTGG'
                    'CTGTTATGGCCAGTAACTTTAGCTTGTTTTGTGCTTGCTGCTGTTTACAGAATAA'
                    'ATTGGATCACCGGTGGAATTGCTATCGCAATGGCTTGTCTTGTAGGCTTGATGTG'
                    'GCTCAGCTACTTCATTGCTTCTTTCAGACTGTTTGCGCGTACGCGTTCCATGTGG'
                    'TCATTCAATCCAGAAACTAACATTCTTCTCAACGTGCCACTCCATGGCACTATTC'
                    'TGACCAGACCGCTTCTAGAAAGTGAACTCGTAATCGGAGCTGTGATCCTTCGTGG'
                    'ACATCTTCGTATTGCTGGACACCATCTAGGACGCTGTGACATCAAGGACCTGCCT'
                    'AAAGAAATCACTGTTGCTACATCACGAACGCTTTCTTATTACAAATTGGGAGCTT'
                    'CGCAGCGTGTAGCAGGTGACTCAGGTTTTGCTGCATACAGTCGCTACAGGATTGG'
                    'CAACTATAAATTAAACACAGACCATTCCAGTAGCAGTGACAATATTGCTTTGCTT'
                    'GTACAGTAA'),
                'translation': (
                    'MADSNGTITVEELKKLLEQWNLVIGFLFLTWICLLQFAYANRNRFLYIIKLIFLW'
                    'LLWPVTLACFVLAAVYRINWITGGIAIAMACLVGLMWLSYFIASFRLFARTRSMW'
                    'SFNPETNILLNVPLHGTILTRPLLESELVIGAVILRGHLRIAGHHLGRCDIKDLP'
                    'KEITVATSRTLSYYKLGASQRVAGDSGFAAYSRYRIGNYKLNTDHSSSSDNIALL'
                    'VQ'),
            })

    def testCanonicalName(self):
        """
        Abbreviated names must be found.
        """
        self.assertEqual(_FEATURES.canonicalName('m'), 'membrane glycoprotein')
        self.assertEqual(_FEATURES.canonicalName('e'), 'envelope protein')
