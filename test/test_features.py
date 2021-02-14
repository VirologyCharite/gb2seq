from unittest import TestCase

from dark.reads import DNARead

from sars2seq.features import Features

_FEATURES = Features()


class TestFeatures(TestCase):
    """
    Test the Features class.
    """
    def testGetFeatures(self):
        """
        The getitem method must return a dict.
        """
        self.assertIsInstance(_FEATURES['spike'], dict)

    def testUnknownFeature(self):
        """
        If an unknown feature is asked for, a KeyError must be raised.
        """
        self.assertRaisesRegex(KeyError, "^'xx'$", _FEATURES.__getitem__, 'xx')

    def testPassingDict(self):
        """
        It must be possible to initialize a Features instance via a dict.
        """
        value = {
            'name': 'spike',
            'sequence': 'ATTC',
            'start': 0,
            'stop': 4,
        }
        features = Features({'spike': value})

        self.assertIn('spike', features)
        self.assertEqual(value, features['spike'])

    def testPassingRefence(self):
        """
        It must be possible to pass a reference
        """
        reference = DNARead('refId', 'ATTC')
        features = Features({}, reference)
        self.assertIs(reference, features.reference)

    def testAliasLookup(self):
        """
        A getitem request that results in a alias lookup must return
        the expected dict.
        """
        self.assertNotIn('membrane', _FEATURES)
        self.assertEqual(_FEATURES['membrane'],
                         _FEATURES['membrane glycoprotein'])

    def testMembrane(self):
        """
        Test the return value when requesting the membrane protein.
        """
        self.assertEqual(
            _FEATURES['membrane'],
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
                    'VQ*'),
            })

    def testCanonicalName(self):
        """
        Converting abbreviated names into canonical names must work.
        """
        self.assertEqual(_FEATURES.canonicalName('m'), 'membrane glycoprotein')
        self.assertEqual(_FEATURES.canonicalName('e'), 'envelope protein')

    def testExpectedNames(self):
        """
        Test the full set of expected coronavirus feature names.
        """
        expected = set((
            "2'-O-ribose methyltransferase",
            "3'-to-5' exonuclease",
            "3'UTR",
            "3C-like proteinase",
            "5'UTR",
            "ORF10 protein",
            "ORF1a polyprotein",
            "ORF1ab polyprotein",
            "ORF3a protein",
            "ORF6 protein",
            "ORF7a protein",
            "ORF7b",
            "ORF8 protein",
            "RNA-dependent RNA polymerase",
            "endoRNAse",
            "envelope protein",
            "helicase",
            "leader protein",
            "membrane glycoprotein",
            "nsp10",
            "nsp11",
            "nsp2",
            "nsp3",
            "nsp4",
            "nsp6",
            "nsp7",
            "nsp8",
            "nsp9",
            "nucleocapsid phosphoprotein",
            "stem loop 1",
            "stem loop 2",
            "stem loop 3",
            "stem loop 4",
            "stem loop 5",
            "surface glycoprotein",
        ))
        self.assertEqual(expected, set(_FEATURES))
