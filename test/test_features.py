from unittest import TestCase

from dark.reads import DNARead

from sars2seq.features import Features, ALIASES

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
        Test the returned dictionary when requesting the membrane protein.
        """
        expected = {
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
        }

        for name in 'membrane glycoprotein', 'membrane', 'm', 'orf5':
            self.assertEqual(expected, _FEATURES[name])
            self.assertEqual(expected, _FEATURES[name.upper()])

    def testOffsetInUnknownFeature(self):
        """
        If a genome offset for an unknown feature is requested, a KeyError
        must be raised.
        """
        self.assertRaisesRegex(KeyError, "^'xx'$",
                               _FEATURES.referenceOffset, 'xx', 10)

    def testAaOffsetInMembrane(self):
        """
        Test we get the correct genome offset given an AA offset in the
        membrane protein.
        """
        offset = 5
        # 26522 is from the membrane test above.
        self.assertEqual(
            26522 + 3 * offset,
            _FEATURES.referenceOffset('membrane', offset, aa=True))

    def testNtOffsetInMembraneDefaultIsNt(self):
        """
        Test we get the correct genome offset given an offset in the
        membrane protein when we don't specify whether the feature offset is
        aa or nt (i.e., check that nt is the default).
        """
        offset = 5
        # 26522 is from the membrane test above.
        self.assertEqual(
            26522 + offset,
            _FEATURES.referenceOffset('membrane', offset))

    def testNtOffsetInMembrane(self):
        """
        Test we get the correct genome offset given a nucleotide offset in the
        membrane protein when we explicitly pass aa=False.
        """
        offset = 5
        # 26522 is from the membrane test above.
        self.assertEqual(
            26522 + offset,
            _FEATURES.referenceOffset('membrane', offset, aa=False))

    def testFeaturesAtMembraneOffset(self):
        """
        Test we get the membrane protein back if we ask what features are at
        an offset it contains.
        """
        # 26522 is from the membrane test above.
        self.assertEqual(
            {
                'membrane glycoprotein'
            },
            _FEATURES.featuresAt(26522))

    def testFeaturesAtTooHighOffset(self):
        """
        Test we get nothing back if we ask what features are at an offset that
        is much bigger than the genome length.
        """
        self.assertEqual(set(), _FEATURES.featuresAt(1E9))

    def testFeaturesAtNegativeOffset(self):
        """
        Test we get nothing back if we ask what features are at a negative
        offset.
        """
        self.assertEqual(set(), _FEATURES.featuresAt(-1))

    def testFeaturesAtZeroOffset(self):
        """
        Test we get nothing back if we ask what features are at offset zero
        and do not ask for untranslated features.
        """
        self.assertEqual(set(), _FEATURES.featuresAt(0))

    def testFeaturesAtZeroOffsetIncludeUntranslated(self):
        """
        Test we get the 5'UTR back if we ask what features are at offset zero
        if we ask for untranslated features to also be returned.
        """
        self.assertEqual(
            {
                "5'UTR",
            },
            _FEATURES.featuresAt(0, includeUntranslated=True))

    def testFeaturesAtOrf1abOffset(self):
        """
        Test we get the expected result if we ask what features are at the
        first offset of Orf1ab.
        """
        self.assertEqual(
            {
                'leader protein',
                'ORF1ab polyprotein',
                'ORF1a polyprotein',
            },
            _FEATURES.featuresAt(265))

    def testFeaturesAtNsp2Offset(self):
        """
        Test we get the expected result if we ask what features are at an
        offset of NSP2, but we don't ask for features without a translation.
        """
        self.assertEqual(
            {
                'nsp2',
                'ORF1ab polyprotein',
                'ORF1a polyprotein',
            },
            _FEATURES.featuresAt(2700))

    def testFeaturesAtRdRPOffsetWithStemLoops(self):
        """
        Test we get the expected result if we ask what features are at an
        offset of the RdRP that is also in two stem loops, but we don't request
        that untranslated features are included.
        """
        self.assertEqual(
            {
                'ORF1ab polyprotein',
                'RNA-dependent RNA polymerase',
            },
            _FEATURES.featuresAt(13500))

    def testFeaturesAtRdRPOffsetWithStemLoopsIncludeUntranslated(self):
        """
        Test we get the expected result if we ask what features, including
        untranslated ones, are at an offset of the RdRP that is also in two
        stem loops.
        """
        self.assertEqual(
            {
                'ORF1ab polyprotein',
                'RNA-dependent RNA polymerase',
                'stem loop 1',
                'stem loop 2',
            },
            _FEATURES.featuresAt(13500, includeUntranslated=True))

    def testFeaturesAtRdRPOffset(self):
        """
        Test we get the expected result if we ask what features are at an
        offset of the RdRP but we don't ask for untranslated features.
        """
        self.assertEqual(
            {
                'ORF1ab polyprotein',
                'RNA-dependent RNA polymerase',
            },
            _FEATURES.featuresAt(13550))

    def testCanonicalName(self):
        """
        Converting abbreviated names into canonical names must work.
        """
        self.assertEqual(_FEATURES.canonicalName('s'), 'surface glycoprotein')
        self.assertEqual(_FEATURES.canonicalName('e'), 'envelope protein')
        self.assertEqual(_FEATURES.canonicalName('m'), 'membrane glycoprotein')
        self.assertEqual(_FEATURES.canonicalName('n'),
                         'nucleocapsid phosphoprotein')

    def testExpectedNames(self):
        """
        Test the full set of expected coronavirus feature names.
        """
        expected = {
            "2'-O-ribose methyltransferase",
            "3'-to-5' exonuclease",
            "3'UTR",
            "3C-like proteinase",
            "5'UTR",
            "endoRNAse",
            "envelope protein",
            "helicase",
            "leader protein",
            "membrane glycoprotein",
            "nsp2",
            "nsp3",
            "nsp4",
            "nsp6",
            "nsp7",
            "nsp8",
            "nsp9",
            "nsp10",
            "nsp11",
            "nucleocapsid phosphoprotein",
            "ORF10 protein",
            "ORF1a polyprotein",
            "ORF1ab polyprotein",
            "ORF3a protein",
            "ORF6 protein",
            "ORF7a protein",
            "ORF7b",
            "ORF8 protein",
            "RNA-dependent RNA polymerase",
            "stem loop 1",
            "stem loop 2",
            "stem loop 3",
            "stem loop 4",
            "stem loop 5",
            "surface glycoprotein",
        }
        self.assertEqual(expected, set(_FEATURES))

    def testAliasKeysLowerCase(self):
        """
        Alphanumeric alias keys must be lower case in order to be found.
        """
        self.assertTrue(all(key.islower() for key in ALIASES if key.isalpha()))
