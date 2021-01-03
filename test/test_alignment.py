from unittest import TestCase

from dark.reads import DNARead

from sars2seq.alignment import Alignment


class TestAlignment(TestCase):
    """
    Test the Alignment class.
    """
    def testNtSequences(self):
        """
        It must be possible to retrieve aligned nucleotide sequences.
        """
        feature = {
            'name': 'spike',
            'sequence': 'ATTC',
            'start': 0,
            'stop': 3,
        }
        alignment = Alignment(DNARead('genome', 'GGATTCGG'), 'referenceId',
                              feature)

        genomeNt, referenceNt = alignment.ntSequences()

        self.assertEqual('ATTC', genomeNt.sequence)
        self.assertEqual('genome (spike)', genomeNt.id)

        self.assertEqual('ATTC', referenceNt.sequence)
        self.assertEqual('referenceId (spike)', referenceNt.id)

    def testNtSequencesChangesString(self):
        """
        It must be possible to retrieve aligned nucleotide sequences
        and check on changes using a string specification.
        """
        feature = {
            'name': 'spike',
            'sequence': 'ATTC',
            'start': 0,
            'stop': 3,
        }
        alignment = Alignment(DNARead('genome', 'GGATTCGG'), 'referenceId',
                              feature)

        # Note: 1-based locations.
        testCount, errorCount, result = alignment.check(
            'A1A T2A A3T T4T', True)

        self.assertEqual(4, testCount)
        self.assertEqual(3, errorCount)
        self.assertEqual((True, True), result['A1A'])
        self.assertEqual((True, False), result['T2A'])
        self.assertEqual((False, True), result['A3T'])
        self.assertEqual((False, False), result['T4T'])

    def testNtSequencesChangesTuple(self):
        """
        It must be possible to retrieve aligned nucleotide sequences
        and check on changes using a tuple specification.
        """
        feature = {
            'name': 'spike',
            'sequence': 'ATTC',
            'start': 0,
            'stop': 3,
        }
        alignment = Alignment(DNARead('genome', 'GGATTCGG'), 'referenceId',
                              feature)

        # Note: 0-based offsets.
        testCount, errorCount, result = alignment.check(
            (('A', 0, 'A'), ('T', 1, 'A'), ('A', 2, 'T'), ('T', 3, 'T')), True)

        self.assertEqual(4, testCount)
        self.assertEqual(3, errorCount)
        self.assertEqual((True, True), result[('A', 0, 'A')])
        self.assertEqual((True, False), result[('T', 1, 'A')])
        self.assertEqual((False, True), result[('A', 2, 'T')])
        self.assertEqual((False, False), result[('T', 3, 'T')])

    def testNtSequencesGenomeSNP(self):
        """
        The genome must be able to have a SNP relative to the reference.
        """
        referenceSequence = 'TGGCGTGGA' + ('T' * 20) + 'CAAATCGG'
        genomeFeature = 'TGGCGTGGA' + ('T' * 9) + 'A' + ('T' * 10) + 'CAAATCGG'
        genomeSequence = 'CCCGG' + genomeFeature + 'CCCCCCC'

        feature = {
            'name': 'spike',
            'sequence': referenceSequence,
            'start': 5,
            'stop': 5 + len(referenceSequence),
        }
        alignment = Alignment(DNARead('genome', genomeSequence),
                              'referenceId', feature)

        genomeNt, referenceNt = alignment.ntSequences()

        expected = 'TGGCGTGGA' + ('T' * 9) + 'A' + ('T' * 10) + 'CAAATCGG'
        self.assertEqual(expected, genomeNt.sequence)
        self.assertEqual('genome (spike)', genomeNt.id)

        self.assertEqual(referenceSequence, referenceNt.sequence)
        self.assertEqual('referenceId (spike)', referenceNt.id)

        testCount, errorCount, result = alignment.check('T19A', True)

        self.assertEqual(1, testCount)
        self.assertEqual(0, errorCount)
        self.assertEqual((True, True), result['T19A'])

    def testNtSequencesGenomeGap(self):
        """
        The genome must be able to have a gap relative to the reference.
        """
        referenceSequence = 'TGGCGTGGA' + ('T' * 20) + 'CAAATCGG'
        genomeFeature = 'TGGCGTGGA' + ('T' * 19) + 'CAAATCGG'
        genomeSequence = 'CCCGG' + genomeFeature + 'CCCCCCC'

        feature = {
            'name': 'spike',
            'sequence': referenceSequence,
            'start': 5,
            'stop': 5 + len(referenceSequence),
        }
        alignment = Alignment(DNARead('genome', genomeSequence),
                              'referenceId', feature)

        genomeNt, referenceNt = alignment.ntSequences()

        expected = 'TGGCGTGGA-' + ('T' * 19) + 'CAAATCGG'
        self.assertEqual(expected, genomeNt.sequence)
        self.assertEqual('genome (spike)', genomeNt.id)

        self.assertEqual(referenceSequence, referenceNt.sequence)
        self.assertEqual('referenceId (spike)', referenceNt.id)

        testCount, errorCount, result = alignment.check('T10-', True)

        self.assertEqual(1, testCount)
        self.assertEqual(0, errorCount)
        self.assertEqual((True, True), result['T10-'])
