"""
Tests for the SARS2Genome class and functions.

See test_genomes.py for tests of the specific genomes in ../data
"""

from unittest import TestCase
from io import StringIO

from dark.reads import DNARead

from sars2seq.features import Features
from sars2seq.genome import SARS2Genome, getNonGapOffsets, alignmentEnd
from sars2seq.translate import NoSlipperySequenceError


class TestSARS2Genome(TestCase):
    """
    Test the SARS2Genome class.
    """
    def testNtSequences(self):
        """
        It must be possible to retrieve aligned nucleotide sequences.
        """
        features = Features(
            {
                'spike': {
                    'name': 'spike',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        referenceNt, genomeNt = genome.ntSequences('spike')

        self.assertEqual('ATTC', genomeNt.sequence)
        self.assertEqual('genId (spike)', genomeNt.id)

        self.assertEqual('ATTC', referenceNt.sequence)
        self.assertEqual('refId (spike)', referenceNt.id)

    def testNtSequencesChangesString(self):
        """
        It must be possible to retrieve aligned nucleotide sequences
        and check on changes using a string specification.
        """
        features = Features(
            {
                'spike': {
                    'name': 'spike',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        # Note: 1-based locations.
        testCount, errorCount, result = genome.checkFeature(
            'spike', 'A1A T2A A3T T4T', True)

        self.assertEqual(4, testCount)
        self.assertEqual(3, errorCount)
        self.assertEqual((True, 'A', True, 'A'), result['A1A'])
        self.assertEqual((True, 'T', False, 'T'), result['T2A'])
        self.assertEqual((False, 'T', True, 'T'), result['A3T'])
        self.assertEqual((False, 'C', False, 'C'), result['T4T'])

    def testNtSequencesChangesIndexErrorRaise(self):
        """
        If we check on nucleotide sequences with an out-of-range
        check, an IndexError must be raised.
        """
        features = Features(
            {
                'spike': {
                    'name': 'spike',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        error = (r"^Index 99999 out of range trying to access feature "
                 r"'spike' of length 4 sequence 'refId \(spike\)' via "
                 r"expected change specification 'A100000A'\.$")
        self.assertRaisesRegex(IndexError, error, genome.checkFeature,
                               'spike', 'A100000A', True)

    def testNtSequencesChangesIndexErrorPrint(self):
        """
        If we check on nucleotide sequences with an out-of-range
        check, an error must be printed if we pass onError='print'
        and the expected error result must be returned.
        """
        features = Features(
            {
                'spike': {
                    'name': 'spike',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        err = StringIO()

        # Two lines of error output are printed.
        error = (
            r"Index 99999 out of range trying to access feature "
            r"'spike' of length 4 sequence 'refId (spike)' via "
            r"expected change specification 'A100000A'."
            "\n"
            r"Index 99999 out of range trying to access feature "
            r"'spike' of length 4 sequence 'genId (spike)' via "
            r"expected change specification 'A100000A'."
            "\n"
        )
        testCount, errorCount, result = genome.checkFeature(
            'spike', 'A100000A', nt=True, onError='print', errFp=err)
        self.assertEqual(error, err.getvalue())

        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual((False, None, False, None), result['A100000A'])

    def testNtSequencesChangesIndexErrorIgnore(self):
        """
        If we check on nucleotide sequences with an out-of-range
        check, no error should be printed if we pass onError='ignore'
        and the expected error result must be returned.
        """
        features = Features(
            {
                'spike': {
                    'name': 'spike',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        err = StringIO()
        testCount, errorCount, result = genome.checkFeature(
            'spike', 'A100000A', nt=True, onError='ignore', errFp=err)
        self.assertEqual('', err.getvalue())

        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual((False, None, False, None), result['A100000A'])

    def testNtSequencesChangesTuple(self):
        """
        It must be possible to retrieve aligned nucleotide sequences
        and check on changes using a tuple specification.
        """
        features = Features(
            {
                'spike': {
                    'name': 'spike',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        # Note: 0-based offsets.
        testCount, errorCount, result = genome.checkFeature(
            'spike',
            (('A', 0, 'A'), ('T', 1, 'A'), ('A', 2, 'T'), ('T', 3, 'T')), True)

        self.assertEqual(4, testCount)
        self.assertEqual(3, errorCount)
        self.assertEqual((True, 'A', True, 'A'), result[('A', 0, 'A')])
        self.assertEqual((True, 'T', False, 'T'), result[('T', 1, 'A')])
        self.assertEqual((False, 'T', True, 'T'), result[('A', 2, 'T')])
        self.assertEqual((False, 'C', False, 'C'), result[('T', 3, 'T')])

    def testNtSequencesGenomeSNP(self):
        """
        The genome must be able to have a SNP relative to the reference.
        """
        referenceSequence = 'TGGCGTGGA' + ('T' * 20) + 'CAAATCGG'
        genomeFeature = 'TGGCGTGGA' + ('T' * 9) + 'A' + ('T' * 10) + 'CAAATCGG'
        genomeSequence = 'CCCGG' + genomeFeature + 'CCCCCCC'

        features = Features(
            {
                'spike': {
                    'name': 'spike',
                    'sequence': referenceSequence,
                    'start': 0,
                    'stop': len(referenceSequence),
                },
            },
            DNARead('refId', referenceSequence))

        genome = SARS2Genome(DNARead('genId', genomeSequence), features)

        referenceNt, genomeNt = genome.ntSequences('spike')

        expected = 'TGGCGTGGA' + ('T' * 9) + 'A' + ('T' * 10) + 'CAAATCGG'
        self.assertEqual(expected, genomeNt.sequence)
        self.assertEqual('genId (spike)', genomeNt.id)

        self.assertEqual(referenceSequence, referenceNt.sequence)
        self.assertEqual('refId (spike)', referenceNt.id)

        testCount, errorCount, result = genome.checkFeature(
            'spike', 'T19A', True)

        self.assertEqual(1, testCount)
        self.assertEqual(0, errorCount)
        self.assertEqual((True, 'T', True, 'A'), result['T19A'])

    def testNtSequencesGenomeGap(self):
        """
        The genome must be able to have a gap relative to the reference.
        """
        referenceSequence = 'TGGCGTGGA' + ('T' * 20) + 'CAAATCGG'
        genomeFeature = 'TGGA' + ('T' * 19) + 'CAAATCGG'
        genomeSequence = 'CCCGGTGGCG' + genomeFeature + 'CCCCCCC'

        features = Features(
            {
                'spike': {
                    'name': 'spike',
                    'sequence': referenceSequence,
                    'start': 5,
                    'stop': len(referenceSequence),
                },
            },
            DNARead('refId', referenceSequence))

        genome = SARS2Genome(DNARead('genId', genomeSequence), features)

        # The genome offset is initialized to None and isn't set until
        # after ntSequences is called.
        # self.assertEqual(None, alignment.genomeOffset)

        referenceNt, genomeNt = genome.ntSequences('spike')

        # self.assertEqual(5, alignment.genomeOffset)

        self.assertEqual(referenceSequence[5:], referenceNt.sequence)
        self.assertEqual('refId (spike)', referenceNt.id)

        expected = 'TGGA-' + ('T' * 19) + 'CAAATCGG'
        self.assertEqual(expected, genomeNt.sequence)
        self.assertEqual('genId (spike)', genomeNt.id)

        testCount, errorCount, result = genome.checkFeature(
            'spike', 'T5-', True)

        self.assertEqual(1, testCount)
        self.assertEqual(0, errorCount)
        self.assertEqual((True, 'T', True, '-'), result['T5-'])

    def testAaSequencesChangesTranslationErrorRaise(self):
        """
        Check that a TranslationError is raised when checking AA
        sequences.
        """
        features = Features(
            {
                'orf1ab': {
                    'name': 'ORF1ab polyprotein',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        error = r"^No slippery sequence found\.$"
        self.assertRaisesRegex(
            NoSlipperySequenceError, error, genome.checkFeature,
            'orf1ab', 'A100000A', False)

    def testAaSequencesChangesTranslationErrorPrint(self):
        """
        Check that a TranslationError is printed when checking AA
        sequences and onError='print' and that the expected result
        is returned.
        """
        features = Features(
            {
                'orf1ab': {
                    'name': 'ORF1ab polyprotein',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        err = StringIO()
        error = 'No slippery sequence found.\n'

        testCount, errorCount, result = genome.checkFeature(
            'orf1ab', 'A100000A', nt=False, onError='print', errFp=err)
        self.assertEqual(error, err.getvalue())

        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual((False, None, False, None), result['A100000A'])

    def testAaSequencesChangesTranslationErrorIgnore(self):
        """
        Check that no error is printed when checking AA sequences and
        onError='ignore' and that the expected result is returned.
        """
        features = Features(
            {
                'orf1ab': {
                    'name': 'ORF1ab polyprotein',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        err = StringIO()

        testCount, errorCount, result = genome.checkFeature(
            'orf1ab', 'A100000A', nt=False, onError='ignore', errFp=err)
        self.assertEqual('', err.getvalue())

        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual((False, None, False, None), result['A100000A'])

    def testAaSequencesTranslationNoSlipperySequenceRaise(self):
        """
        The aaSequences function must raise if it can't translate an
        'ORF1ab polyprotein' sequence due to a missing slippery sequence.
        """
        features = Features(
            {
                'ORF1ab polyprotein': {
                    'name': 'ORF1ab polyprotein',
                    'sequence': 'ATTC',
                    'start': 0,
                    'stop': 4,
                },
            },
            DNARead('refId', 'ATTC'))

        genome = SARS2Genome(DNARead('genId', 'GGATTCGG'), features)

        error = r'^No slippery sequence found\.$'
        self.assertRaisesRegex(NoSlipperySequenceError, error,
                               genome.aaSequences, 'ORF1ab polyprotein')


class TestAlignmentEnd(TestCase):
    """
    Test the alignmentEnd function.
    """
    def testOffsetIndexError(self):
        """
        If the start index is greater than the length of the passed string, an
        IndexError should be raised.
        """
        error = r'^string index out of range$'
        self.assertRaisesRegex(IndexError, error, alignmentEnd, '', 4, 1)

    def testLengthTooLarge(self):
        """
        If the initial offset plus the length is greater than the length of the
        (non-gaps) in the passed string, an IndexError should be raised.
        """
        error = r'^string index out of range$'
        self.assertRaisesRegex(IndexError, error, alignmentEnd, 'ACCG', 2, 5)

    def testSequenceTooShort(self):
        """
        If the passed sequence is long enough (with respect to the passed
        offset and length) but doesn't have enough non-gap characters, an
        IndexError should be raised.
        """
        error = r'^string index out of range$'
        self.assertRaisesRegex(
            IndexError, error, alignmentEnd, 'CC--------T-', 2, 5)

    def testEmptyString(self):
        """
        Looking for a zero length section of a zero length string starting
        at offset zero should get a result of zero.
        """
        self.assertEqual(0, alignmentEnd('', 0, 0))

    def testNonZeroNoGaps(self):
        """
        Passing a non-zero start offset and a string with no gaps should work.
        """
        self.assertEqual(3, alignmentEnd('ACCTA', 1, 2))

    def testZeroWithOneGap(self):
        """
        Passing a zero start offset and a string with one gap should work.
        """
        self.assertEqual(3, alignmentEnd('A-CCTA', 0, 2))

    def testZeroWithTwoGaps(self):
        """
        Passing a zero start offset and a string with two gaps should work.
        """
        self.assertEqual(4, alignmentEnd('A--CCTA', 0, 2))

    def testZeroWithTwoGapsNonContiguous(self):
        """
        Passing a zero start offset and a string with two gaps that are not
        contiguous should work.
        """
        self.assertEqual(5, alignmentEnd('A-C-CTA', 0, 3))

    def testNonZeroWithTwoGapsNonContiguous(self):
        """
        Passing a non-zero start offset and a string with two gaps that are not
        contiguous should work.
        """
        self.assertEqual(7, alignmentEnd('TTA-C-CTA', 2, 3))


class TestGetNonGapOffsets(TestCase):
    """
    Test the getNonGapOffsets function.
    """
    def testEmpty(self):
        """
        An empty string should get back an empty dictionary.
        """
        self.assertEqual({}, getNonGapOffsets(''))

    def testOnlyGaps(self):
        """
        An string of gaps should get back an empty dictionary.
        """
        self.assertEqual({}, getNonGapOffsets('---'))

    def testGapsBefore(self):
        """
        If there are gaps before the bases, the offsets must be correct.
        """
        self.assertEqual({0: 2, 1: 3}, getNonGapOffsets('--CC'))
