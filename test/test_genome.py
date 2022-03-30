"""
Tests for the SARS2Genome class and functions.

See test_genomes.py for tests of the specific genomes in ../data
"""

from unittest import TestCase
from io import StringIO

from dark.reads import DNARead

from sars2seq import DATA_DIR
from sars2seq.change import splitChange
from sars2seq.features import (
    Features, AmbiguousFeatureError, MissingFeatureError)
from sars2seq.genome import SARS2Genome, getGappedOffsets, alignmentEnd
from sars2seq.translate import NoSlipperySequenceError

from .fasta import getSequence


ALPHA_ID = 'EPI_ISL_601443 hCoV-19/England/MILK-9E05B3/2020'


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

    def testNtSequencesIndexErrorRaise(self):
        """
        If we check on nucleotide sequences with an out-of-range
        check, an IndexError must be raised.
        """
        features = Features(
            {
                'spike': {
                    'name': 'spike',
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

    def testNtSequencesIndexErrorPrint(self):
        """
        If we check on nucleotide sequences with an out-of-range
        check, an error must be printed if we pass onError='print'
        and the expected error result must be returned.
        """
        features = Features(
            {
                'spike': {
                    'name': 'spike',
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

    def testNtSequencesIndexErrorIgnore(self):
        """
        If we check on nucleotide sequences with an out-of-range
        check, no error should be printed if we pass onError='ignore'
        and the expected error result must be returned.
        """
        features = Features(
            {
                'spike': {
                    'name': 'spike',
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
                    'start': 5,
                    'stop': len(referenceSequence),
                },
            },
            DNARead('refId', referenceSequence))

        genome = SARS2Genome(DNARead('genId', genomeSequence), features)
        referenceNt, genomeNt = genome.ntSequences('spike')

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

    def testAaSequencesTranslationErrorRaise(self):
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

    def testAaSequencesTranslationErrorPrint(self):
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

    def testAaSequencesTranslationErrorIgnore(self):
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


class TestGetGappedOffsets(TestCase):
    """
    Test the getGappedOffsets function.
    """
    def testEmpty(self):
        """
        An empty string should get back an empty dictionary.
        """
        self.assertEqual({}, getGappedOffsets(''))

    def testOnlyGaps(self):
        """
        An string of gaps should get back an empty dictionary.
        """
        self.assertEqual({}, getGappedOffsets('---'))

    def testGapsBefore(self):
        """
        If there are gaps before the bases, the offsets must be correct.
        """
        self.assertEqual({0: 2, 1: 3}, getGappedOffsets('--CC'))


class TestOffsetInfo(TestCase):
    """
    Test the SARS2Genome.offsetInfo method.
    """
    def testRelativeOffsetWithNoFeatureName(self):
        """
        If relativeToFeature is passed as True, but a feature name is not
        given, a RuntimeError must be raised.
        """
        features = Features({}, DNARead('refId', 'AA'))
        genome = SARS2Genome(DNARead('genId', 'AA'), features=features)
        error = (
            r'^If relativeToFeature is True, a feature name must be given\.$')
        self.assertRaisesRegex(RuntimeError, error, genome.offsetInfo, 0,
                               relativeToFeature=True)

    def testNotARelativeOffsetButAaIsTrue(self):
        """
        If relativeToFeature is passed as False, a RuntimeError must be raised
        if aa is passed as True.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 0,
                    'stop': 6,
                },
            },
            DNARead('refId', 'AA'))
        genome = SARS2Genome(DNARead('genId', 'AA'), features=features)
        error = (r'^You cannot pass aa=True unless the offset you pass is '
                 r'relative to the feature\.$')
        self.assertRaisesRegex(RuntimeError, error, genome.offsetInfo, 0,
                               featureName='surface glycoprotein', aa=True)

    def testNoFeaturesAtOffset(self):
        """
        If there are no features at the given offset, a MissingFeatureError
        must be raised. This must be true regardless of the passed value of
        relativeToFeature.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 0,
                    'stop': 6,
                },
            },
            DNARead('refId', 'AA'))
        genome = SARS2Genome(DNARead('genId', 'AA'), features=features)
        error = (r"^Feature 'surface glycoprotein' \(located at sites 1-6\) "
                 r"does not overlap site 101\. There are no features at "
                 r"that site\.$")
        for relativeToFeature in False, True:
            self.assertRaisesRegex(MissingFeatureError, error,
                                   genome.offsetInfo, 100,
                                   featureName='surface glycoprotein',
                                   relativeToFeature=relativeToFeature)

    def testFeatureNotFoundAtOffsetButOneOtherFeatureIs(self):
        """
        If the requested feature is not found at the given offset, a
        MissingFeatureError must be raised and the error should include
        the one other feature found at that offset. This must be true
        regardless of the passed value of relativeToFeature.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 0,
                    'stop': 6,
                },
                'nsp10': {
                    'name': 'nsp10',
                    'start': 10,
                    'stop': 16,
                },
            },
            DNARead('refId', 'AA'))
        genome = SARS2Genome(DNARead('genId', 'AA'), features=features)
        error = (r"^Requested feature 'surface glycoprotein' \(located at "
                 r"sites 1-6\) does not overlap site 11. The feature\(s\) "
                 r"at that site are: 'nsp10'\.$")
        for relativeToFeature in False, True:
            self.assertRaisesRegex(MissingFeatureError, error,
                                   genome.offsetInfo, 10,
                                   featureName='surface glycoprotein',
                                   relativeToFeature=relativeToFeature)

    def testFeatureNotFoundAtOffsetButTwoOtherFeaturesAre(self):
        """
        If the requested feature is not found at the given offset, a
        MissingFeatureError must be raised and the error should include
        the one other feature found at that offset. This must be true
        regardless of the passed value of relativeToFeature.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 0,
                    'stop': 6,
                },
                'nsp10': {
                    'name': 'nsp10',
                    'start': 10,
                    'stop': 16,
                },
                'nsp11': {
                    'name': 'nsp11',
                    'start': 10,
                    'stop': 16,
                },
            },
            DNARead('refId', 'AA'))
        genome = SARS2Genome(DNARead('genId', 'AA'), features=features)
        error = (r"^Requested feature 'surface glycoprotein' \(located at "
                 r"sites 1-6\) does not overlap site 11. The feature\(s\) "
                 r"at that site are: 'nsp10', 'nsp11'\.$")
        for relativeToFeature in False, True:
            self.assertRaisesRegex(MissingFeatureError, error,
                                   genome.offsetInfo, 10,
                                   featureName='surface glycoprotein',
                                   relativeToFeature=relativeToFeature)

    def testMultipleFeaturesAtOffsetButNoFeatureRequested(self):
        """
        If multiple features are found at an offset but no feature name is
        passed, an AmbiguousFeatureError must be raised.
        """
        features = Features(
            {
                'nsp10': {
                    'name': 'nsp10',
                    'start': 10,
                    'stop': 16,
                },
                'nsp11': {
                    'name': 'nsp11',
                    'start': 10,
                    'stop': 16,
                },
            },
            DNARead('refId', 'AA'))
        genome = SARS2Genome(DNARead('genId', 'AA'), features=features)
        error = (r"^There are multiple features at site 13: 'nsp10', "
                 r"'nsp11'. Pass a feature name to specify which one you "
                 r"want\.$")

        self.assertRaisesRegex(AmbiguousFeatureError, error,
                               genome.offsetInfo, 12)

    def testNoFeaturesAtOffsetZero(self):
        """
        If there are no features at offset zero, the return information
        should indicate that, but a codon and amino acid are still returned.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 1,
                    'stop': 6,
                },
            },
            DNARead('refId', 'GTTCCCG'))

        genome = SARS2Genome(DNARead('genId', 'ATTCCCG'), features)

        self.assertEqual(
            {
                'alignmentOffset': 0,
                'featureName': None,
                'featureNames': set(),
                'reference': {
                    'aa': 'V',
                    'codon': 'GTT',
                    'frame': 0,
                    'id': 'refId',
                    'aaOffset': 0,
                    'ntOffset': 0,
                },
                'genome': {
                    'aa': 'I',
                    'codon': 'ATT',
                    'frame': 0,
                    'id': 'genId',
                    'aaOffset': 0,
                    'ntOffset': 0,
                }
            },
            genome.offsetInfo(0))

    def testOneFeatureAtOffsetZero(self):
        """
        If there is one feature at offset zero, the return information
        should indicate that.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 0,
                    'stop': 6,
                },
            },
            DNARead('refId', 'GTTCCC'))

        genome = SARS2Genome(DNARead('genId', 'ATTCCC'), features)

        self.assertEqual(
            {
                'alignmentOffset': 0,
                'featureName': 'surface glycoprotein',
                'featureNames': {'surface glycoprotein'},
                'reference': {
                    'aa': 'V',
                    'codon': 'GTT',
                    'frame': 0,
                    'id': 'refId',
                    'aaOffset': 0,
                    'ntOffset': 0,
                },
                'genome': {
                    'aa': 'I',
                    'codon': 'ATT',
                    'frame': 0,
                    'id': 'genId',
                    'aaOffset': 0,
                    'ntOffset': 0,
                }
            },
            genome.offsetInfo(0))

    def testOffsetAtLastNucleotideOfLastCodonOfFeature(self):
        """
        If the requested offset is the final nucleotide of the last codon of a
        feature, the return information should indicate that.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 0,
                    'stop': 6,
                },
            },
            DNARead('refId', 'GTTCCC'))

        genome = SARS2Genome(DNARead('genId', 'ATTCCC'), features)

        self.assertEqual(
            {
                'alignmentOffset': 5,
                'featureName': 'surface glycoprotein',
                'featureNames': {'surface glycoprotein'},
                'reference': {
                    'aa': 'P',
                    'codon': 'CCC',
                    'frame': 2,
                    'id': 'refId',
                    'aaOffset': 1,
                    'ntOffset': 5,
                },
                'genome': {
                    'aa': 'P',
                    'codon': 'CCC',
                    'frame': 2,
                    'id': 'genId',
                    'aaOffset': 1,
                    'ntOffset': 5,
                }
            },
            genome.offsetInfo(5))

    def testOffsetAtFirstNucleotideOfIncompleteFinalCodonOfFeature(self):
        """
        If the requested offset is the first nucleotide of an incomplete final
        codon of the genome, the return information should indicate that.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 0,
                    'stop': 6,
                },
            },
            DNARead('refId', 'GTTCCCA'))

        genome = SARS2Genome(DNARead('genId', 'ATTCCCA'), features)

        self.assertEqual(
            {
                'alignmentOffset': 6,
                'featureName': None,
                'featureNames': set(),
                'reference': {
                    'aa': 'X',
                    'codon': 'A',
                    'frame': 0,
                    'id': 'refId',
                    'aaOffset': 2,
                    'ntOffset': 6,
                },
                'genome': {
                    'aa': 'X',
                    'codon': 'A',
                    'frame': 0,
                    'id': 'genId',
                    'aaOffset': 2,
                    'ntOffset': 6,
                }
            },
            genome.offsetInfo(6))

    def testOffsetAtSecondNucleotideOfIncompleteFinalCodonOfFeature(self):
        """
        If the requested offset is the second nucleotide of an incomplete final
        codon of the genome, the return information should indicate that.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 0,
                    'stop': 6,
                },
            },
            DNARead('refId', 'GTTCCCAT'))

        genome = SARS2Genome(DNARead('genId', 'ATTCCCAT'), features)

        self.assertEqual(
            {
                'alignmentOffset': 7,
                'featureName': None,
                'featureNames': set(),
                'reference': {
                    'aa': 'X',
                    'codon': 'AT',
                    'frame': 1,
                    'id': 'refId',
                    'aaOffset': 2,
                    'ntOffset': 7,
                },
                'genome': {
                    'aa': 'X',
                    'codon': 'AT',
                    'frame': 1,
                    'id': 'genId',
                    'aaOffset': 2,
                    'ntOffset': 7,
                }
            },
            genome.offsetInfo(7))

    def testRelativeToFeature(self):
        """
        It must be possible to request information about an offset that is
        relative to the start of the feature.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 2,
                    'stop': 11,
                },
            },
            DNARead('refId', 'AATGTTCCCTTTAAA'))

        genome = SARS2Genome(DNARead('genId', 'AATGTACGCTTTAAA'), features)

        self.assertEqual(
            {
                'alignmentOffset': 5,
                'featureName': 'surface glycoprotein',
                'featureNames': {'surface glycoprotein'},
                'reference': {
                    'aa': 'S',
                    'codon': 'TCC',
                    'frame': 0,
                    'id': 'refId',
                    'aaOffset': 1,
                    'ntOffset': 5,
                },
                'genome': {
                    'aa': 'T',
                    'codon': 'ACG',
                    'frame': 0,
                    'id': 'genId',
                    'aaOffset': 1,
                    'ntOffset': 5,
                }
            },
            genome.offsetInfo(3, relativeToFeature=True,
                              featureName='surface glycoprotein'))

    def testRelativeToFeatureWithAaOffset(self):
        """
        It must be possible to request information about an offset that is
        relative to the start of the feature when the offset is given in terms
        of amino acids.
        """
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 2,
                    'stop': 11,
                },
            },
            DNARead('refId', 'AATGTTCCCTTTAAA'))

        genome = SARS2Genome(DNARead('genId', 'AATGTACGCTTTAAA'), features)

        self.assertEqual(
            {
                'alignmentOffset': 5,
                'featureName': 'surface glycoprotein',
                'featureNames': {'surface glycoprotein'},
                'reference': {
                    'aa': 'S',
                    'codon': 'TCC',
                    'frame': 0,
                    'id': 'refId',
                    'aaOffset': 1,
                    'ntOffset': 5,
                },
                'genome': {
                    'aa': 'T',
                    'codon': 'ACG',
                    'frame': 0,
                    'id': 'genId',
                    'aaOffset': 1,
                    'ntOffset': 5,
                }
            },
            genome.offsetInfo(1, aa=True, relativeToFeature=True,
                              featureName='surface glycoprotein'))

    def testInitialGapInAlignedGenomeOffsetZero(self):
        """
        If the alignment results in a gap character at the start of the
        genome and we request offset zero, we should see the gap in the
        codon and get back a '-' amino acid.
        """
        features = Features({}, DNARead('refId', 'GTTCCCAAATTGCTACTTTGATTGAG'))

        genome = SARS2Genome(
            DNARead('genId', 'TTCCCAAATTGCTACTTTGATTGAG'),
            features=features)

        print(genome.referenceAligned.sequence)
        print(genome.genomeAligned.sequence)

        self.assertEqual(
            {
                'alignmentOffset': 0,
                'featureName': None,
                'featureNames': set(),
                'reference': {
                    'aa': 'V',
                    'codon': 'GTT',
                    'frame': 0,
                    'id': 'refId',
                    'aaOffset': 0,
                    'ntOffset': 0,
                },
                'genome': {
                    'aa': '-',
                    'codon': '-TT',
                    'frame': 0,
                    'id': 'genId',
                    'aaOffset': 0,
                    'ntOffset': 0,
                }
            },
            genome.offsetInfo(0))

    def testInitialGapInAlignedGenomeOffsetOne(self):
        """
        If the alignment results in a gap character at the start of the
        genome and we request offset 1, we should get the first three
        nucleotides of the reference back as its codon, and the same
        from the genome. But the frames will differ because the genome
        does not have the first nucleotide that's in the reference.
        """
        seq = 'GTTCCCAAATTGCTACTTTGATTGAG'
        features = Features({}, DNARead('refId', seq))
        genome = SARS2Genome(DNARead('genId', seq[1:]), features=features)

        # Check the alignment is as expected.
        self.assertEqual(len(seq), len(genome.referenceAligned.sequence))
        self.assertEqual(0, genome.referenceAligned.sequence.count('-'))
        self.assertEqual(1, genome.genomeAligned.sequence.count('-'))
        self.assertEqual('-', genome.genomeAligned.sequence[0])

        self.assertEqual(
            {
                'alignmentOffset': 1,
                'featureName': None,
                'featureNames': set(),
                'reference': {
                    'aa': 'V',
                    'codon': 'GTT',
                    'frame': 1,
                    'id': 'refId',
                    'aaOffset': 0,
                    'ntOffset': 1,
                },
                'genome': {
                    'aa': 'F',
                    'codon': 'TTC',
                    'frame': 0,
                    'id': 'genId',
                    'aaOffset': 0,
                    'ntOffset': 0,
                }
            },
            genome.offsetInfo(1))

    def testInitialGapInAlignedGenomeWithFeatureAaOffset(self):
        """
        If the alignment results in a gap character at the start of the
        genome and we request an amino acid offset relative to a feature,
        we should get the expected identical result from the reference and
        the genome.
        """
        seq = 'GTTCCCAAATTGCTACTTTGATTGAG'
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 2,
                    'stop': 11,
                },
            },
            DNARead('refId', seq))

        genome = SARS2Genome(DNARead('genId', seq[1:]), features=features)

        # Check the alignment is as expected.
        self.assertEqual(len(seq), len(genome.referenceAligned.sequence))
        self.assertEqual(0, genome.referenceAligned.sequence.count('-'))
        self.assertEqual(1, genome.genomeAligned.sequence.count('-'))
        self.assertEqual('-', genome.genomeAligned.sequence[0])

        self.assertEqual(
            {
                'alignmentOffset': 8,
                'featureName': 'surface glycoprotein',
                'featureNames': {'surface glycoprotein'},
                'reference': {
                    'aa': 'I',
                    'codon': 'ATT',
                    'frame': 0,
                    'id': 'refId',
                    'aaOffset': 2,
                    'ntOffset': 8,
                },
                'genome': {
                    'aa': 'I',
                    'codon': 'ATT',
                    'frame': 0,
                    'id': 'genId',
                    'aaOffset': 2,
                    'ntOffset': 7,
                }
            },
            genome.offsetInfo(2, featureName='surface glycoprotein', aa=True,
                              relativeToFeature=True))

    def testInitialGapInAlignedGenomeWithFeatureNtOffsetAllFrames(self):
        """
        If the alignment results in a gap character at the start of the
        genome and we request a nucleotide offset all frames in a feature, we
        should get the expected identical result from the reference and the
        genome.
        """
        seq = 'GTTCCCAAATTGCTACTTTGATTGAG'
        features = Features(
            {
                'surface glycoprotein': {
                    'name': 'surface glycoprotein',
                    'start': 2,
                    'stop': 11,
                },
            },
            DNARead('refId', seq))

        genome = SARS2Genome(DNARead('genId', seq[1:]), features=features)

        # Check the alignment is as expected.
        self.assertEqual(len(seq), len(genome.referenceAligned.sequence))
        self.assertEqual(0, genome.referenceAligned.sequence.count('-'))
        self.assertEqual(1, genome.genomeAligned.sequence.count('-'))
        self.assertEqual('-', genome.genomeAligned.sequence[0])

        for frame in range(3):
            self.assertEqual(
                {
                    'alignmentOffset': 8 + frame,
                    'featureName': 'surface glycoprotein',
                    'featureNames': {'surface glycoprotein'},
                    'reference': {
                        'aa': 'I',
                        'codon': 'ATT',
                        'frame': frame,
                        'id': 'refId',
                        'aaOffset': 2,
                        'ntOffset': 8 + frame,
                    },
                    'genome': {
                        'aa': 'I',
                        'codon': 'ATT',
                        'frame': frame,
                        'id': 'genId',
                        'aaOffset': 2,
                        'ntOffset': 7 + frame,
                    }
                },
                genome.offsetInfo(6 + frame,
                                  featureName='surface glycoprotein',
                                  relativeToFeature=True))

    def testAlphaN501YByNucleotideOffsetRelativeToFeatureAllFrames(self):
        """
        We must be able to see the N501Y change in the spike of Alpha
        relative to the feature when using any frame and a nucleotide offset.
        """
        genomeRead = getSequence(DATA_DIR / 'EPI_ISL_601443.fasta')
        genome = SARS2Genome(genomeRead)
        start = genome.features['surface glycoprotein']['start']
        for frame in range(3):
            # The 500 comes from the N501Y.
            offset = 500 * 3 + frame

            # The offset in the genome should be 22990 + frame (this is based
            # on manually looking at the MAFFT alignment). We can calculate the
            # offset as the offset in the aligned genome (that corresponds to
            # the offset in the unaligned reference), minus the number of gaps
            # in the aligned genome up to that point. That produces the offset
            # in the original unaligned genome (assuming it didn't contain '-'
            # characters when we received it).
            expectedGenomeOffset = (
                start + offset - genome.genomeAligned.sequence[
                    :genome.gappedOffsets[start + offset]].count('-'))
            self.assertEqual(22990 + frame, expectedGenomeOffset)

            self.assertEqual(
                {
                    'alignmentOffset': start + offset,
                    'featureName': 'surface glycoprotein',
                    'featureNames': {'surface glycoprotein'},
                    'reference': {
                        'aa': 'N',
                        'codon': 'AAT',
                        'frame': frame,
                        'id': 'NC_045512.2',
                        'aaOffset': 500,
                        'ntOffset': start + offset,
                    },
                    'genome': {
                        'aa': 'Y',
                        'codon': 'TAT',
                        'frame': frame,
                        'id': ALPHA_ID,
                        # Alpha has 3 deletions before site 501 (at 69, 70,
                        # and 144).
                        'aaOffset': 497,
                        'ntOffset': expectedGenomeOffset,
                    }
                },
                genome.offsetInfo(offset, relativeToFeature=True,
                                  featureName='spike'))

    def testAlphaN501YWithNucleotideOffset(self):
        """
        We must be able to see the N501Y change in the spike of Alpha when
        the offset is given in nucleotides in any frame.
        """
        genomeRead = getSequence(DATA_DIR / 'EPI_ISL_601443.fasta')
        genome = SARS2Genome(genomeRead)
        start = genome.features['surface glycoprotein']['start']
        for frame in range(3):
            offset = start + 1500 + frame

            # See comment in test above re the 22990 test below.
            expectedGenomeOffset = (
                offset - genome.genomeAligned.sequence[
                    :genome.gappedOffsets[offset]].count('-'))
            self.assertEqual(22990 + frame, expectedGenomeOffset)
            self.assertEqual(
                {
                    'alignmentOffset': offset,
                    'featureName': 'surface glycoprotein',
                    'featureNames': {'surface glycoprotein'},
                    'reference': {
                        'aa': 'N',
                        'codon': 'AAT',
                        'frame': frame,
                        'id': 'NC_045512.2',
                        'aaOffset': 500,
                        'ntOffset': offset,
                    },
                    'genome': {
                        'aa': 'Y',
                        'codon': 'TAT',
                        'frame': frame,
                        'id': ALPHA_ID,
                        # Alpha has 3 deletions before site 501 (at 69, 70,
                        # and 144).
                        'aaOffset': 497,
                        'ntOffset': expectedGenomeOffset,
                    }
                },
                genome.offsetInfo(offset, featureName='spike'))

    def testAlphaN501YWithAaOffset(self):
        """
        We must be able to see the N501Y change in the spike of Alpha
        using an amino acid offset.
        """
        genomeRead = getSequence(DATA_DIR / 'EPI_ISL_601443.fasta')
        genome = SARS2Genome(genomeRead)
        start = genome.features['surface glycoprotein']['start']
        offset = 500

        # See comment in test above re the 22990 test below.
        expectedGenomeOffset = (
            start + offset * 3 - genome.genomeAligned.sequence[
                :genome.gappedOffsets[start + offset * 3]].count('-'))
        self.assertEqual(22990, expectedGenomeOffset)

        self.assertEqual(
            {
                'alignmentOffset': start + offset * 3,
                'featureName': 'surface glycoprotein',
                'featureNames': {'surface glycoprotein'},
                'reference': {
                    'aa': 'N',
                    'codon': 'AAT',
                    'frame': 0,
                    'id': 'NC_045512.2',
                    'aaOffset': offset,
                    'ntOffset': start + offset * 3,
                },
                'genome': {
                    'aa': 'Y',
                    'codon': 'TAT',
                    'frame': 0,
                    'id': ALPHA_ID,
                    # Alpha has 3 deletions before site 501 (at 69, 70,
                    # and 144).
                    'aaOffset': offset - 3,
                    'ntOffset': expectedGenomeOffset,
                }
            },
            genome.offsetInfo(offset, aa=True, relativeToFeature=True,
                              featureName='spike'))

    def testAlphaSpikeSubstitutions(self):
        """
        We must be able to see all the substitutions in the spike of Alpha.
        """
        genomeRead = getSequence(DATA_DIR / 'EPI_ISL_601443.fasta')
        genome = SARS2Genome(genomeRead)

        for change in 'N501Y A570D D614G P681H T716I S982A D1118H'.split():
            referenceAa, offset, genomeAa = splitChange(change)
            offsetInfo = genome.offsetInfo(
                offset, aa=True, relativeToFeature=True, featureName='spike')
            self.assertEqual(offsetInfo['reference']['aa'], referenceAa)
            self.assertEqual(offsetInfo['genome']['aa'], genomeAa)

    def testAlphaSpikeDeletions(self):
        """
        We must be able to see all the deletions in the spike of Alpha.
        """
        genomeRead = getSequence(DATA_DIR / 'EPI_ISL_601443.fasta')
        genome = SARS2Genome(genomeRead)

        for change in 'H69- V70- Y144-'.split():
            referenceAa, offset, genomeAa = splitChange(change)
            offsetInfo = genome.offsetInfo(
                offset, aa=True, relativeToFeature=True, featureName='spike')
            self.assertEqual(offsetInfo['reference']['aa'], referenceAa)
            self.assertEqual(offsetInfo['genome']['aa'], genomeAa)
