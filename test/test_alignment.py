"""
Tests for the Gb2Alignment class and functions.

See test_genomes.py for tests of the specific genomes in ../data
"""

from unittest import TestCase
from io import StringIO

from dark.reads import DNARead

from gb2seq import DATA_DIR
from gb2seq.alignment import (
    Gb2Alignment,
    getGappedOffsets,
    alignmentEnd,
    offsetInfoMultipleGenomes,
    AlignmentError,
)
from gb2seq.change import splitChange
from gb2seq.features import Features, AmbiguousFeatureError, MissingFeatureError
from gb2seq.translate import NoSlipperySequenceError

from .fasta import getSequence


ALPHA_ID = "EPI_ISL_601443 hCoV-19/England/MILK-9E05B3/2020"


class TestGb2Alignment(TestCase):
    """
    Test the Gb2Alignment class.
    """

    def testUnequalPreAlignedSequences(self):
        """
        If pre-aligned sequences are passed that are not of equal length, a
        AlignmentError must be raised.
        """
        error = (
            r"^The length of the given pre-aligned reference \(5\) does "
            r"not match the length of the given pre-aligned genome "
            r"\(4\)\.$"
        )
        self.assertRaisesRegex(
            AlignmentError,
            error,
            Gb2Alignment,
            DNARead("ref", "ATGC"),
            features={},
            referenceAligned=DNARead("ref-aln", "AT-GC"),
            genomeAligned=DNARead("gen-aln", "AT-G"),
        )

    def testPreAlignedSequencesBothStartWithAGap(self):
        """
        If pre-aligned sequences are passed and both start with a gap, an
        AlignmentError must be raised.
        """
        error = (
            r"^The reference and genome alignment sequences both start "
            r'with a "-" character\.$'
        )
        self.assertRaisesRegex(
            AlignmentError,
            error,
            Gb2Alignment,
            DNARead("ref", "ATGC"),
            features={},
            referenceAligned=DNARead("ref-aln", "-ATGC"),
            genomeAligned=DNARead("gen-aln", "-ATGC"),
        )

    def testPreAlignedSequencesBothEndWithAGap(self):
        """
        If pre-aligned sequences are passed and both end with a gap, an
        AlignmentError must be raised.
        """
        error = (
            r"^The reference and genome alignment sequences both end "
            r'with a "-" character\.$'
        )
        self.assertRaisesRegex(
            AlignmentError,
            error,
            Gb2Alignment,
            DNARead("ref", "ATGC"),
            features={},
            referenceAligned=DNARead("ref-aln", "ATGC-"),
            genomeAligned=DNARead("gen-aln", "ATGC-"),
        )

    def testNtSequences(self):
        """
        It must be possible to retrieve aligned nucleotide sequences.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        referenceNt, genomeNt = alignment.ntSequences("spike")

        self.assertEqual("ATTC", genomeNt.sequence)
        self.assertEqual("genId (spike)", genomeNt.id)

        self.assertEqual("ATTC", referenceNt.sequence)
        self.assertEqual("refId (spike)", referenceNt.id)

    def testNtSequencesChangesString(self):
        """
        It must be possible to retrieve aligned nucleotide sequences
        and check on changes using a string specification.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        # Note: 1-based locations.
        testCount, errorCount, result = alignment.checkFeature(
            "spike", "A1A T2A A3T T4T", False
        )

        self.assertEqual(4, testCount)
        self.assertEqual(3, errorCount)
        self.assertEqual((True, "A", True, "A"), result["A1A"])
        self.assertEqual((True, "T", False, "T"), result["T2A"])
        self.assertEqual((False, "T", True, "T"), result["A3T"])
        self.assertEqual((False, "C", False, "C"), result["T4T"])

    def testNtSequencesIndexErrorRaise(self):
        """
        If we check on nucleotide sequences with an out-of-range
        check, an IndexError must be raised.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        error = (
            r"^Index 99999 out of range trying to access feature "
            r"'spike' of length 4 sequence 'refId \(spike\)' via "
            r"expected change specification 'A100000A'\.$"
        )
        self.assertRaisesRegex(
            IndexError, error, alignment.checkFeature, "spike", "A100000A", False
        )

    def testNtSequencesIndexErrorPrint(self):
        """
        If we check on nucleotide sequences with an out-of-range
        check, an error must be printed if we pass onError='print'
        and the expected error result must be returned.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

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
        testCount, errorCount, result = alignment.checkFeature(
            "spike", "A100000A", aa=False, onError="print", errFp=err
        )
        self.assertEqual(error, err.getvalue())

        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual((False, None, False, None), result["A100000A"])

    def testNtSequencesIndexErrorIgnore(self):
        """
        If we check on nucleotide sequences with an out-of-range
        check, no error should be printed if we pass onError='ignore'
        and the expected error result must be returned.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        err = StringIO()
        testCount, errorCount, result = alignment.checkFeature(
            "spike", "A100000A", aa=False, onError="ignore", errFp=err
        )
        self.assertEqual("", err.getvalue())

        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual((False, None, False, None), result["A100000A"])

    def testNtSequencesChangesTuple(self):
        """
        It must be possible to retrieve aligned nucleotide sequences
        and check on changes using a tuple specification.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        # Note: 0-based offsets.
        testCount, errorCount, result = alignment.checkFeature(
            "spike",
            (("A", 0, "A"), ("T", 1, "A"), ("A", 2, "T"), ("T", 3, "T")),
            aa=False,
        )

        self.assertEqual(4, testCount)
        self.assertEqual(3, errorCount)
        self.assertEqual((True, "A", True, "A"), result[("A", 0, "A")])
        self.assertEqual((True, "T", False, "T"), result[("T", 1, "A")])
        self.assertEqual((False, "T", True, "T"), result[("A", 2, "T")])
        self.assertEqual((False, "C", False, "C"), result[("T", 3, "T")])

    def testNtSequencesGenomeSNP(self):
        """
        The genome must be able to have a SNP relative to the reference.
        """
        referenceSequence = "TGGCGTGGA" + ("T" * 20) + "CAAATCGG"
        genomeFeature = "TGGCGTGGA" + ("T" * 9) + "A" + ("T" * 10) + "CAAATCGG"
        genomeSequence = "CCCGG" + genomeFeature + "CCCCCCC"

        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": len(referenceSequence),
                },
            },
            DNARead("refId", referenceSequence),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", genomeSequence), features)

        referenceNt, genomeNt = alignment.ntSequences("spike")

        expected = "TGGCGTGGA" + ("T" * 9) + "A" + ("T" * 10) + "CAAATCGG"
        self.assertEqual(expected, genomeNt.sequence)
        self.assertEqual("genId (spike)", genomeNt.id)

        self.assertEqual(referenceSequence, referenceNt.sequence)
        self.assertEqual("refId (spike)", referenceNt.id)

        testCount, errorCount, result = alignment.checkFeature("spike", "T19A", False)

        self.assertEqual(1, testCount)
        self.assertEqual(0, errorCount)
        self.assertEqual((True, "T", True, "A"), result["T19A"])

    def testNtSequencesGenomeGap(self):
        """
        The genome must be able to have a gap relative to the reference.
        """
        referenceSequence = "TGGCGTGGA" + ("T" * 20) + "CAAATCGG"
        genomeFeature = "TGGA" + ("T" * 19) + "CAAATCGG"
        genomeSequence = "CCCGGTGGCG" + genomeFeature + "CCCCCCC"

        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 5,
                    "stop": len(referenceSequence),
                },
            },
            DNARead("refId", referenceSequence),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", genomeSequence), features)
        referenceNt, genomeNt = alignment.ntSequences("spike")

        self.assertEqual(referenceSequence[5:], referenceNt.sequence)
        self.assertEqual("refId (spike)", referenceNt.id)

        expected = "TGGA-" + ("T" * 19) + "CAAATCGG"
        self.assertEqual(expected, genomeNt.sequence)
        self.assertEqual("genId (spike)", genomeNt.id)

        testCount, errorCount, result = alignment.checkFeature("spike", "T5-", False)

        self.assertEqual(1, testCount)
        self.assertEqual(0, errorCount)
        self.assertEqual((True, "T", True, "-"), result["T5-"])

    def testAaSequencesTranslationErrorRaise(self):
        """
        Check that a TranslationError is raised when checking AA sequences.
        """
        features = Features(
            {
                "orf1ab": {
                    "name": "ORF1ab polyprotein",
                    "sequence": "ATTC",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        error = r"^No slippery sequence found\.$"
        self.assertRaisesRegex(
            NoSlipperySequenceError,
            error,
            alignment.checkFeature,
            "orf1ab",
            "A100000A",
            True,
        )

    def testAaSequencesTranslationErrorPrint(self):
        """
        Check that a TranslationError is printed when checking AA
        sequences and onError='print' and that the expected result
        is returned.
        """
        features = Features(
            {
                "orf1ab": {
                    "name": "ORF1ab polyprotein",
                    "sequence": "ATTC",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        err = StringIO()
        error = "No slippery sequence found.\n"

        testCount, errorCount, result = alignment.checkFeature(
            "orf1ab", "A100000A", aa=True, onError="print", errFp=err
        )
        self.assertEqual(error, err.getvalue())

        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual((False, None, False, None), result["A100000A"])

    def testAaSequencesTranslationErrorIgnore(self):
        """
        Check that no error is printed when checking AA sequences and
        onError='ignore' and that the expected result is returned.
        """
        features = Features(
            {
                "orf1ab": {
                    "name": "ORF1ab polyprotein",
                    "sequence": "ATTC",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        err = StringIO()

        testCount, errorCount, result = alignment.checkFeature(
            "orf1ab", "A100000A", aa=True, onError="ignore", errFp=err
        )
        self.assertEqual("", err.getvalue())

        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual((False, None, False, None), result["A100000A"])

    def testAaSequencesTranslationNoSlipperySequenceRaise(self):
        """
        The aaSequences function must raise if it can't translate an
        'ORF1ab polyprotein' sequence due to a missing slippery sequence.
        """
        features = Features(
            {
                "ORF1ab polyprotein": {
                    "name": "ORF1ab polyprotein",
                    "sequence": "ATTC",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        error = r"^No slippery sequence found\.$"
        self.assertRaisesRegex(
            NoSlipperySequenceError, error, alignment.aaSequences, "ORF1ab polyprotein"
        )


class TestAlignmentEnd(TestCase):
    """
    Test the alignmentEnd function.
    """

    def testOffsetIndexError(self):
        """
        If the start index is greater than the length of the passed string, an
        IndexError should be raised.
        """
        error = r"^string index out of range$"
        self.assertRaisesRegex(IndexError, error, alignmentEnd, "", 4, 1)

    def testLengthTooLarge(self):
        """
        If the initial offset plus the length is greater than the length of the
        (non-gaps) in the passed string, an IndexError should be raised.
        """
        error = r"^string index out of range$"
        self.assertRaisesRegex(IndexError, error, alignmentEnd, "ACCG", 2, 5)

    def testSequenceTooShort(self):
        """
        If the passed sequence is long enough (with respect to the passed
        offset and length) but doesn't have enough non-gap characters, an
        IndexError should be raised.
        """
        error = r"^string index out of range$"
        self.assertRaisesRegex(IndexError, error, alignmentEnd, "CC--------T-", 2, 5)

    def testEmptyString(self):
        """
        Looking for a zero length section of a zero length string starting
        at offset zero should get a result of zero.
        """
        self.assertEqual(0, alignmentEnd("", 0, 0))

    def testNonZeroNoGaps(self):
        """
        Passing a non-zero start offset and a string with no gaps should work.
        """
        self.assertEqual(3, alignmentEnd("ACCTA", 1, 2))

    def testZeroWithOneGap(self):
        """
        Passing a zero start offset and a string with one gap should work.
        """
        self.assertEqual(3, alignmentEnd("A-CCTA", 0, 2))

    def testZeroWithTwoGaps(self):
        """
        Passing a zero start offset and a string with two gaps should work.
        """
        self.assertEqual(4, alignmentEnd("A--CCTA", 0, 2))

    def testZeroWithTwoGapsNonContiguous(self):
        """
        Passing a zero start offset and a string with two gaps that are not
        contiguous should work.
        """
        self.assertEqual(5, alignmentEnd("A-C-CTA", 0, 3))

    def testNonZeroWithTwoGapsNonContiguous(self):
        """
        Passing a non-zero start offset and a string with two gaps that are not
        contiguous should work.
        """
        self.assertEqual(7, alignmentEnd("TTA-C-CTA", 2, 3))


class TestGetGappedOffsets(TestCase):
    """
    Test the getGappedOffsets function.
    """

    def testEmpty(self):
        """
        An empty string should get back an empty dictionary.
        """
        self.assertEqual({}, getGappedOffsets(""))

    def testOnlyGaps(self):
        """
        An string of gaps should get back an empty dictionary.
        """
        self.assertEqual({}, getGappedOffsets("---"))

    def testGapsBefore(self):
        """
        If there are gaps before the bases, the offsets must be correct.
        """
        self.assertEqual({0: 2, 1: 3, 2: 4}, getGappedOffsets("--CC"))


class TestOffsetInfo(TestCase):
    """
    Test the Gb2Alignment.offsetInfo method.
    """

    def testRelativeOffsetWithNoFeatureName(self):
        """
        If relativeToFeature is passed as True, but a feature name is not
        given, a ValueError must be raised.
        """
        features = Features({}, DNARead("refId", "AA"), sars2=True)
        alignment = Gb2Alignment(DNARead("genId", "AA"), features=features)
        error = r"^If relativeToFeature is True, a feature name must be given\.$"
        self.assertRaisesRegex(
            ValueError, error, alignment.offsetInfo, 0, relativeToFeature=True
        )

    def testNotARelativeOffsetButAaIsTrue(self):
        """
        If relativeToFeature is passed as False, a ValueError must be raised
        if aa is passed as True.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 0,
                    "stop": 6,
                },
            },
            DNARead("refId", "AA"),
            sars2=True,
        )
        alignment = Gb2Alignment(DNARead("genId", "AA"), features=features)
        error = (
            r"^You cannot pass aa=True unless the offset you pass is "
            r"relative to the feature\.$"
        )
        self.assertRaisesRegex(
            ValueError,
            error,
            alignment.offsetInfo,
            0,
            featureName="surface glycoprotein",
            aa=True,
        )

    def testOffsetTooBig(self):
        """
        If the requested offset is beyond the end of the reference genome, an
        IndexError must be raised.
        """
        features = Features(
            {
                "nsp10": {
                    "name": "nsp10",
                    "start": 10,
                    "stop": 16,
                },
                "nsp11": {
                    "name": "nsp11",
                    "start": 10,
                    "stop": 16,
                },
            },
            DNARead("refId", "AA"),
            sars2=True,
        )
        alignment = Gb2Alignment(DNARead("genId", "AA"), features=features)
        error = r"^Request for offset 12 in a reference genome of length 2\.$"

        self.assertRaisesRegex(
            IndexError,
            error,
            alignment.offsetInfo,
            12,
        )

    def testNoFeaturesAtOffset(self):
        """
        If there are no features at the given offset, a MissingFeatureError
        must be raised. This must be true regardless of the passed value of
        relativeToFeature.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 0,
                    "stop": 6,
                },
            },
            DNARead("refId", "A" * 200),
            sars2=True,
        )
        alignment = Gb2Alignment(DNARead("genId", "AA"), features=features)
        error = (
            r"^Feature 'surface glycoprotein' \(located at sites 1-6\) "
            r"does not overlap site 101\. There are no features at "
            r"that site\.$"
        )
        for relativeToFeature in False, True:
            self.assertRaisesRegex(
                MissingFeatureError,
                error,
                alignment.offsetInfo,
                100,
                featureName="surface glycoprotein",
                relativeToFeature=relativeToFeature,
            )

    def testFeatureNotFoundAtOffsetButOneOtherFeatureIs(self):
        """
        If the requested feature is not found at the given offset, a
        MissingFeatureError must be raised and the error should include
        the one other feature found at that offset. This must be true
        regardless of the passed value of relativeToFeature.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 0,
                    "stop": 6,
                },
                "nsp10": {
                    "name": "nsp10",
                    "start": 10,
                    "stop": 16,
                },
            },
            DNARead("refId", "A" * 25),
            sars2=True,
        )
        alignment = Gb2Alignment(DNARead("genId", "AA"), features=features)
        error = (
            r"^Requested feature 'surface glycoprotein' \(located at "
            r"sites 1-6\) does not overlap site 11. The feature\(s\) "
            r"at that site are: 'nsp10' \(11 - 16\)\.$"
        )
        for relativeToFeature in False, True:
            self.assertRaisesRegex(
                MissingFeatureError,
                error,
                alignment.offsetInfo,
                10,
                featureName="surface glycoprotein",
                relativeToFeature=relativeToFeature,
            )

    def testFeatureNotFoundAtOffsetButTwoOtherFeaturesAre(self):
        """
        If the requested feature is not found at the given offset, a
        MissingFeatureError must be raised and the error should include
        the one other feature found at that offset. This must be true
        regardless of the passed value of relativeToFeature.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 0,
                    "stop": 6,
                },
                "nsp10": {
                    "name": "nsp10",
                    "start": 10,
                    "stop": 16,
                },
                "nsp11": {
                    "name": "nsp11",
                    "start": 10,
                    "stop": 16,
                },
            },
            DNARead("refId", "A" * 25),
            sars2=True,
        )
        alignment = Gb2Alignment(DNARead("genId", "AA"), features=features)
        error = (
            r"^Requested feature 'surface glycoprotein' \(located at "
            r"sites 1-6\) does not overlap site 11. The feature\(s\) "
            r"at that site are: 'nsp10' \(11 - 16\), 'nsp11' \(11 - 16\)\.$"
        )
        for relativeToFeature in False, True:
            self.assertRaisesRegex(
                MissingFeatureError,
                error,
                alignment.offsetInfo,
                10,
                featureName="surface glycoprotein",
                relativeToFeature=relativeToFeature,
            )

    def testMultipleFeaturesAtOffsetButNoFeatureRequestedNoAmbiguous(self):
        """
        If multiple features are found at an offset but no feature name is
        passed and allowAmbiguous is False, an AmbiguousFeatureError must be
        raised.
        """
        features = Features(
            {
                "nsp10": {
                    "name": "nsp10",
                    "start": 10,
                    "stop": 16,
                },
                "nsp11": {
                    "name": "nsp11",
                    "start": 10,
                    "stop": 16,
                },
            },
            DNARead("refId", "AAAAAAAAAAAAAA"),
            sars2=True,
        )
        alignment = Gb2Alignment(DNARead("genId", "AA"), features=features)
        error = (
            r"^There are multiple features at site 13: 'nsp10' \(11 - 16\), 'nsp11' "
            r"\(11 - 16\). Pass a feature name to specify which one you want\.$"
        )

        self.assertRaisesRegex(
            AmbiguousFeatureError, error, alignment.offsetInfo, 12, allowAmbiguous=False
        )

    def testMultipleFeaturesAtOffsetButNoFeatureRequestedAllowAmbiguous(self):
        """
        If multiple features are found at an offset but no feature name is
        passed and allowAmbiguous is true,...
        """
        features = Features(
            {
                "nsp10": {
                    "name": "nsp10",
                    "start": 10,
                    "stop": 16,
                },
                "nsp11": {
                    "name": "nsp11",
                    "start": 10,
                    "stop": 16,
                },
            },
            DNARead("refId", "AAAAAAAAAAAAAA"),
            sars2=True,
        )
        alignment = Gb2Alignment(DNARead("genId", "AA"), features=features)

        self.assertEqual(
            {
                "alignmentOffset": 11,
                "featureName": "nsp10",
                "featureNames": {"nsp10", "nsp11"},
                "reference": {
                    "aa": "K",
                    "codon": "AAA",
                    "frame": 1,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 11,
                },
                "genome": {
                    "aa": "-",
                    "codon": "AA-",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
            },
            alignment.offsetInfo(11, allowAmbiguous=True),
        )

    def testNoFeaturesAtOffsetZero(self):
        """
        If there are no features at offset zero, the return information
        should indicate that, but a codon and amino acid are still returned.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 1,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCCG"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "ATTCCCG"), features)

        self.assertEqual(
            {
                "alignmentOffset": 0,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "V",
                    "codon": "GTT",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
                "genome": {
                    "aa": "I",
                    "codon": "ATT",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
            },
            alignment.offsetInfo(0),
        )

    def testEdlibNoAmbiguous(self):
        """
        If the edlib aligner is used, but without ambiguous nucleotide codes,
        it should align completely different sequences.
        """
        features = Features({}, DNARead("refId", "CGTTCCCG"), sars2=True)

        alignment = Gb2Alignment(
            DNARead("genId", "RWWMMMR"), features, aligner="edlib", matchAmbiguous=False
        )

        self.assertEqual(
            {
                "alignmentOffset": 0,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "R",
                    "codon": "CGT",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
                "genome": {
                    "aa": "X",
                    "codon": "RWW",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
            },
            alignment.offsetInfo(0),
        )

        self.assertEqual("-", alignment.genomeAligned.sequence[-1])

    def testEdlibAmbiguous(self):
        """
        If the edlib aligner is used with ambiguous nucleotide codes,
        the alignment should match properly (in this case the first nt
        of the aligned genome will be set to '-' because the rest of the
        genome matches perfectly (they are all ambiguous codes that mach the
        reference).
        """
        features = Features({}, DNARead("refId", "CGTTCCCG"), sars2=True)

        alignment = Gb2Alignment(
            DNARead("genId", "RWWMMMR"), features, aligner="edlib", matchAmbiguous=True
        )

        self.assertEqual(
            {
                "alignmentOffset": 0,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "R",
                    "codon": "CGT",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
                "genome": {
                    "aa": "-",
                    "codon": "-RW",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
            },
            alignment.offsetInfo(0),
        )

        self.assertEqual("-", alignment.genomeAligned.sequence[0])

    def testLowCoverageGenome(self):
        """
        If the genome has insufficient coverage of the reference, None
        must be returned.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 1,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCCG"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GTT"), features)

        self.assertIsNone(alignment.offsetInfo(0, minReferenceCoverage=0.9))

    def testOneFeatureAtOffsetZero(self):
        """
        If there is one feature at offset zero, the return information
        should indicate that.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 0,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "ATTCCC"), features)

        self.assertEqual(
            {
                "alignmentOffset": 0,
                "featureName": "surface glycoprotein",
                "featureNames": {"surface glycoprotein"},
                "reference": {
                    "aa": "V",
                    "codon": "GTT",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
                "genome": {
                    "aa": "I",
                    "codon": "ATT",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
            },
            alignment.offsetInfo(0),
        )

    def testOffsetAtLastNucleotideOfLastCodonOfFeature(self):
        """
        If the requested offset is the final nucleotide of the last codon of a
        feature, the return information should indicate that.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 0,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "ATTCCC"), features)

        self.assertEqual(
            {
                "alignmentOffset": 5,
                "featureName": "surface glycoprotein",
                "featureNames": {"surface glycoprotein"},
                "reference": {
                    "aa": "P",
                    "codon": "CCC",
                    "frame": 2,
                    "id": "refId",
                    "aaOffset": 1,
                    "ntOffset": 5,
                },
                "genome": {
                    "aa": "P",
                    "codon": "CCC",
                    "frame": 2,
                    "id": "genId",
                    "aaOffset": 1,
                    "ntOffset": 5,
                },
            },
            alignment.offsetInfo(5),
        )

    def testOffsetAtFirstNucleotideOfIncompleteFinalCodonOfFeature(self):
        """
        If the requested offset is the first nucleotide of an incomplete final
        codon of the genome, the return information should indicate that.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 0,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCCA"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "ATTCCCA"), features)

        self.assertEqual(
            {
                "alignmentOffset": 6,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "X",
                    "codon": "A",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 2,
                    "ntOffset": 6,
                },
                "genome": {
                    "aa": "X",
                    "codon": "A",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 2,
                    "ntOffset": 6,
                },
            },
            alignment.offsetInfo(6),
        )

    def testOffsetAtSecondNucleotideOfIncompleteFinalCodonOfFeature(self):
        """
        If the requested offset is the second nucleotide of an incomplete final
        codon of the genome, the return information should indicate that.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 0,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCCAT"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "ATTCCCAT"), features)

        self.assertEqual(
            {
                "alignmentOffset": 7,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "X",
                    "codon": "AT",
                    "frame": 1,
                    "id": "refId",
                    "aaOffset": 2,
                    "ntOffset": 7,
                },
                "genome": {
                    "aa": "X",
                    "codon": "AT",
                    "frame": 1,
                    "id": "genId",
                    "aaOffset": 2,
                    "ntOffset": 7,
                },
            },
            alignment.offsetInfo(7),
        )

    def testRelativeToFeature(self):
        """
        It must be possible to request information about an offset that is
        relative to the start of the feature.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 2,
                    "stop": 11,
                },
            },
            DNARead("refId", "AATGTTCCCTTTAAA"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "AATGTACGCTTTAAA"), features)

        self.assertEqual(
            {
                "alignmentOffset": 5,
                "featureName": "surface glycoprotein",
                "featureNames": {"surface glycoprotein"},
                "reference": {
                    "aa": "S",
                    "codon": "TCC",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 1,
                    "ntOffset": 5,
                },
                "genome": {
                    "aa": "T",
                    "codon": "ACG",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 1,
                    "ntOffset": 5,
                },
            },
            alignment.offsetInfo(
                3, relativeToFeature=True, featureName="surface glycoprotein"
            ),
        )

    def testRelativeToFeatureWithAaOffset(self):
        """
        It must be possible to request information about an offset that is
        relative to the start of the feature when the offset is given in terms
        of amino acids.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 2,
                    "stop": 11,
                },
            },
            DNARead("refId", "AATGTTCCCTTTAAA"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "AATGTACGCTTTAAA"), features)

        self.assertEqual(
            {
                "alignmentOffset": 5,
                "featureName": "surface glycoprotein",
                "featureNames": {"surface glycoprotein"},
                "reference": {
                    "aa": "S",
                    "codon": "TCC",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 1,
                    "ntOffset": 5,
                },
                "genome": {
                    "aa": "T",
                    "codon": "ACG",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 1,
                    "ntOffset": 5,
                },
            },
            alignment.offsetInfo(
                1, aa=True, relativeToFeature=True, featureName="surface glycoprotein"
            ),
        )

    def testInitialGapInAlignedGenomeOffsetZero(self):
        """
        If the alignment results in a gap character at the start of the
        genome and we request offset zero, we should see the gap in the
        codon and get back a '-' amino acid.
        """
        features = Features(
            {}, DNARead("refId", "GTTCCCAAATTGCTACTTTGATTGAG"), sars2=True
        )

        alignment = Gb2Alignment(
            DNARead("genId", "TTCCCAAATTGCTACTTTGATTGAG"), features=features
        )

        self.assertEqual(
            {
                "alignmentOffset": 0,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "V",
                    "codon": "GTT",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
                "genome": {
                    "aa": "-",
                    "codon": "-TT",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
            },
            alignment.offsetInfo(0),
        )

    def testInitialGapInAlignedGenomeOffsetOne(self):
        """
        If the alignment results in a gap character at the start of the
        genome and we request offset 1, we should get the first three
        nucleotides of the reference back as its codon, and the same
        from the genome. But the frames will differ because the genome
        does not have the first nucleotide that's in the reference.
        """
        seq = "GTTCCCAAATTGCTACTTTGATTGAG"
        features = Features({}, DNARead("refId", seq), sars2=True)
        alignment = Gb2Alignment(DNARead("genId", seq[1:]), features)

        # Check the alignment is as expected.
        self.assertEqual(len(seq), len(alignment.referenceAligned.sequence))
        self.assertEqual(0, alignment.referenceAligned.sequence.count("-"))
        self.assertEqual(1, alignment.genomeAligned.sequence.count("-"))
        self.assertEqual("-", alignment.genomeAligned.sequence[0])

        self.assertEqual(
            {
                "alignmentOffset": 1,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "V",
                    "codon": "GTT",
                    "frame": 1,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 1,
                },
                "genome": {
                    "aa": "F",
                    "codon": "TTC",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
            },
            alignment.offsetInfo(1),
        )

    def testInitialGapInAlignedGenomeWithFeatureAaOffset(self):
        """
        If the alignment results in a gap character at the start of the
        genome and we request an amino acid offset relative to a feature,
        we should get the expected identical result from the reference and
        the genome.
        """
        seq = "GTTCCCAAATTGCTACTTTGATTGAG"
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 2,
                    "stop": 11,
                },
            },
            DNARead("refId", seq),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", seq[1:]), features)

        # Check the alignment is as expected.
        self.assertEqual(len(seq), len(alignment.referenceAligned.sequence))
        self.assertEqual(0, alignment.referenceAligned.sequence.count("-"))
        self.assertEqual(1, alignment.genomeAligned.sequence.count("-"))
        self.assertEqual("-", alignment.genomeAligned.sequence[0])

        self.assertEqual(
            {
                "alignmentOffset": 8,
                "featureName": "surface glycoprotein",
                "featureNames": {"surface glycoprotein"},
                "reference": {
                    "aa": "I",
                    "codon": "ATT",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 2,
                    "ntOffset": 8,
                },
                "genome": {
                    "aa": "I",
                    "codon": "ATT",
                    "frame": 0,
                    "id": "genId",
                    "aaOffset": 2,
                    "ntOffset": 7,
                },
            },
            alignment.offsetInfo(
                2, featureName="surface glycoprotein", aa=True, relativeToFeature=True
            ),
        )

    def testInitialGapInAlignedGenomeWithFeatureNtOffsetAllFrames(self):
        """
        If the alignment results in a gap character at the start of the
        genome and we request a nucleotide offset all frames in a feature, we
        should get the expected identical result from the reference and the
        genome.
        """
        seq = "GTTCCCAAATTGCTACTTTGATTGAG"
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 2,
                    "stop": 11,
                },
            },
            DNARead("refId", seq),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", seq[1:]), features)

        # Check the alignment is as expected.
        self.assertEqual(len(seq), len(alignment.referenceAligned.sequence))
        self.assertEqual(0, alignment.referenceAligned.sequence.count("-"))
        self.assertEqual(1, alignment.genomeAligned.sequence.count("-"))
        self.assertEqual("-", alignment.genomeAligned.sequence[0])

        for frame in range(3):
            self.assertEqual(
                {
                    "alignmentOffset": 8 + frame,
                    "featureName": "surface glycoprotein",
                    "featureNames": {"surface glycoprotein"},
                    "reference": {
                        "aa": "I",
                        "codon": "ATT",
                        "frame": frame,
                        "id": "refId",
                        "aaOffset": 2,
                        "ntOffset": 8 + frame,
                    },
                    "genome": {
                        "aa": "I",
                        "codon": "ATT",
                        "frame": frame,
                        "id": "genId",
                        "aaOffset": 2,
                        "ntOffset": 7 + frame,
                    },
                },
                alignment.offsetInfo(
                    6 + frame,
                    featureName="surface glycoprotein",
                    relativeToFeature=True,
                ),
            )

    def testAlphaN501YByNucleotideOffsetRelativeToFeatureAllFrames(self):
        """
        We must be able to see the N501Y change in the spike of Alpha
        relative to the feature when using any frame and a nucleotide offset.
        """
        genomeRead = getSequence(DATA_DIR / "EPI_ISL_601443.fasta")
        alignment = Gb2Alignment(genomeRead, features=Features(sars2=True))
        start = alignment.features["surface glycoprotein"]["start"]
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
                start
                + offset
                - alignment.genomeAligned.sequence[
                    : alignment.gappedOffsets[start + offset]
                ].count("-")
            )
            self.assertEqual(22990 + frame, expectedGenomeOffset)

            self.assertEqual(
                {
                    "alignmentOffset": start + offset,
                    "featureName": "surface glycoprotein",
                    "featureNames": {"surface glycoprotein"},
                    "reference": {
                        "aa": "N",
                        "codon": "AAT",
                        "frame": frame,
                        "id": "NC_045512.2",
                        "aaOffset": 500,
                        "ntOffset": start + offset,
                    },
                    "genome": {
                        "aa": "Y",
                        "codon": "TAT",
                        "frame": frame,
                        "id": ALPHA_ID,
                        # Alpha has 3 deletions before site 501 (at 69, 70,
                        # and 144).
                        "aaOffset": 497,
                        "ntOffset": expectedGenomeOffset,
                    },
                },
                alignment.offsetInfo(
                    offset, relativeToFeature=True, featureName="spike"
                ),
            )

    def testAlphaN501YWithNucleotideOffset(self):
        """
        We must be able to see the N501Y change in the spike of Alpha when
        the offset is given in nucleotides in any frame.
        """
        genomeRead = getSequence(DATA_DIR / "EPI_ISL_601443.fasta")
        alignment = Gb2Alignment(genomeRead, features=Features(sars2=True))
        start = alignment.features["surface glycoprotein"]["start"]
        for frame in range(3):
            offset = start + 1500 + frame

            # See comment in test above re the 22990 test below.
            expectedGenomeOffset = offset - alignment.genomeAligned.sequence[
                : alignment.gappedOffsets[offset]
            ].count("-")
            self.assertEqual(22990 + frame, expectedGenomeOffset)
            self.assertEqual(
                {
                    "alignmentOffset": offset,
                    "featureName": "surface glycoprotein",
                    "featureNames": {"surface glycoprotein"},
                    "reference": {
                        "aa": "N",
                        "codon": "AAT",
                        "frame": frame,
                        "id": "NC_045512.2",
                        "aaOffset": 500,
                        "ntOffset": offset,
                    },
                    "genome": {
                        "aa": "Y",
                        "codon": "TAT",
                        "frame": frame,
                        "id": ALPHA_ID,
                        # Alpha has 3 deletions before site 501 (at 69, 70,
                        # and 144).
                        "aaOffset": 497,
                        "ntOffset": expectedGenomeOffset,
                    },
                },
                alignment.offsetInfo(offset, featureName="spike"),
            )

    def testAlphaN501YWithAaOffset(self):
        """
        We must be able to see the N501Y change in the spike of Alpha
        using an amino acid offset.
        """
        genomeRead = getSequence(DATA_DIR / "EPI_ISL_601443.fasta")
        alignment = Gb2Alignment(genomeRead, features=Features(sars2=True))
        start = alignment.features["surface glycoprotein"]["start"]
        offset = 500

        # See comment in test above re the 22990 test below.
        expectedGenomeOffset = (
            start
            + offset * 3
            - alignment.genomeAligned.sequence[
                : alignment.gappedOffsets[start + offset * 3]
            ].count("-")
        )
        self.assertEqual(22990, expectedGenomeOffset)

        self.assertEqual(
            {
                "alignmentOffset": start + offset * 3,
                "featureName": "surface glycoprotein",
                "featureNames": {"surface glycoprotein"},
                "reference": {
                    "aa": "N",
                    "codon": "AAT",
                    "frame": 0,
                    "id": "NC_045512.2",
                    "aaOffset": offset,
                    "ntOffset": start + offset * 3,
                },
                "genome": {
                    "aa": "Y",
                    "codon": "TAT",
                    "frame": 0,
                    "id": ALPHA_ID,
                    # Alpha has 3 deletions before site 501 (at 69, 70,
                    # and 144).
                    "aaOffset": offset - 3,
                    "ntOffset": expectedGenomeOffset,
                },
            },
            alignment.offsetInfo(
                offset, aa=True, relativeToFeature=True, featureName="spike"
            ),
        )

    def testAlphaSpikeSubstitutions(self):
        """
        We must be able to see all the substitutions in the spike of Alpha.
        """
        genomeRead = getSequence(DATA_DIR / "EPI_ISL_601443.fasta")
        alignment = Gb2Alignment(genomeRead, features=Features(sars2=True))

        for change in "N501Y A570D D614G P681H T716I S982A D1118H".split():
            referenceAa, offset, genomeAa = splitChange(change)
            offsetInfo = alignment.offsetInfo(
                offset, aa=True, relativeToFeature=True, featureName="spike"
            )
            self.assertEqual(offsetInfo["reference"]["aa"], referenceAa)
            self.assertEqual(offsetInfo["genome"]["aa"], genomeAa)

    def testAlphaSpikeDeletions(self):
        """
        We must be able to see all the deletions in the spike of Alpha.
        """
        genomeRead = getSequence(DATA_DIR / "EPI_ISL_601443.fasta")
        alignment = Gb2Alignment(genomeRead, features=Features(sars2=True))

        for change in "H69- V70- Y144-".split():
            referenceAa, offset, genomeAa = splitChange(change)
            offsetInfo = alignment.offsetInfo(
                offset, aa=True, relativeToFeature=True, featureName="spike"
            )
            self.assertEqual(offsetInfo["reference"]["aa"], referenceAa)
            self.assertEqual(offsetInfo["genome"]["aa"], genomeAa)


class TestCoverage(TestCase):
    """
    Test the coverage function.
    """

    def testFullyCoveredGenome(self):
        """
        Get the coverage of a genome that is fully covered.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "ATTC"), features)

        self.assertEqual((4, 4), alignment.coverage())

    def testFullyCoveredFeature(self):
        """
        Get the coverage of a feature that is fully covered.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "GGATTCGG"), features)

        self.assertEqual((4, 4), alignment.coverage("spike"))

    def testHalfCoveredGenome(self):
        """
        Get the coverage of a genome that is half covered.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "AT"), features)

        self.assertEqual((2, 4), alignment.coverage())

    def testQuarterCoveredFeature(self):
        """
        Get the coverage of a feature that is one-quarter covered.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "A"), features)

        self.assertEqual((1, 4), alignment.coverage("spike"))

    def testNInGenome(self):
        """
        Get the coverage of a genome with an N in it.
        """
        features = Features(
            {
                "spike": {
                    "name": "spike",
                    "start": 0,
                    "stop": 4,
                },
            },
            DNARead("refId", "ATTC"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "ATTN"), features)

        self.assertEqual((3, 4), alignment.coverage())


class TestOffsetInfoMultipleGenomes(TestCase):
    """
    Test the offsetInfoMultipleGenomes function.
    """

    def testNoReferences(self):
        """
        Passing no Gb2Alignment instances must result in a ValueError.
        """
        error = r"^No Gb2Alignment instances given\.$"
        self.assertRaisesRegex(ValueError, error, offsetInfoMultipleGenomes, [], 0)

    def testDifferingReferences(self):
        """
        Passing Gb2Alignment instances with different reference ids must
        result in a ValueError.
        """
        features1 = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 1,
                    "stop": 6,
                },
            },
            DNARead("refId1", "GTTCCCG"),
            sars2=True,
        )

        features2 = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 1,
                    "stop": 6,
                },
            },
            DNARead("refId2", "GTTCCCG"),
            sars2=True,
        )

        alignment1 = Gb2Alignment(DNARead("genId", "ATTCCCG"), features1)
        alignment2 = Gb2Alignment(DNARead("genId", "ATTCCCG"), features2)

        error = (
            r"^Gb2Alignment instances with differing reference ids "
            r"passed to offsetInfoMultipleGenomes: refId1, refId2\.$"
        )
        self.assertRaisesRegex(
            ValueError, error, offsetInfoMultipleGenomes, [alignment1, alignment2], 0
        )

    def testOneGenome(self):
        """
        Passing in one Gb2Alignment should get the expected result and its
        components should be the same result as a single call to offsetInfo.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 1,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCCG"),
            sars2=True,
        )

        alignment = Gb2Alignment(DNARead("genId", "ATTCCCG"), features)
        multipleResult = offsetInfoMultipleGenomes([alignment], 0)

        self.assertEqual(
            {
                "alignmentOffset": 0,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "V",
                    "codon": "GTT",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
                "genomes": [
                    {
                        "aa": "I",
                        "codon": "ATT",
                        "frame": 0,
                        "id": "genId",
                        "aaOffset": 0,
                        "ntOffset": 0,
                    },
                ],
            },
            multipleResult,
        )

        oneResult = alignment.offsetInfo(0)
        self.assertEqual(oneResult["reference"], multipleResult["reference"])
        self.assertEqual(oneResult["genome"], multipleResult["genomes"][0])

        for key in "alignmentOffset", "featureName", "featureNames":
            self.assertEqual(oneResult[key], multipleResult[key])

    def testTwoGenomes(self):
        """
        Passing in two Gb2Alignment should get the expected result and its
        genome components should be the same as come back from single calls to
        offsetInfo.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 1,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCCG"),
            sars2=True,
        )

        alignment1 = Gb2Alignment(DNARead("genId1", "ATTCCCG"), features)
        alignment2 = Gb2Alignment(DNARead("genId2", "ATTCCCG"), features)
        multipleResult = offsetInfoMultipleGenomes([alignment1, alignment2], 0)

        self.assertEqual(
            {
                "alignmentOffset": 0,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "V",
                    "codon": "GTT",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
                "genomes": [
                    {
                        "aa": "I",
                        "codon": "ATT",
                        "frame": 0,
                        "id": "genId1",
                        "aaOffset": 0,
                        "ntOffset": 0,
                    },
                    {
                        "aa": "I",
                        "codon": "ATT",
                        "frame": 0,
                        "id": "genId2",
                        "aaOffset": 0,
                        "ntOffset": 0,
                    },
                ],
            },
            multipleResult,
        )

        oneResult = alignment1.offsetInfo(0)
        self.assertEqual(oneResult["reference"], multipleResult["reference"])
        self.assertEqual(oneResult["genome"], multipleResult["genomes"][0])

        oneResult = alignment2.offsetInfo(0)
        self.assertEqual(oneResult["reference"], multipleResult["reference"])
        self.assertEqual(oneResult["genome"], multipleResult["genomes"][1])

        for key in "alignmentOffset", "featureName", "featureNames":
            self.assertEqual(oneResult[key], multipleResult[key])

    def testThreeGenomesOneInsufficientlyCovered(self):
        """
        Passing in three Gb2Alignment should get the expected result,
        including when one insufficiently covers the reference and its genome
        components should be the same as come back from single calls to
        offsetInfo.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 1,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCCG"),
            sars2=True,
        )

        alignment1 = Gb2Alignment(DNARead("genId1", "ATTCCCG"), features)
        alignment2 = Gb2Alignment(DNARead("genId2", "ATT"), features)
        alignment3 = Gb2Alignment(DNARead("genId3", "ATTCCCG"), features)
        multipleResult = offsetInfoMultipleGenomes(
            [alignment1, alignment2, alignment3], 0, minReferenceCoverage=0.9
        )

        self.assertEqual(
            {
                "alignmentOffset": 0,
                "featureName": None,
                "featureNames": set(),
                "reference": {
                    "aa": "V",
                    "codon": "GTT",
                    "frame": 0,
                    "id": "refId",
                    "aaOffset": 0,
                    "ntOffset": 0,
                },
                "genomes": [
                    {
                        "aa": "I",
                        "codon": "ATT",
                        "frame": 0,
                        "id": "genId1",
                        "aaOffset": 0,
                        "ntOffset": 0,
                    },
                    None,
                    {
                        "aa": "I",
                        "codon": "ATT",
                        "frame": 0,
                        "id": "genId3",
                        "aaOffset": 0,
                        "ntOffset": 0,
                    },
                ],
            },
            multipleResult,
        )

        oneResult = alignment1.offsetInfo(0)
        self.assertEqual(oneResult["reference"], multipleResult["reference"])
        self.assertEqual(oneResult["genome"], multipleResult["genomes"][0])

        oneResult = alignment3.offsetInfo(0)
        self.assertEqual(oneResult["reference"], multipleResult["reference"])
        self.assertEqual(oneResult["genome"], multipleResult["genomes"][2])

        for key in "alignmentOffset", "featureName", "featureNames":
            self.assertEqual(oneResult[key], multipleResult[key])

    def testThreeGenomesAllInsufficientlyCovered(self):
        """
        Passing in three Gb2Alignment instances, all with insufficient cover,
        should result in a None value.
        """
        features = Features(
            {
                "surface glycoprotein": {
                    "name": "surface glycoprotein",
                    "start": 1,
                    "stop": 6,
                },
            },
            DNARead("refId", "GTTCCCG"),
            sars2=True,
        )

        alignment1 = Gb2Alignment(DNARead("genId1", "ATT"), features)
        alignment2 = Gb2Alignment(DNARead("genId2", "ATT"), features)
        alignment3 = Gb2Alignment(DNARead("genId3", "ATT"), features)

        self.assertIsNone(
            offsetInfoMultipleGenomes(
                [alignment1, alignment2, alignment3], 0, minReferenceCoverage=0.5
            )
        )
