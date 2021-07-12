from Bio.Seq import Seq
from unittest import TestCase

from sars2seq.translate import (
    translate, NoSlipperySequenceError, NoStopCodonError,
    StopCodonTooDistantError, SLIPPERY_SEQUENCE, translateSpike,
    TranslatedSequenceLengthError, KNOWN_INSERTIONS)


class TestTranslate(TestCase):
    """
    Test the translate function.
    """
    def testNoSlipperySequencs(self):
        """
        An ORF1ab polyprotein sequence must have a slippery sequence.
        """
        error = r'^No slippery sequence found\.$'
        self.assertRaisesRegex(NoSlipperySequenceError, error, translate,
                               'AAATTT', 'ORF1ab polyprotein')

    def testNoStopCodonFollowingTheSlipperySequence(self):
        """
        An ORF1ab polyprotein sequence must have a stop codon after the
        slippery sequence.
        """
        error = (r'^Could not find a stop codon downstream from the start of '
                 r'the slippery sequence at location 13001\.$')
        sequence = 'A' * 13000 + SLIPPERY_SEQUENCE
        self.assertRaisesRegex(NoStopCodonError, error, translate,
                               sequence, 'ORF1ab polyprotein')

    def testDistantStopCodonFollowingTheSlipperySequence(self):
        """
        An ORF1ab polyprotein sequence must have a stop codon not too far
        downstream of the slippery sequence.
        """
        error = (r'The stop codon was too far \(107 nucleotides\) downstream '
                 r'\(max allowed distance is 20\) from the start of the '
                 r'slippery sequence at location 13001\.$')
        sequence = 'A' * 13000 + SLIPPERY_SEQUENCE + 'A' * 100 + 'TAA'
        self.assertRaisesRegex(StopCodonTooDistantError, error, translate,
                               sequence, 'ORF1ab polyprotein')

    def testEmpty(self):
        """
        An empty nt sequence must translate to an empty aa sequence.
        """
        self.assertEqual('', translate(''))

    def testIncomplete(self):
        """
        An incomplete nt codon must translate to an X.
        """
        self.assertEqual('X', translate('AA'))

    def testAAA(self):
        """
        An AAA codon must translate to a Lysine (K).
        """
        self.assertEqual('K', translate('AAA'))

    def testNameWithAAA(self):
        """
        An AAA codon must translate to a Lysine (K) when a name other than
        'ORF1ab polyprotein' is passed.
        """
        self.assertEqual('K', translate('AAA', 'name'))

    def testAAATT(self):
        """
        An AAATT sequence must translate to a KX.
        """
        self.assertEqual('KX', translate('AAATT'))

    def testAAATTT(self):
        """
        An AAATTT sequence must translate to a KF.
        """
        self.assertEqual('KF', translate('AAATTT'))

    def testORF1abPolyprotein(self):
        """
        Test an ORF1ab polyprotein. The translation goes all the way through
        the end of the slippery sequence, then continues starting at the final
        nucleotide of the slippery sequence.
        """
        slipperySeq = 'TTTAAAC'
        repeats = int(15000 / 3)
        seq = 'AA' + ('AAA' * repeats) + slipperySeq + 'CCCTAAAA'
        # The sequence that gets translated is:
        # AAA 'repeats' times, then AA TTTAAAC C CCCTAAAA
        # Regrouping, we have:
        # AAA 'repeats' times, then AAT TTA AAC CCC CTA AAA
        # K   'repeats' times, then  N   L   N   P   L   K
        expected = 'K' * repeats + 'NLNPLK'
        self.assertEqual(expected, translate(seq, 'ORF1ab polyprotein'))


class TestTranslateSpike(TestCase):
    """
    Tests for the translate.translateSpike function.
    """
    def testNoGapsCorrectSequence(self):
        """
        A sequence with no gaps must be translated correctly.
        """
        seq = 'TTGGTTGTTTATTACCAC'
        self.assertEqual(Seq(seq).translate(), translateSpike(seq))

    def testNoGapsCorrectSequenceNotMultipleOfThree(self):
        """
        A sequence with no gaps that is not a multiple of three must
        raise an AssertionError.
        """
        seq = 'TTGGTTGTTTATTACCA'
        error = (r'^The length of a sequence to be translated must '
                 r'be a multiple of 3 but is 17\.$')
        self.assertRaisesRegex(TranslatedSequenceLengthError, error,
                               translateSpike, seq)

    def testInFrameGapCorrectLength(self):
        """
        A sequence with an in frame gap must have the correct length.
        """
        seq = 'TTG---GTTTATTACCAC'
        self.assertEqual(len(seq) / 3, len(translateSpike(seq)))

    def testPlus1GapCorrectLength(self):
        """
        A sequence with an out of frame gap (AA-) must have the correct length.
        """
        seq = 'TT---GGTTTATTACCAC'
        self.assertEqual(len(seq) / 3, len(translateSpike(seq)))

    def testPlus2GapCorrectLength(self):
        """
        A sequence with an out of frame gap (A--) must have the correct length.
        """
        seq = 'TTGG---TTTATTACCAC'
        self.assertEqual(len(seq) / 3, len(translateSpike(seq)))

    def testInFrameGapCorrectLocation(self):
        """
        A sequence with an in frame gap must be in the correct location.
        """
        seq = 'TTG---GTTTATTACCAC'
        self.assertEqual('L-VYYH', translateSpike(seq))

    def testInFrameGapAmbiguousCorrectLocation(self):
        """
        A sequence with an in frame gap and an ambiguity must be translated
        correcty.
        """
        seq = 'TTG---GTTTANTACCAC'
        self.assertEqual('L-VXYH', translateSpike(seq))

    def testPlus1GapCorrectLocation(self):
        """
        A sequence with an out of frame gap (TT-) must be in the correct
        location.
        """
        seq = 'TT---GGTTTATTACCAC'
        self.assertEqual('L-VYYH', translateSpike(seq))

    def testPlus1GapAdjacentAmbiguityCorrectLocation(self):
        """
        A sequence with an out of frame gap (TT-) must be translated correctly.
        """
        seq = 'TT---NGTTTATTACCAC'
        self.assertEqual('X-VYYH', translateSpike(seq))

    def testPlus2GapCorrectLocation(self):
        """
        A sequence with an out of frame gap (G--) must be in the correct
        location.
        """
        seq = 'TTGG---TTTATTACCAC'
        self.assertEqual('LV-YYH', translateSpike(seq))

    def test6970S71FCorrectLocation(self):
        """
        The 69-70 deletion with a substitution leading to S71F must be
        aligned correctly.
        """
        seq = 'CATGCTAT------CTTTGGGACC'
        self.assertEqual('HAI--FGT', translateSpike(seq))

    def test6970G72VCorrectLocation(self):
        """
        The 69-70 deletion with a substitution leading to G72V must be
        aligned correctly.
        """
        seq = 'CATGCTAT------CTCTGTGACC'
        self.assertEqual('HAI--SVT', translateSpike(seq))

    def testB16172_156_157GapCorrectLocation(self):
        """
        The gap at 156/157 in B.1.617.2 must be in the correct location.
        """
        seq = 'GAAAGTG------GAGTTTATTCTAGT'
        self.assertEqual('ESG--VYSS', translateSpike(seq))


class TestKnownInsertions(TestCase):
    """
    Tests for KNOWN_INSERTIONS.
    """
    def testNoDuplicates(self):
        """
        There must be no duplicated sequences of insertions.
        """
        seqs = [knownInsertion[0] for knownInsertion in KNOWN_INSERTIONS]
        self.assertEqual(len(seqs), len(set(seqs)))

    def testSliceStartLowerSliceStop(self):
        """
        sliceStart must be smaller than sliceStop.
        """
        for t, findStart, findStop, sliceStart, sliceStop in KNOWN_INSERTIONS:
            self.assertTrue(sliceStart < sliceStop)

    def testFindStartLowerFindStop(self):
        """
        findStart must be smaller than findStop.
        """
        for t, findStart, findStop, sliceStart, sliceStop in KNOWN_INSERTIONS:
            self.assertTrue(findStart < findStop)

    def testSliceStopLowerFindStop(self):
        """
        sliceStop must be smaller than findStop.
        """
        for t, findStart, findStop, sliceStart, sliceStop in KNOWN_INSERTIONS:
            self.assertTrue(sliceStop < findStop)

    def testSliceStartHigherFindStart(self):
        """
        sliceStart must be higher than findStart.
        """
        for t, findStart, findStop, sliceStart, sliceStop in KNOWN_INSERTIONS:
            self.assertTrue(findStart < sliceStart)
