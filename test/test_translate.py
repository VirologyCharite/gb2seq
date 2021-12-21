from Bio.Seq import Seq
from unittest import TestCase

from dark.reads import AARead
from sars2seq.translate import (
    translate, NoSlipperySequenceError, NoStopCodonError,
    StopCodonTooDistantError, TranslatedReferenceAndGenomeLengthError,
    TranslatedSequenceLengthError, SLIPPERY_SEQUENCE, translateSpike,
    KNOWN_INSERTIONS, getSubstitutionsString)


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


class TestGetSubstitutionsString(TestCase):
    """
    Test the getSubstitutionsString function.
    """
    def testUnequalLengths(self):
        """
        If sequences with unequal lengths are passed, a
        TranslatedReferenceAndGenomeLengthError must be raised.
        """
        reference = AARead('id', 'A')
        genome = AARead('id', 'MK')
        error = r'^Reference and genome lengths unequal \(1 != 2\)\.$'
        self.assertRaisesRegex(
            TranslatedReferenceAndGenomeLengthError,
            error, getSubstitutionsString, reference, genome)

    def testEmpty(self):
        """
        If the empty string is passed for both reference and genome, the empty
        string must be returned.
        """
        reference = AARead('id', '')
        genome = AARead('id', '')
        self.assertEqual('', getSubstitutionsString(reference, genome))

    def testOneLetterIdentical(self):
        """
        If two identical one-AA sequences are passed, the empty string must
        be returned.
        """
        reference = AARead('id', 'K')
        genome = AARead('id', 'K')
        self.assertEqual('', getSubstitutionsString(reference, genome))

    def testOneLetterDifferent(self):
        """
        If two different one-AA sequences are passed, a string showing the
        change at position 1 must be retuned.
        """
        reference = AARead('id', 'S')
        genome = AARead('id', 'K')
        self.assertEqual('S1K', getSubstitutionsString(reference, genome))

    def testOneLetterReferenceGap(self):
        """
        If two different one-AA sequences are passed, with a reference gap,
        a string showing the change at position 1 must be retuned.
        """
        reference = AARead('id', '-')
        genome = AARead('id', 'K')
        self.assertEqual('-1K', getSubstitutionsString(reference, genome))

    def testOneLetterGenomeGap(self):
        """
        If two different one-AA sequences are passed, with a genome gap,
        a string showing the change at position 1 must be retuned.
        """
        reference = AARead('id', 'S')
        genome = AARead('id', '-')
        self.assertEqual('S1-', getSubstitutionsString(reference, genome))

    def testOneLetterGenomeX(self):
        """
        If two different one-AA sequences are passed, with an X in the genome
        gap, a string showing the change at position 1 must be retuned.
        """
        reference = AARead('id', 'S')
        genome = AARead('id', 'X')
        self.assertEqual('no coverage 1-1',
                         getSubstitutionsString(reference, genome))

    def testTwoLettersBothDifferent(self):
        """
        If two different two-AA sequences are passed, a string showing the
        change at positions 1 and 2 must be retuned.
        """
        reference = AARead('id', 'SP')
        genome = AARead('id', 'KL')
        self.assertEqual('S1K; P2L', getSubstitutionsString(reference, genome))

    def testInitialStringOfXs(self):
        """
        If the genome starts with Xs, they must be summarized correctly.
        """
        reference = AARead('id', 'TRSP')
        genome = AARead('id', 'XXXL')
        self.assertEqual('no coverage 1-3; P4L',
                         getSubstitutionsString(reference, genome))

    def testFinalStringOfXs(self):
        """
        If the genome ends with Xs, they must be summarized correctly.
        """
        reference = AARead('id', 'TRSP')
        genome = AARead('id', 'LXXX')
        self.assertEqual('T1L; no coverage 2-4',
                         getSubstitutionsString(reference, genome))

    def testStringOfXs(self):
        """
        If the genome has a string of Xs, they must be summarized correctly.
        """
        reference = AARead('id', 'STRSP')
        genome = AARead('id', 'KXXXL')
        self.assertEqual('S1K; no coverage 2-4; P5L',
                         getSubstitutionsString(reference, genome))

    def testTwoStringsOfXs(self):
        """
        If the genome has two strings of Xs, they must be summarized correctly.
        """
        reference = AARead('id', 'STRSPFFFFFA')
        genome = AARead('id', 'KXXXLXXXXXT')
        self.assertEqual('S1K; no coverage 2-4; P5L; no coverage 6-10; A11T',
                         getSubstitutionsString(reference, genome))

    def testUnreportedXsIssue21(self):
        """
        Test with the protein sequences that caused the issue in
        https://github.com/VirologyCharite/sars2seq/issues/21
        The output should include the fact that site 417 is not covered
        in the genome.
        """
        reference = AARead(
            'id',
            'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTW'
            'FHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVI'
            'KVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREF'
            'VFKNIDGYFKIYSKHTPINLVR---DLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGD'
            'SSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSN'
            'FRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYG'
            'VSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKV'
            'GGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYR'
            'VVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTT'
            'DAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYS'
            'TGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGA'
            'ENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNR'
            'ALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLA'
            'DAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAA'
            'LQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQ'
            'ALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRAS'
            'ANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDG'
            'KAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFK'
            'EELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWP'
            'WYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*')

        genome = AARead(
            'id',
            'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTW'
            'FHVI--SGTNGTKRFDNPVLPFNDGVYFASIEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVI'
            'KVCEFQFCNDPFLXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXPFLMDLEGKQGNFKNLREF'
            'VFKNIDGYFKIYSKHTPII-VREPEDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGD'
            'SSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSN'
            'FRVQPTESIVRFPNITNLCPFDEVFNATRFASVYAXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXYKLPDDFTGCVIAWNSNKLDSKV'
            'SGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGFNCYFPLRSYSFRPTYGVGHQPYR'
            'VVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLKGTGVLTESNKKFLPFQQFGRDIADTT'
            'DAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYS'
            'TGSNVFQTRAGCLIGAEYVNNSYECDIPIGAGICASYQTQTKSHRRARSVASQSIIAYTMSLGV'
            'ENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNXXLQYGSFCTQLKR'
            'ALTGIAVEQDKNTQEVFAQVKQIYKTPPIKYFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLA'
            'DAGFIKQYGDCLGDIAARDLICAQKFKGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAA'
            'LQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNHNAQ'
            'ALNTLVKQLSSKFGAISSVLNDIFSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRAS'
            'ANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDG'
            'KAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFK'
            'EELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWP'
            'WYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*')

        self.assertEqual('',
                         getSubstitutionsString(reference, genome))

    def testUnreportedXsIssue21Simple(self):
        """
        Small test to trigger the issue in
        https://github.com/VirologyCharite/sars2seq/issues/21
        """
        reference = AARead('id', 'C-LABF')
        genome = AARead('id', 'CSMXXF')

        # Note that the code is currently returning '-2S; L2M; no coverage 3'
        # which has several problems.
        self.assertEqual('-2S; L2M; no coverage 4-5',
                         getSubstitutionsString(reference, genome))
