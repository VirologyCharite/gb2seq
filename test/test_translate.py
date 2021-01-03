from unittest import TestCase

from sars2seq.translate import translate


class TestTranslate(TestCase):
    """
    Test the translate function.
    """
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
        the end of the slippery sequence, then continues at the beginning of
        the slippery sequence.
        """
        slipperySeq = 'TTTAAAC'
        repeats = int(15000 / 3)
        seq = 'AA' + ('AAA' * repeats) + slipperySeq + 'CCCCCCTAAAA'
        # The sequence that gets translated is:
        # AAA 'repeats' times, then AA TTTAAAC TTTAAAC CCCCCCTAAAA
        # Regrouping:
        # AAA 'repeats' times, then AAT TTA AAC TTT AAA CCC CCC CTA AAA
        # K   'repeats' times, then  N   L   N   F   K   P   P   L   K
        expected = 'K' * repeats + 'NLNFKPPLK'
        self.assertEqual(expected, translate(seq, 'ORF1ab polyprotein'))
