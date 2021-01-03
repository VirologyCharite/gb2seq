from unittest import TestCase

from sars2seq.change import splitChange


class TestSplitChange(TestCase):
    """
    Test the splitChange function.
    """
    def testNoBases(self):
        """
        If no base are given, a ValueError must be raised.
        """
        error = (r"^Change string '300' does not include a reference or "
                 r"genome base\.$")
        self.assertRaisesRegex(ValueError, error, splitChange, '300')

    def testNoReferenceBase(self):
        """
        If no reference base is given, the result must be as expected.
        """
        self.assertEqual((None, 299, 'D'), splitChange('300D'))

    def testReferenceBaseIsGap(self):
        """
        If the reference base is a gap, the result must be as expected.
        """
        self.assertEqual(('-', 299, 'D'), splitChange('-300D'))

    def testNoGenomeBase(self):
        """
        If no genome base is given, the result must be as expected.
        """
        self.assertEqual(('D', 299, None), splitChange('D300'))

    def testGenomeBaseIsGap(self):
        """
        If the genome base is a gap, the result must be as expected.
        """
        self.assertEqual((None, 299, '-'), splitChange('300-'))

    def testBothBases(self):
        """
        If a reference and a genome base are given, the result must be as
        expected.
        """
        self.assertEqual(('D', 299, 'E'), splitChange('D300E'))

    def testBothBasesAreGaps(self):
        """
        If the reference and a genome base are both gaps, the result must be as
        expected.
        """
        self.assertEqual(('-', 299, '-'), splitChange('-300-'))
