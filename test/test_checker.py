from unittest import TestCase

from .fasta import getSequence

from sars2seq import DATA_DIR
from sars2seq.alignment import SARS2Alignment
from sars2seq.checker import Checker, AAChecker, NTChecker
from sars2seq.features import Features

REF_GB = DATA_DIR / "NC_045512.2.gb"
FEATURES = Features(REF_GB)


class Test_EPI_ISL_601443(TestCase):
    """
    Test the EPI_ISL_601433 sequence. This is the variant of concern
    (VOC 202012/01) referred to in https://www.gov.uk/government/publications/
    investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201
    """

    genomeRead = getSequence(DATA_DIR / "EPI_ISL_601443.fasta")
    alignment = SARS2Alignment(genomeRead, FEATURES)

    def testIndexError(self):
        """
        If an check on a non-existent index is attempted, an IndexError must
        be raised.
        """
        checker = Checker("spike", "N500001Y", aa=True)
        error = (
            r"^Index 500000 out of range trying to access feature "
            r"'spike' of length 1274 sequence 'NC_045512.2 \(surface "
            r"glycoprotein\)' via expected change specification "
            r"'N500001Y'\."
        )
        self.assertRaisesRegex(IndexError, error, checker, self.alignment)

    def testN501Y(self):
        """
        The variant has the N501Y change.
        """
        checker = Checker("spike", "N501Y", aa=True)
        self.assertTrue(checker(self.alignment))

    def testN501YandA570D(self):
        """
        The variant has the A570D change. Check with the base Checker class.
        """
        checker = Checker("spike", "N501Y", aa=True) and Checker(
            "spike", "A570D", aa=True
        )
        self.assertTrue(checker(self.alignment))

    def testN501YAA(self):
        """
        Check if the variant has the NY501 change.
        """
        checker = AAChecker("spike", "N501Y")
        self.assertTrue(checker(self.alignment))

    def testN501YandA570DAA(self):
        """
        Check if the variant has the NY501 and A570D changes.
        """
        checker = AAChecker("spike", "N501Y") & AAChecker("spike", "A570D")
        self.assertTrue(checker(self.alignment))

    def testN501YorA570DAA(self):
        """
        Check if the variant has the NY501 or the A570D change (both are
        present).
        """
        checker = AAChecker("spike", "N501Y") | AAChecker("spike", "A570D")
        self.assertTrue(checker(self.alignment))

    def testN501YandA571DAA(self):
        """
        Check if the variant has the NY501 and a A571D change. The latter is
        not the case, so this returns False.
        """
        checker = AAChecker("spike", "N501Y") & AAChecker("spike", "A571D")
        self.assertFalse(checker(self.alignment))

    def testN501YorA571DAA(self):
        """
        Check if the variant has the NY501 or a A571D change. The latter is not
        the case.
        """
        checker = AAChecker("spike", "N501Y") | AAChecker("spike", "A571D")
        self.assertTrue(checker(self.alignment))

    def testA571DN501YorAA(self):
        """
        Check if the variant has the NY501 or a A571D change. The former is not
        the case.
        """
        checker = AAChecker("spike", "A571D") | AAChecker("spike", "N501Y")
        self.assertTrue(checker(self.alignment))

    def testSpikeAAandNucleocapsidNtMutations(self):
        """
        The spike should have the N501Y change, several deletions, and the
        nucleocapsid genome should have the expected changes.
        """
        checker = (
            NTChecker("N", "G7C A8T T9A G608A G609A G610C C704T")
            & AAChecker("spike", "N501Y")
            & AAChecker("spike", "H69-")
            & AAChecker("spike", "V70-")
            & AAChecker("spike", "Y144-")
        )

        self.assertTrue(checker(self.alignment))

    def testSpikeAAandNucleocapsidNtMutationsCombined(self):
        """
        The spike should have the N501Y change, several deletions, and the
        nucleocapsid genome should have the expected changes.
        """
        checker = NTChecker("N", "G7C A8T T9A G608A G609A G610C C704T") & AAChecker(
            "spike", "N501Y H69- V70- Y144-"
        )

        self.assertTrue(checker(self.alignment))

    def testAndOperator(self):
        """
        The Checker object should stay unchanged after using the & operator
        """
        checker1 = AAChecker("spike", "N501Y")
        checker1_func_id = id(checker1._func)
        checker2 = AAChecker("spike", "H69-")

        checker3 = checker1 & checker2

        self.assertNotEqual(checker1, checker3)
        self.assertIsNot(checker1, checker3)
        self.assertEqual(id(checker1._func), checker1_func_id)

    def testOrOperator(self):
        """
        The Checker object should stay unchanged after using the | operator
        """
        checker1 = AAChecker("spike", "N501Y")
        checker1_func_id = id(checker1._func)
        checker2 = AAChecker("spike", "H69-")

        checker3 = checker1 | checker2

        self.assertNotEqual(checker1, checker3)
        self.assertIsNot(checker1, checker3)
        self.assertEqual(id(checker1._func), checker1_func_id)
