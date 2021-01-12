from unittest import TestCase

from os.path import dirname, join

import sars2seq
from sars2seq.checker import Checker, AAChecker, NTChecker
from sars2seq.fasta import getSequence
from sars2seq.features import Features
from sars2seq.genome import SARS2Genome

DATA_DIR = join(dirname(dirname(sars2seq.__file__)), 'data')
REF_GB = join(DATA_DIR, 'NC_045512.2.gb')
FEATURES = Features(REF_GB)


class Test_EPI_ISL_601443(TestCase):
    """
    Test the EPI_ISL_601433 sequence. This is the variant of concern
    (VOC 202012/01) referred to in https://www.gov.uk/government/publications/
    investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201
    """
    genomeRead = getSequence(join(DATA_DIR, 'EPI_ISL_601443.fasta'))
    genome = SARS2Genome(genomeRead, FEATURES)

    def testN501Y(self):
        """
        The variant has the N501Y change.
        """
        checker = Checker('spike', 'N501Y', False)
        self.assertTrue(checker(self.genome))

    def testN501YandA570D(self):
        """
        The variant has the A570D change. Check with the base Checker class.
        """
        checker = (Checker('spike', 'N501Y', False) and
                   Checker('spike', 'A570D', False))
        self.assertTrue(checker(self.genome))

    def testN501YAA(self):
        """
        Check if the variant has the NY501 change.
        """
        checker = AAChecker('spike', 'N501Y')
        self.assertTrue(checker(self.genome))

    def testN501YandA570DAA(self):
        """
        Check if the variant has the NY501 and A570D changes.
        """
        checker = AAChecker('spike', 'N501Y') & AAChecker('spike', 'A570D')
        self.assertTrue(checker(self.genome))

    def testN501YorA570DAA(self):
        """
        Check if the variant has the NY501 or the A570D change (both are
        present).
        """
        checker = AAChecker('spike', 'N501Y') | AAChecker('spike', 'A570D')
        self.assertTrue(checker(self.genome))

    def testN501YandA571DAA(self):
        """
        Check if the variant has the NY501 and a A571D change. The latter is
        not the case, so this returns False.
        """
        checker = AAChecker('spike', 'N501Y') & AAChecker('spike', 'A571D')
        self.assertFalse(checker(self.genome))

    def testN501YorA571DAA(self):
        """
        Check if the variant has the NY501 or a A571D change. The latter is not
        the case.
        """
        checker = AAChecker('spike', 'N501Y') | AAChecker('spike', 'A571D')
        self.assertTrue(checker(self.genome))

    def testA571DN501YorAA(self):
        """
        Check if the variant has the NY501 or a A571D change. The former is not
        the case.
        """
        checker = AAChecker('spike', 'A571D') | AAChecker('spike', 'N501Y')
        self.assertTrue(checker(self.genome))

    def testSpikeAAandNucleocapsidNtMutations(self):
        """
        The spike should have the N501Y change, several deletions, and the
        nucleocapsid genome should have the expected changes.
        """
        checker = (NTChecker('N', 'G7C A8T T9A G608A G609A G610C C704T') &
                   AAChecker('spike', 'N501Y') &
                   AAChecker('spike', 'H69-') &
                   AAChecker('spike', 'V70-') &
                   AAChecker('spike', 'Y144-'))

        self.assertTrue(checker(self.genome))

    def testSpikeAAandNucleocapsidNtMutationsCombined(self):
        """
        The spike should have the N501Y change, several deletions, and the
        nucleocapsid genome should have the expected changes.
        """
        checker = (NTChecker('N', 'G7C A8T T9A G608A G609A G610C C704T') &
                   AAChecker('spike', 'N501Y H69- V70- Y144-'))

        self.assertTrue(checker(self.genome))
