from unittest import TestCase

from os.path import dirname, join

import sars2seq
from sars2seq.fasta import getSequence
from sars2seq.features import Features
from sars2seq.sequence import SARS2Sequence


DATA_DIR = join(dirname(dirname(sars2seq.__file__)), 'data')
REF_GB = join(DATA_DIR, 'NC_045512.2.gb')
FEATURES = Features(REF_GB)


class _TestMixin:
    """
    Test the SARS2Sequence class.
    """
    @classmethod
    def setup_class(cls):
        """
        Set up the class.
        """
        cls.seq = SARS2Sequence(cls.read, FEATURES)

    def testLength(self):
        self.assertGreater(len(self.read), 28000)

    def checkLocation(self, read, location, expected):
        """
        Check that a 1-based sequence location has the expected value.
        """
        self.assertEqual(expected, read.sequence[location - 1])

    def checkLocationIdentical(self, read1, read2, location):
        """
        Check that a 1-based sequence location has the same value in two reads.
        """
        self.assertEqual(read1.sequence[location - 1],
                         read2.sequence[location - 1])


class TestEPI_ISL_402125(TestCase):
    """
    Test the EPI_ISL_402125 sequence. This should be the same as the NCBI
    reference.
    """
    def testIdentical(self):
        """
        The EPI_ISL_402125 sequence from GISAID should be the same as the
        NC_045512.2 reference sequence from NCBI.
        """
        self.assertEqual(
            getSequence(join(DATA_DIR, 'EPI_ISL_402125.fasta')).sequence,
            getSequence(join(DATA_DIR, 'NC_045512.2.fasta')).sequence)


class TestEPI_ISL_601443(TestCase, _TestMixin):
    """
    Test the EPI_ISL_601433 sequence. This is the variant of concern referred
    to in https://www.gov.uk/government/publications/
    investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201
    """
    read = getSequence(join(DATA_DIR, 'EPI_ISL_601443.fasta'))

    def testSpikeDeletions(self):
        """
        The spike protein should have the three deletions.
        """
        spike = self.seq.feature('S')
        sequenceAa, _ = spike.aaSequences()
        for location in 69, 70, 144:
            self.checkLocation(sequenceAa, location, '-')

    def testSpikeMutationsAa(self):
        """
        The spike protein should have the expected amino acid changes. Note
        that the UK report does not include mention of D614G.
        """
        spike = self.seq.feature('S')
        sequenceAa, _ = spike.aaSequences()
        for location, aa in ((501, 'Y'), (570, 'D'), (614, 'G'), (681, 'H'),
                             (716, 'I'), (982, 'A'), (1118, 'H')):
            self.checkLocation(sequenceAa, location, aa)

    def testNucloecapsidMutationsNt(self):
        """
        The nucleocapsid genome should have the expected changes.
        Note that the UK report does not include mention of R203K or G204R.
        """
        spike = self.seq.feature('N')
        sequenceNt, _ = spike.ntSequences()
        # The 704 below is due to a change at the 1-based amino acid
        # location 235. That's a 0-based offset of 234 = 702 in the
        # genome. The change is an S -> via a TCT -> TTT mutation in the
        # middle position, or 703 in 0-based, and 704 in 1-based.
        for location, nt in (7, 'C'), (8, 'T'), (9, 'A'), (704, 'T'):
            self.checkLocation(sequenceNt, location, nt)

    def testNucloecapsidMutationsAa(self):
        """
        The nucleocapsid protein should have the expected amino acid changes.
        Note that the UK report does not include mention of R203K or G204R.
        """
        spike = self.seq.feature('N')
        sequenceAa, _ = spike.aaSequences()
        for location, aa in (3, 'L'), (203, 'K'), (204, 'R'), (235, 'F'):
            self.checkLocation(sequenceAa, location, aa)

    def testORF8MutationsAa(self):
        """
        The ORF8 protein should have the expected amino acid changes.
        """
        orf8 = self.seq.feature('orf8')
        sequenceAa, _ = orf8.aaSequences()
        for location, aa in (27, '*'), (52, 'I'), (73, 'C'):
            self.checkLocation(sequenceAa, location, aa)


class TestNC_045512(TestCase, _TestMixin):
    """
    Test the NC_045512.2 sequence, which should test as equal seeing as it is
    the default feature reference.
    """
    read = getSequence(join(DATA_DIR, 'NC_045512.2.fasta'))

    def testSpikeIdenticalNt(self):
        """
        The spike nucleotides should be identical.
        """
        spike = self.seq.feature('S')
        sequenceNt, referenceNt = spike.ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

    def testSpikeIdentical(self):
        """
        The spike protein should be identical.
        """
        spike = self.seq.feature('S')
        sequenceAa, referenceAa = spike.aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)
