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
    Mixin for SARS2Sequence class tests.
    """
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

    def checkChanges(self, changes, sequence, reference):
        """
        Check that a set of changes all happened as expected.

        @param changes: A C{str} specification in the form of space-separated
            RNS strings, where R is a reference base, N is an integer offset,
            and S is a sequence base. So, e.g., 'L28S P1003Q' indicates that
            we expected a change from 'L' to 'S' at offset 28 and from 'P' to
            'Q' at offset 1003.
        """
        for change in changes.split():
            location = int(change[1:-1])
            self.checkLocation(reference, location, change[0])
            self.checkLocation(sequence, location, change[-1])


class Test_EPI_ISL_402125(TestCase):
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


class Test_EPI_ISL_601443(TestCase, _TestMixin):
    """
    Test the EPI_ISL_601433 sequence. This is the variant of concern
    (VOC 202012/01) referred to in https://www.gov.uk/government/publications/
    investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201
    """
    read = getSequence(join(DATA_DIR, 'EPI_ISL_601443.fasta'))
    seq = SARS2Sequence(read, FEATURES)

    def testSpikeDeletionsAa(self):
        """
        The spike protein should have the three deletions.
        """
        self.checkChanges('H69- V70- Y144-',
                          *self.seq.feature('spike').aaSequences())

    def testSpikeMutationsAa(self):
        """
        The spike protein should have the expected amino acid changes. Note
        that the UK report does not include mention of D614G (but they also
        don't say what reference their SNPs are relative to)
        """
        self.checkChanges('N501Y A570D D614G P681H T716I S982A D1118H',
                          *self.seq.feature('spike').aaSequences())

    def testORF1aDeletionsAa(self):
        """
        The ORF1a protein should have the expected deletions.
        """
        self.checkChanges('S3675- G3676- F3677-',
                          *self.seq.feature('orf1a').aaSequences())

    def testORF1aMutationsAa(self):
        """
        The ORF1a protein should have the expected amino acid changes.
        """
        self.checkChanges('T1001I A1708D I2230T',
                          *self.seq.feature('orf1a').aaSequences())

    def testORF1abDeletionsAa(self):
        """
        The ORF1ab protein should have the expected deletions.
        """
        self.checkChanges('S3675- G3676- F3677-',
                          *self.seq.feature('orf1ab').aaSequences())

    # TODO: Is this actually correct???
    def testORF1abInsertionsAa(self):
        """
        The ORF1ab protein should have the expected insertions.
        """
        self.checkChanges('-4402F -4403K',
                          *self.seq.feature('orf1ab').aaSequences())

    def testORF1abMutationsAa(self):
        """
        The ORF1ab protein should have the expected amino acid changes.
        """
        self.checkChanges('T1001I A1708D I2230T P4717L',
                          *self.seq.feature('orf1ab').aaSequences())

    def testNucleocapsidMutationsNt(self):
        """
        The nucleocapsid genome should have the expected changes.
        """
        # The 704 below is due to a change at the 1-based amino acid
        # location 235. That's a 0-based offset of 234 = 702 in the
        # genome. The change is an S -> via a TCT -> TTT mutation in the
        # middle position, or 703 in 0-based, and 704 in 1-based.
        self.checkChanges('G7C A8T T9A G608A G609A G610C C704T',
                          *self.seq.feature('N').ntSequences())

    def testNucleocapsidMutationsAa(self):
        """
        The nucleocapsid protein should have the expected amino acid changes.
        Note that the UK report does not include mention of R203K or G204R.
        """
        self.checkChanges('D3L R203K G204R S235F',
                          *self.seq.feature('N').aaSequences())

    def testORF8MutationsAa(self):
        """
        The ORF8 protein should have the expected amino acid changes.
        """
        self.checkChanges('Q27* R52I Y73C',
                          *self.seq.feature('orf8').aaSequences())

    def testSNPs(self):
        """
        The SNPs in Table 1 of the report mentioned the docstring of this class
        should all be present.
        """
        # Table 1 seems incorrect. Instead of C5388A, the sequence has a
        # G. Instead of C23271A the sequence has a G. Instead of C23604A
        # the sequence has a T. That's as far as I checked. I checked at
        # the given offset and also +/- 18 nucleotides due to the
        # deletions. They don't say what their SNPs/deletions are relative
        # to or how they got their offsets.
        for location, nt in ((3267, 'T'), (5388, 'G'), (6954, 'C'),
                             (23063, 'T'), (23271, 'G'), (23604, 'T')):
            self.checkLocation(self.read, location, nt)


class Test_BavPat2(TestCase, _TestMixin):
    """
    Test the BavPat2 sequence. This is Bavarian patient #2.
    """
    read = getSequence(join(DATA_DIR, 'BavPat2.fasta'))
    seq = SARS2Sequence(read, FEATURES)

    def testSpikeMutationsNt(self):
        """
        The spike genome should have the expected change.
        """
        self.checkChanges('A1841G',
                          *self.seq.feature('spike').ntSequences())

    def testSpikeMutationsAa(self):
        """
        The spike protein should have the expected amino acid change.
        """
        self.checkChanges('D614G',
                          *self.seq.feature('spike').aaSequences())

    def testORF1aMutationsNt(self):
        """
        The ORF1a genome should have the expected change.
        """
        self.checkChanges('C2772T',
                          *self.seq.feature('orf1a').ntSequences())

    def testORF1abInsertionsAa(self):
        """
        The ORF1ab protein should have the expected insertions.
        """
        self.checkChanges('-4402F -4403K',
                          *self.seq.feature('orf1ab').aaSequences())

    def testNucleocapsidIdentical(self):
        """
        The nucleocapsid genome should be identical to the reference.
        """
        sequenceNt, referenceNt = self.seq.feature('N').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

    def testORF8Identical(self):
        """
        The ORF8 genome should be identical to the reference.
        """
        sequenceNt, referenceNt = self.seq.feature('orf8').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

        sequenceAa, referenceAa = self.seq.feature('orf8').aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)

    def testEnvelopeIdentical(self):
        """
        The envelope should be identical to the reference.
        """
        sequenceNt, referenceNt = self.seq.feature('E').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

        sequenceAa, referenceAa = self.seq.feature('E').aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)

    def testMembraneIdentical(self):
        """
        The membrane should be identical to the reference.
        """
        sequenceNt, referenceNt = self.seq.feature('M').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

        sequenceAa, referenceAa = self.seq.feature('M').aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)

    def testRdRpIdentical(self):
        """
        The polymerase should be identical to the reference.
        """
        sequenceNt, referenceNt = self.seq.feature('rdrp').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

        sequenceAa, referenceAa = self.seq.feature('rdrp').aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)


class Test_NC_045512(TestCase, _TestMixin):
    """
    Test the NC_045512.2 sequence, which should test as equal seeing as it is
    the default feature reference.
    """
    read = getSequence(join(DATA_DIR, 'NC_045512.2.fasta'))
    seq = SARS2Sequence(read, FEATURES)

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
