from unittest import TestCase

from os.path import dirname, join

import sars2seq
from sars2seq.fasta import getSequence
from sars2seq.features import Features
from sars2seq.genome import SARS2Genome


DATA_DIR = join(dirname(dirname(sars2seq.__file__)), 'data')
REF_GB = join(DATA_DIR, 'NC_045512.2.gb')
FEATURES = Features(REF_GB)


class _Mixin:
    """
    Mixin for SARS2Genome class tests.
    """
    def testLength(self):
        self.assertGreater(len(self.genomeRead), 28000)

    def check(self, name, changes, nt):
        """
        Check that a set of changes all happened as expected.

        @param name: The C{str} name of the feature to check (e.g., 'nsp2').
        @param changes: A C{str} specification in the form of space-separated
            RNS strings, where R is a reference base, N is an integer offset,
            and S is a sequence base. So, e.g., 'L28S P1003Q' indicates that
            we expected a change from 'L' to 'S' at offset 28 and from 'P' to
            'Q' at offset 1003.
        @param nt: If C{True} check nucleotide sequences. Else protein.
        """
        feature = self.genome.feature(name)
        _, errorCount, result = feature.check(changes, nt)
        if errorCount:
            for change, (referenceOK, genomeOK) in result.items():
                if not referenceOK:
                    self.fail(f'Reference base check failed on {change!r}')
                if not genomeOK:
                    self.fail(f'Genome base check failed on {change!r}')


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


class Test_EPI_ISL_601443(TestCase, _Mixin):
    """
    Test the EPI_ISL_601433 sequence. This is the variant of concern
    (VOC 202012/01) referred to in https://www.gov.uk/government/publications/
    investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201
    """
    genomeRead = getSequence(join(DATA_DIR, 'EPI_ISL_601443.fasta'))
    genome = SARS2Genome(genomeRead, FEATURES)

    def testSpikeDeletionsAa(self):
        """
        The spike protein should have the three deletions.
        """
        self.check('spike', 'H69- V70- Y144-', nt=False)

    def testSpikeMutationsAa(self):
        """
        The spike protein should have the expected amino acid changes. Note
        that the UK report does not include mention of D614G (but they also
        don't say what reference their SNPs are relative to)
        """
        self.check('spike', 'N501Y A570D D614G P681H T716I S982A D1118H',
                   False)

    def testSpikeDeletionVariant(self):
        """
        The genome must fulfil all the requirements of a spike deletion
        variant.
        """
        testCount, errorCount, result = self.genome.checkVariant(
            'spikeDeletion')
        self.assertEqual(2, testCount)
        self.assertEqual(0, errorCount)
        self.assertEqual(
            {
                'spike': {
                    'aa': {
                        '69-': (True, True),
                        '70-': (True, True),
                    },
                },
            },
            result
        )

    def testN501YVariant(self):
        """
        The genome must be an N501Y variant.
        """
        testCount, errorCount, result = self.genome.checkVariant('N501Y')
        self.assertEqual(1, testCount)
        self.assertEqual(0, errorCount)
        self.assertEqual(
            {
                'spike': {
                    'aa': {
                        'N501Y': (True, True),
                    },
                },
            },
            result
        )

    def testVariantOfConcern20201201(self):
        """
        The genome must fulfil all the requirements of the UK variant of
        concern 202012/01.
        """
        testCount, errorCount, _ = self.genome.checkVariant('VOC_20201201_UK')
        self.assertEqual(20, testCount)
        self.assertEqual(0, errorCount)

    def testORF1aDeletionsAa(self):
        """
        The ORF1a protein should have the expected deletions.
        """
        self.check('orf1a', 'S3675- G3676- F3677-', False)

    def testORF1aMutationsAa(self):
        """
        The ORF1a protein should have the expected amino acid changes.
        """
        self.check('orf1a', 'T1001I A1708D I2230T', False)

    def testORF1abDeletionsAa(self):
        """
        The ORF1ab protein should have the expected deletions.
        """
        self.check('orf1ab', 'S3675- G3676- F3677-', False)

    # TODO: Is this actually correct???
    def testORF1abInsertionsAa(self):
        """
        The ORF1ab protein should have the expected insertions.
        """
        self.check('orf1ab', '-4402F -4403K', False)

    def testORF1abMutationsAa(self):
        """
        The ORF1ab protein should have the expected amino acid changes.
        """
        self.check('orf1ab', 'T1001I A1708D I2230T P4717L', False)

    def testNucleocapsidMutationsNt(self):
        """
        The nucleocapsid genome should have the expected changes.
        """
        # The 704 below is due to a change at the 1-based amino acid
        # location 235. That's a 0-based offset of 234 = 702 in the
        # genome. The change is an S -> via a TCT -> TTT mutation in the
        # middle position, or 703 in 0-based, and 704 in 1-based.
        self.check('N', 'G7C A8T T9A G608A G609A G610C C704T', True)

    def testNucleocapsidMutationsAa(self):
        """
        The nucleocapsid protein should have the expected amino acid changes.
        Note that the UK report does not include mention of R203K or G204R.
        """
        self.check('N', 'D3L R203K G204R S235F', False)

    def testORF8MutationsAa(self):
        """
        The ORF8 protein should have the expected amino acid changes.
        """
        self.check('orf8', 'Q27* R52I Y73C', False)

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
        reference = self.genomeRead.sequence
        for location, nt in ((3267, 'T'), (5388, 'G'), (6954, 'C'),
                             (23063, 'T'), (23271, 'G'), (23604, 'T')):
            self.assertEqual(nt, reference[location - 1])


class Test_BavPat2(TestCase, _Mixin):
    """
    Test the BavPat2 sequence. This is Bavarian patient #2.
    """
    genomeRead = getSequence(join(DATA_DIR, 'BavPat2.fasta'))
    genome = SARS2Genome(genomeRead, FEATURES)

    def testSpikeMutationsNt(self):
        """
        The spike genome should have the expected change.
        """
        self.check('spike', 'A1841G', True)

    def testSpikeMutationsAa(self):
        """
        The spike protein should have the expected amino acid change.
        """
        self.check('spike', 'D614G', False)

    def testORF1aMutationsNt(self):
        """
        The ORF1a genome should have the expected change.
        """
        self.check('orf1a', 'C2772T', True)

    def testORF1abInsertionsAa(self):
        """
        The ORF1ab protein should have the expected insertions.
        """
        self.check('orf1ab', '-4402F -4403K', False)

    def testSpikeDeletionVariant(self):
        """
        The genome is not a spike deletion variant.
        """
        testCount, errorCount, result = self.genome.checkVariant(
            'spikeDeletion')
        self.assertEqual(2, testCount)
        self.assertEqual(2, errorCount)
        self.assertEqual(
            {
                'spike': {
                    'aa': {
                        '69-': (True, False),
                        '70-': (True, False),
                    },
                },
            },
            result
        )

    def testN501YVariant(self):
        """
        The genome must not be an N501Y variant.
        """
        testCount, errorCount, result = self.genome.checkVariant('N501Y')
        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual(
            {
                'spike': {
                    'aa': {
                        'N501Y': (True, False),
                    },
                },
            },
            result
        )

    def testVariantOfConcern20201201(self):
        """
        The genome must not have any of the UK variant of concern 202012/01
        changes.
        """
        testCount, errorCount, _ = self.genome.checkVariant('VOC_20201201_UK')
        self.assertEqual(20, testCount)
        self.assertEqual(20, errorCount)

    def testNucleocapsidIdentical(self):
        """
        The nucleocapsid genome should be identical to the reference.
        """
        sequenceNt, referenceNt = self.genome.feature('N').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

    def testORF8Identical(self):
        """
        The ORF8 genome should be identical to the reference.
        """
        sequenceNt, referenceNt = self.genome.feature('orf8').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

        sequenceAa, referenceAa = self.genome.feature('orf8').aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)

    def testEnvelopeIdentical(self):
        """
        The envelope should be identical to the reference.
        """
        sequenceNt, referenceNt = self.genome.feature('E').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

        sequenceAa, referenceAa = self.genome.feature('E').aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)

    def testMembraneIdentical(self):
        """
        The membrane should be identical to the reference.
        """
        sequenceNt, referenceNt = self.genome.feature('M').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

        sequenceAa, referenceAa = self.genome.feature('M').aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)

    def testRdRpIdentical(self):
        """
        The polymerase should be identical to the reference.
        """
        sequenceNt, referenceNt = self.genome.feature('rdrp').ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

        sequenceAa, referenceAa = self.genome.feature('rdrp').aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)


class Test_NC_045512(TestCase, _Mixin):
    """
    Test the NC_045512.2 sequence, which should test as equal seeing as it is
    the default feature reference.
    """
    genomeRead = getSequence(join(DATA_DIR, 'NC_045512.2.fasta'))
    genome = SARS2Genome(genomeRead, FEATURES)

    def testSpikeIdenticalNt(self):
        """
        The spike nucleotides should be identical.
        """
        spike = self.genome.feature('S')
        sequenceNt, referenceNt = spike.ntSequences()
        self.assertEqual(sequenceNt.sequence, referenceNt.sequence)

    def testSpikeIdentical(self):
        """
        The spike protein should be identical.
        """
        spike = self.genome.feature('S')
        sequenceAa, referenceAa = spike.aaSequences()
        self.assertEqual(sequenceAa.sequence, referenceAa.sequence)

    def testSpikeDeletionVariant(self):
        """
        The genome is not a spike deletion variant.
        """
        testCount, errorCount, result = self.genome.checkVariant(
            'spikeDeletion')
        self.assertEqual(2, testCount)
        self.assertEqual(2, errorCount)
        self.assertEqual(
            {
                'spike': {
                    'aa': {
                        '69-': (True, False),
                        '70-': (True, False),
                    },
                },
            },
            result
        )

    def testN501YVariant(self):
        """
        The genome must not be an N501Y variant.
        """
        testCount, errorCount, result = self.genome.checkVariant('N501Y')
        self.assertEqual(1, testCount)
        self.assertEqual(1, errorCount)
        self.assertEqual(
            {
                'spike': {
                    'aa': {
                        'N501Y': (True, False),
                    },
                },
            },
            result
        )

    def testVariantOfConcern20201201(self):
        """
        The genome must not have any of the UK variant of concern 202012/01
        changes.
        """
        testCount, errorCount, _ = self.genome.checkVariant('VOC_20201201_UK')
        self.assertEqual(20, testCount)
        self.assertEqual(20, errorCount)
