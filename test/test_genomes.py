"""
Tests for specific genomes in ../data

See test_genome.py for tests of more basic SARS2Genome functionality.
"""

from unittest import TestCase

from .fasta import getSequence

from sars2seq.features import Features, DATA_DIR
from sars2seq.genome import SARS2Genome
from sars2seq.variants import VARIANTS

REF_GB = DATA_DIR / 'NC_045512.2.gb'
FEATURES = Features(REF_GB)


class _Mixin:
    """
    Mixin for SARS2Genome class tests.
    """
    def testLength(self):
        self.assertGreater(len(self.genomeRead), 28000)

    def check(self, featureName, changes, nt):
        """
        Check that a set of changes all happened as expected.

        @param featureName: The C{str} name of the feature to check (e.g.,
            'nsp2').
        @param changes: A C{str} specification in the form of space-separated
            RNS strings, where R is a reference base, N is an integer offset,
            and S is a sequence base. So, e.g., 'L28S P1003Q' indicates that
            we expected a change from 'L' to 'S' at offset 28 and from 'P' to
            'Q' at offset 1003.
        @param nt: If C{True} check nucleotide sequences. Else protein.
        """
        _, errorCount, result = self.genome.checkFeature(
            featureName, changes, nt)
        if errorCount:
            for change, (referenceOK, _, genomeOK, _) in result.items():
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
            getSequence(DATA_DIR / 'EPI_ISL_402125.fasta').sequence,
            getSequence(DATA_DIR / 'NC_045512.2.fasta').sequence)


class Test_EPI_ISL_601443(TestCase, _Mixin):
    """
    Test the EPI_ISL_601433 sequence. This is the variant of concern
    (VOC 202012/01) referred to in https://www.gov.uk/government/publications/
    investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201
    """
    genomeRead = getSequence(DATA_DIR / 'EPI_ISL_601443.fasta')
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
                        '69-': (True, 'H', True, '-'),
                        '70-': (True, 'V', True, '-'),
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
                        'N501Y': (True, 'N', True, 'Y'),
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

    def testVariantOfConcern20201201ByDict(self):
        """
        The genome must fulfil all the requirements of the UK variant of
        concern 202012/01 when a dictionary describing the variant is passed.
        """
        changes = VARIANTS['VOC_20201201_UK']['changes']
        testCount, errorCount, _ = self.genome.checkVariant(changes)
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

    def testORF1abTranslationBug(self):
        """
        The ORF1ab protein should not have insertions at 4402 or 4403. This
        ensures that https://github.com/VirologyCharite/sars2seq/issues/9 is
        fixed and does not revert.
        """
        _, errorCount, _ = self.genome.checkFeature(
            'orf1ab', '-4402F -4403K', False)
        self.assertEqual(2, errorCount)

    def testORF1abDeletionsAa(self):
        """
        The ORF1ab protein should have the expected deletions.
        """
        self.check('orf1ab', 'S3675- G3676- F3677-', False)

    def testORF1abMutationsAa(self):
        """
        The ORF1ab protein should have the expected amino acid changes.
        """
        self.check('orf1ab', 'T1001I A1708D I2230T', False)

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
        # The locations checked here are not the same as those given in
        # Table 1 because the VOC starts at an offset of 54 bases into the
        # Wuhan reference and there are deletions (rising to a total of 18
        # in the higher offsets) that have to be adjusted for.
        reference = self.genomeRead.sequence
        for location, nt in (
                (3267, 'T'), (5388, 'A'), (6954, 'C'), (23045, 'T'),
                (23253, 'A'), (23686, 'A'), (23691, 'T'), (24488, 'G'),
                (24896, 'C'), (27954, 'T'), (28030, 'T'), (28093, 'G'),
                (28262, 'C'), (28263, 'T'), (28264, 'A'), (28959, 'T')):
            self.assertEqual(nt, reference[location - 54 - 1],
                             f'Failed on {location}{nt}, saw '
                             f'{self.genomeRead.sequence[location - 54 - 1]}')


class Test_BavPat2(TestCase, _Mixin):
    """
    Test the BavPat2 sequence. This is Bavarian patient #2.
    """
    genomeRead = getSequence(DATA_DIR / 'BavPat2.fasta')
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
                        '69-': (True, 'H', False, 'H'),
                        '70-': (True, 'V', False, 'V'),
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
                        'N501Y': (True, 'N', False, 'N'),
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
        genomeNt, referenceNt = self.genome.ntSequences('N')
        self.assertEqual(genomeNt.sequence, referenceNt.sequence)

    def testORF8Identical(self):
        """
        The ORF8 genome should be identical to the reference.
        """
        genomeNt, referenceNt = self.genome.ntSequences('orf8')
        self.assertEqual(genomeNt.sequence, referenceNt.sequence)

        genomeAa, referenceAa = self.genome.aaSequences('orf8')
        self.assertEqual(genomeAa.sequence, referenceAa.sequence)

    def testEnvelopeIdentical(self):
        """
        The envelope should be identical to the reference.
        """
        genomeNt, referenceNt = self.genome.ntSequences('E')
        self.assertEqual(genomeNt.sequence, referenceNt.sequence)

        genomeAa, referenceAa = self.genome.aaSequences('E')
        self.assertEqual(genomeAa.sequence, referenceAa.sequence)

    def testMembraneIdentical(self):
        """
        The membrane should be identical to the reference.
        """
        genomeNt, referenceNt = self.genome.ntSequences('M')
        self.assertEqual(genomeNt.sequence, referenceNt.sequence)

        genomeAa, referenceAa = self.genome.aaSequences('M')
        self.assertEqual(genomeAa.sequence, referenceAa.sequence)

    def testRdRpIdentical(self):
        """
        The polymerase should be identical to the reference.
        """
        genomeNt, referenceNt = self.genome.ntSequences('rdrp')
        self.assertEqual(genomeNt.sequence, referenceNt.sequence)

        genomeAa, referenceAa = self.genome.aaSequences('rdrp')
        self.assertEqual(genomeAa.sequence, referenceAa.sequence)


class Test_NC_045512(TestCase, _Mixin):
    """
    Test the NC_045512.2 sequence, which should test as equal seeing as it is
    the default feature reference.
    """
    genomeRead = getSequence(DATA_DIR / 'NC_045512.2.fasta')
    genome = SARS2Genome(genomeRead, FEATURES)

    def testSpikeIdenticalNt(self):
        """
        The spike nucleotides should be identical.
        """
        genomeNt, referenceNt = self.genome.ntSequences('S')
        self.assertEqual(genomeNt.sequence, referenceNt.sequence)

    def testSpikeIdentical(self):
        """
        The spike protein should be identical.
        """
        genomeAa, referenceAa = self.genome.aaSequences('S')
        self.assertEqual(genomeAa.sequence, referenceAa.sequence)

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
                        '69-': (True, 'H', False, 'H'),
                        '70-': (True, 'V', False, 'V'),
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
                        'N501Y': (True, 'N', False, 'N'),
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


class Test_EPI_ISL_678597(TestCase, _Mixin):
    """
    Test the EPI_ISL_678597 sequence. This is the South African variant of
    concern.
    """
    genomeRead = getSequence(DATA_DIR / 'EPI_ISL_678597.fasta')
    genome = SARS2Genome(genomeRead, FEATURES)

    def testSpikeDeletionsAa(self):
        """
        The spike protein should have the three deletions.
        """
        self.check('spike', 'L241- L242- A243-', nt=False)

    def testSpikeMutationsAa(self):
        """
        The spike protein should have the expected amino acid changes.
        """
        self.check('spike', 'D80A D215G K417N E484K N501Y D614G A701V',
                   False)

    def testSpikeDeletionVariant(self):
        """
        The genome is not a spike deletion variant.
        """
        testCount, errorCount, result = self.genome.checkVariant(
            'spikeDeletion')
        self.assertEqual(2, errorCount)

    def testN501YVariant(self):
        """
        The genome is an N501Y variant.
        """
        testCount, errorCount, result = self.genome.checkVariant('N501Y')
        self.assertEqual(1, testCount)
        self.assertEqual(0, errorCount)
        self.assertEqual(
            {
                'spike': {
                    'aa': {
                        'N501Y': (True, 'N', True, 'Y'),
                    },
                },
            },
            result
        )

    def test501YV2Variant(self):
        """
        The genome must be a 501Y.V2 variant.
        """
        testCount, errorCount, _ = self.genome.checkVariant('501Y.V2')
        self.assertEqual(8, testCount)
        self.assertEqual(0, errorCount)

    def testNotVariantOfConcern20201201(self):
        """
        The genome is not a UK variant of concern 202012/01.
        """
        testCount, errorCount, _ = self.genome.checkVariant('VOC_20201201_UK')
        self.assertEqual(20, testCount)
        self.assertEqual(16, errorCount)
        self.check('spike', 'N501Y', False)
        self.check('orf1ab', 'S3675- G3676- F3677-', False)

    def testORF1aDeletionsAa(self):
        """
        The ORF1a protein should have the expected deletions.
        """
        self.check('orf1a', 'S3675- G3676- F3677-', False)

    def testORF1aMutationsAa(self):
        """
        The ORF1a protein should have the expected amino acid changes.
        """
        self.check('orf1a', 'T265I K1655N K3353R', False)

    def testNucleocapsidMutationsNt(self):
        """
        The nucleocapsid genome should have the expected changes.
        """
        self.check('N', 'C614T', True)

    def testNucleocapsidMutationsAa(self):
        """
        The nucleocapsid protein should have the expected amino acid changes.
        """
        self.check('N', 'T205I', False)

    def testORF3aMutationsAa(self):
        """
        The ORF3a protein should have the expected amino acid changes.
        """
        self.check('orf3a', 'Q57H S171L', False)
