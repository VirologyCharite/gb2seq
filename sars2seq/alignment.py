from Bio.Seq import Seq

from dark.aligners import mafft
from dark.reads import AARead, DNARead, Reads


DEBUG = False
SLICE = slice(300)


def translate(seq, name):
    """
    Translate a sequence.

    @param seq: A C{str} nucelotide sequence.
    @param name: A C{str} feature name.
    @return: A translated C{str} amino acid sequence.
    """
    if name == 'ORF1ab polyprotein':
        # See Fields Virology (figure 10.6a on page 421, 7th edition or
        # figure 28.7a on page 836, 6th edition) plus
        # https://www.ncbi.nlm.nih.gov/nuccore/NC_045512 for details of
        # what happens below. Note that the nucelotide sequence we are
        # passed is the one that's made from the alignment with the
        # reference ORF1ab nucleotide sequence (in sequence.py) and so is
        # just that ORF and does not include the leading ~265 nucleotides
        # of the 5' UTR. As a result, the offset used to begin the search
        # for the slippery sequence is 13000, which is chosen to be a bit
        # before 13468 - 265.  There are various occurrences of the
        # slippery sequence in the reference genome (and hence probably in
        # other CoV genomes), but only one in this region and with a stop
        # codon shortly (<20 nt) downstream.
        slipperySeq = 'TTTAAAC'
        slipperyLen = len(slipperySeq)
        offset = seq.find(slipperySeq, 13000)
        stop = seq.find('TAA', offset + slipperyLen)
        if DEBUG:
            print(f'LEN: {len(seq)}, OFFSET: {offset}, STOP: {stop}')
        assert offset > -1 and stop > -1 and stop - offset < 20
        seq = seq[:offset + slipperyLen] + seq[offset:]

    # Pad with 'N' to avoid a 'BiopythonWarning: Partial codon' warning.
    remainder = len(seq) % 3
    seq += 'N' * (3 - remainder if remainder else 0)

    return Seq(seq).translate()


class Alignment:
    """
    Hold information about a matched feature.

    @param genome: A C{dark.reads.DNARead} instance holding the SARS-CoV-2
        genome of interest.
    @param referenceId: A C{str} with the id of the reference.
    @param feature: A C{dict} with information about the reference feature.
    @param window: The C{int} number of nucleotides to include in the region
        of the genome that is searched for the reference sequence. This is
        in addition to the range of the feature in the reference. So e.g., if
        the feature is in positions 1000-1300 in the reference and the window
        is 500, then the alignment call that tries to find the reference in the
        genome will look in the region 500-1800 of the genome. This is to
        allow for indels in the genome that cause its feature to be located
        at quite different offsets than the location of the feature in the
        reference.  It greatly speeds up the alignment and also prevents cases
        where a short prefix of the reference matches very early in the
        genome, followed by a huge gap and then the rest of the reference
        matching the genome. The default value below was arrived at by
        experimentation.
    """
    def __init__(self, genome, referenceId, feature, window=500):
        self.genome = genome
        self.referenceId = referenceId
        self.feature = feature
        self.window = window
        self.genomeOffset = None
        self.genomeNt = self.genomeAa = None
        self.referenceNt = self.referenceAa = None

    def ntSequences(self):
        """
        Get the aligned nucelotide sequences.

        @return: A 2-C{tuple} of C{dark.reads.DNARead} instances, holding
            the nucleotides for the feature as located in the genome being
            examined and the corresponding nucleotides from the reference
            genome.
        """
        if self.genomeNt:
            assert self.referenceNt
            return self.genomeNt, self.referenceNt

        name = self.feature['name']

        subsequence = self.genome.sequence[
            max(0, self.feature['start'] - self.window):
            self.feature['stop'] + self.window]

        genomeRead = DNARead(
            self.genome.id + f' ({name})', subsequence)
        referenceRead = DNARead(
            self.referenceId + f' ({name})', self.feature['sequence'])

        assert len(genomeRead) >= len(referenceRead)

        if DEBUG:
            print('NT MATCH:')
            print('gen  nt:', genomeRead.sequence[SLICE])
            print('ref  nt:', referenceRead.sequence[SLICE])

        self.genomeNt, self.referenceNt = self._align(genomeRead,
                                                      referenceRead)

        # Figure out where the found feature falls in the genome. This
        # assumes the first occurrence of the feature string in the genome
        # is the right one.
        self.genomeOffset = self.genome.sequence.find(
            self.genomeNt.sequence.replace('-', ''))
        assert self.genomeOffset > -1

        return self.genomeNt, self.referenceNt

    def aaSequences(self):
        """
        Match the genome and the reference at the amino acid level.

        @return: A 2-C{tuple} of C{dark.reads.AARead} instances, holding
            the amino acids of the genome that is being examined and those
            of the reference sequence.
        """
        if self.genomeAa:
            assert self.referenceAa
            return self.genomeAa, self.referenceAa

        # We can't do the protein work unless we've already figured out
        # where the feature is at the nucelotide level.
        if not self.genomeNt:
            self.genomeNt, self.referenceNt = self.ntSequences()

        name = self.feature['name']

        genomeRead = AARead(
            self.genome.id + f' ({name})',
            translate(self.genomeNt.sequence.replace('-', ''), name))

        referenceTranslation = self.feature.get(
            'translation', translate(self.feature['sequence'], name))
        referenceRead = AARead(self.referenceId + f' ({name})',
                               referenceTranslation)

        if DEBUG:
            print(f'AA MATCH {name}:')
            print('gen  nt:', self.genome.sequence[SLICE])
            print('ref  nt:', self.reference.sequence[SLICE])
            print('gen  aa:', genomeRead.sequence[SLICE])
            print('ref  aa:', referenceRead.sequence[SLICE])

        self.genomeAa, self.referenceAa = self._align(
            genomeRead, referenceRead, nt=False)

        return self.genomeAa, self.referenceAa

    def _align(self, genomeRead, referenceRead, nt=True):
        """
        Align two sequences to find where the reference read fits in the
        genome read.

        @param genomeRead: A C{dark.reads.Read} instance.
        @param referenceRead: A C{dark.reads.Read} instance.
        @param nt: If C{True} the sequences are nucleotide. Else protein.
        @return: A 2-C{tuple} of aligned C{dark.reads.Read} instances, with the
            genome and the reference sequences. Both may contain gaps ('-').
        """
        alignment = mafft(
            Reads([genomeRead, referenceRead]),
            options='--anysymbol --preservecase' + (
                ' --nuc' if nt else ' --amino'))

        genomeResult, referenceResult = list(alignment)
        if DEBUG:
            print('process alignment')

        # One but not both of the genome and reference can begin with a
        # gap, and the same goes for the sequence ends. The reason is that
        # an aligner has no reason to put gaps at the extremes of both
        # sequences at once.
        assert not (genomeResult.sequence.startswith('-') and
                    referenceResult.sequence.startswith('-'))
        assert not (genomeResult.sequence.endswith('-') and
                    referenceResult.sequence.endswith('-'))

        # In the alignment, the reference may have leading and/or trailing '-'
        # chars. The length of these gives us the offsets in the original
        # full-length genome where the match is.

        if DEBUG:
            print('gen   al', genomeResult.sequence[SLICE])
            print('ref   al', referenceResult.sequence[SLICE])
            with open('/tmp/align-%s.fasta' %
                      ('nt' if nt else 'aa'), 'w') as fp:
                print('seq', genomeResult.sequence, file=fp)
                print('res', referenceResult.sequence, file=fp)

        if referenceResult.sequence.startswith('-'):
            offset = (len(referenceResult) -
                      len(referenceResult.sequence.lstrip('-')))
            genomeResult = genomeResult[offset:]
            referenceResult = referenceResult[offset:]

            if DEBUG:
                print(f'CLIPPED {offset} on left')
                print('gen   al', genomeResult.sequence[SLICE])
                print('ref   al', referenceResult.sequence[SLICE])

        if referenceResult.sequence.endswith('-'):
            offset = (len(referenceResult) -
                      len(referenceResult.sequence.rstrip('-')))

            if DEBUG:
                print(f'CLIPPING {offset} on right')
                print('gen   al', genomeResult.sequence)
                print('ref   al', referenceResult.sequence)

            genomeResult = genomeResult[:len(referenceResult) - offset]
            referenceResult = referenceResult[:len(referenceResult) - offset]

            if DEBUG:
                print(f'CLIPPED {offset} on right')
                print('gen   al', genomeResult.sequence[SLICE])
                print('ref   al', referenceResult.sequence[SLICE])

        if DEBUG:
            print('gen  end', genomeResult.sequence[SLICE])
            print('ref  end', referenceResult.sequence[SLICE])
            print('END process alignment')

        return genomeResult, referenceResult
