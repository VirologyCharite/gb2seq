from dark.aligners import mafft
from dark.reads import AARead, DNARead, Reads

from sars2seq.change import splitChange
from sars2seq.translate import translate

DEBUG = False
SLICE = slice(300)


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

    def check(self, changes, nt):
        """
        Check that a set of changes all happened as expected.

        @param changes: Either a C{str} or a iterable of 3-C{tuple}s. If a
            C{str}, gives a specification in the form of space-separated
            RNG strings, where R is a reference base, N is a 1-based location,
            and G is a sequence base. So, e.g., 'L28S P1003Q' indicates that
            we expected a change from 'L' to 'S' at offset 28 and from 'P' in
            the reference to 'Q' in the genome we're examining at offset 1003.
            The reference or genome base (but not both) may be absent.
            If an iterable of 3-C{tuple}s, each tuple should have a C{str}
            expected reference base, a 0-based offset, and a C{str} expected
            genome base.  Note that the string format (meant for humans) uses
            1-based locations whereas the tuple format uses 0-based offsets.
        @param nt: If C{True} check nucleotide sequences. Else protein.
        @raise ValueError: If a change string cannot be parsed.
        @return: A 3-C{tuple} with the number of checks done, the number of
            errors, and a C{dict} keyed by changes in C{changes}, with values
            a 2-C{tuple} of Booleans to indicate success or failure of the
            check for the reference and the genome respectively.
        """
        genome, reference = self.ntSequences() if nt else self.aaSequences()
        genome, reference = genome.sequence, reference.sequence
        result = {}
        testCount = errorCount = 0

        if isinstance(changes, str):
            for change in changes.split():
                refBase, offset, genBase = splitChange(change)
                result[change] = (
                    refBase is None or reference[offset] == refBase,
                    genBase is None or genome[offset] == genBase)
        else:
            for change in changes:
                refBase, offset, genBase = change
                result[change] = (reference[offset] == refBase,
                                  genome[offset] == genBase)

        errorCount = len([v for v in result.values() if v != (True, True)])
        testCount = len(result)

        return testCount, errorCount, result
