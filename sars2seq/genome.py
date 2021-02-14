from dark.aligners import mafft
from dark.reads import AARead, DNARead, Reads

from sars2seq.change import splitChange
from sars2seq.translate import translate
from sars2seq.variants import VARIANTS

DEBUG = False
SLICE = slice(300)


def getNonGapOffsets(s):
    """
    Make a dictionary mapping non-gap offsets to the equivalent offset in a
    gapped sequence.

    @param s: A C{str} sequence, possibly with gap ('-') characters.
    @return: A C{dict} mapping C{int} offsets to equivalent C{int} offsets
        in C{s}, where gaps are ignored.
    """
    result = {}
    gapCount = index = 0

    for base in s:
        if base == '-':
            gapCount += 1
        else:
            result[index] = index + gapCount
            index += 1

    assert set(range(len(s) - s.count('-'))) == set(result)

    return result


def alignmentEnd(s, startOffset, length):
    """
    @param s: A C{str} sequence, possibly with gap ('-') characters.
    """
    found = 0
    index = startOffset

    while found < length:
        if s[index] != '-':
            found += 1
        index += 1

    return index


class SARS2Genome:
    """
    Methods for working with SARS-CoV-2 genomes.

    @param genome: A C{dark.reads.Read} instance.
    @param features: An C{Features} instance.
    @raise ReferenceWithGapsError: If the reference has gaps.
    """
    def __init__(self, genome, features):
        # Geneious uses ? to indicate unknown nucleotides, at least when it
        # exports a consensus / alignment. Replace with N.
        self.genome = DNARead(genome.id, genome.sequence.replace('?', 'N'))
        self.features = features
        self.getAlignment()
        self._cache = {'aa': {}, 'nt': {}}

    def getAlignment(self):
        """
        Align the reference and the genome.
        """
        self.genomeAligned, self.referenceAligned = mafft(
            Reads([self.genome, self.features.reference]),
            options='--anysymbol --preservecase --nuc')

        # One but not both of the aligned genome and reference can begin
        # with a gap, and the same goes for the sequence ends. The reason
        # is that an aligner has no reason to put gaps at the extremes of
        # both sequences at once.
        assert not (self.genomeAligned.sequence.startswith('-') and
                    self.referenceAligned.sequence.startswith('-'))
        assert not (self.genomeAligned.sequence.endswith('-') and
                    self.referenceAligned.sequence.endswith('-'))

        self.offsetMap = getNonGapOffsets(self.referenceAligned.sequence)

    def ntSequences(self, featureName):
        """
        Get the aligned nucelotide sequences.

        @param featureName: A C{str} feature name.
        @return: A 2-C{tuple} of C{dark.reads.DNARead} instances, holding
            the nucleotides for the feature as located in the genome being
            examined and the corresponding nucleotides from the reference
            genome.
        """
        try:
            return self._cache['nt'][featureName]
        except KeyError:
            pass

        feature = self.features[featureName]
        name = feature['name']
        length = feature['stop'] - feature['start']
        offset = self.offsetMap[feature['start']]
        end = alignmentEnd(self.referenceAligned.sequence, offset, length)

        genomeRead = DNARead(self.genome.id + f' ({name})',
                             self.genomeAligned.sequence[offset:end])

        referenceRead = DNARead(self.features.reference.id + f' ({name})',
                                self.referenceAligned.sequence[offset:end])

        if DEBUG:
            print('NT MATCH:')
            print('gen  nt:', genomeRead.sequence[SLICE])
            print('ref  nt:', referenceRead.sequence[SLICE])

        self._cache['nt'][featureName] = genomeRead, referenceRead
        return genomeRead, referenceRead

    def aaSequences(self, featureName):
        """
        Match the genome and the reference at the amino acid level.

        @param featureName: A C{str} feature name.
        @return: A 2-C{tuple} of C{dark.reads.AARead} instances, holding
            the amino acids of the genome that is being examined and those
            of the reference sequence.
        """
        try:
            return self._cache['aa'][featureName]
        except KeyError:
            pass

        genomeRead, referenceRead = self.ntSequences(featureName)

        feature = self.features[featureName]
        name = feature['name']

        genomeRead = AARead(
            self.genome.id + f' ({name})',
            translate(genomeRead.sequence.replace('-', ''), name))

        referenceTranslation = feature.get(
            'translation', translate(feature['sequence'], name))
        referenceRead = AARead(self.features.reference.id + f' ({name})',
                               referenceTranslation)

        if DEBUG:
            print(f'AA MATCH {name}:')
            print('gen  nt:', self.genome.sequence[SLICE])
            print('ref  nt:', self.features.reference.sequence[SLICE])
            print('gen  aa:', genomeRead.sequence[SLICE])
            print('ref  aa:', referenceRead.sequence[SLICE])

        genomeAa, referenceAa = mafft(
            Reads([genomeRead, referenceRead]),
            options='--anysymbol --preservecase --amino')

        self._cache['aa'][featureName] = genomeAa, referenceAa
        return genomeAa, referenceAa

    def _checkChange(self, base, offset, read, change, featureName):
        """
        Check that a base occurs at an offset.

        @param base: A C{str} nucleotide base (or amino acid), or C{None}
            in which case there was no expectation that the read had any
            particular value and C{True} is returned.
        @param offset: A 0-based C{int} offset into C{sequence}.
        @param read: A C{dark.reads.Read} instance.
        @param change: The C{str} or C{tuple} change specification.
        @param featureName: A C{str} feature name.
        @raise IndexError: If the offset is out of range.
        @return: A C{list} containing a C{bool} indicating whether the read
            sequence has C{base} at C{offset} (or C{True} if C{base} is
            C{None}) and the base found at the offset.
        """
        try:
            actual = read.sequence[offset]
        except IndexError:
            raise IndexError(f'Index {offset} out of range trying to access '
                             f'feature {featureName!r} of length '
                             f'{len(read)} sequence {read.id!r} via '
                             f'expected change specification {change!r}.')
        else:
            return [(base is None or actual == base), actual]

    def checkFeature(self, featureName, changes, nt):
        """
        Check that a set of changes all happened as expected.

        @param featureName: A C{str} feature name.
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
        @raise IndexError: If the offset of a change exceeds the length of the
            sequence being checked.
        @return: A 3-C{tuple} with the number of checks done, the number of
            errors, and a C{dict} keyed by changes in C{changes}, with values
            a 2-C{tuple} of Booleans to indicate success or failure of the
            check for the reference and the genome respectively.
        """
        genome, reference = (self.ntSequences(featureName) if nt else
                             self.aaSequences(featureName))
        result = {}
        testCount = errorCount = 0

        if isinstance(changes, str):
            for change in changes.split():
                referenceBase, offset, genomeBase = splitChange(change)
                result[change] = tuple(
                    self._checkChange(referenceBase, offset, reference, change,
                                      featureName) +
                    self._checkChange(genomeBase, offset, genome, change,
                                      featureName))
        else:
            for change in changes:
                referenceBase, offset, genomeBase = change
                result[change] = tuple(
                    self._checkChange(referenceBase, offset, reference, change,
                                      featureName) +
                    self._checkChange(genomeBase, offset, genome, change,
                                      featureName))

        errorCount = sum((v[0] is False or v[2] is False)
                         for v in result.values())
        testCount = len(result)

        return testCount, errorCount, result

    def checkVariant(self, variant):
        """
        Check that a set of changes in different features all happened as
        expected.

        @param variant: The C{str} name of a key in the VARIANTS C{dict}.
        @return: A 3-C{tuple} with the number of checks done, the number of
            errors, and a C{dict} keyed by changes in C{changes}, with values
            a 2-C{tuple} of Booleans to indicate success or failure of the
            check for the reference and the genome respectively.
        """
        result = {}
        testCountTotal = errorCountTotal = 0
        changeDict = VARIANTS[variant]['changes']

        for featureName in changeDict:
            result[featureName] = {}
            for what in 'aa', 'nt':
                try:
                    changes = changeDict[featureName][what]
                except KeyError:
                    pass
                else:
                    result[featureName][what] = {}
                    testCount, errorCount, result_ = self.checkFeature(
                        featureName, changes, what == 'nt')
                    testCountTotal += testCount
                    errorCountTotal += errorCount
                    for change, values in result_.items():
                        result[featureName][what][change] = values

        return testCountTotal, errorCountTotal, result
