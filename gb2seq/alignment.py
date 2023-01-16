from typing import Dict, Iterable, Iterator, Optional, TextIO, Tuple, Union

import argparse
from Bio.Seq import Seq

from dark.aligners import edlibAlign, mafft
from dark.reads import AARead, DNARead, Reads

from gb2seq import Gb2SeqError
from gb2seq.change import splitChange
from gb2seq.features import Features
from gb2seq.translate import (
    translate,
    translateSARS2Spike,
    TranslationError,
    TranslatedReferenceAndGenomeLengthError,
)
from gb2seq.variants import VARIANTS

DEBUG = False
SLICE = slice(300)

MAFFT_OPTIONS = "--anysymbol --preservecase --retree 1 --reorder"

# Supported alignment methods.
ALIGNERS = ("edlib", "mafft")
DEFAULT_ALIGNER = "mafft"


class ReferenceInsertionError(Gb2SeqError):
    "A genome resulted in MAFFT suggesting a reference insertion."


class AlignmentError(Gb2SeqError):
    "There is an unexpected problem in the alignment."


def addAlignerOption(parser: argparse.ArgumentParser) -> None:
    """
    Add a command line option for specifying an aligner.

    @param parser: An argparse argument parser.
    """
    parser.add_argument(
        "--aligner",
        default=DEFAULT_ALIGNER,
        choices=ALIGNERS,
        help="The alignment method.",
    )


def getGappedOffsets(s: str) -> dict:
    """
    Make a dictionary mapping offsets in a sequence with no gaps to the
    equivalent offset in a gapped sequence.

    In more detail, we are given C{s}, which is an aligned sequence
    (potentially) with gaps in it and we want a mapping that will allow us to
    go from an offset in the original string C{s} (i.e., before it was aligned)
    to the equivalent offset in the aligned (with gaps) string. This function
    returns such a mapping (a C{dict})

    @param s: A C{str} sequence, possibly with gap ('-') characters.
    @return: A C{dict} mapping C{int} offsets (in the sequence without its
        gaps) to equivalent C{int} offsets in C{s}, when gaps are ignored.

    """
    result = {}
    gapCount = index = 0

    # To build the mapping, we walk through s. If we see a gap we just
    # count it. Non-gap offsets go into the map with their offset plus the
    # number of gaps that have been encountered before them.
    for base in s:
        if base == "-":
            gapCount += 1
        else:
            result[index] = index + gapCount
            index += 1

    # Sanity check that the result mapping has an entry for all the offsets
    # in the original string (before it was aligned).  This check can be
    # removed one day, seeing as there are tests for this function.
    assert set(result) == set(range(len(s) - s.count("-")))

    # If we have a non-empty result, add a final offset corresponding to
    # the length of the original string. This allows us to use Python
    # indexing in the case wanting to convert an offset that is one beyond
    # the end of the string.
    if result:
        assert index == len(s) - s.count("-")
        result[index] = index + gapCount

    return result


def alignmentEnd(s: str, startOffset: int, length: int) -> int:
    """
    Find the offset where an aligned sequence (i.e., potentially with gaps in
    it) of a given length ends.

    In other words, starting at C{startOffset}, look through C{s} until we have
    seen C{length} non-gap chars and return the corresponding index in C{s}.

    @param s: A C{str} sequence, possibly with gap ('-') characters.
    @param startOffset: The C{int} offset to start scanning from.
    @param length: The C{int} length of the (ungapped) original sequence
        that was aligned (to produce C{s}, with the alignment gaps in it).
    @return: The C{int} offset in C{s} (the alignment) where the original
        sequence ends, taking the gaps in C{s} into account.
    """
    nonGapCount = 0
    index = startOffset

    while nonGapCount < length:
        nonGapCount += s[index] != "-"
        index += 1

    return index


class Gb2Alignment:
    """
    Methods for working with a genome aligned to a GenBank reference.

    @param genome: A C{dark.reads.Read} instance.
    @param features: An C{Features} instance.
    @param referenceAligned: A C{dark.reads.Read} instance with an aligned
        reference sequence, or C{None} if the alignment should be done here.
        If not C{None} then C{genomeAligned} must also be given.
    @param genomeAligned: A C{dark.reads.Read} instance with an aligned
        genome sequence, or C{None} if the alignment should be done here.
        If not C{None} then C{referenceAligned} must also be given.
    @param aligner: A C{str} specifying the alignment algorithm to use.
        Either 'mafft' (the default) or 'edlib' (experimental).
    @param matchAmbiguous: If C{True}, and the edlib aligner is in use, count
        ambiguous nucleotides that are possibly correct as actually being
        correct. Otherwise, the edlib aligner is strict and nucleotide codes
        only match themselves.
    @param untranslatable: A C{dict} with C{str} keys and values. If any of
        the keys appears in a codon, the corresponding value is added to the
        translation. This can be used e.g., to make occurrences of '?'
        translate into '-' or 'X'.
    @raise AlignmentError: If the genome and reference sequences both start
        or end with '-' characters, or if pre-aligned sequences are not of
        the same length.
    @raise ReferenceWithGapsError: If the reference has gaps.
    """

    def __init__(
        self,
        genome: DNARead,
        features: Features,
        referenceAligned: Optional[DNARead] = None,
        genomeAligned: Optional[DNARead] = None,
        aligner: str = DEFAULT_ALIGNER,
        matchAmbiguous: bool = True,
        untranslatable: Optional[Dict[str, str]] = None,
    ):
        if aligner not in ALIGNERS:
            raise ValueError(f"Unknown aligner {aligner!r}.")

        self.genome = genome
        self.features = Features() if features is None else features
        self._matchAmbiguous = matchAmbiguous
        self._aligner = aligner
        self._untranslatable = untranslatable
        self._getAlignment(referenceAligned, genomeAligned)
        self.gappedOffsets = getGappedOffsets(self.referenceAligned.sequence)
        self._cache: dict = {"aa": {}, "nt": {}}

    def _getAlignment(
        self,
        referenceAligned: Optional[DNARead] = None,
        genomeAligned: Optional[DNARead] = None,
    ) -> None:
        """
        Align the reference and the genome.

        @param referenceAligned: A C{dark.reads.Read} instance with an aligned
            reference sequence, or C{None} if the alignment should be done here.
            If not C{None} then C{genomeAligned} must also be given.
        @param genomeAligned: A C{dark.reads.Read} instance with an aligned
            genome sequence, or C{None} if the alignment should be done here.
            If not C{None} then C{referenceAligned} must also be given.
        @raise ValueError: If an already aligned genome is given but no
            aligned reference is also given, or vice versa.
        @raise AlignmentError: If the genome and reference sequences both start
            or end with '-' characters, or if pre-aligned sequences are not of
            the same length.
        """
        if bool(referenceAligned) + bool(genomeAligned) == 1:
            raise ValueError(
                "Either both or neither of referenceAligned and "
                "genomeAligned can be given, not a mix."
            )

        preAligned = referenceAligned is not None

        if not preAligned:
            if self.features.reference.sequence == self.genome.sequence:
                # No need to align if the sequences are identical.
                referenceAligned = self.features.reference
                genomeAligned = self.genome
            else:
                reads = Reads([self.features.reference, self.genome])
                if self._aligner == "mafft":
                    referenceAligned, genomeAligned = mafft(
                        reads, options=MAFFT_OPTIONS
                    )
                else:
                    assert self._aligner == "edlib"
                    referenceAligned, genomeAligned = edlibAlign(
                        reads, matchAmbiguous=self._matchAmbiguous
                    )

            if DEBUG:
                print("ALIGNING")
                print(
                    f"ref {self.features.reference.id:10s}: "
                    f"{self.features.reference.sequence}"
                )
                print(f"gen {self.genome.id:10s}: {self.genome.sequence}")
                print("ALIGNED")
                print(f"ref {referenceAligned.id:10s}: " f"{referenceAligned.sequence}")
                print(f"gen {genomeAligned.id:10s}: " f"{genomeAligned.sequence}")

        if len(referenceAligned) != len(genomeAligned):
            if preAligned:
                error = (
                    f"The length of the given pre-aligned reference "
                    f"({len(referenceAligned)}) does not match the "
                    f"length of the given pre-aligned genome "
                    f"({len(genomeAligned)})."
                )
            else:
                error = (
                    f"The length of the aligned reference produced "
                    f"by {self._aligner} ({len(referenceAligned)}) does "
                    f"not match the length of the aligned genome "
                    f"({len(genomeAligned)})."
                )
            raise AlignmentError(error)

        # One but not both of the aligned genome and reference can begin
        # with a gap, and the same goes for the sequence ends. The reason
        # is that an aligner has no reason to put gaps at the extremes of
        # both sequences at once.  Note that this could happen when making
        # a multiple sequence alignment involving more than just these two
        # sequences, but we're not (yet?) trying to handle that case
        # (wherein, I think, it would be necessary to walk both sequences
        # and elide all sites where both have a gap, but I'm not 100% sure
        # that doing that would be optimal since the alignment between the
        # two sequences of specific interest would have been calculated in
        # the presence of other sequences, so may not be pairwise optimal).
        if referenceAligned.sequence.startswith(
            "-"
        ) and genomeAligned.sequence.startswith("-"):
            raise AlignmentError(
                "The reference and genome alignment sequences "
                'both start with a "-" character.'
            )

        if referenceAligned.sequence.endswith("-") and genomeAligned.sequence.endswith(
            "-"
        ):
            raise AlignmentError(
                "The reference and genome alignment sequences "
                'both end with a "-" character.'
            )

        self.referenceAligned: DNARead = referenceAligned
        self.genomeAligned: DNARead = genomeAligned

    def ntSequences(
        self,
        featureName: str,
        reverseComplement: bool = False,
        addFeatureToId: bool = True,
        raiseOnReferenceGaps: bool = True,
    ) -> Tuple[DNARead, DNARead]:
        """
        Get the aligned nucelotide sequences.

        @param featureName: A C{str} feature name.
        @param reverseComplement: If C{True} and the feature does not have a
            forward orientation in the genome, its reverse complement will be
            returned.
        @param addFeatureToId: If C{True}, the feature name (in parentheses,
            preceded by a space) will be appended to the id of the returned
            sequences.
        @param raiseOnReferenceGaps: A C{bool}. If C{True} and the aligner
            suggests a gap in the reference, raise a ReferenceInsertionError.
        @raise AssertionError: If the aligned reference and genome sequences
            are not the same length.
        @raise ReferenceInsertionError: (as above for raiseOnReferenceGaps).
        @return: A 2-C{tuple} of C{dark.reads.DNARead} instances, with 1) the
            nucleotides for the feature as located in the reference genome and
            2) the corresponding nucleotides from the genome being examined.
        """
        cacheKey = featureName, reverseComplement, addFeatureToId

        try:
            return self._cache["nt"][cacheKey]
        except KeyError:
            pass

        feature = self.features[featureName]
        name = feature["name"]
        length = feature["stop"] - feature["start"]
        # 'offset' is the offset in the (possibly gapped) alignment.
        offset = self.gappedOffsets[feature["start"]]
        end = alignmentEnd(self.referenceAligned.sequence, offset, length)
        idSuffix = f" ({name})" if addFeatureToId else ""

        referenceNt = DNARead(
            self.features.reference.id + idSuffix,
            self.referenceAligned.sequence[offset:end],
        )

        genomeNt = DNARead(
            self.genome.id + idSuffix, self.genomeAligned.sequence[offset:end]
        )

        # Use True as a default for 'forward' in the following, as some
        # features may be unannotated regions, and these do not have a
        # 'forward' key.
        if reverseComplement and not feature.get("forward", True):
            referenceNt = referenceNt.reverseComplement()
            genomeNt = genomeNt.reverseComplement()

        assert len(referenceNt) == len(genomeNt)

        # In general, there should not be insertions to the reference.
        # There are SARS-CoV-2 lineages with insertions in the Spike
        # (e.g. B.1.214.2) that we can correct for in the downstream
        # processing, therefore the error is not raised for the SARS-CoV-2
        # Spike.
        if (
            raiseOnReferenceGaps
            and "-" in referenceNt.sequence
            and (not self.features.sars2 or name != "surface glycoprotein")
        ):
            raise ReferenceInsertionError(
                f"{self._aligner} aligner suggested a reference insertion "
                f"into {featureName!r}."
            )

        if DEBUG:
            print("NT MATCH:")
            print("ref  nt:", referenceNt.sequence[SLICE])
            print("gen  nt:", genomeNt.sequence[SLICE])

        self._cache["nt"][cacheKey] = referenceNt, genomeNt

        return referenceNt, genomeNt

    def aaSequences(
        self,
        featureName: str,
        addFeatureToId: bool = True,
        raiseOnReferenceGaps: bool = True,
    ) -> Tuple[DNARead, DNARead]:
        """
        Match the genome and the reference at the amino acid level.

        @param featureName: A C{str} feature name.
        @param addFeatureToId: If C{True}, the feature name (in parentheses,
            preceded by a space) will be appended to the id of the returned
            sequences.
        @param raiseOnReferenceGaps: A C{bool}. If C{True} and the aligner
            suggests a gap in the reference, raise a ReferenceInsertionError.
        @raise ReferenceInsertionError: (as above for raiseOnReferenceGaps).
        @raise TranslationError: or one of its sub-classes (see translate.py)
            if a feature nucleotide sequence cannot be translated.
        @return: A 2-C{tuple} of C{dark.reads.AARead} instances, with 1) the
            amino acids for the feature as located in the reference genome and
            2) the corresponding amino acids from the genome being examined.
        """
        cacheKey = featureName, addFeatureToId

        try:
            return self._cache["aa"][cacheKey]
        except KeyError:
            pass

        referenceNt, genomeNt = self.ntSequences(
            featureName,
            raiseOnReferenceGaps=raiseOnReferenceGaps,
            reverseComplement=True,
        )
        feature = self.features[featureName]
        name = feature["name"]
        idSuffix = f" ({name})" if addFeatureToId else ""

        gapCount = genomeNt.sequence.count("-")
        if (
            self.features.sars2
            and name == "surface glycoprotein"
            and gapCount
            and gapCount % 3 == 0
        ):
            referenceAa = referenceAaAligned = AARead(
                self.features.reference.id + idSuffix,
                translateSARS2Spike(referenceNt.sequence),
            )

            genomeAa = genomeAaAligned = AARead(
                self.genome.id + idSuffix, translateSARS2Spike(genomeNt.sequence)
            )

            if len(referenceAaAligned) != len(genomeAaAligned):
                raise TranslatedReferenceAndGenomeLengthError(
                    f"Genome and reference AA sequences lengths differ "
                    f"({len(genomeAaAligned)} != {len(referenceAaAligned)})."
                )
        else:
            if "-" in referenceNt.sequence:
                if raiseOnReferenceGaps:
                    raise ReferenceInsertionError(
                        f"{self._aligner} aligner suggested a reference insertion "
                        f"into {featureName!r}."
                    )
                referenceNtNoGaps = DNARead(
                    referenceNt.id, referenceNt.sequence.replace("-", "")
                )
            else:
                referenceNtNoGaps = referenceNt

            referenceAa = AARead(
                self.features.reference.id + idSuffix,
                feature.get(
                    "translation",
                    translate(
                        referenceNtNoGaps.sequence,
                        name,
                        sars2=self.features.sars2,
                        untranslatable=self._untranslatable,
                    ),
                ),
            )

            genomeNtNoGaps = DNARead(genomeNt.id, genomeNt.sequence.replace("-", ""))

            if len(genomeNtNoGaps) == 0:
                # This feature has been deleted in the genome. Give back an
                # "aligned" genome that is all gaps. This matches what
                # ntSequences will return.
                referenceAaAligned = referenceAa
                genomeAaAligned = AARead(
                    self.features.reference.id + idSuffix, "-" * len(referenceAaAligned)
                )
            else:

                genomeAa = AARead(
                    self.genome.id + idSuffix,
                    translate(
                        genomeNtNoGaps.sequence,
                        name,
                        sars2=self.features.sars2,
                        untranslatable=self._untranslatable,
                    ),
                )

                if self._aligner == "mafft":
                    referenceAaAligned, genomeAaAligned = mafft(
                        Reads([referenceAa, genomeAa]), options=MAFFT_OPTIONS
                    )
                else:
                    assert self._aligner == "edlib"
                    # This is an AA alignment, so set matchAmbiguous=False
                    # unconditionally since the edlib aligner only knows about
                    # nucleotide ambiguities.
                    try:
                        referenceAaAligned, genomeAaAligned = edlibAlign(
                            Reads([referenceAa, genomeAa]), matchAmbiguous=False
                        )
                    except Exception as e:
                        raise AlignmentError(
                            f"Could not edlib align: {e}\n"
                            f"  ref nt: {referenceNt.sequence}\n"
                            f"  gen nt: {genomeNt.sequence}"
                            f"  ref aa: {referenceAa.sequence}\n"
                            f"  gen aa: {genomeAa.sequence}"
                        )

        if DEBUG:
            print(f"AA MATCH {name}:")

            print(
                f"ref nt aligned {len(referenceNt.sequence)}:",
                referenceNt.sequence[SLICE],
            )
            print(f"gen nt aligned {len(genomeNt.sequence)}:", genomeNt.sequence[SLICE])

            print(
                f"ref aa        {len(referenceAa.sequence)}:",
                referenceAa.sequence[SLICE],
            )
            print(f"gen aa        {len(genomeAa.sequence)}:", genomeAa.sequence[SLICE])

            print(
                f"ref aa aligned {len(referenceAaAligned.sequence)}:",
                referenceAaAligned.sequence[SLICE],
            )
            print(
                f"gen aa aligned {len(genomeAaAligned.sequence)}:",
                genomeAaAligned.sequence[SLICE],
            )

        self._cache["aa"][cacheKey] = referenceAaAligned, genomeAaAligned
        return referenceAaAligned, genomeAaAligned

    def _checkChange(
        self,
        base: str,
        offset: int,
        read: DNARead,
        change: Union[str, Tuple[Optional[str], int, Optional[str]]],
        featureName: str,
        onError,
        errFp,
    ) -> [bool, Optional[str]]:
        """
        Check that a base occurs at an offset.

        @param base: A C{str} nucleotide base (or amino acid), or C{None}
            in which case there was no expectation that the read had any
            particular value and C{True} is returned.
        @param offset: A 0-based C{int} offset into C{sequence}.
        @param read: A C{dark.reads.Read} instance.
        @param change: The C{str} or C{tuple} change specification.
        @param featureName: A C{str} feature name.
        @param onError: A C{str} indicating what to do if an error is
            encountered. Must be one of 'raise', 'ignore', or 'print' (in which
            case an error message will be printed to C{errFp}.
        @param errFp: An open file pointer to write error messages to, if any.
            Only used if C{onError} is 'print'.
        @raise IndexError: If the offset is out of range.
        @return: A C{list} containing a C{bool} indicating whether the read
            sequence has C{base} at C{offset} (or C{True} if C{base} is
            C{None}) and the C{str} base found at the offset. If there is an
            error (and C{onError} is not 'raise'), return [False, None].
        """
        try:
            actual = read.sequence[offset]
        except IndexError:
            mesg = (
                f"Index {offset} out of range trying to access feature "
                f"{featureName!r} of length {len(read)} sequence "
                f"{read.id!r} via expected change specification "
                f"{change!r}."
            )
            if onError == "raise":
                raise IndexError(mesg)
            elif onError == "print":
                print(mesg, file=errFp)
                return [False, None]
            else:
                assert onError == "ignore"
                return [False, None]
        else:
            return [(base is None or actual == base), actual]

    def checkFeature(
        self,
        featureName: str,
        changes: Union[str, Iterable[Tuple[str, int, str]]],
        aa: bool = False,
        onError: str = "raise",
        errFp: Optional[TextIO] = None,
        raiseOnReferenceGaps: bool = True,
    ) -> Tuple[
        int,
        int,
        Dict[
            Union[str, Tuple[str, int, str]],
            Tuple[bool, Optional[str], bool, Optional[str]],
        ],
    ]:
        """Check that a set of changes all happened as expected.

        @param featureName: A C{str} feature name.
        @param changes: Either a C{str} or a iterable of 3-C{tuple}s. If a
            C{str}, gives a specification in the form of space-separated
            RNG strings, where R is a reference base, N is a 1-based site,
            and G is a sequence base. So, e.g., 'L28S P1003Q' indicates that
            we expected a change from 'L' to 'S' at site 28 and from 'P' in
            the reference to 'Q' in the genome we're examining at site 1003.
            The reference or genome base (but not both) may be absent.
            If an iterable of 3-C{tuple}s, each tuple should have a C{str}
            expected reference base, a 0-based offset, and a C{str} expected
            genome base. Note that the string format (meant for humans) uses
            1-based locations whereas the tuple format uses 0-based offsets.
        @param aa: If C{True} check amino acid sequences. Else nucleotide
        @param onError: A C{str} indicating what to do if an error is
            encountered. Must be one of 'raise', 'ignore', or 'print' (in which
            case an error message will be printed to C{errFp}.
        @param errFp: An open file pointer to write error messages to, if any.
            Only used if C{onError} is 'print'.
        @param raiseOnReferenceGaps: A C{bool}. If C{True} and the aligner
            suggests a gap in the reference, raise a ReferenceInsertionError.
        @raise ValueError: If a change string cannot be parsed.
        @raise IndexError: If the offset of a change exceeds the length of the
            sequence being checked.
        @return: A 3-C{tuple} with the C{int} number of checks done, the
            C{int} number of errors, and a C{dict} keyed by changes in
            C{changes}. The C{dict} values are 4-C{tuple}s of (Boolean,
            reference base, Boolean, genome base) to indicate success or
            failure of the check for the reference and the reference base
            found, then the same thing for the genome. E.g., a tuple of
            (True, 'A', False, 'T') would indicate that the expected reference
            base was found, that it was an 'A', but that the expected genome
            base was not found and that instead a 'T' was found. If C{aa} is
            C{True} and there is a translation error, the 4-tuple will contain
            (False, None, False, None).
        """

        def _getChanges(
            changes: Union[str, Iterable[Tuple[str, int, str]]]
        ) -> Iterator[Tuple[Union[str, Tuple[str, int, str]], str, int, str]]:
            if isinstance(changes, str):
                for change in changes.split():
                    referenceBase, offset, genomeBase = splitChange(change)
                    yield change, referenceBase, offset, genomeBase
            else:
                for change in changes:
                    referenceBase, offset, genomeBase = change
                    yield change, referenceBase, offset, genomeBase

        result: Dict[
            Union[str, Tuple[Optional[str], int, Optional[str]]],
        ] = {}
        testCount = errorCount = 0

        try:
            reference, genome = (
                self.aaSequences(featureName, raiseOnReferenceGaps=raiseOnReferenceGaps)
                if aa
                else self.ntSequences(featureName)
            )
        except TranslationError as e:
            if onError == "raise":
                raise
            elif onError == "print":
                print(e, file=errFp)

            translationError = True
        else:
            translationError = False

        for change, referenceBase, offset, genomeBase in _getChanges(changes):
            if translationError:
                result[change] = (False, None, False, None)
            else:
                result[change] = tuple(
                    self._checkChange(
                        referenceBase,
                        offset,
                        reference,
                        change,
                        featureName,
                        onError,
                        errFp,
                    )
                    + self._checkChange(
                        genomeBase, offset, genome, change, featureName, onError, errFp
                    )
                )

        errorCount = sum((v[0] is False or v[2] is False) for v in result.values())
        testCount = len(result)

        return testCount, errorCount, result

    def checkVariant(
        self,
        variant: Union[str, dict],
        onError: str = "raise",
        errFp: Optional[TextIO] = None,
        raiseOnReferenceGaps: bool = True,
    ) -> Tuple[int, int, Dict[str, Tuple[bool, Optional[str], bool, Optional[str]]]]:
        """
        Check that a set of changes in different features all happened as
        expected.

        @param variant: Either the C{str} name of a key in the known VARIANTS
            C{dict} or else a C{dict} with the same structure as the 'changes'
            value C{dict} in the known variants (see variants.py).
        @param onError: A C{str} indicating what to do if an error is
            encountered. Must be one of 'raise', 'ignore', or 'print' (in which
            case an error message will be printed to C{errFp}.
        @param errFp: An open file pointer to write error messages to, if any.
            Only used if C{onError} is 'print'.
        @param raiseOnReferenceGaps: A C{bool}. If C{True} and the aligner
            suggests a gap in the reference, raise a ReferenceInsertionError.
        @raise ValueError: If incorrect arguments are passed (see below).
        @return: A 3-C{tuple} with the number of checks done, the number of
            errors, and a C{dict} keyed by C{str} changes in C{changes}, with values
            a 2-C{tuple} of Booleans to indicate success or failure of the
            check for the reference and the genome respectively. There is no
            specific indication of any TranslationError when checking amino
            acid sequences.
        """
        if onError == "print" and not errFp:
            raise ValueError(
                'If you pass onError="print" you must also '
                "give a file descriptor via errFp."
            )

        result = {}
        testCountTotal = errorCountTotal = 0

        if isinstance(variant, str):
            changeDict = VARIANTS[variant]["changes"]
        else:
            assert isinstance(variant, dict)
            changeDict = variant

        for featureName in changeDict:
            result[featureName] = {}
            for what in "aa", "nt":
                try:
                    changes = changeDict[featureName][what]
                except KeyError:
                    pass
                else:
                    testCount, errorCount, result_ = self.checkFeature(
                        featureName,
                        changes,
                        aa=(what == "aa"),
                        onError=onError,
                        errFp=errFp,
                        raiseOnReferenceGaps=raiseOnReferenceGaps,
                    )
                    testCountTotal += testCount
                    errorCountTotal += errorCount
                    result[featureName][what] = result_

        return testCountTotal, errorCountTotal, result

    def coverage(self, featureName: Optional[str] = None) -> Tuple[int, int]:
        """
        Get the coverage of a feature or the whole genome.

        @param featureName: A C{str} feature name, or C{None} to examine the
            whole genome.
        @return: A 2-C{tuple} with the C{int} number of aligned genome sites
            for which there is coverage of the reference, and the C{int} number
            of reference sites. The quotient of the two gives the coverage
            fraction.
        """
        referenceNt, genomeNt = (
            self.ntSequences(featureName)
            if featureName
            else (self.referenceAligned, self.genomeAligned)
        )

        referenceCount = coveredCount = 0
        for referenceBase, genomeBase in zip(referenceNt.sequence, genomeNt.sequence):
            if referenceBase != "-":
                referenceCount += 1
                coveredCount += genomeBase not in "N-"

        return coveredCount, referenceCount

    def offsetInfo(
        self,
        offset: int,
        relativeToFeature: bool = False,
        aa: bool = False,
        featureName: Optional[str] = None,
        includeUntranslated: bool = False,
        minReferenceCoverage: Optional[float] = None,
        allowAmbiguous: bool = True,
    ) -> Union[dict, None]:
        """
        Get information about genome features at an offset.

        @param offset: An C{int} offset.
        @param relativeToFeature: If C{True}, the offset is relative to the
            start of the feature that occurs at this offset.
        @param aa: If C{True}, the offset is a number of amino acids, else a
            number of nucleotides.
        @param featureName: If not C{None} and multiple features occur at this
            offset, use the feature with this name.
        @param includeUntranslated: If C{True}, also return features that are
            not translated.
        @param minReferenceCoverage: The C{float} fraction of non-N bases
            required in the genome (or feature, if one is given) in order
            for it to be processed. If the required coveragel is not met,
            C{None} is returned.
        @allowAmbiguous: If C{True}, do not raise an error if multiple features
            are found for the offset and no feature name is given to
            disambiguate. Instead, use the first feature name returned by
            C{getFeatureNames}.
        @raise KeyError: If the feature name is unknown.
        @raise ValueError: If incorrect arguments are passed (see below).
        @raise AmbiguousFeatureError: If multiple features occur at the offset
            and C{featureName} does not indicate the one to use.
        @return: A C{dict} with information about what is found in the
            reference and the genome at the offset. See the C{result}
            dictionary below. If C{minReferenceCoverage} is not C{None} and
            the coverage of the genome (or feature, if given) is lower,
            C{None} is returned.
        """
        if minReferenceCoverage is not None:
            genomeCount, referenceCount = self.coverage(featureName)
            if genomeCount / referenceCount < minReferenceCoverage:
                return None

        if relativeToFeature:
            if featureName is None:
                raise ValueError(
                    "If relativeToFeature is True, a feature name must be given."
                )
            referenceOffset = self.features.referenceOffset(featureName, offset, aa)
        else:
            if aa:
                raise ValueError(
                    "You cannot pass aa=True unless the offset "
                    "you pass is relative to the feature."
                )
            referenceOffset = offset

        feature, features = self.features.getFeature(
            referenceOffset, featureName, includeUntranslated, allowAmbiguous
        )

        # The 'None' values in the following are filled in below. They're set
        # up-front here to explicitly show what's in the returned dictionary.
        #
        # The 'aa' value will be None if the reference or genome sequence ends
        # within three nucleotides of the given offset (i.e., the sequence is
        # too short to get a 3-nucleotide codon). The 'aa' value in the
        # 'genome' dict will be '-' if the aligned genome contains a gap at any
        # of the codon locations in the reference and 'X' if the offset falls
        # in the last two nucleotides of the genome (and so the codon is
        # too short for translation).
        #
        # The 'frame' value will be 0, 1, or 2.
        result = {
            "featureName": feature["name"] if feature else None,
            "featureNames": features,
            "alignmentOffset": None,
            "reference": {
                "aa": None,
                "codon": None,
                "frame": None,
                "id": self.features.reference.id,
                "aaOffset": None,
                "ntOffset": None,
            },
            "genome": {
                "aa": None,
                "codon": None,
                "frame": None,
                "id": self.genome.id,
                "aaOffset": None,
                "ntOffset": None,
            },
        }

        # Reference.
        if relativeToFeature:
            referenceCodonAaOffset = offset if aa else offset // 3
            referenceFrame = 0 if aa else offset % 3
        else:
            referenceCodonAaOffset, referenceFrame = divmod(
                referenceOffset - (feature["start"] if feature else 0), 3
            )

        codonOffset = referenceOffset - referenceFrame
        codon = self.features.reference.sequence[codonOffset : codonOffset + 3]
        result["reference"].update(
            {
                "aa": str(Seq(codon).translate()) if len(codon) == 3 else "X",
                "codon": codon,
                "frame": referenceFrame,
                "aaOffset": referenceCodonAaOffset,
                "ntOffset": referenceOffset,
            }
        )

        # Genome.
        gappedOffset = self.gappedOffsets[referenceOffset]
        if relativeToFeature or feature:
            genomeGappedStart = self.gappedOffsets[feature["start"]]
            assert genomeGappedStart <= gappedOffset
            gapCount = self.genomeAligned.sequence[
                genomeGappedStart:gappedOffset
            ].count("-")
        else:
            genomeGappedStart = 0
            gapCount = self.genomeAligned.sequence[:gappedOffset].count("-")

        genomeCodonAaOffset, genomeFrame = divmod(
            gappedOffset - genomeGappedStart - gapCount, 3
        )

        codonOffset = gappedOffset - genomeFrame
        codon = self.genomeAligned.sequence[codonOffset : codonOffset + 3]
        result["genome"].update(
            {
                "aa": (
                    "-"
                    if "-" in codon
                    else (str(Seq(codon).translate()) if len(codon) == 3 else "X")
                ),
                "codon": codon,
                "frame": genomeFrame,
                "aaOffset": genomeCodonAaOffset,
                "ntOffset": gappedOffset
                - (self.genomeAligned.sequence[:gappedOffset].count("-")),
            }
        )

        result["alignmentOffset"] = gappedOffset

        return result


def offsetInfoMultipleGenomes(
    genomes: Iterable[Gb2Alignment],
    offset: int,
    relativeToFeature: bool = False,
    aa: bool = False,
    featureName: Optional[str] = None,
    includeUntranslated: bool = False,
    minReferenceCoverage: Optional[float] = None,
) -> Optional[dict]:
    """
    Get information about an offset for multiple Gb2Alignment instances.

    @param genomes: An interable of C{Gb2Alignment} instances.
    @param offset: An C{int} offset.
    @param relativeToFeature: If C{True}, the offset is relative to the
        start of the feature that occurs at this offset.
    @param aa: If C{True}, the offset is a number of amino acids, else a
        number of nucleotides.
    @param featureName: If not C{None} and multiple features occur at this
        offset, use the feature with this name.
    @param includeUntranslated: If C{True}, also return features that are
        not translated.
    @param minReferenceCoverage: The C{float} fraction of non-N bases required
        in the genome (or feature, if one is given) in order for it to be
        processed. If the required coveragel is not met, C{None} is returned.
    @raise KeyError: If the feature name is unknown.
    @raise ValueError: If an incorrect arguments is passed
        (see Gb2Alignment.offsetInfo) or if all genomes were not aligned
        against the same reference id.
    @raise AmbiguousFeatureError: If multiple features occur at the offset
        and C{featureName} does not indicate the one to use.
    @return: A C{dict} with information about what is found in the
        reference and the genomes at the offset. This is identical to the
        result dictionary returned by C{Gb2Alignment.offsetInfo}, but the
        'genome' key is replaced with a 'genomes' key that holds a C{list}
        of genome results (each is a C{dict}, as returned by
        C{Gb2Alignment.offsetInfo}). The order of genome info in the list
        matches the order in the passed C{genomes}.  If no genomes are passed,
        or all passed genomes have insufficient coverage, C{None} is returned.
    """
    ids = set(genome.features.reference.id for genome in genomes)

    if len(ids) == 0:
        raise ValueError("No Gb2Alignment instances given.")
    elif len(ids) != 1:
        raise ValueError(
            f"Gb2Alignment instances with differing reference ids "
            f'passed to offsetInfoMultipleGenomes: {", ".join(sorted(ids))}.'
        )

    first = True
    result = None
    genomeResults = []

    for genome in genomes:
        thisResult = genome.offsetInfo(
            offset,
            relativeToFeature=relativeToFeature,
            aa=aa,
            featureName=featureName,
            includeUntranslated=includeUntranslated,
            minReferenceCoverage=minReferenceCoverage,
        )

        if first and thisResult:
            result = thisResult.copy()
            # Make sure nothing crucial has changed in the returned dict
            # from Gb2Alignment.offsetInfo
            assert "genome" in result
            assert "genomes" not in result
            del result["genome"]
            first = False

        genomeResults.append(thisResult["genome"] if thisResult else None)

    if result is None:
        return None
    else:
        result["genomes"] = genomeResults
        return result
