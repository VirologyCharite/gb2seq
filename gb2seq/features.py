import sys
from os import environ
import itertools
from collections import UserDict
from pathlib import Path
from typing import Dict, Optional, Set, Union

try:
    from importlib.resources import files, as_file
except ImportError:
    from importlib_resources import files, as_file


# from warnings import warn
import json
import argparse

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from dark.aaVars import STOP_CODONS
from dark.genbank import GenomeRanges
from dark.reads import DNARead

from gb2seq import Gb2SeqError
from gb2seq.sars2 import SARS_COV_2_ALIASES, SARS_COV_2_TRANSLATED

# Set ENTREZ_EMAIL in your environment to have your requests to NCBI Entez
# be accompanied by your address. If you don't do this you'll see warning
# messages and be limited to a lower rate of querying.
Entrez.email = environ.get("ENTREZ_EMAIL")


class ReferenceWithGapError(Gb2SeqError):
    "A GenBank reference sequence had a gap."


class MissingFeatureError(Gb2SeqError):
    "A feature expected at an offset is not present."


class AmbiguousFeatureError(Gb2SeqError):
    "More than one feature is referred to by an offset."


class UnknownFeatureNameError(Gb2SeqError):
    "A feature name is unknown."


class Features(UserDict):
    """
    Manage sequence features.

    @param referenceSpecification: Either:
        * A C{str} name or C{Path} of a GenBank file containing the features.
        * A C{str} GenBank accession id.
        * A C{dict} of pre-prepared features, in which case C{reference}
              must not be C{None}. Passing a C{dict} is provided for testing.
        * C{None}, in which case the default reference, NC_045512.2.gb, is
              loaded.
    @param reference: A C{dark.reads.DNARead} instance if C{spec} is a C{dict},
        else C{None}.
    @param sars2: A C{bool} indicating whether we are dealing with SARS-CoV-2
        features (in which case some defaults can be set).
    @param alsoInclude: A C{set} of feature types to also return data for
        (beyond those normally returned), or C{None}.
    @raise ValueError: If a reference is passed with a string or Path
        specification.
    @raise ReferenceWithGapError: If the reference or one of its features has a
        gap in its nucleotide sequence.
    """

    def __init__(
        self,
        referenceSpecification: Union[str, Dict, Path] = None,
        reference: Optional[DNARead] = None,
        sars2: bool = False,
        translated: Optional[Set[str]] = None,
        aliases: Optional[Dict[str, str]] = None,
        addUnannotatedRegions: bool = False,
        alsoInclude: Optional[set[str]] = None,
    ) -> None:
        super().__init__()
        self.sars2 = sars2
        self.translatedNames: Set[str] = set()

        if referenceSpecification is None:
            if sars2:
                _gb = files("gb2seq").joinpath("data").joinpath("NC_045512.2.gb")
                with as_file(_gb) as fp:
                    try:
                        record = SeqIO.read(fp, "genbank")
                    except Exception as e:
                        print(
                            "Could not parse default SARS-CoV-2 GenBank record",
                            e,
                            file=sys.stderr,
                        )
                        sys.exit(1)
                    else:
                        self._initializeFromGenBankRecord(record, alsoInclude)
            else:
                raise ValueError(
                    "A reference specification must be provided for non-SARS-CoV-2 "
                    "features."
                )
        elif isinstance(referenceSpecification, SeqRecord):
            record = referenceSpecification
            self._initializeFromGenBankRecord(record, alsoInclude)
        elif isinstance(referenceSpecification, str):
            if reference is not None:
                raise ValueError(
                    "A reference cannot be passed with a string specification."
                )
            path = Path(referenceSpecification)
            if path.exists():
                # A file argument can either be in GenBank format or contain a JSON
                # object (the saved output of annotate-genome.py).
                jsonError = seqError = None
                with open(path) as fp:
                    try:
                        record = json.load(fp)
                    except json.decoder.JSONDecodeError as e:
                        jsonError = e
                    else:
                        self.data.update(record["features"])
                        self.reference = DNARead(record["id"], record["sequence"])

                if jsonError:
                    with open(path) as fp:
                        try:
                            record = SeqIO.read(fp, "genbank")
                        except Exception as e:
                            seqError = e
                        else:
                            self._initializeFromGenBankRecord(record, alsoInclude)

                if jsonError and seqError:
                    print(
                        f"Could not read {referenceSpecification!r} as a JSON or "
                        f"GenBank file. Here are the parsing errors.\nJSON: "
                        f"{jsonError}\nGenBank: {seqError}",
                        file=sys.stderr,
                    )
                    sys.exit(1)
            else:
                print(
                    f"Fetching GenBank record for {referenceSpecification!r}.",
                    file=sys.stderr,
                )
                try:
                    client = Entrez.efetch(
                        db="nucleotide",
                        rettype="gb",
                        retmode="text",
                        id=referenceSpecification,
                    )
                    try:
                        record = SeqIO.read(client, "genbank")
                    except Exception as e:
                        print(
                            "Could not parse fetched GenBank record:",
                            e,
                            file=sys.stderr,
                        )
                        sys.exit(1)
                    else:
                        self._initializeFromGenBankRecord(record, alsoInclude)
                    finally:
                        client.close()
                except Exception as e:
                    print("Could not fetch GenBank record:", e, file=sys.stderr)
                    sys.exit(1)
        elif isinstance(referenceSpecification, Path):
            if reference is not None:
                raise ValueError(
                    "A reference cannot be passed with a Path specification."
                )
            with open(referenceSpecification) as fp:
                record = SeqIO.read(fp, "genbank")
            self._initializeFromGenBankRecord(record, alsoInclude)
        elif isinstance(referenceSpecification, dict):
            self.data.update(referenceSpecification)
            self.reference = reference
        else:
            raise ValueError(f"Unrecognized specification {referenceSpecification!r}.")

        if sars2:
            self.translatedNames = (
                SARS_COV_2_TRANSLATED if translated is None else translated
            )
            self.aliasDict = SARS_COV_2_ALIASES if aliases is None else aliases
        else:
            if translated is None:
                self.translatedNames = set(
                    name for (name, value) in self.items() if "translation" in value
                )
            else:
                self.translatedNames = translated
            self.aliasDict = {} if aliases is None else aliases

        if addUnannotatedRegions:
            self._addUnannotatedRegions()

    def _initializeFromGenBankRecord(
        self, record: SeqIO.SeqRecord, alsoInclude: Optional[set[str]]
    ) -> None:
        """
        Initialize from a GenBank record.

        @param record: A BioPython C{SeqRecord} sequence record.
        @param alsoInclude: A C{set} of feature types to also return data for
            (beyond those normally returned) or C{None}.
        @raise ReferenceWithGapError: If the reference sequence or the
            sequence of any of its features has a gap.
        """
        self.reference = DNARead(record.id, str(record.seq))
        alsoInclude = alsoInclude or set()

        for feature in record.features:
            type_ = feature.type
            value = {
                "type": type_,
            }

            if type_ == "3'UTR" or type_ == "5'UTR":
                name = type_

            elif type_ == "stem_loop":
                for n in itertools.count(1):
                    name = f"stem loop {n}"
                    if name not in self:
                        break
                try:
                    value["function"] = feature.qualifiers["function"][0]
                except KeyError:
                    pass

            elif type_ == "repeat_region":
                for n in itertools.count(1):
                    name = f"repeat region {n}"
                    if name not in self:
                        break
                value["type"] = feature.qualifiers["rpt_type"][0]

            else:
                if "product" in feature.qualifiers:
                    # This covers "CDS", "mat_peptide", "ncRNA".
                    name = feature.qualifiers["product"][0]
                    value["product"] = name
                else:
                    # This covers "source", "gap", "gene", "misc_feature"
                    # and any other feature that is not annotated as having
                    # a product.
                    if type_ not in alsoInclude:
                        continue

            start = int(feature.location.start)
            stop = int(feature.location.end)
            genomeRanges = GenomeRanges(str(feature.location))

            # We can only handle a single range at the moment.
            if len(genomeRanges.ranges) == 1:
                assert start == genomeRanges.ranges[0][0]
                assert stop == genomeRanges.ranges[0][1]
                forward = genomeRanges.ranges[0][2]
            elif self.sars2 and name == "ORF1ab polyprotein":
                assert len(genomeRanges.ranges) == 2
                assert start == genomeRanges.ranges[0][0]
                assert stop == genomeRanges.ranges[1][1]
                forward = True
            else:
                if not self.sars2:
                    # At some point (soon) we should emit a warning. But let's first try
                    # to fix things so we can translate anything.
                    #
                    # warn(
                    #     f"Multiple reference genome ranges {genomeRanges} found "
                    #     f"for feature {name!r} will not be translated reliably."
                    # )
                    pass

            sequence = str(record.seq)[start:stop]

            value.update(
                {
                    "forward": forward,
                    "name": name,
                    "sequence": sequence,
                    "start": start,
                    "stop": stop,
                }
            )

            for optional in "translation", "note":
                try:
                    value[optional] = feature.qualifiers[optional][0]
                except KeyError:
                    pass

            # If there is a translation, add an amino acid '*' stop
            # indicator if there is not one already and the sequence ends
            # with a stop codon.
            try:
                translation = value["translation"]
            except KeyError:
                pass
            else:
                if not translation.endswith("*"):
                    if forward:
                        codon = value["sequence"][-3:]
                    else:
                        codon = (
                            DNARead("id", sequence).reverseComplement().sequence[-3:]
                        )
                    if codon.upper() in STOP_CODONS:
                        value["translation"] += "*"

            # The following elaborate dance makes sure that we use a name
            # for the feature that is unique, including when case is
            # ignored.  This makes it possible to use the feature names in
            # self as file names even on a case-insensitive filesystem.
            # Otherwise you will certainly end up with a mess, due to
            # genome annotations that have multiple features with identical
            # names (e.g., "hypothetical protein" several times, along with
            # "Hypothetical protein").  E.g., download the GenBank file for
            # NC_003310.1 (https://www.ncbi.nlm.nih.gov/nuccore/NC_003310.1/)
            # and run this:
            #
            # $ grep product NC_003310.1.gb | grep -i hypothetical
            existingName = self.getKeyIgnoringCase(name)
            if existingName:
                if self[existingName] != value:
                    lowercaseNames = set(map(str.lower, self))
                    for n in itertools.count(2):
                        adjustedName = f"{name} {n}"
                        if adjustedName.lower() not in lowercaseNames:
                            name = adjustedName
                            value["name"] = adjustedName
                            self[name] = value
                            break
            else:
                self[name] = value

        self._checkForGaps()

    def _addUnannotatedRegions(self) -> None:
        """
        Find unannotated regions and add them as features.
        """
        annotatedOffsets: Set[int] = set()
        for feature in self.data.values():
            annotatedOffsets.update(range(feature["start"], feature["stop"]))
            feature["annotated"] = True

        start = None
        unannotatedRegionCount = 0

        def _addNew(stop: int) -> None:
            nonlocal unannotatedRegionCount
            unannotatedRegionCount += 1
            name = f"unannotated region {unannotatedRegionCount}"
            assert name not in self.data
            self.data[name] = {
                "name": name,
                "annotated": False,
                "start": start,
                "stop": offset,
                "sequence": self.reference.sequence[start:stop],
            }

        for offset in range(len(self.reference)):
            if offset in annotatedOffsets:
                if start is not None:
                    _addNew(offset)
                    start = None
            else:
                if start is None:
                    start = offset

        if start is not None:
            _addNew(offset)

    def getKeyIgnoringCase(self, name: str) -> Union[str, None]:
        """
        Find a key in self.data, ignoring case.

        @param name: A C{str} name to look up.
        @return: An existing C{str} key that matches C{name} when case is
            ignored, or C{None} if no such key exists.
        """
        name = name.lower()
        for thisName in self.data:
            if thisName.lower() == name:
                return thisName
        else:
            return None

    def __getitem__(self, name: str) -> Union[dict, None]:
        """
        Find a feature by name. This produces a dictionary with the following
        keys and values for the feature:

            {
                'function': A string description of the feature's function,
                    if one was present in the GenBank file.
                'product': A string description of the product of the feature.
                    if one was present in the GenBank file.
                'name': The string name.
                'note': A string note, if one was is present in the GenBank
                    file.
                'sequence': The string nucleotide sequence.
                'start': A 0-based integer offset of the first nucleotide of
                    the feature in the complete genome.
                'stop': A 0-based integer offset of the nucleotide after the
                    last nucleotide in the feature in the complete genome.
                    Thus start and stop can be used in the regular Python
                    way of slicing out substrings.
                'translation': The amino acid translation of the feature, if
                    one is provided by the GenBank record. Note that
                    translations may not be provided for all features! If you
                    need to know what is actually translated, use the
                    TRANSLATED set defined at the top of this file.
            }


        @param name: A C{str} feature name to look up.
        @raise UnknownFeatureNameError: If the feature name is unknown.
        @return: A C{dict} for the feature, as above.
        """
        return self.data[self.canonicalName(name)]

    def canonicalName(self, name: str) -> str:
        """
        Get the canonical name for a feature.

        @param name: A C{str} feature name to look up.
        @raise UnknownFeatureNameError: If C{name} is not a known feature name.
        @return: A C{str} canonical feature name.
        """
        if name in self:
            return name

        nameLower = name.lower()
        for featureName in self:
            if nameLower == featureName.lower():
                return featureName

        alias = self.aliasDict.get(nameLower)
        if alias is None:
            raise UnknownFeatureNameError(name)
        else:
            assert alias in self
            return alias

    def translated(self, name: str) -> bool:
        """
        Is a feature translated.

        @param name: A C{str} feature name.
        @return: A C{bool} to indicate whether the feature is translated.
        """
        return name in self.translatedNames

    def aliases(self, name: str) -> set:
        """
        Get all aliases for a name.

        @param name: A C{str} feature name.
        @return: A C{set} of C{str} canonical names.
        """
        try:
            canonicalName = self.canonicalName(name)
        except UnknownFeatureNameError:
            return set()

        result = {canonicalName}

        for alias, canonical in self.aliasDict.items():
            if canonical == canonicalName:
                result.add(alias)

        return result

    def _checkForGaps(self) -> None:
        """
        Check there are no gaps in the reference or any feature sequence.

        @raise ReferenceWithGapError: If the reference sequence or the
            sequence of any of its features has a gap.
        """
        referenceId = self.reference.id
        if self.reference.sequence.find("-") > -1:
            raise ReferenceWithGapError(
                f"Reference sequence {referenceId!r} has a gap!"
            )

        for featureInfo in self.values():
            if featureInfo["sequence"].find("-") > -1:
                raise ReferenceWithGapError(
                    f'Feature {featureInfo["name"]!r} sequence in '
                    f"{referenceId!r} has a gap!"
                )

    def referenceOffset(self, name: str, offset: int, aa: bool = False) -> int:
        """
        Get the (nucleotide) offset in the reference, given an offset in a
        feature.

        @param name: A C{str} feature name.
        @param offset: An C{int} offset into the feature.
        @param aa: If C{True}, the offset is a number of amino acids, else a
            number of nucleotides.
        @raise UnknownFeatureNameError: If the name is unknown.
        @return: An C{int} nucleotide offset into the reference genome.
        """
        return self[name]["start"] + offset * (3 if aa else 1)

    def getFeatureNames(self, offset: int, includeUntranslated: bool = False) -> set:
        """
        Get the names of all features that overlap a given offset.

        @param offset: An C{int} offset into the genome.
        @param includeUntranslated: If C{True}, also return features that are
            not translated.
        @return: A C{set} of C{str} feature names.
        """
        result = set()

        for name, feature in self.items():
            if feature["start"] <= offset < feature["stop"] and (
                includeUntranslated or self.translated(name)
            ):
                result.add(name)

        return result

    def getFeature(
        self,
        offset: int,
        featureName: Optional[str] = None,
        includeUntranslated: bool = False,
        allowAmbiguous: bool = True,
    ):
        """
        Find a single feature at an offset.

        @param offset: An C{int} offset into the genome.
        @param featureName: A C{str} feature name. Used for disambiguation when
            multiple features are present at the offset. If C{None}, a feature
            will be returned if there is only one at the given offset.
        @param includeUntranslated: If C{True}, also consider features that are
            not translated.
        @allowAmbiguous: If C{True}, do not raise an error if multiple features
            are found for the offset and no feature name is given to
            disambiguate. Instead, use the alphabetically first feature name
            returned by C{getFeatureNames}.
        @raise MissingFeatureError: if the requested feature does not overlap
            the given C{offset} or if there are no features at the offset.
        @raise AmbiguousFeatureError: if there are multiple features at the
            offset and a feature name is not given.
        @raise UnknownFeatureNameError: If the feature name is unknown (can be raised
            via the call to self.canonicalName)
        @return: A 2-tuple, containing 1) a features C{dict} (as returned by
            __getitem__), if the requested feature is among those present at
            the offset, and 2) the set of all C{str} feature names present at
            the offset.
        """
        features = self.getFeatureNames(offset, includeUntranslated)

        if featureName is None:
            if features:
                # There are some features here, but we weren't told which
                # one to use. Only proceed if there's just one or if we've
                # been told to allow ambiguity (in which case we use the
                # alphabetically first feature name).
                if len(features) == 1 or allowAmbiguous:
                    featureName = list(sorted(features))[0]
                    feature = self[featureName]
                else:
                    present = ", ".join(
                        f"{f!r} ({self[f]['start'] + 1} - {self[f]['stop']})"
                        for f in sorted(features)
                    )
                    raise AmbiguousFeatureError(
                        f"There are multiple features at site {offset + 1}: "
                        f"{present}. Pass a feature name to specify which "
                        f"one you want."
                    )
            else:
                # There were no features at this offset.
                feature = None
        else:
            canonicalName = self.canonicalName(featureName)
            feature = self[canonicalName]
            if canonicalName not in features:
                if features:
                    present = ", ".join(
                        f"{f!r} ({self[f]['start'] + 1} - {self[f]['stop']})"
                        for f in sorted(features)
                    )
                    raise MissingFeatureError(
                        f"Requested feature {featureName!r} (located at sites "
                        f'{feature["start"] + 1}-{feature["stop"]}) does not '
                        f"overlap site {offset + 1}. The feature(s) at "
                        f"that site are: {present}."
                    )
                else:
                    raise MissingFeatureError(
                        f"Feature {featureName!r} (located at sites "
                        f'{feature["start"] + 1}-{feature["stop"]}) '
                        f"does not overlap site {offset + 1}. There are no "
                        f"features at that site."
                    )

        return feature, features

    def toString(
        self,
        name: str,
        maxSequenceLength: int = 80,
        oneBased: bool = True,
        prefix: str = "",
    ) -> str:
        """
        Get a string representation of a feature suitable for printing.

        @param name: A C{str} feature name.
        @param maxSequenceLength: The C{int} maximum sequence (prefix) length to
            include. Pass -1 to specify the full sequence or 0 to exclude sequences
            from the result.
        @param oneBased: A C{bool}. If C{True} the feature location should be printed
            1-based instead of 0-based.
        @param prefix: A C{str} to put at the start of each line.
        @return: A C{str}.
        """
        canonicalName = self.canonicalName(name)
        feature = self[canonicalName]
        sequence = feature["sequence"]
        start = feature["start"]
        stop = feature["stop"]
        result = [
            f"{name}:",
            f"  start: {start + bool(oneBased)}",
            f"  stop: {feature['stop']}",
            f"  length (nt): {stop - start}",
        ]

        for name in "product", "note", "function":
            try:
                result.append(f"  {name}: {feature[name]}")
            except KeyError:
                pass

        if "translation" in feature:
            if "forward" in feature:
                if feature["forward"]:
                    result.append("  feature is translated left-to-right.")
                else:
                    result.append(
                        "  feature is translated right-to-left from the "
                        "reverse complement."
                    )
                    rc = DNARead("id", sequence).reverseComplement().sequence
                    result.append(
                        "  reverse complement: "
                        + (
                            (rc[:maxSequenceLength] + "...")
                            if maxSequenceLength > 0 and len(rc) > maxSequenceLength
                            else rc
                        )
                    )

        if maxSequenceLength:
            result.append(
                "  sequence: "
                + (
                    (sequence[:maxSequenceLength] + "...")
                    if maxSequenceLength > 0 and len(sequence) > maxSequenceLength
                    else sequence
                )
            )

            try:
                translation = feature["translation"]

                result.extend(
                    [
                        f"  length (aa): {len(translation)}",
                        "  translation: "
                        + (
                            (translation[:maxSequenceLength] + "...")
                            if maxSequenceLength > 0
                            and len(translation) > maxSequenceLength
                            else translation
                        ),
                    ]
                )
            except KeyError:
                result.append("  region is not translated.")

        return "\n".join(prefix + s for s in result)


def addFeatureOptions(
    parser: argparse.ArgumentParser, referenceHelpInfo: str = ""
) -> None:
    """
    Add standard command-line options that can then be passed to the Feature
    constructor.

    @param parser: An argparse parser to add options to.
    @param referenceHelpInfo: A C{str} to append to the usage information for the
        --reference option. Note that you need to put a space at the start, if
        you want one.
    """
    parser.add_argument(
        "--reference",
        metavar="file.gb",
        help=f"The GenBank file to read for features and sequences.{referenceHelpInfo}",
    )

    parser.add_argument(
        "--sars2",
        action="store_true",
        help="The sequence is from a SARS-CoV-2 virus.",
    )

    parser.add_argument(
        "--addUnannotatedRegions",
        action="store_true",
        help=(
            "Also add unannotated regions (i.e., genome regions that have "
            'no features). These will be named "unannotated region 1", '
            '"unannotated region 2", and so on, as needed.'
        ),
    )
