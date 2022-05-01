import sys
from os import environ
import itertools
from collections import UserDict
from pathlib import Path

from Bio import Entrez, SeqIO

from dark.aa import STOP_CODONS
from dark.reads import DNARead

from sars2seq import Sars2SeqError, DATA_DIR

# Set ENTREZ_EMAIL in your environment to have your requests to NCBI Entez
# be accompanied by your address. If you don't do this you'll see warning
# messages and be limited to a lower rate of querying.
Entrez.email = environ.get("ENTREZ_EMAIL")


class ReferenceWithGapError(Sars2SeqError):
    "A GenBank reference sequence had a gap."


class MissingFeatureError(Sars2SeqError):
    "A feature expected at an offset is not present."


class AmbiguousFeatureError(Sars2SeqError):
    "More than one feature is referred to by an offset."


# Provide convenient aliases for feature names. The alias is the key, the
# canonical name (as found in the GenBank file) is the value.
#
# Alphanumeric feature aliases must have lower case keys. If not they will not
# be detected (and the test suite will fail).
#
# At some point we might want to make it possible to pass a custom set of
# aliases to the Features constructor.

ALIASES = {
    "2": "2'-O-ribose methyltransferase",
    "3clpro": "3C-like proteinase",
    "3utr": "3'UTR",
    "5utr": "5'UTR",
    "e": "envelope protein",
    "endornase": "endoRNAse",
    "envelope": "envelope protein",
    "exon": "3'-to-5' exonuclease",
    "exonuclease": "3'-to-5' exonuclease",
    "leader": "leader protein",
    "m": "membrane glycoprotein",
    "membrane": "membrane glycoprotein",
    "mpro": "3C-like proteinase",
    "n": "nucleocapsid phosphoprotein",
    "nsp1": "leader protein",
    "nsp5": "3C-like proteinase",
    "nsp12": "RNA-dependent RNA polymerase",
    "nsp13": "helicase",
    "nsp14": "3'-to-5' exonuclease",
    "nsp15": "endoRNAse",
    "orf4": "envelope protein",
    "orf5": "membrane glycoprotein",
    "orf1a": "ORF1a polyprotein",
    "orf1ab": "ORF1ab polyprotein",
    "orf3a": "ORF3a protein",
    "orf6": "ORF6 protein",
    "orf7a": "ORF7a protein",
    "orf7b": "ORF7b",
    "orf8": "ORF8 protein",
    "orf9": "nucleocapsid phosphoprotein",
    "orf10": "ORF10 protein",
    "rdrp": "RNA-dependent RNA polymerase",
    "s": "surface glycoprotein",
    "sl1": "stem loop 1",
    "sl2": "stem loop 2",
    "sl3": "stem loop 3",
    "sl4": "stem loop 4",
    "sl5": "stem loop 5",
    "s": "surface glycoprotein",
    "spike": "surface glycoprotein",
}

# Name of translated features, with (case sensitive!) names matching those in
# the GenBank file ../data/NC_045512.2.gb
#
# At some point we might want to make it possible to pass a custom set of
# translated names to the Features constructor.

TRANSLATED = {
    "3'-to-5' exonuclease",  # nsp14
    "3C-like proteinase",  # nsp5
    "endoRNAse",  # nsp15
    "envelope protein",  # ORF4
    "helicase",  # nsp13
    "leader protein",  # nsp1
    "membrane glycoprotein",  # ORF5
    "nsp2",
    "nsp3",
    "nsp4",
    "nsp6",
    "nsp7",
    "nsp8",
    "nsp9",
    "nsp10",
    "nsp11",
    "nucleocapsid phosphoprotein",  # ORF9
    "ORF1a polyprotein",
    "ORF1ab polyprotein",
    "ORF3a protein",
    "ORF6 protein",
    "ORF7a protein",
    "ORF7b",
    "ORF8 protein",
    "ORF10 protein",
    "RNA-dependent RNA polymerase",  # nsp12
    "surface glycoprotein",
}


class Features(UserDict):
    """
    Manage sequence features from the information in C{gbFile}.

    @param spec: Either:
        * A C{str} name or C{Path} of a GenBank file containing the features.
        * A C{str} GenBank accession id.
        * A C{dict} of pre-prepared features, in which case C{reference}
              must not be C{None}. Passing a C{dict} is provided for testing.
        * C{None}, in which case the default reference, NC_045512.2.gb, is
              loaded.
    @param reference: A C{dark.reads.DNARead} instance if C{spec} is a C{dict},
        else C{None}.
    @raise ValueError: If a reference is passed with a string or Path
        specification.
    @raise ReferenceWithGapError: If the reference or one of its features has a
        gap in its nucleotide sequence.
    """

    REF_GB = DATA_DIR / "NC_045512.2.gb"

    def __init__(self, spec=None, reference=None):
        super().__init__()
        spec = self.REF_GB if spec is None else spec

        if isinstance(spec, str):
            if reference is not None:
                raise ValueError(
                    "A reference cannot be passed with a string " "specification."
                )
            path = Path(spec)
            if path.exists():
                with open(path) as fp:
                    record = SeqIO.read(fp, "genbank")
            else:
                print(f"Fetching GenBank record for {spec!r}.", file=sys.stderr)
                client = Entrez.efetch(
                    db="nucleotide", rettype="gb", retmode="text", id=spec
                )
                record = SeqIO.read(client, "gb")
                client.close()
            self._initializeFromGenBankRecord(record)
        elif isinstance(spec, Path):
            if reference is not None:
                raise ValueError(
                    "A reference cannot be passed with a Path " "specification."
                )
            with open(spec) as fp:
                record = SeqIO.read(fp, "genbank")
            self._initializeFromGenBankRecord(record)
        elif isinstance(spec, dict):
            self.update(spec)
            self.reference = reference
        else:
            raise ValueError(f"Unrecognized specification {spec!r}.")

    def _initializeFromGenBankRecord(self, record):
        """
        Initialize from a GenBank record.

        @param record: A BioPython C{SeqRecord} sequence record.
        @raise ReferenceWithGapError: If the reference sequence or the
            sequence of any of its features has a gap.
        """
        self.reference = DNARead(record.id, str(record.seq))

        for feature in record.features:
            type_ = feature.type
            start = int(feature.location.start)
            stop = int(feature.location.end)
            value = {}

            if type_ == "3'UTR" or type_ == "5'UTR":
                name = type_

            elif type_ == "stem_loop":
                for n in itertools.count(1):
                    name = f"stem loop {n}"
                    if name not in self:
                        break
                value["function"] = feature.qualifiers["function"][0]

            elif type_ in {"CDS", "mat_peptide"}:
                name = feature.qualifiers["product"][0]
                value["product"] = name

            elif type_ in {"source", "gene"}:
                continue

            else:
                raise ValueError(f"Unknown feature type {type_!r}.")

            value.update(
                {
                    "name": name,
                    "sequence": str(record.seq)[start:stop],
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
                if (
                    not translation.endswith("*")
                    and value["sequence"][-3:].upper() in STOP_CODONS
                ):
                    value["translation"] += "*"

            if name in self:
                assert self[name] == value
            else:
                self[name] = value

        self._checkForGaps()

    def __getitem__(self, name):
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
        @raise KeyError: If the name is unknown.
        @return: A C{dict} for the feature, as above.
        """
        return self.data[self.canonicalName(name)]

    def canonicalName(self, name):
        """
        Get the canonical name for a feature.

        @param name: A C{str} feature name to look up.
        @raise KeyError: If the name is unknown.
        @return: A C{str} canonical name.
        """
        if name in self:
            return name

        nameLower = name.lower()
        for featureName in self:
            if nameLower == featureName.lower():
                return featureName

        alias = ALIASES.get(nameLower)
        if alias is None:
            raise KeyError(name)
        else:
            assert alias in self
            return alias

    def aliases(self, name):
        """
        Get all aliases for a name.

        @param name: A C{str} feature name.
        @return: A C{set} of C{str} canonical names.
        """
        canonicalName = self.canonicalName(name)
        result = {canonicalName}

        for alias, canonical in ALIASES.items():
            if canonical == canonicalName:
                result.add(alias)

        return result

    def _checkForGaps(self):
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

    def referenceOffset(self, name, offset, aa=False):
        """
        Get the (nucleotide) offset in the reference, given an offset in a
        feature.

        @param name: A C{str} feature name.
        @param offset: An C{int} offset into the feature.
        @param aa: If C{True}, the offset is a number of amino acids, else a
            number of nucleotides.
        @raise KeyError: If the name is unknown.
        @return: An C{int} nucleotide offset into the reference genome.
        """
        return self[name]["start"] + offset * (3 if aa else 1)

    def getFeatureNames(self, offset, includeUntranslated=False):
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
                name in TRANSLATED or includeUntranslated
            ):
                result.add(name)

        return result

    def getFeature(self, offset, featureName=None, includeUntranslated=False):
        """
        Find a single feature at an offset.

        @param offset: An C{int} offset into the genome.
        @param featureName: A C{str} feature name. Used for disambiguation when
            multiple features are present at the offset. If C{None}, a feature
            will be returned if there is only one at the given offset.
        @param includeUntranslated: If C{True}, also consider features that are
            not translated.
        @raise MissingFeatureError: if the requested feature does not overlap
            the given C{offset} or if there are no features at the offset.
        @raise AmbiguousFeatureError: if there are multiple features at the
            offset and a feature name is not given.
        @return: A 2-tuple, containing 1) a features C{dict} (as returned by
            __getitem__), if the requested feature is among those present at
            the offset, and 2) the set of all C{str} feature names present at
            the offset.
        """
        features = self.getFeatureNames(offset, includeUntranslated)

        if featureName is None:
            if features:
                # There are some features here, but we weren't told which one
                # to use. Only proceed if there's just one.
                if len(features) == 1:
                    featureName = list(features)[0]
                    feature = self[featureName]
                else:
                    present = ", ".join(f"{f!r}" for f in sorted(features))
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
                    present = ", ".join(f"{f!r}" for f in sorted(features))
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
