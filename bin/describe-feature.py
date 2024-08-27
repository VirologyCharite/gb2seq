#!/usr/bin/env python

import sys
import os
import argparse
from collections import defaultdict

from dark.fasta import FastaReads

from gb2seq.alignment import Gb2Alignment, addAlignerOption
from gb2seq.features import Features, addFeatureOptions, UnknownFeatureNameError
from gb2seq.sars2 import SARS_COV_2_ALIASES


def printNames(features, sars2):
    """
    Print feature names, each with all known aliases (if any).

    @param features: A C{Features} instance.
    @param sars2: C{True} if we have a SARS-CoV-2 GenBank record.
    """

    if sars2:

        def key(name):
            """
            Make a sort key for SARS-CoV-2 feature names.

            @param name: A C{str} feature name.
            @return: A C{str} C{int} 2-tuple for sorting feature names.
            """
            if name.startswith("nsp"):
                return "nsp", int(name[3:])
            elif name.startswith("ORF"):
                return "orf", int(name[3:].split()[0].rstrip("ab"))
            else:
                return name.lower(), 0

    else:

        def key(name):
            return name.lower()

    featureNames = sorted(features, key=key)

    if sars2:
        aka = defaultdict(set)
        for alias, name in SARS_COV_2_ALIASES.items():
            aka[name].add(alias)

        for featureName in featureNames:
            try:
                akas = " " + ", ".join(sorted(aka[featureName]))
            except KeyError:
                akas = ""

            print(f"{featureName}:{akas}")
    else:
        print("\n".join(featureNames))


def reportGenomeFeature(features, name, alignment, maxSequenceLength, oneBased):
    """
    Print details of a feature in a passed genome.

    @param features: A C{Features} instance.
    @param name: A C{str} feature name.
    @param alignment: A C{Gb2seq} instance with the aligned genome whose feature
        should be reported.
    @param maxSequenceLength: The maximum sequence length to print. Longer sequences
        will be truncated. Use 0 or C{None} to skip printing sequences.
    @param oneBased: If true, print one-based sites instead of zero-based offsets.

    """
    print(f"  Genome {alignment.genome.id}:")
    feature = features[name]
    try:
        alignedStart = alignment.gappedOffsets[feature["start"]]
        alignedStop = alignment.gappedOffsets[feature["stop"]]
    except KeyError:
        print(
            f"Offset {feature['stop']} not found in gappedOffsets "
            f"(len {len(alignment.gappedOffsets)})."
        )
        sys.exit(1)

    gappedSequence = alignment.genomeAligned.sequence
    sequence = gappedSequence[alignedStart:alignedStop].replace("-", "")
    absoluteStart = len(gappedSequence[:alignedStart].replace("-", ""))
    absoluteStop = len(gappedSequence[:alignedStop].replace("-", ""))
    _, genomeNt = alignment.ntSequences(name, raiseOnReferenceGaps=False)

    print(f"    start: {absoluteStart + bool(oneBased)}")
    print(f"    stop: {absoluteStop}")
    print(f"    length (nt): {len(genomeNt.sequence)}")
    print(f"    aligned (to ref) start: {alignedStart + bool(oneBased)}")
    print(f"    aligned (to ref) stop: {alignedStop}")

    if maxSequenceLength:
        print(
            "    sequence: " + (genomeNt.sequence[:maxSequenceLength] + "...")
            if maxSequenceLength > 0 and len(sequence) > maxSequenceLength
            else genomeNt.sequence
        )

    if "translation" in feature:
        _, genomeAa = alignment.aaSequences(name, raiseOnReferenceGaps=False)

        if maxSequenceLength:
            aaSequence = genomeAa.sequence.replace("-", "")
            print(
                "    translation: " + (aaSequence[:maxSequenceLength] + "...")
                if maxSequenceLength > 0 and len(aaSequence) > maxSequenceLength
                else aaSequence
            )


def main(args):
    """
    Describe a genome feature or features.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    features = Features(
        args.reference,
        sars2=args.sars2,
        addUnannotatedRegions=args.addUnannotatedRegions,
    )

    if args.names:
        printNames(features, args.sars2)
        return

    wantedNames = []

    if args.name:
        for name in args.name:
            try:
                wantedNames.append(features.canonicalName(name))
            except UnknownFeatureNameError:
                forgot = (
                    " It looks like you forgot to use --sars2."
                    if name.lower() in SARS_COV_2_ALIASES
                    else ""
                )
                print(f"Feature {name!r} is not known.{forgot}", file=sys.stderr)
                sys.exit(1)
    else:
        wantedNames = list(features)
        # All features are wanted, so print a little title.
        print(f"Features for {features.reference.id}:")

    if args.genome or not os.isatty(0):
        fp = open(args.genome) if args.genome else sys.stdin
        alignments = [
            Gb2Alignment(read, features, aligner=args.aligner)
            for read in FastaReads(fp)
        ]
        if args.genome:
            fp.close()
    else:
        alignments = []

    if args.sortBy == "name":
        wantedNames.sort()
    else:
        assert args.sortBy == "site"

        def key(name):
            return features[name]["start"]

        wantedNames.sort(key=key)

    for name in wantedNames:
        print("Reference:")
        print(
            features.toString(
                name,
                maxSequenceLength=args.maxLen,
                oneBased=args.oneBased,
                prefix="  ",
            )
        )

        for alignment in alignments:
            reportGenomeFeature(features, name, alignment, args.maxLen, args.oneBased)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Describe a genome feature in the reference and "
            "(optionally) an unannotated genome."
        ),
    )

    parser.add_argument(
        "--genome",
        metavar="genome.fasta",
        help="The FASTA file containing the genome to examine for features.",
    )

    parser.add_argument(
        "--name",
        metavar="NAME",
        action="append",
        help=(
            "The feature to print information for (all features are "
            "printed if not specified). May be repeated."
        ),
    )

    parser.add_argument(
        "--maxLen",
        type=int,
        default=80,
        help=(
            "The maximum sequence length to print. Longer sequences will "
            "be truncated."
        ),
    )

    parser.add_argument(
        "--names",
        action="store_true",
        help="Only print feature names and aliases (if any).",
    )

    parser.add_argument(
        "--zeroBased",
        dest="oneBased",
        action="store_false",
        help="Print zero-based offsets instead of one-based sites.",
    )

    parser.add_argument(
        "--sortBy",
        choices=("name", "site"),
        default="name",
        help=(
            "The order in which to print features (default is the order given "
            "on the command line)."
        ),
    )

    addAlignerOption(parser)
    addFeatureOptions(parser)

    args = parser.parse_args()

    main(args)
