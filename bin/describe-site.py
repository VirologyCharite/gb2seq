#!/usr/bin/env python

import sys
import os
from json import dumps
import argparse

from dark.fasta import FastaReads

from gb2seq import Gb2SeqError
from gb2seq.alignment import Gb2Alignment, addAlignerOption
from gb2seq.features import Features, addFeatureOptions


def report(genome, args, includeGenome=True):
    """
    Report what's found at a site for a given genome (or report insufficient
    coverage to standard error).

    @param genome: A C{Gb2Alignment} instance.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    @param includeGenome: If C{True}, include information about the genome
        (not just the reference).
    """
    try:
        offsetInfo = genome.offsetInfo(
            args.site - 1,
            relativeToFeature=args.relativeToFeature,
            aa=args.aa,
            featureName=args.feature,
            includeUntranslated=args.includeUntranslated,
            minReferenceCoverage=args.minReferenceCoverage,
        )
    except Gb2SeqError as e:
        print(e, file=sys.stderr)
        sys.exit(1)

    if offsetInfo is None:
        what = args.featureName or "genome"
        print(f"Insufficient {what} coverage", file=sys.stderr)
        return

    if args.genomeAaOnly:
        print(offsetInfo["genome"]["aa"])
    else:
        offsetInc = int(bool(args.oneBased))

        offsetInfo["alignmentOffset"] += offsetInc
        offsetInfo["reference"]["ntOffset"] += offsetInc

        if includeGenome:
            offsetInfo["genome"]["ntOffset"] += offsetInc
        else:
            del offsetInfo["genome"]

        if args.includeFeature:
            featureName = offsetInfo["featureName"]
            if featureName:
                assert "feature" not in offsetInfo
                offsetInfo["feature"] = genome.features[featureName]

        # TODO: what should we print if the user doesn't want JSON? Some kind
        # of textual summary, I guess. When that's implemented, remove the
        # "or True" below.
        if args.json or True:
            # Make the featureNames into a sorted list (it is by default a
            # set), so it can be printed as JSON.
            offsetInfo["featureNames"] = sorted(offsetInfo["featureNames"])
            print(dumps(offsetInfo, indent=4, sort_keys=True))
        else:
            print(offsetInfo)


def main(args):
    """
    Describe a site in a genome or genomes.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    @return: An C{int} exit status.
    """
    features = Features(
        args.reference,
        sars2=args.sars2,
        addUnannotatedRegions=args.addUnannotatedRegions,
    )

    if args.genome is None and os.isatty(0):
        alignment = Gb2Alignment(features.reference, features, aligner=args.aligner)
        report(alignment, args, False)
    else:
        count = 0
        fp = open(args.genome) if args.genome else sys.stdin
        for count, read in enumerate(FastaReads(fp), start=1):
            alignment = Gb2Alignment(read, features, aligner=args.aligner)
            report(alignment, args)

        if args.verbose:
            print(f"Examined {count} genomes.", file=sys.stderr)

        if args.genome:
            fp.close()

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Describe a site of a genome(s).",
    )

    parser.add_argument(
        "--genome",
        metavar="file.fasta",
        help="The FASTA file containing the genome(s) to examine.",
    )

    parser.add_argument(
        "--site",
        metavar="N",
        type=int,
        required=True,
        help="The (1-based) site to find information for.",
    )

    parser.add_argument(
        "--feature",
        metavar="FEATURE",
        help=(
            "The feature to examine (e.g., nsp2). This is required if you "
            "use --aa or --relativeToFeature"
        ),
    )

    parser.add_argument(
        "--includeUntranslated",
        action="store_true",
        help=(
            "Include untranslated features (if no feature name is given and "
            "it is necessary to identify the intended feature just based on "
            "offset)."
        ),
    )

    parser.add_argument(
        "--aa",
        action="store_true",
        help=("The given site is an amino acid count (the default is " "nucleotides)."),
    )

    parser.add_argument(
        "--relativeToFeature",
        action="store_true",
        help="The given site is relative to the start of the feature.",
    )

    parser.add_argument("--json", action="store_true", help="Print the result as JSON.")

    parser.add_argument(
        "--genomeAaOnly",
        action="store_true",
        help="Only print the amino acid from the genome.",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print information about proceesing to standard error.",
    )

    parser.add_argument(
        "--includeFeature",
        action="store_true",
        help="Also print information about the feature at the site.",
    )

    parser.add_argument(
        "--minReferenceCoverage",
        metavar="coverage",
        type=float,
        help=(
            "The fraction of non-N bases required in the genome(s) in order "
            "for them to be processed. Genomes with lower coverage will be "
            "ignored, with a message printed to standard error. Note that "
            "the denominator used to compute the coverage fraction is the "
            "length of the reference. I.e., coverage is computed as number "
            "of non-N bases in the genome divided by the length of the "
            "reference."
        ),
    )

    parser.add_argument(
        "--zeroBased",
        dest="oneBased",
        action="store_false",
        help="Print zero-based offsets instead of one-based sites.",
    )

    addFeatureOptions(parser)
    addAlignerOption(parser)

    args = parser.parse_args()

    sys.exit(main(args))
