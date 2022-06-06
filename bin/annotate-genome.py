#!/usr/bin/env python

import sys
import argparse
from json import dumps

from dark.fasta import FastaReads

from sars2seq.alignment import addAlignerOption
from sars2seq.annotate import annotateGenome, summarizeDifferences
from sars2seq.features import Features, addFeatureOptions


def main(args: argparse.Namespace) -> int:
    """
    Create, optionally summarize (to standard error), and print (as JSON, to
    standard out) annotations for a genome.

    @param args: An C{argparse.Namespace} instance with command-line values.
    """
    features = Features(args.reference, sars2=args.sars2)
    genome = list(FastaReads(args.genome))[0]

    try:
        annotations = annotateGenome(features, genome, args.aligner)
    except Exception as e:
        print(
            f"Failed to annotate {genome.id!r} from {args.reference!r}: {e}",
            file=sys.stderr,
        )
        return 1
    else:
        print(dumps(annotations, sort_keys=True, indent=4))

        if args.summarize:
            print(summarizeDifferences(annotations), file=sys.stderr)

        return 0


def makeParser() -> argparse.ArgumentParser:
    """
    Make the command-line parser.

    @return: An C{argparse.ArgumentParser}.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Print JSON annotations for a genome (or genomes).",
    )

    parser.add_argument(
        "--genome",
        metavar="file.fasta",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="The FASTA file containing the SARS-CoV-2 genome(s) to examine.",
    )

    parser.add_argument(
        "--summarize",
        action="store_true",
        help="Write a summary of differences to standard error.",
    )

    addAlignerOption(parser)
    addFeatureOptions(parser)

    return parser


if __name__ == "__main__":
    sys.exit(main(makeParser().parse_args()))
