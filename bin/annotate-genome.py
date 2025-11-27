#!/usr/bin/env python

import sys
import argparse
from json import dumps

from dark.fasta import FastaReads

from gb2seq.alignment import addAlignerOption
from gb2seq.annotate import annotateGenome, summarizeDifferences
from gb2seq.features import Features, addFeatureOptions


def main(args: argparse.Namespace) -> int:
    """
    Create, optionally summarize (to standard error), and print (as JSON, to
    standard out) annotations for a genome.

    @param args: An C{argparse.Namespace} instance with command-line values.
    """
    features = Features(args.reference, sars2=args.sars2, alsoInclude={"gene"})
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
        if args.format == "json":
            print(dumps(annotations, sort_keys=True, indent=4))
        else:
            assert args.format == "feature-table"
            print(f">Feature {genome.id}")

            for feature in annotations["features"].values():
                start = feature["genome"]["start"] + 1
                stop = feature["genome"]["stop"]

                reference = feature["reference"]
                type_ = reference["type"]
                note = reference.get("note")
                product = reference.get("product") or reference.get("name")

                print(f"{start}\t{stop}\t{type_}")

                if type_ == "gene":
                    print(f"\t\t\tgene\t{product}")
                else:
                    for what, value in (
                        ("note", note),
                        ("product", product),
                    ):
                        if value:
                            print(f"\t\t\t{what}\t{value}")

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
        type=argparse.FileType("rb"),
        default=sys.stdin,
        help="The FASTA file containing the SARS-CoV-2 genome(s) to examine.",
    )

    parser.add_argument(
        "--summarize",
        action="store_true",
        help="Write a summary of differences to standard error.",
    )

    parser.add_argument(
        "--format",
        choices=("json", "feature-table"),
        default="json",
        help=(
            "The output format. The feature-table format is described at "
            "https://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html"
        ),
    )

    addAlignerOption(parser)
    addFeatureOptions(parser)

    return parser


if __name__ == "__main__":
    sys.exit(main(makeParser().parse_args()))
