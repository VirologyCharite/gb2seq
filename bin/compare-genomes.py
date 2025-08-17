#!/usr/bin/env python

import sys
import argparse
import rich
from typing import TextIO

from dark.aligners import align
from dark.dna import compareDNAReads, matchToString
from dark.reads import (
    Read,
    Reads,
    addFASTACommandLineOptions,
    parseFASTACommandLineOptions,
)

from gb2seq.alignment import Gb2Alignment, addAlignerOption
from gb2seq.compare import compare, print_results
from gb2seq.features import Features, addFeatureOptions


def parse_args() -> argparse.Namespace:
    """
    Make an argument parser and use it to parse the command line.
    """
    parser = argparse.ArgumentParser(
        description="Compare genome sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--alignment-file",
        help=(
            "The (optional) filename to save the aligned genomes to as FASTA. "
            "The file will contain 1) the reference; 2) the sequence from the "
            "genome-1 file; 3) the reference as aligned to the sequence from "
            "the genome-1 file; 4) the sequence from the genome-1 file as "
            "aligned to the reference; and then, for each sequence in the "
            "genome-2 file, the sequence and then the sequence aligned to the "
            "reference. Note that the sequence identifiers will be annotated with "
            "prefixes to help make it clear what's what in the file."
        ),
    )

    parser.add_argument(
        "--rich",
        action="store_true",
        help="Use rich to produce colored output.",
    )

    parser.add_argument(
        "--force-terminal",
        "-R",
        action="store_true",
        help=(
            "Force the 'rich' library to produce colored output to override its "
            "logic about when to color output, e.g., for piping output into "
            "'less -R' (implies --rich)."
        ),
    )

    parser.add_argument(
        "--simple-compare",
        action="store_true",
        help=(
            "As well as doing the substitution analysis, print a high-level sequence "
            "comparision (e.g., percent identity, number of mismatches, gaps, etc.)."
        ),
    )

    parser.add_argument(
        "--sars2-overlapping-features",
        action="store_false",
        dest="sars2_simplify_features",
        help=(
            "When processing SARS-CoV-2 genomes, it's usually not useful to print "
            "identical substitution details for features that overlap one another "
            "(e.g., for NSP proteins found in ORF1a/b or the stem loops found in the "
            "3'UTR). The default behavior is to omit output for the longer feature "
            "so only the more specific information is printed. Use this option if you "
            "are processing SARS-CoV-2 genomes and you want all overlapping features "
            "to be printed."
        ),
    )

    parser.add_argument(
        "--verbose",
        type=int,
        default=0,
        help=(
            "Print extra information. Use 1 for a summary of sites, features, "
            "and synonymous/non-synonymous changes, 2 to also print information "
            "about offsets, frames, and bases involved in substitutions."
        ),
    )

    addAlignerOption(parser)
    addFeatureOptions(parser)
    addFASTACommandLineOptions(parser)

    return parser.parse_args()


def simple_compare(
    a_read: Read,
    b_read: Read,
    aligner: str,
) -> None:
    """
    Align two sequences and print a simple identity comparison.
    """
    a, b = list(align(Reads((a_read, b_read)), aligner=aligner))
    print(matchToString(compareDNAReads(a, b), a, b))


def save_alignment(alignment: Gb2Alignment, count: int, fp: TextIO) -> None:
    """
    Save a reference and genome from before and after an alignment, adding
    prefixes to the sequence IDs to make it clear what's what.
    """
    reads = Reads(
        [
            Read(
                f"REFERENCE: {alignment.features.reference.id}",
                alignment.features.reference.sequence,
            ),
            Read(f"GENOME {count}: {alignment.genome.id}", alignment.genome.sequence),
            Read(
                f"REFERENCE ALIGNED: {alignment.referenceAligned.id}",
                alignment.referenceAligned.sequence,
            ),
            Read(
                f"GENOME {count} ALIGNED: {alignment.genomeAligned.id}",
                alignment.genomeAligned.sequence,
            ),
        ]
    )
    reads.save(fp)


def print_title(a_read: Read, b_read: Read, use_rich: bool):
    if use_rich:
        color = "bold red"
        rich.print(
            f"Summary of changes from [{color}]{a_read.id!r}[/] to "
            f"[{color}]{b_read.id!r}[/]"
        )
    else:
        print(f"Summary of changes from {a_read.id!r} to {b_read.id!r}")


def main():
    args = parse_args()

    if not (args.reference or args.sars2):
        sys.exit(
            "You must either supply a GenBank reference file (via --reference) or "
            "use --sars2 for matching against the original SARS-CoV-2 Wuhan "
            "reference.",
        )

    if not args.sars2_simplify_features and not args.sars2:
        args.sars2 = True
        print(
            "Assuming --sars2 because you used --sars2-overlapping-features. "
            "Also specify --sars2 to silence this message.",
            file=sys.stderr,
        )

    if args.force_terminal:
        args.rich = True
        rich.reconfigure(force_terminal=True)

    features = Features(
        referenceSpecification=args.reference,
        sars2=args.sars2,
        addUnannotatedRegions=args.addUnannotatedRegions,
        alsoInclude=args.alsoInclude,
    )
    assert features.reference

    reads_and_alignments = [
        (read, Gb2Alignment(read, features, aligner=args.aligner))
        for read in parseFASTACommandLineOptions(args)
    ]

    n = len(reads_and_alignments)

    if n == 0:
        sys.exit("Your input had no FASTA sequences!")
    elif n == 1:
        sys.exit("Your input had only one FASTA sequence!")

    if args.alignment_file:
        with open(args.alignment_file, "w") as fp:
            for index, (_, alignment) in enumerate(reads_and_alignments):
                save_alignment(alignment, index + 1, fp)

    for i in range(n):
        a_read, a_alignment = reads_and_alignments[i]

        for j in range(i + 1, n):
            b_read, b_alignment = reads_and_alignments[j]

            print_title(a_read, b_read, args.rich)

            if args.simple_compare:
                simple_compare(a_read, b_read, args.aligner)

            changes = compare(
                a_read,
                b_read,
                a_alignment,
                b_alignment,
                features,
                args.verbose,
                args.sars2,
                args.sars2_simplify_features,
            )

            print_results(changes, len(features.reference), args.verbose, args.rich)


if __name__ == "__main__":
    main()
