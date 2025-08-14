#!/usr/bin/env python

import sys
import argparse
import rich

from dark.aligners import align
from dark.dna import compareDNAReads, matchToString
from dark.fasta import FastaReads
from dark.reads import Read, Reads

from gb2seq.alignment import Gb2Alignment, addAlignerOption
from gb2seq.compare import compare, print_results
from gb2seq.features import Features, addFeatureOptions


def parse_args() -> argparse.Namespace:
    """
    Make an argument parser and use it to parse the command line.
    """
    parser = argparse.ArgumentParser(
        description="Compare two genome sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--genome-1",
        "-1",
        required=True,
        help="The FASTA file containing the first genome sequence.",
    )

    parser.add_argument(
        "--genome-2",
        "-2",
        required=True,
        help="The FASTA file containing the second genome sequence.",
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
            "reference."
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

    return parser.parse_args()


def simple_compare(
    a_read: Read,
    b_read: Read,
    aligner: str,
) -> None:
    a, b = list(align(Reads((a_read, b_read)), aligner=aligner))
    match = compareDNAReads(a, b)
    print(matchToString(match, a, b))


def main() -> int:
    args = parse_args()

    if not (args.reference or args.sars2):
        print(
            "You must either supply a GenBank reference file (via --reference) or "
            "use --sars2 for matching against the original SARS-CoV-2 Wuhan "
            "reference).",
            file=sys.stderr,
        )

        return 1

    if args.force_terminal:
        args.rich = True
        rich.reconfigure(force_terminal=True)

    features = Features(
        sars2=args.sars2, addUnannotatedRegions=args.addUnannotatedRegions
    )
    assert features.reference
    (a_read,) = list(FastaReads(args.genome_1))

    a_alignment = Gb2Alignment(a_read, features, aligner=args.aligner)

    if args.alignment_file:
        reads = Reads([features.reference, a_read, a_alignment.genomeAligned])
        reads.save(args.alignment_file)

    for b_read in FastaReads(args.genome_2):
        b_alignment = Gb2Alignment(b_read, features, aligner=args.aligner)

        if args.simple_compare:
            simple_compare(a_read, b_read, args.aligner)

        if args.alignment_file:
            reads = Reads([b_read, b_alignment.genomeAligned])
            reads.save(args.alignment_file, mode="a")

        changes = compare(
            a_read,
            b_read,
            a_alignment,
            b_alignment,
            features,
            args.verbose,
            True,
            True,
        )

        print_results(changes, len(features.reference), args.verbose, args.rich)

    return 0


if __name__ == "__main__":
    sys.exit(main())
