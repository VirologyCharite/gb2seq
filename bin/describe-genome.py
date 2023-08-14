#!/usr/bin/env python

import sys
import os
import argparse
from json import load
from math import log10
from os.path import exists, join
from contextlib import contextmanager

from dark.fasta import FastaReads
from dark.aa import compareAaReads, matchToString as aaMatchToString
from dark.dna import compareDNAReads, matchToString as dnaMatchToString
from dark.reads import Read, Reads

from gb2seq.alignment import Gb2Alignment, addAlignerOption
from gb2seq.features import Features, addFeatureOptions, UnknownFeatureNameError
from gb2seq.translate import TranslationError
from gb2seq.variants import VARIANTS


@contextmanager
def genomeFilePointer(read, args, suffix):
    """
    Open a file whose name is derived from the reference read, the feature
    being examined and whether nucelotides or amino acids are involved.

    @param read: The C{dark.reads.Read} that was read from the input FASTA file
        (this is the overall genome from which the feature was obtained).
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    @param suffix: The C{str} suffix for the filename.
    @return: This is a context manager decorator, so it yields a file pointer
        and then closes it.
    """
    if args.outDir:
        prefix = read.id.split()[0].replace("/", "_")
        filename = join(args.outDir, f"{prefix}{suffix}")
        with open(filename, "w") as fp:
            yield fp
    else:
        yield sys.stdout


@contextmanager
def featureFilePointers(read, feature, args=None):
    """
    Return a dictionary of file pointers for output streams on a per-read
    (i.e., per input genome) basis.

    @param read: The C{dark.reads.Read} that was read from the input FASTA file
        (this is the overall genome from which the feature was obtained).
    @param feature: The C{str} name of the feature.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """

    def fp(suffix, nt):
        """
        Get a file pointer.

        @param suffix: A C{str} file name suffix.
        @param nt: If C{True} the sequences are nucelotide, else protein.
        @return: An file pointer open for writing.
        """
        if args.outDir:
            prefix = (
                (read.id.split()[0] + "-" + feature)
                .replace("/", args.slashReplacement)
                .replace(" ", args.spaceReplacement)
            )
            filename = join(args.outDir, f"{prefix}{('-nt' if nt else '-aa')}{suffix}")
            return open(filename, "w")
        else:
            return sys.stdout

    fps = {}

    try:
        if args.printNtMatch:
            fps["nt-match"] = fp("-match.txt", True)
        if args.printAaMatch:
            fps["aa-match"] = fp("-match.txt", False)

        if args.printNtSequence:
            fps["nt-sequence"] = fp("-sequence.fasta", True)
        if args.printAaSequence:
            fps["aa-sequence"] = fp("-sequence.fasta", False)

        if args.printNtAlignment:
            fps["nt-align"] = fp("-align.fasta", True)
        if args.printAaAlignment:
            fps["aa-align"] = fp("-align.fasta", False)

        yield fps

    finally:
        if args.outDir:
            for fp in fps.values():
                fp.close()


def printDiffs(read1, read2, nt, referenceOffset, fp, indent=""):
    """
    Print differences between sequences.

    @param read1: A C{dark.reads.Read} instance.
    @param read2: A C{dark.reads.Read} instance.
    @param nt: If C{True} the sequences are nucelotide, else protein.
    @param referenceOffset: The C{int} 0-based offset of the feature in the
        reference.
    @param indent: A C{str} prefix for each output line.
    """
    len1, len2 = len(read1), len(read2)
    width = int(log10(max(len1, len2))) + 1
    headerPrinted = False
    multiplier = 1 if nt else 3
    what = "nt" if nt else "aa"
    header = "%sDifferences: site, %s1, %s2, ref nt %s" % (
        indent,
        what,
        what,
        "site" if nt else "codon start",
    )

    for site, (a, b) in enumerate(zip(read1.sequence, read2.sequence)):
        if a != b:
            if not headerPrinted:
                print(header, file=fp)
                headerPrinted = True
            print(
                "%s  %*d %s %s %5d"
                % (
                    indent,
                    width,
                    site + 1,
                    a,
                    b,
                    referenceOffset + (multiplier * site) + 1,
                ),
                file=fp,
            )


def printVariantSummary(genome, fp, args):
    """
    Print a summary of whether the genome fulfils the various variant
    properties.

    @param genome: A C{Gb2Alignment} instance.
    @param fp: An open file pointer to write to.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    print("Variant summary:", file=fp)
    for variant in args.checkVariant:
        testCount, errorCount, tests = genome.checkVariant(
            variant, args.onError, sys.stderr
        )
        successCount = testCount - errorCount
        print(f'  {VARIANTS[variant]["description"]}:', file=fp)
        print(f"  {testCount} checks, {successCount} passed.", file=fp)
        for feature in tests:
            for type_ in tests[feature]:
                passed = set()
                failed = set()
                for change, (_, _, genOK, _) in tests[feature][type_].items():
                    if genOK:
                        passed.add(change)
                    else:
                        failed.add(change)
                print(f"    {feature} {type_}:", file=fp, end="")
                if passed:
                    print(" PASS:", ", ".join(sorted(passed)), file=fp, end="")
                if failed:
                    print(" FAIL:", ", ".join(sorted(failed)), file=fp, end="")
                print(file=fp)


def processFeature(featureName, genome, fps, featureNumber, args):
    """
    Process a feature from a genome.

    @param featureName: A C{str} feature name.
    @param genome: A C{Gb2Alignment} instance.
    @param fps: A C{dict} of file pointers for the various output streams.
    @param featureNumber: The C{int} 0-based count of the features requested.
        This will be zero for the first feature, 1 for the second, etc.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    @raises UnknownFeatureNameError: If the feature is unknown.
    """
    referenceNt, genomeNt = genome.ntSequences(featureName)
    feature = genome.features[featureName]

    if args.printAaMatch or args.printAaSequence or args.printAaAlignment:
        try:
            referenceAa, genomeAa = genome.aaSequences(featureName)
        except TranslationError as e:
            if args.onError == "raise":
                raise
            elif args.onError == "print":
                print(
                    f"Could not translate feature {featureName} in genome "
                    f"{genome.genome.id}: {e}",
                    file=sys.stderr,
                )
            referenceAa = genomeAa = None

    newlineNeeded = False

    if args.printNtMatch:
        fp = fps["nt-match"]
        if featureNumber:
            print(file=fp)
        print(f"Feature: {featureName} nucleotide match", file=fp)
        print(f'  Reference nt location {feature["start"] + 1}', file=fp)
        match = compareDNAReads(referenceNt, genomeNt)
        print(
            dnaMatchToString(
                match, referenceNt, genomeNt, matchAmbiguous=False, indent="  "
            ),
            file=fp,
        )
        printDiffs(referenceNt, genomeNt, True, feature["start"], fp, indent="    ")
        newlineNeeded = True

    if args.printAaMatch and genomeAa:
        fp = fps["aa-match"]
        if newlineNeeded or featureNumber:
            print(file=fp)
        print(f"Feature: {featureName} amino acid match", file=fp)
        match = compareAaReads(referenceAa, genomeAa)
        print(aaMatchToString(match, referenceAa, genomeAa, indent="  "), file=fp)
        printDiffs(referenceAa, genomeAa, False, feature["start"], fp, indent="    ")

    if args.printNtSequence:
        noGaps = Read(genomeNt.id, genomeNt.sequence.replace("-", ""))
        Reads([noGaps]).save(fps["nt-sequence"])

    if args.printAaSequence and genomeAa:
        noGaps = Read(genomeAa.id, genomeAa.sequence.replace("-", ""))
        Reads([noGaps]).save(fps["aa-sequence"])

    if args.printNtAlignment:
        Reads([genomeNt, referenceNt]).save(fps["nt-align"])

    if args.printAaAlignment and genomeAa:
        Reads([genomeAa, referenceAa]).save(fps["aa-align"])


def main(args):
    """
    Describe a genome or genomes.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    @return: An C{int} exit status.
    """
    outDir = args.outDir
    if outDir:
        if not exists(outDir):
            os.makedirs(outDir)

    features = Features(
        args.reference,
        sars2=args.sars2,
        addUnannotatedRegions=args.addUnannotatedRegions,
    )

    if args.feature:
        unknownFeatureNames = set()
        for featureName in args.feature:
            try:
                features[featureName]
            except UnknownFeatureNameError:
                unknownFeatureNames.add(featureName)

        if unknownFeatureNames:
            if len(unknownFeatureNames) == 1:
                print(
                    f"Feature name {unknownFeatureNames.pop()!r} is unknown.",
                    file=sys.stderr,
                )
            else:
                unknownStr = ", ".join(
                    f"{featureName!r}" for featureName in sorted(unknownFeatureNames)
                )
                print(
                    f"The following {len(unknownFeatureNames)} features are "
                    f"unknown: {unknownStr}.",
                    file=sys.stderr,
                )
            sys.exit(1)
        else:
            wantedFeatures = args.feature
    else:
        if args.noFeatures:
            wantedFeatures = []
        else:
            wantedFeatures = sorted(features)

    if not (args.checkVariant or wantedFeatures):
        print("No action specified - I have nothing to do!", file=sys.stderr)
        return 1

    if args.canonicalNames:
        wantedFeatures = map(features.canonicalName, wantedFeatures)

    if args.variantFile:
        try:
            VARIANTS.update(load(open(args.variantFile, encoding="utf-8")))
        except Exception as e:
            print(
                f"Could not parse variant JSON in {args.variantFile!r}: {e}",
                file=sys.stderr,
            )
            sys.exit(1)

    count = ignoredDueToCoverageCount = 0

    for count, read in enumerate(FastaReads(args.genome), start=1):
        if args.minReferenceCoverage is not None:
            coverage = (len(read) - read.sequence.upper().count("N")) / len(
                features.reference
            )
            if coverage < args.minReferenceCoverage:
                ignoredDueToCoverageCount += 1
                print(
                    f"Genome {read.id!r} ignored due to low "
                    f"({coverage * 100.0:.2f}%) coverage of the reference.",
                    file=sys.stderr,
                )
                continue

        alignment = Gb2Alignment(read, features, aligner=args.aligner)

        if args.checkVariant:
            with genomeFilePointer(read, args, "-variant-summary.txt") as fp:
                print(read.id, file=fp)
                printVariantSummary(alignment, fp, args)

        for i, featureName in enumerate(wantedFeatures):
            with featureFilePointers(read, featureName, args) as fps:
                processFeature(featureName, alignment, fps, i, args)

    if not args.quiet:
        print(f'Examined {count} genome{"" if count == 1 else "s"}.', file=sys.stderr)

        if args.minReferenceCoverage is not None:
            print(
                f"Ignored {ignoredDueToCoverageCount} genomes due to low " f"coverage.",
                file=sys.stderr,
            )

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Describe a genome (or genomes).",
    )

    parser.add_argument(
        "--genome",
        metavar="file.fasta",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="The FASTA file containing the genome(s) to examine.",
    )

    parser.add_argument(
        "--feature",
        action="append",
        metavar="FEATURE",
        help="The feature to describe (e.g., nsp2). May be repeated.",
    )

    parser.add_argument(
        "--outDir",
        metavar="DIR",
        help=(
            "The directory to write alignments and sequences to. If not "
            "specified, standard output is used."
        ),
    )

    parser.add_argument(
        "--checkVariant",
        action="append",
        help=(
            f"Check whether the genome matches the changes in a known "
            f"variant. The checked variant(s) must either be found in the "
            f'known variants (currently {", ".join(sorted(VARIANTS))}) '
            f"or else be given in a JSON file using --variantFile. "
            f"In case of conflict, names in any given --variantFile have "
            f"precedence over the predefined names. May be repeated."
        ),
    )

    parser.add_argument(
        "--variantFile",
        metavar="VARIANT-FILE.json",
        help=(
            "A JSON file of variant information. See gb2seq/variants.py "
            "for the required format."
        ),
    )

    parser.add_argument(
        "--printNtSequence",
        "--printNTSequence",
        action="store_true",
        help="Print the nucleotide sequence.",
    )

    parser.add_argument(
        "--printAaSequence",
        "--printAASequence",
        action="store_true",
        help="Print the amino acid sequence.",
    )

    parser.add_argument(
        "--printNtMatch",
        "--printNTMatch",
        action="store_true",
        help="Print details of the nucleotide match with the reference.",
    )

    parser.add_argument(
        "--printAaMatch",
        "--printAAMatch",
        action="store_true",
        help="Print details of the amino acid match with the reference.",
    )

    parser.add_argument(
        "--printNtAlignment",
        "--printNTAlignment",
        action="store_true",
        help="Print the nucleotide alignment with the reference.",
    )

    parser.add_argument(
        "--printAaAlignment",
        "--printAAAlignment",
        action="store_true",
        help="Print the amino acid alignment with the reference.",
    )

    parser.add_argument(
        "--canonicalNames",
        action="store_true",
        help=(
            "Use canonical feature names for output files, as oppposed to "
            "aliases that might be given on the command line. This can be "
            "used to ensure that output files have predictable names."
        ),
    )

    parser.add_argument(
        "--noFeatures",
        action="store_true",
        help="Do not look up any features by default.",
    )

    parser.add_argument(
        "--slashReplacement",
        default="_",
        help=(
            "The character(s) used to replace slashes in filenames (if "
            "slashes are found in a sequence id or feature name)."
        ),
    )

    parser.add_argument(
        "--spaceReplacement",
        default=" ",
        help=(
            "The character(s) used to replace spaces in filenames (if "
            "spaces are found in a sequence id or feature name)."
        ),
    )

    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print a summary of what was done.",
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
        "--onError",
        choices=("ignore", "print", "raise"),
        default="print",
        help=(
            "What to do if an error occurs (e.g., due to translating or an "
            "index out of range."
        ),
    )

    addFeatureOptions(parser)
    addAlignerOption(parser)

    args = parser.parse_args()

    sys.exit(main(args))
