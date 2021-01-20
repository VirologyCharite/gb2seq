#!/usr/bin/env python

import sys
import os
import argparse
from math import log10
from os.path import exists, join
from contextlib import contextmanager

from dark.fasta import FastaReads
from dark.aa import compareAaReads, matchToString as aaMatchToString
from dark.dna import compareDNAReads, matchToString as dnaMatchToString
from dark.reads import Read, Reads

from sars2seq.features import Features
from sars2seq.genome import SARS2Genome
from sars2seq.translate import TranslationError
from sars2seq.variants import VARIANTS


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
        prefix = read.id.split()[0].replace('/', '_')
        filename = join(args.outDir, f"{prefix}{suffix}")
        with open(filename, 'w') as fp:
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
            prefix = read.id.split()[0].replace('/', '_')
            filename = join(
                args.outDir,
                f"{prefix}-{feature}{('-nt' if nt else '-aa')}{suffix}")
            return open(filename, 'w')
        else:
            return sys.stdout

    fps = {}

    try:
        if args.printNtMatch:
            fps['nt-match'] = fp('-match.txt', True)
        if args.printAaMatch:
            fps['aa-match'] = fp('-match.txt', False)

        if args.printNtSequence:
            fps['nt-sequence'] = fp('-sequence.fasta', True)
        if args.printAaSequence:
            fps['aa-sequence'] = fp('-sequence.fasta', False)

        if args.printNtAlignment:
            fps['nt-align'] = fp('-align.fasta', True)
        if args.printAaAlignment:
            fps['aa-align'] = fp('-align.fasta', False)

        yield fps

    finally:
        if args.outDir:
            for fp in fps.values():
                fp.close()


def printDiffs(read1, read2, nt, referenceOffset, fp, indent=''):
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
    what = 'nt' if nt else 'aa'
    header = '%sDifferences: site, %s1, %s2, ref nt %s' % (
        indent, what, what, 'site' if nt else 'codon start')

    for site, (a, b) in enumerate(zip(read1.sequence, read2.sequence)):
        if a != b:
            if not headerPrinted:
                print(header, file=fp)
                headerPrinted = True
            print('%s  %*d %s %s %5d' % (
                indent, width, site + 1, a, b,
                referenceOffset + (multiplier * site) + 1), file=fp)


def printVariantSummary(genome, fp, args):
    """
    Print a summary of whether the genome fulfils the various
    variant properties.

    @param genome: A C{SARS2Genome} instance.
    @param fp: An open file pointer to write to.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    print('Variant summary:', file=fp)
    for variant in args.checkVariant:
        testCount, errorCount, tests = genome.checkVariant(variant,
                                                           args.window)
        successCount = testCount - errorCount
        print(f'  {VARIANTS[variant]["description"]}:', file=fp)
        print(f'  {testCount} checks, {successCount} passed.',
              file=fp)
        if successCount == 0:
            continue
        for feature in tests:
            for type_ in tests[feature]:
                found = set()
                for change, (_, _, genOK, _) in tests[feature][type_].items():
                    if genOK:
                        found.add(change)
                if found:
                    print(f'    {feature} {type_}: ',
                          ', '.join(sorted(found)), file=fp)


def processFeature(featureName, features, genome, fps, featureNumber, args):
    """
    Process a feature from a genome.

    @param featureName: A C{str} feature name.
    @param features: A C{Features} instance.
    @param genome: A C{SARS2Genome} instance.
    @param fps: A C{dict} of file pointers for the various output streams.
    @param featureNumber: The C{int} 0-based count of the features requested.
        This will be zero for the first feature, 1 for the second, etc.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    result = genome.feature(featureName, args.window)
    feature = features.getFeature(featureName)
    genomeNt, referenceNt = result.ntSequences()

    if args.printAaMatch or args.printAaSequence or args.printAaAlignment:
        try:
            genomeAa, referenceAa = result.aaSequences()
        except TranslationError as e:
            print(f'Could not translate feature {featureName} in genome '
                  f'{genome.genome.id}: {e}', file=sys.stderr)
            genomeAa = referenceAa = None

    newlineNeeded = False

    if args.printNtMatch:
        fp = fps['nt-match']
        if featureNumber:
            print(file=fp)
        print(f'Feature: {featureName} nucleotide match', file=fp)
        print(f'  Reference nt location {feature["start"] + 1}, genome nt '
              f'location {result.genomeOffset + 1}', file=fp)
        match = compareDNAReads(referenceNt, genomeNt)
        print(dnaMatchToString(match, referenceNt, genomeNt,
                               matchAmbiguous=False, indent='  '), file=fp)
        printDiffs(referenceNt, genomeNt, True, feature['start'], fp,
                   indent='    ')
        newlineNeeded = True

    if args.printAaMatch and genomeAa:
        fp = fps['aa-match']
        if newlineNeeded or featureNumber:
            print(file=fp)
        print(f'Feature: {featureName} amino acid match', file=fp)
        match = compareAaReads(referenceAa, genomeAa)
        print(aaMatchToString(match, referenceAa, genomeAa, indent='  '),
              file=fp)
        printDiffs(referenceAa, genomeAa, False, feature['start'], fp,
                   indent='    ')

    if args.printNtSequence:
        noGaps = Read(genomeNt.id, genomeNt.sequence.replace('-', ''))
        Reads([noGaps]).save(fps['nt-sequence'])

    if args.printAaSequence and genomeAa:
        noGaps = Read(genomeAa.id, genomeAa.sequence.replace('-', ''))
        Reads([noGaps]).save(fps['aa-sequence'])

    if args.printNtAlignment:
        Reads([genomeNt, referenceNt]).save(fps['nt-align'])

    if args.printAaAlignment and genomeAa:
        Reads([genomeAa, referenceAa]).save(fps['aa-align'])


def main(args):
    """
    Describe a SARS-CoV-2 genome.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    outDir = args.outDir
    if outDir:
        if not exists(outDir):
            os.makedirs(outDir)

    features = Features(args.gbFile)

    if args.feature:
        if args.canonicalNames:
            wantedFeatures = map(features.canonicalName, args.feature)
        else:
            wantedFeatures = args.feature
    else:
        if args.noFeatures:
            wantedFeatures = []
        else:
            wantedFeatures = sorted(features.featuresDict())

    for read in FastaReads(args.genome):
        genome = SARS2Genome(read, features)

        if args.checkVariant:
            with genomeFilePointer(read, args, '-variant-summary.txt') as fp:
                print(read.id, file=fp)
                printVariantSummary(genome, fp, args)

        for i, featureName in enumerate(wantedFeatures):
            with featureFilePointers(read, featureName, args) as fps:
                processFeature(featureName, features, genome, fps, i, args)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Describe a SARS-CoV-2 genome (or genomes).')

    parser.add_argument(
        '--genome', metavar='file.fasta',
        help='The FASTA file containing the SARS-CoV-2 genome(s) to examine.')

    parser.add_argument(
        '--feature', action='append', metavar='FEATURE',
        help='The feature to describe (e.g., nsp2). May be repeated.')

    parser.add_argument(
        '--outDir', metavar='DIR',
        help=('The directory to write alignments and sequences to. If not '
              'specified, standard output is used.'))

    parser.add_argument(
        '--checkVariant', action='append', choices=sorted(VARIANTS),
        help='Check whether the genome fulfils a known variant.')

    parser.add_argument(
        '--printNtSequence', default=False, action='store_true',
        help='Print the nucleotide sequence.')

    parser.add_argument(
        '--printAaSequence', default=False, action='store_true',
        help='Print the amino acid sequence.')

    parser.add_argument(
        '--printNtMatch', default=False, action='store_true',
        help='Print details of the nucleotide match with the reference.')

    parser.add_argument(
        '--printAaMatch', default=False, action='store_true',
        help='Print details of the amino acid match with the reference.')

    parser.add_argument(
        '--printNtAlignment', '--printNTAlignment', default=False,
        action='store_true',
        help='Print the nucleotide alignment with the reference.')

    parser.add_argument(
        '--printAaAlignment', '--printAAAlignment', default=False,
        action='store_true',
        help='Print the amino acid alignment with the reference.')

    parser.add_argument(
        '--canonicalNames', default=False, action='store_true',
        help=('Use canonical feature names for output files, as oppposed to '
              'aliases that might be given on the command line. This can be '
              'used to ensure that output files have predictable names.'))

    parser.add_argument(
        '--noFeatures', default=False, action='store_true',
        help='Do not look up any features by default.')

    parser.add_argument(
        '--gbFile', metavar='file.gb', default=Features.REF_GB,
        help='The Genbank file to read for SARS-CoV-2 features.')

    parser.add_argument(
        '--window', type=int, default=500,
        help=('The size of the window (of nucleotides) surrounding the '
              'feature (in the reference) to examine in the genome.'))

    args = parser.parse_args()

    main(args)
