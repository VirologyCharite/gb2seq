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


def printDiffs(read1, read2, fp, indent=''):
    """
    Print differences between sequences.

    @param read1: A C{dark.reads.Read} instance.
    @param read2: A C{dark.reads.Read} instance.
    @param indent: A C{str} prefix for each output line.
    """
    len1, len2 = len(read1), len(read2)
    width = int(log10(max(len1, len2))) + 1
    headerPrinted = False
    for site, (a, b) in enumerate(zip(read1.sequence, read2.sequence),
                                  start=1):
        if a != b:
            if not headerPrinted:
                print('%sDifferences (site, %s, %s):' % (
                    indent, read1.id, read2.id))
                headerPrinted = True
            print('%s  %*d %s %s' % (indent, width, site, a, b), file=fp)

    if not headerPrinted:
        print('No sequence differences found.', file=fp)


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
        _, errorCount, _ = genome.checkVariant(variant)
        print(f'  {VARIANTS[variant]["description"]}:',
              'Yes' if errorCount == 0 else 'No', file=fp)


def processFeature(feature, genome, fps, args):
    """
    Process a feature from a genome.

    @param feature: A C{str} feature name.
    @param genome: A C{SARS2Genome} instance.
    @param fps: A C{dict} of file pointers for the various output streams.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    result = genome.feature(feature)
    genomeNt, referenceNt = result.ntSequences()
    genomeAa, referenceAa = result.aaSequences()

    if args.printNtMatch:
        fp = fps['nt-match']
        match = compareDNAReads(referenceNt, genomeNt)
        print(dnaMatchToString(match, referenceNt, genomeNt,
                               matchAmbiguous=False), file=fp)
        printDiffs(referenceNt, genomeNt, fp, indent='  ')

    if args.printAaMatch:
        fp = fps['aa-match']
        match = compareAaReads(referenceAa, genomeAa)
        print(aaMatchToString(match, referenceAa, genomeAa), file=fp)
        printDiffs(referenceAa, genomeAa, fp, indent='  ')

    if args.printNtSequence:
        noGaps = Read(genomeNt.id, genomeNt.sequence.replace('-', ''))
        Reads([noGaps]).save(fps['nt-sequence'])

    if args.printAaSequence:
        noGaps = Read(genomeAa.id, genomeAa.sequence.replace('-', ''))
        Reads([noGaps]).save(fps['aa-sequence'])

    if args.printNtAlignment:
        Reads([genomeNt, referenceNt]).save(fps['nt-align'])

    if args.printAaAlignment:
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
        wantedFeatures = sorted(features.featuresDict())

    for read in FastaReads(args.genome):
        genome = SARS2Genome(read, features)

        if args.checkVariant:
            with genomeFilePointer(read, args, '-variant-summary.txt') as fp:
                printVariantSummary(genome, fp, args)

        for feature in wantedFeatures:
            with featureFilePointers(read, feature, args) as fps:
                processFeature(feature, genome, fps, args)


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
        '--printNtAlignment', default=False, action='store_true',
        help='Print the nucleotide alignment with the reference.')

    parser.add_argument(
        '--printAaAlignment', default=False, action='store_true',
        help='Print the amino acid alignment with the reference.')

    parser.add_argument(
        '--canonicalNames', default=False, action='store_true',
        help=('Use canonical feature names for output files, as oppposed to '
              'aliases that might be given on the command line. This can be '
              'used to ensure that output files have predictable names.'))

    parser.add_argument(
        '--gbFile', metavar='file.gb', default=Features.REF_GB,
        help='The Genbank file to read for SARS-CoV-2 features.')

    args = parser.parse_args()

    main(args)
