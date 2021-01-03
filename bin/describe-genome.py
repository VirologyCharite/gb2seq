#!/usr/bin/env python

import os
import argparse
from math import log10
from os.path import exists, join

from dark.fasta import FastaReads
from dark.aa import compareAaReads, matchToString as aaMatchToString
from dark.dna import compareDNAReads, matchToString as dnaMatchToString
from dark.reads import Read, Reads

from sars2seq.features import Features
from sars2seq.genome import SARS2Genome
from sars2seq.variants import spikeDeletion, VOC_20201201_UK, N501Y


def save(genome, reference, read, feature, outDir, nt):
    """
    Save the sequences and alignment.

    @param genome: A C{dark.reads.Read} instance (aligned, so with gaps).
    @param reference: A C{dark.reads.Read} instance (aligned, so with gaps).
    @param read: The C{dark.reads.Read} that was read from the input FASTA file
        (this is the overall genome from which the feature was obtained).
    @param feature: The C{str} name of the feature.
    @param outDir: The C{str} output directory (must already exist).
    @param nt: If C{True} the sequences are nucelotide, else protein.
    """
    filenameBase = join(
        outDir,
        read.id.split()[0] + '-' + feature + ('-nt' if nt else '-aa'))

    Reads([genome, reference]).save(filenameBase + '-align.fasta')

    genomeNoGaps = Read(genome.id, genome.sequence.replace('-', ''))
    referenceNoGaps = Read(reference.id, reference.sequence.replace('-', ''))
    Reads([genomeNoGaps, referenceNoGaps]).save(filenameBase + '.fasta')


def printDiffs(read1, read2, indent=''):
    """
    Print differences between sequences.

    @param read1: A C{dark.reads.Read} instance.
    @param read2: A C{dark.reads.Read} instance.
    @param indent: A C{str} prefix for each output line.
    """
    # Print all sites where the sequences differ.
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
            print('%s  %*d %s %s' % (indent, width, site, a, b))

    if not headerPrinted:
        print('No sequence differences found.')


def printVariantSummary(genome):
    """
    Print a summary of whether the genome fulfils the various
    variant properties.

    @param genome: A C{SARS2Genome} instance.
    """
    print('Variant summary:')
    for variant, desc in ((spikeDeletion, 'Spike deletion'),
                          (VOC_20201201_UK, 'UK VOC202012/01'),
                          (N501Y, 'N501K change')):
        _, errorCount, _ = genome.checkVariant(variant)
        print(f'  {desc}:', 'Yes' if errorCount == 0 else 'No')


def main(args):
    """
    Describe a SARS-CoV-2 genome.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    features = Features(args.gbFile)
    for read in FastaReads(args.genome):
        genome = SARS2Genome(read, features)

        if args.summarizeVariants:
            printVariantSummary(genome)

        wantedFeatures = args.feature or list(features.featuresDict())

        for feature in wantedFeatures:
            print(f'Summary for {feature!r}')
            result = genome.feature(feature)

            genomeNt, referenceNt = result.ntSequences()
            match = compareDNAReads(referenceNt, genomeNt)
            print('  DNA:')
            print(dnaMatchToString(match, referenceNt, genomeNt,
                                   matchAmbiguous=False, indent='   '))
            if match['match']['nonGapMismatchCount']:
                printDiffs(referenceNt, genomeNt, indent='     ')

            genomeAa, referenceAa = result.aaSequences()
            match = compareAaReads(referenceAa, genomeAa)
            print('  AA:')
            print(aaMatchToString(match, referenceAa, genomeAa,
                                  indent='   '))
            if (match['match']['nonGapMismatchCount'] or
                    match['match']['gapMismatchCount']):
                printDiffs(referenceAa, genomeAa, indent='     ')

            outDir = args.outDir
            if outDir:
                if not exists(outDir):
                    os.makedirs(outDir)
                save(genomeNt, referenceNt, read, feature, outDir, True)
                save(genomeAa, referenceAa, read, feature, outDir, False)


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
        help='The directory to write alignments and sequences to.')

    parser.add_argument(
        '--summarizeVariants', default=False, action='store_true',
        help=('Summarize whether the genome fulfils any of the known variant '
              'profiles.'))

    parser.add_argument(
        '--gbFile', metavar='file.gb', default=Features.REF_GB,
        help='The Genbank file to read for SARS-CoV-2 features.')

    args = parser.parse_args()

    main(args)
