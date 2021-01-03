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
from sars2seq.sequence import SARS2Sequence


def save(sequence, reference, read, feature, outDir, nt):
    """
    Save the sequences and alignment.

    @param sequence: A C{dark.reads.Read} instance (aligned, so with gaps).
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

    Reads([sequence, reference]).save(filenameBase + '-align.fasta')

    sequenceNoGaps = Read(sequence.id, sequence.sequence.replace('-', ''))
    referenceNoGaps = Read(reference.id, reference.sequence.replace('-', ''))
    Reads([sequenceNoGaps, referenceNoGaps]).save(filenameBase + '.fasta')


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


def main(args):
    """
    Describe a SARS-CoV-2 sequence.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    for read in FastaReads(args.fastaFile):
        seq = SARS2Sequence(read, Features(args.gbFile))

        features = args.feature or seq.featureNames()

        for feature in features:
            print(f'Matching {feature}')
            result = seq.feature(feature)

            sequenceNt, referenceNt = result.ntSequences()
            match = compareDNAReads(referenceNt, sequenceNt)
            print('  DNA:')
            print(dnaMatchToString(match, referenceNt, sequenceNt,
                                   matchAmbiguous=False, indent='   '))
            if match['match']['nonGapMismatchCount']:
                printDiffs(referenceNt, sequenceNt, indent='     ')

            sequenceAa, referenceAa = result.aaSequences()
            match = compareAaReads(referenceAa, sequenceAa)
            print('  AA:')
            print(aaMatchToString(match, referenceAa, sequenceAa,
                                  indent='   '))
            if match['match']['nonGapMismatchCount']:
                printDiffs(referenceAa, sequenceAa, indent='     ')

            outDir = args.outDir
            if outDir:
                if not exists(outDir):
                    os.makedirs(outDir)
                save(sequenceNt, referenceNt, read, feature, outDir, True)
                save(sequenceAa, referenceAa, read, feature, outDir, False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Describe a SARS-CoV-2 sequence.')

    parser.add_argument(
        'fastaFile', metavar='file.fasta',
        help='The FASTA file to examine.')

    parser.add_argument(
        '--feature', action='append', metavar='FEATURE',
        help='The feature to describe (e.g., nsp2). May be repeated.')

    parser.add_argument(
        '--outDir', metavar='DIR',
        help='The directory to write alignments and sequences to.')

    parser.add_argument(
        '--gbFile', metavar='file.gb', default=Features.REF_GB,
        help='The Genbank file to read for features.')

    args = parser.parse_args()

    main(args)
