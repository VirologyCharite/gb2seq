#!/usr/bin/env python

import argparse
from math import log10

from dark.fasta import FastaReads
from dark.aa import compareAaReads, matchToString as aaMatchToString
from dark.dna import compareDNAReads, matchToString as dnaMatchToString

from sars2seq.features import Features
from sars2seq.sequence import SARS2Sequence


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
    features = Features(args.gbFile)

    for read in FastaReads(args.fastaFile):
        seq = SARS2Sequence(read, features)

        if args.feature:
            features = [args.feature]
        else:
            features = seq.featureNames()

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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Describe a SARS-CoV-2 sequence.')

    parser.add_argument(
        'fastaFile', metavar='file.fasta',
        help='The FASTA file to examine.')

    parser.add_argument(
        '--feature',
        help='The feature to describe (e.g., nsp2).')

    parser.add_argument(
        '--gbFile', metavar='file.gb', default=Features.REF_GB,
        help='The Genbank file to read for features.')

    args = parser.parse_args()

    main(args)
