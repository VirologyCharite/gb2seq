#!/usr/bin/env python

import sys
import os
from json import dumps
import argparse

from dark.fasta import FastaReads

import sars2seq
from sars2seq.features import Features
from sars2seq.genome import SARS2Genome, addAlignerOption


def report(genome, args, includeGenome=True):
    """
    Report what's found at a site for a given genome.

    @param genome: A C{SARS2Genome} instance.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    @param includeGenome: If C{True}, include information about the genome
        (not just the reference).
    """
    try:
        offsetInfo = genome.offsetInfo(
            args.site - 1, relativeToFeature=args.relativeToFeature,
            aa=args.aa, featureName=args.feature,
            includeUntranslated=args.includeUntranslated,
            minReferenceCoverage=args.minReferenceCoverage)
    except sars2seq.Sars2SeqError as e:
        print(e, file=sys.stderr)
        sys.exit(1)

    if args.genomeAaOnly:
        print(offsetInfo['genome']['aa'])
    else:
        if not includeGenome:
            del offsetInfo['genome']

        if args.includeFeature:
            featureName = offsetInfo['featureName']
            if featureName:
                assert 'feature' not in offsetInfo
                offsetInfo['feature'] = genome.features[featureName]

        # TODO: what should we print if the user doesn't want JSON? Some kind
        # of textual summary, I guess. When that's implemented, remove the
        # "or True" below.
        if args.json or True:
            # Make the featureNames into a sorted list (it is by default a
            # set), so it can be printed as JSON.
            offsetInfo['featureNames'] = sorted(offsetInfo['featureNames'])
            print(dumps(offsetInfo, indent=4, sort_keys=True))
        else:
            print(offsetInfo)


def main(args):
    """
    Describe a site in a SARS-CoV-2 genome or genomes.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    @return: An C{int} exit status.
    """
    features = Features(args.gbFile)

    count = ignoredDueToCoverageCount = 0

    if args.genome is None and os.isatty(0):
        genome = SARS2Genome(features.reference, features,
                             aligner=args.aligner)
        report(genome, args, False)
    else:
        fp = open(args.genome) if args.genome else sys.stdin
        for count, read in enumerate(FastaReads(fp), start=1):
            if args.minReferenceCoverage is not None:
                coverage = ((len(read) - read.sequence.upper().count('N')) /
                            len(features.reference))
                if coverage < args.minReferenceCoverage:
                    ignoredDueToCoverageCount += 1
                    if args.verbose:
                        print(f'Genome {read.id!r} ignored due to low '
                              f'({coverage * 100.0:.2f}%) coverage of the '
                              f'reference.', file=sys.stderr)
                    continue

            genome = SARS2Genome(read, features, aligner=args.aligner)
            report(genome, args)

        if args.verbose:
            print(f'Examined {count} genomes.', file=sys.stderr)

            if args.minReferenceCoverage is not None:
                print(f'Ignored {ignoredDueToCoverageCount} genomes due to '
                      f'low coverage.', file=sys.stderr)
        if args.genome:
            fp.close()

    return 0


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Describe a site of a SARS-CoV-2 genome(s).')

    parser.add_argument(
        '--genome', metavar='file.fasta',
        help='The FASTA file containing the SARS-CoV-2 genome(s) to examine.')

    parser.add_argument(
        '--site', metavar='N', type=int, required=True,
        help='The (1-based) site to find information for.')

    parser.add_argument(
        '--feature', metavar='FEATURE',
        help=('The feature to examine (e.g., nsp2). This is required if you '
              'use --aa or --relativeToFeature'))

    parser.add_argument(
        '--includeUntranslated', action='store_true',
        help=('Include untranslated features (if no feature name is given and '
              'it is necessary to identify the intended feature just based on '
              'offset).'))

    parser.add_argument(
        '--aa', action='store_true',
        help=('The given site is an amino acid count (the default is '
              'nucleotides).'))

    parser.add_argument(
        '--relativeToFeature', action='store_true',
        help='The given site is relative to the start of the feature.')

    parser.add_argument(
        '--json', action='store_true',
        help='Print the result as JSON.')

    parser.add_argument(
        '--genomeAaOnly', action='store_true',
        help='Only print the amino acid from the genome.')

    parser.add_argument(
        '--verbose', action='store_true',
        help='Print information about proceesing to standard error.')

    parser.add_argument(
        '--includeFeature', action='store_true',
        help='Also print information about the feature at the site.')

    parser.add_argument(
        '--minReferenceCoverage', metavar='coverage', type=float,
        help=('The fraction of non-N bases required in the genome(s) in order '
              'for them to be processed. Genomes with lower coverage will be '
              'ignored, with a message printed to standard error. Note that '
              'the denominator used to compute the coverage fraction is the '
              'length of the reference. I.e., coverage is computed as number '
              'of non-N bases in the genome divided by the length of the '
              'reference.'))

    parser.add_argument(
        '--gbFile', metavar='file.gb', default=Features.REF_GB,
        help='The Genbank file to read for SARS-CoV-2 features.')

    addAlignerOption(parser)

    args = parser.parse_args()

    sys.exit(main(args))
