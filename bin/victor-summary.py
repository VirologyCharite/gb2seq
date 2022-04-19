#!/usr/bin/env python

import sys
import os
import argparse
from math import log10
from os.path import exists, join
from contextlib import contextmanager
from collections import defaultdict

from dark.fasta import FastaReads
from dark.aa import compareAaReads, matchToString as aaMatchToString
from dark.dna import compareDNAReads, matchToString as dnaMatchToString
from dark.reads import Read, Reads

from sars2seq.features import Features
from sars2seq.genome import SARS2Genome, addAlignerOption
from sars2seq.variants import VARIANTS

CHANGE_SETS = {
    'UK': {'H69-', 'V70-', 'N501Y', 'D614G', 'P681H'},
    'ZA': {'K417N', 'E484K', 'N501Y', 'D614G'},
    'Japan': {'K417T', 'E484K', 'N501Y', 'D614G'},
    'Mink': {'H69-', 'V70-', 'Y453F', 'D614G'},
}


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

    print('referenceOffset', referenceOffset)
    for site, (a, b) in enumerate(zip(read1.sequence, read2.sequence)):
        if a != b:
            if not headerPrinted:
                print(header)
                headerPrinted = True
            print('%s  %*d %s %s %5d' % (
                indent, width, site + 1, a, b,
                referenceOffset + (multiplier * site) + 1), file=fp)


def key(s):
    return int(s.split('(')[0][1:-1])


def printVariantSummary(genome, fp, args):
    """
    Print a summary of whether the genome fulfils the various
    variant properties.

    @param genome: A C{SARS2Genome} instance.
    @param fp: An open file pointer to write to.
    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    namedMatches = defaultdict(list)
    foundSets = defaultdict(list)
    shortId = genome.genome.id.split()[0]

    for variant in args.checkVariant:
        testCount, errorCount, tests = genome.checkVariant(variant)
        successCount = testCount - errorCount

        for feature in tests:
            for type_ in tests[feature]:
                found = set()
                nonRef = set()
                notCovered = set()
                for change, (refOK, refBase, genOK, genBase) in tests[
                        feature][type_].items():
                    if not refOK:
                        print(f'  Ref mismatch for {change}: {refBase}',
                              file=fp)
                    if genOK is True:
                        found.add(change)
                    else:
                        if genBase != refBase:
                            if (type_, genBase) in (('aa', 'X'), ('nt', 'N')):
                                notCovered.add(change)
                            else:
                                # Note that the '(' in the string made here
                                # will be used for splitting in the sort
                                # 'key' function above.
                                nonRef.add(f'{change}({genBase})')

                if found:
                    assert successCount == len(found)
                    print(f'  {successCount}/{testCount} changes found:',
                          ', '.join(sorted(found, key=key)), file=fp)
                    foundSets[tuple(sorted(found))].append(shortId)
                else:
                    print(f'  0/{testCount} changes found.', file=fp)

                if nonRef:
                    print('  Unexpected changes:', ', '.join(
                        sorted(nonRef, key=key)), file=fp)

                if notCovered:
                    print('  No coverage:',
                          ', '.join(sorted(notCovered, key=key)), file=fp)

                if type_ == 'aa':
                    matched = set()
                    for changeSet, changes in CHANGE_SETS.items():
                        if not (changes - found):
                            extras = found - changes
                            if extras:
                                matched.add(changeSet + ' + ' +
                                            ', '.join(sorted(extras, key=key)))
                            else:
                                matched.add(changeSet)

                    if matched:
                        for match in matched:
                            namedMatches[match].append(shortId)
                        print('  Matched:',
                              ', '.join(sorted(matched)), file=fp)
                    else:
                        if len(found - {'D614G'}):
                            print('  Unnamed combination of changes.', file=fp)

    return namedMatches, foundSets


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
    result = genome.feature(featureName)
    feature = features.getFeature(featureName)
    referenceNt, genomeNt = result.ntSequences()
    referenceAa, genomeAa = result.aaSequences()

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

    if args.printAaMatch:
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
        if args.noFeatures:
            wantedFeatures = []
        else:
            wantedFeatures = sorted(features)

    namedMatches = defaultdict(list)
    foundSets = defaultdict(list)

    reads = list(FastaReads(args.genome))

    print('SEQUENCE SHORT NAMES\n')
    maxLen = 0
    nameSummary = []
    for read in reads:
        shortId = read.id.split()[0]
        if len(shortId) > maxLen:
            maxLen = len(shortId)
        nameSummary.append((shortId, read.id))
        read.id = shortId

    for shortId, longId in nameSummary:
        print(f'{shortId:{maxLen}s} = {longId}')

    print('\nPER-SEQUENCE RESULTS\n')

    for read in reads:
        genome = SARS2Genome(read, features, aligner=args.aligner)

        if args.checkVariant:
            with genomeFilePointer(read, args, '-variant-summary.txt') as fp:
                nCount = genome.genome.sequence.count('N')
                genomeLen = len(genome.genome)
                nonNCount = genomeLen - nCount
                coverage = nonNCount / genomeLen
                print(f'{read.id} (coverage {nonNCount}/{genomeLen} = '
                      f'{coverage * 100.0:.2f} %)', file=fp)

                theseNamedMatches, theseFoundSets = printVariantSummary(
                    genome, fp, args)

                for match, ids in theseNamedMatches.items():
                    namedMatches[match].extend(ids)

                for match, ids in theseFoundSets.items():
                    foundSets[match].extend(ids)

                print(file=fp)

        for i, featureName in enumerate(wantedFeatures):
            with featureFilePointers(read, featureName, args) as fps:
                processFeature(featureName, features, genome, fps, i, args)

    print('\nSUMMARY\n')

    if namedMatches:
        print('Named change sets:')
        for changeSet in sorted(CHANGE_SETS):
            desc = ', '.join(sorted(CHANGE_SETS[changeSet], key=key))
            print(f'  {changeSet}: {desc}')
        print()

        print('Known variant combinations matched (count):')
        for match in sorted(namedMatches):
            print(f'  {match} ({len(namedMatches[match])}):')
            for name in sorted(namedMatches[match]):
                print(f'    {name}')
        if foundSets:
            print()

    if foundSets:
        print('Sets of changes found (count):')
        for match in sorted(foundSets):
            desc = ', '.join(sorted(match, key=key))
            print(f'  {desc} ({len(foundSets[match])}):')
            for name in sorted(foundSets[match]):
                print(f'    {name}')


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
        '--noFeatures', default=False, action='store_true',
        help='Do not look up any features by default.')

    parser.add_argument(
        '--gbFile', metavar='file.gb', default=Features.REF_GB,
        help='The GenBank file to read for SARS-CoV-2 features.')

    addAlignerOption(parser)

    args = parser.parse_args()

    main(args)
