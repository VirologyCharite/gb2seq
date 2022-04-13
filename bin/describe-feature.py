#!/usr/bin/env python

import argparse
from collections import defaultdict

from sars2seq.features import Features, ALIASES


def printNames(features):
    """
    Print feature names, each with all known aliases (if any).

    @param features: A C{Features} instance.
    """

    def key(name):
        """
        Make a sort key for feature names.

        @param name: A C{str} feature name.
        @return: A C{str} C{int} 2-tuple for sorting feature names.
        """
        if name.startswith('nsp'):
            return 'nsp', int(name[3:])
        elif name.startswith('ORF'):
            return 'orf', int(name[3:].split()[0].rstrip('ab'))
        else:
            return name.lower(), 0

    featureNames = sorted(features, key=key)
    aka = defaultdict(set)
    for alias, name in ALIASES.items():
        aka[name].add(alias)

    for featureName in featureNames:
        try:
            akas = ' ' + ', '.join(sorted(aka[featureName]))
        except KeyError:
            akas = ''

        print(f'{featureName}:{akas}')


def main(args):
    """
    Describe SARS-CoV-2 annotations.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    features = Features(args.gbFile)

    if args.names:
        printNames(features)
        return

    if args.name:
        wantedName = features.canonicalName(args.name)
    else:
        wantedName = None
        print(f'Features for {features.reference.id}:')

    for featureName, feature in sorted(features.items()):
        if wantedName and featureName != wantedName:
            continue
        print(f'{featureName}:')
        print('  start:', feature['start'])
        print('  stop:', feature['stop'])
        print('  length:', feature['stop'] - feature['start'])
        try:
            print('  product:', feature['product'])
        except KeyError:
            pass
        try:
            print('  function:', feature['function'])
        except KeyError:
            pass

        sequence = feature['sequence']
        print(f'  sequence    (len {len(sequence):5d} nt):',
              (sequence[:args.maxLen] + '...') if len(sequence) > args.maxLen
              else sequence)

        try:
            translation = feature['translation']
        except KeyError:
            # Some features (e.g., UTR, stem loops) do not have a translation.
            pass
        else:
            print(f'  translation (len {len(translation):5d} aa):',
                  (translation[:args.maxLen] + '...')
                  if len(translation) > args.maxLen else translation)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Describe a SARS-CoV-2 sequence.')

    parser.add_argument(
        '--gbFile', metavar='file.gb', default=Features.REF_GB,
        help='The Genbank file to examine.')

    parser.add_argument(
        '--name', metavar='NAME',
        help=('The feature to print information for (all features are '
              'printed if not specified).'))

    parser.add_argument(
        '--maxLen', type=int, default=80,
        help=('The maximum sequence length to print. Longer sequences will '
              'be truncated.'))

    parser.add_argument(
        '--names', action='store_true',
        help='Only print feature names and aliases (if any).')

    args = parser.parse_args()

    main(args)
