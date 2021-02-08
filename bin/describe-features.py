#!/usr/bin/env python

import argparse

from sars2seq.features import Features


def main(args):
    """
    Describe SARS-CoV-2 annotations.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    features = Features(args.gbFile)

    print(f'Features for {features.reference.id}:')

    for featureName in sorted(features):
        feature = features[featureName]
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
        '--maxLen', type=int, default=80,
        help=('The maximum sequence length to print. Longer sequences will '
              'be truncated.'))

    args = parser.parse_args()

    main(args)
