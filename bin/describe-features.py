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
    featuresDict = features.featuresDict()

    print(f'Features for {features.id}:')

    for name in sorted(featuresDict):
        feature = featuresDict[name]
        print(f'{name}:')
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
              (sequence[:70] + '...') if len(sequence) > 70 else sequence)

        try:
            translation = feature['translation']
        except KeyError:
            # Some features (e.g., UTR, stem loops) do not have a translation.
            pass
        else:
            print(f'  translation (len {len(translation):5d} aa):',
                  (translation[:70] + '...') if len(translation) > 70
                  else translation)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Describe a SARS-CoV-2 sequence.')

    parser.add_argument(
        '--gbFile', metavar='file.gb', default=Features.REF_GB,
        help='The Genbank file to examine.')

    args = parser.parse_args()

    main(args)
