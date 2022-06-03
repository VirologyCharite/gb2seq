#!/usr/bin/env python

import sys
import argparse
from collections import defaultdict

from sars2seq.features import Features
from sars2seq.sars2 import SARS_COV_2_ALIASES


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
        if name.startswith("nsp"):
            return "nsp", int(name[3:])
        elif name.startswith("ORF"):
            return "orf", int(name[3:].split()[0].rstrip("ab"))
        else:
            return name.lower(), 0

    featureNames = sorted(features, key=key)
    aka = defaultdict(set)
    for alias, name in SARS_COV_2_ALIASES.items():
        aka[name].add(alias)

    for featureName in featureNames:
        try:
            akas = " " + ", ".join(sorted(aka[featureName]))
        except KeyError:
            akas = ""

        print(f"{featureName}:{akas}")


def main(args):
    """
    Describe SARS-CoV-2 annotations.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    """
    features = Features(
        args.reference,
        sars2=args.sars2,
        addUnannotatedRegions=args.addUnannotatedRegions,
    )

    if args.names:
        printNames(features)
        return

    if args.name:
        try:
            wantedName = features.canonicalName(args.name)
        except KeyError:
            forgot = (
                " It looks like you forget to use --sars2."
                if args.name.lower() in SARS_COV_2_ALIASES
                else ""
            )
            print(f"Feature {args.name!r} is not known.{forgot}", file=sys.stderr)
            sys.exit(1)
    else:
        wantedName = None
        # All features are wanted, so print a little title.
        print(f"Features for {features.reference.id}:")

    for name in sorted(features):
        if not wantedName or name == wantedName:
            print(features.toString(name, args.maxLen))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Describe a SARS-CoV-2 sequence.",
    )

    parser.add_argument(
        "--reference",
        "--gbFile",
        metavar="file.gb",
        default=Features.REF_GB,
        help="The GenBank file to examine.",
    )

    parser.add_argument(
        "--name",
        metavar="NAME",
        help=(
            "The feature to print information for (all features are "
            "printed if not specified)."
        ),
    )

    parser.add_argument(
        "--maxLen",
        type=int,
        default=80,
        help=(
            "The maximum sequence length to print. Longer sequences will "
            "be truncated."
        ),
    )

    parser.add_argument(
        "--names",
        action="store_true",
        help="Only print feature names and aliases (if any).",
    )

    parser.add_argument(
        "--sars2",
        action="store_true",
        help="The sequence is from SARS-CoV-2.",
    )

    parser.add_argument(
        "--addUnannotatedRegions",
        action="store_true",
        help=(
            "Add unannotated regions (i.e., genome regions that have "
            'no features). These will be named "unannotated region N" '
            "where N is a natural number."
        ),
    )

    args = parser.parse_args()

    main(args)
