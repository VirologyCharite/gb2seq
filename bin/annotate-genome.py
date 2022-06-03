#!/usr/bin/env python

import sys
import argparse
from json import dumps
from operator import itemgetter

from dark.fasta import FastaReads

from sars2seq.alignment import SARS2Alignment, addAlignerOption, alignmentEnd
from sars2seq.features import Features, addFeatureOptions


def main(args):
    """
    Add features to an unannotated genome and print the result as JSON.

    @param args: A C{Namespace} instance as returned by argparse with
        values for command-line options.
    @return: An C{int} exit status.
    """
    features = Features(args.reference, sars2=args.sars2)
    genome = list(FastaReads(args.genome))[0]
    alignment = SARS2Alignment(genome, features, aligner=args.aligner)
    key = itemgetter("start")
    result = {
        "id": genome.id,
        "sequence": genome.sequence,
        "features": {},
    }

    resultFeatures = result["features"]

    for feature in sorted(features.values(), key=key):
        name = feature["name"]
        forward = feature["forward"]
        gappedOffset = alignment.gappedOffsets[feature["start"]]

        referenceAA, genomeAA = alignment.aaSequences(name, raiseOnReferenceGaps=False)

        start = gappedOffset - alignment.genomeAligned.sequence[:gappedOffset].count(
            "-"
        )

        if forward:
            assert 0 <= start < len(genome)
            orf = genome.findORF(start, forward=True, requireStartCodon=False)
        else:
            reverseGappedOffset = alignment.gappedOffsets[feature["stop"]]
            reverseStart = (
                len(alignment.genomeAligned)
                - reverseGappedOffset
                - alignment.genomeAligned.sequence[reverseGappedOffset:].count("-")
            )
            assert 0 <= reverseStart < len(genome)
            orf = genome.findORF(reverseStart, forward=False, requireStartCodon=False)

        # Small sanity check.
        if not orf["foundStartCodon"]:
            assert genomeAA.sequence[0] != "M"

        stop = start + orf["length"] * 3

        resultFeatures[name] = feature.copy()
        resultFeatures[name].update(
            {
                "start": start,
                "stop": stop,
                "sequence": genome.sequence[start:stop],
                "translation": orf["translation"],
            }
        )

        aaDiffs = {}

        if len(orf["translation"]) != len(referenceAA):
            aaDiffs["lengths"] = {
                "genome length": len(orf["translation"]),
                "reference length": len(referenceAA),
                "difference": len(orf["translation"]) - len(referenceAA),
            }

        subs = []
        # Note the the following zip deliberately stops at the shortest AA
        # sequence.
        for site, (aa1, aa2) in enumerate(
            zip(referenceAA.sequence, orf["translation"]), start=1
        ):
            if aa1 != aa2:
                subs.append(f"{aa1}{site}{aa2}")

        if subs:
            aaDiffs["substitutions"] = subs

        if aaDiffs:
            resultFeatures[name]["aa differences"] = aaDiffs

        if not (orf["foundStartCodon"] and orf["foundStopCodon"]):
            aend = alignmentEnd(
                alignment.referenceAligned.sequence,
                gappedOffset,
                feature["stop"] - feature["start"],
            )
            resultFeatures[name]["debug"] = {
                "genome orf": orf,
                "feature alignment": {
                    "reference": alignment.referenceAligned.sequence[gappedOffset:aend],
                    "genome": alignment.genomeAligned.sequence[gappedOffset:aend],
                },
                "reference feature": feature.copy(),
            }

        if args.reportDifferences:
            if not orf["foundStartCodon"]:
                print(
                    f"No START codon found for genome feature " f"{feature['name']!r}.",
                    file=sys.stderr,
                )

            if not orf["foundStopCodon"]:
                print(
                    f"No STOP  codon found for genome feature " f"{feature['name']!r}.",
                    file=sys.stderr,
                )

    # Sanity check that the Features class can load what we just produced.
    Features(result, sars2=args.sars2)
    print(dumps(result, sort_keys=True, indent=4))

    if args.summarize:
        print(f'Summary of changes for {result["id"]!r}', file=sys.stderr)
        count = 0
        for name, info in result["features"].items():
            try:
                aaDiffs = info["aa differences"]["substitutions"]
            except KeyError:
                aaDiff = "-"
            else:
                aaDiff = ", ".join(aaDiffs)

            try:
                n = info["aa differences"]["lengths"]["difference"]
            except KeyError:
                lenDiff = "-"
            else:
                lenDiff = f"{abs(n)} aa " + ("shorter" if n < 0 else "longer")

            if aaDiff or lenDiff:
                count += 1
                print(
                    "\t".join(
                        (str(count), name, aaDiff, lenDiff, info.get("note", ""))
                    ),
                    file=sys.stderr,
                )

    return 0


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Print JSON annotations for a genome (or genomes).",
    )

    parser.add_argument(
        "--genome",
        metavar="file.fasta",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="The FASTA file containing the SARS-CoV-2 genome(s) to examine.",
    )

    parser.add_argument(
        "--summarize",
        action="store_true",
        help="Write a summary of differences to standard error.",
    )

    parser.add_argument(
        "--reportDifferences",
        action="store_true",
        help=(
            "If the genome ORF differs from the reference due to no start "
            "or stop codon, print a message.",
        ),
    )

    addAlignerOption(parser)
    addFeatureOptions(parser)

    args = parser.parse_args()

    sys.exit(main(args))
