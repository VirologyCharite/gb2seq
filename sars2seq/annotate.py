from operator import itemgetter
from typing import Dict, List, Union

from dark.reads import DNARead

from sars2seq.alignment import SARS2Alignment, alignmentEnd
from sars2seq.features import Features


def annotateGenome(features: Features, genome: DNARead, aligner: str) -> dict:
    """
    Find features in an unannotated genome.

    @param features: A C{Features} instance, with the features from
        a sequence that is annotated.
    @param genome: The C{DNARead} of the genome to annotate.
    @aligner: A C{str} aligner name, such as 'edlib' or 'mafft'.
    @raise ValueError: If the result we plan to return cannot be used to
        initialize a new C{Features} instance.
    @return: A C{dict} with C{str} feature names as keys. See the 'result'
        variable below.
    """
    alignment = SARS2Alignment(genome, features, aligner=aligner)

    result = {
        "id": genome.id,
        "sequence": genome.sequence,
        "features": {},
    }

    resultFeatures = result["features"]

    for feature in sorted(features.values(), key=itemgetter("start")):
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

        aaDiffs: Dict[str, Union[dict, List[str]]] = {}

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

    try:
        # Sanity check that the Features class can load what we just produced.
        Features(result, sars2=features.sars2)
    except Exception as e:
        raise ValueError(
            f"Could not use the result to initialize a " f"fresh Features instance: {e}"
        )

    return result


def summarizeDifferences(annotations: dict) -> None:
    """
    Summarize new annotations.

    @param annotations: The C{dict} result of calling C{annotateGenome}.
    @return: A C{str} summary of differences.
    """
    result = [f'Summary of changes for {annotations["id"]!r}']

    count = 0
    for name, info in annotations["features"].items():
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

        if not (aaDiff == lenDiff == "-"):
            count += 1
            result.append(
                "\t".join((str(count), name, aaDiff, lenDiff, info.get("note", "")))
            )

    return "\n".join(result)
