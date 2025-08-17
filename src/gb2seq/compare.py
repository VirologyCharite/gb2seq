from collections import defaultdict
import rich

from gb2seq.alignment import Gb2Alignment
from gb2seq.features import Features

from dark.aaVars import REVERSE_CODONS
from dark.reads import Read
from dark.utils import pct


def _get_codon(
    sequence: str,
    offset: int,
    feature_offset: int,
    which_sequence: str,
    verbose: int,
    summary: list[str],
):
    assert offset >= feature_offset

    frame = (offset - feature_offset) % 3
    codon_offset = offset - frame
    codon = sequence[codon_offset : codon_offset + 3]

    if verbose > 1:
        summary.append(
            f"    Got codon {codon!r} corresponding to nt {sequence[offset]!r} in "
            f"{which_sequence} at offset {offset} in {which_sequence}. Frame: {frame}, "
            f"Codon offset: {codon_offset}. Feature offset: {feature_offset}."
        )

    return codon


def _analyze_mutation(
    a_offset: int,
    b_offset: int,
    a_read: Read,
    b_read: Read,
    a_feature_offset: int,
    b_feature_offset: int,
    verbose: int,
    summary: list[str],
) -> str:
    """
    What can we say about this from_base to to_base mutation at this site?
    """
    from_base: str = a_read.sequence[a_offset]
    to_base: str = b_read.sequence[b_offset]
    assert from_base != to_base

    original_codon = _get_codon(
        a_read.sequence, a_offset, a_feature_offset, "genome A", verbose, summary
    )
    mutated_codon = _get_codon(
        b_read.sequence, b_offset, b_feature_offset, "genome B", verbose, summary
    )

    if mutated_codon == original_codon:
        raise ValueError(
            f"Codon {mutated_codon!r} in {a_read.id!r} at offset {a_offset} is "
            f"unexpectedly identical to the codon at offset {b_offset} in "
            f"{b_read.id!r}."
        )

    if (set(original_codon) | set(mutated_codon)) - set("ACGT"):
        return "ambiguous"

    original_aa = REVERSE_CODONS.get(original_codon, "*")
    mutated_aa = REVERSE_CODONS.get(mutated_codon, "*")

    if original_aa == mutated_aa:
        if verbose:
            summary.append(
                f"    The {from_base} -> {to_base} mutation is synonymous "
                f"({original_aa})."
            )
        result = "synonymous"
    else:
        if verbose:
            summary.append(
                f"    The {from_base} -> {to_base} mutation is non-synonymous "
                f"({original_aa} -> {mutated_aa})."
            )
        result = "non-synonymous"

    summary.append(
        f"    Original codon: {original_codon} ({original_aa}) -> "
        f"new codon: {mutated_codon} ({mutated_aa})"
    )

    return result


def genome_offset_and_nt(
    reference_offset: int, alignment: Gb2Alignment, verbose: int
) -> tuple[int, str]:
    """
    Get the genome offset (if any) and nucleotide for a given reference offset.

    If the offset doesn't exist in the genome, return an offset of -1.
    """
    offset_info = alignment.offsetInfo(reference_offset)
    assert offset_info is not None
    genome_offset = offset_info["genome"]["ntOffset"]

    if genome_offset is None:
        alignment_offset = offset_info["alignmentOffset"]
        assert alignment.gappedOffsets[alignment_offset] >= len(alignment.genome)
        if verbose > 2:
            print(
                f"Reference offset {reference_offset} has alignment offset "
                f"{alignment_offset}, but that does not exist in length "
                f"{len(alignment.genome)} sequence {alignment.genome.id!r}. Skipping."
            )
        return -1, ""

    return genome_offset, alignment.genome.sequence[genome_offset]


def compare(
    a_read: Read,
    b_read: Read,
    a_alignment: Gb2Alignment,
    b_alignment: Gb2Alignment,
    features: Features,
    verbose: int,
    sars2: bool = False,
    sars2_simplify_features: bool = False,
):
    changes = defaultdict(list)

    assert features.reference

    for offset in range(len(features.reference)):
        a_offset, a_nt = genome_offset_and_nt(offset, a_alignment, verbose)
        if a_offset == -1:
            continue

        b_offset, b_nt = genome_offset_and_nt(offset, b_alignment, verbose)
        if b_offset == -1:
            continue

        if verbose > 2:
            print(
                f"Reference offset {offset} has gapped offset "
                f"{a_alignment.gappedOffsets[offset]} in A alignment, and "
                f"{b_alignment.gappedOffsets[offset]} in B alignment. "
                f"A offset {a_offset} has {a_nt!r}, B offset {b_offset} has {b_nt!r}."
            )

        if a_nt == b_nt:
            continue

        summaries = []
        feature_info = []

        changes[a_nt + b_nt].append((offset, feature_info, summaries))

        features_list = sorted(
            features.getFeatureNames(offset, includeUntranslated=True)
        )

        if features_list:
            for feature_name in features_list:
                if (
                    len(features_list) == 1
                    or not sars2
                    or not sars2_simplify_features
                    or (
                        # We are processing SARS-2 genomes and there are multiple
                        # features. If we've been asked to skip the ones that contain
                        # sub-features, process the feature unless it is one we want to
                        # skip.
                        sars2_simplify_features
                        and not (
                            feature_name.endswith(" polyprotein")
                            or feature_name.endswith(" glycoprotein")
                            or feature_name == "3'UTR"
                        )
                    )
                ):
                    summary = []
                    summaries.append(summary)
                    summary.append(
                        f"  Reference offset {offset}: {a_nt}->{b_nt} change:"
                    )

                    feature = features[feature_name]

                    # summary.append(f"  Change in {feature_name}:")
                    feature_start = feature["start"]

                    if features.translated(feature_name):
                        a_feature_info = a_alignment.offsetInfo(feature_start)
                        assert a_feature_info is not None
                        a_feature_offset = a_feature_info["genome"]["ntOffset"]

                        b_feature_info = b_alignment.offsetInfo(feature_start)
                        assert b_feature_info is not None
                        b_feature_offset = b_feature_info["genome"]["ntOffset"]

                        change_type = _analyze_mutation(
                            a_offset,
                            b_offset,
                            a_read,
                            b_read,
                            a_feature_offset,
                            b_feature_offset,
                            verbose,
                            summary,
                        )
                    else:
                        summary.append(f"  {feature_name} is not translated.")
                        change_type = "untranslated"

                    feature_info.append((feature_name, change_type))
        else:
            summaries.append(
                [
                    f"Reference offset {offset}: {a_nt}->{b_nt} but no features "
                    "found here.",
                ]
            )

    return changes


def print_results(
    changes: dict[str, list[tuple[int, list[tuple[str, str]], list[str]]]],
    reference_length: int,
    verbose: int,
    use_rich: bool,
):
    for (a_nt, b_nt), details in sorted(changes.items()):
        change_type_count = defaultdict(int)
        for _, feature_info, _ in details:
            if feature_info:
                for _, change_type in feature_info:
                    change_type_count[change_type] += 1
            else:
                change_type_count["no feature"] += 1

        change_type_summary = ", ".join(
            f"{change_type}: {count}"
            for change_type, count in sorted(change_type_count.items())
        )

        if use_rich:
            rich.print(
                f"[bold yellow]{a_nt}->{b_nt}[/] substitutions in "
                f"{pct(len(details), reference_length)} sites ({change_type_summary})"
            )
        else:
            print(
                f"{a_nt}->{b_nt} substitutions in "
                f"{pct(len(details), reference_length)} sites ({change_type_summary})"
            )

        if verbose:
            for change_count, (offset, feature_info, summaries) in enumerate(
                details, start=1
            ):
                print(f"  {change_count:2d}: site {offset + 1:,}")
                if feature_info:
                    for (feature_name, change_type), summary in zip(
                        feature_info, summaries
                    ):
                        if use_rich:
                            rich.print(
                                f"      [bold green]{feature_name} ({change_type})"
                            )
                        else:
                            print(f"      {feature_name} ({change_type})")
                        if verbose > 1 and summary:
                            print("      " + "\n      ".join(summary))
                else:
                    print("      No features.")
