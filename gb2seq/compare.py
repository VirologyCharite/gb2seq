from collections import defaultdict
import rich

from gb2seq.alignment import Gb2Alignment
from gb2seq.features import Features

from dark.aaVars import REVERSE_CODONS
from dark.reads import Read
from dark.utils import pct


def _get_codon(
    sequence: str,
    site: int,
    start: int,
    which_sequence: str,
    verbose: int,
    summary: list[str],
):
    offset = site - start
    frame = offset % 3

    match frame:
        case 0:
            sites = site, site + 1, site + 2
        case 1:
            sites = site - 1, site, site + 1
        case 2:
            sites = site - 2, site - 1, site
        case _:
            raise ValueError("Impossible!")

    if verbose > 1:
        summary.append(
            f"    Getting codon for {sequence[site]!r} at site {site} in "
            f"{which_sequence}. Frame: {frame}, Codon sites: {sites}. "
            f"Feature offset: {start}."
        )

    return "".join(sequence[s] for s in sites)


def _analyze_mutation(
    a_site: int,
    b_site: int,
    a_feature_offset: int,
    b_feature_offset: int,
    a_sequence: str,
    b_sequence: str,
    verbose: int,
    summary: list[str],
) -> str:
    """
    What can we say about this from_base to to_base mutation at this site?
    """
    from_base: str = a_sequence[a_site]
    to_base: str = b_sequence[b_site]
    assert from_base != to_base

    original_codon = _get_codon(
        a_sequence, a_site, a_feature_offset, "genome A", verbose, summary
    )
    mutated_codon = _get_codon(
        b_sequence, b_site, b_feature_offset, "genome B", verbose, summary
    )

    assert mutated_codon != original_codon

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
        f"{result}: {original_codon} ({original_aa}) -> {mutated_codon} ({mutated_aa})"
    )

    return result


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
        a_offset_info = a_alignment.offsetInfo(offset)
        assert a_offset_info is not None
        a_offset = a_offset_info["genome"]["ntOffset"]
        a_nt = a_read.sequence[a_offset]

        b_offset_info = b_alignment.offsetInfo(offset)
        assert b_offset_info is not None
        b_offset = b_offset_info["genome"]["ntOffset"]
        if b_offset is None:
            if verbose > 2:
                print(
                    f"Offset {offset} has alignment offset "
                    f"{b_offset_info['alignmentOffset']}, but that does "
                    f"not exist in B genome of len {len(b_read)}. Skipping."
                )
            continue
        else:
            b_nt = b_read.sequence[b_offset]

        if a_nt == b_nt:
            continue

        summaries = []
        feature_info = []

        changes[a_nt + b_nt].append((offset, feature_info, summaries))

        features_list = sorted(features.getFeatureNames(offset))

        if features_list:
            summary = []
            summaries.append(summary)
            summary.append(f"0-based offset {offset}: {a_nt}->{b_nt} change:")
            for feature_name in features_list:
                if (
                    len(features_list) == 1
                    or not sars2
                    or (
                        # We are processing SARS-2 genomes and there are multiple
                        # features. If we've been asked to skip the ones that contain
                        # sub-features, process the feature unless it is one we want to
                        # skip.
                        sars2_simplify_features
                        and not feature_name.endswith(" polyprotein")
                        and not feature_name.endswith(" glycoprotein")
                    )
                ):
                    feature = features[feature_name]

                    summary.append(f"  Change in {feature_name}:")
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
                            a_feature_offset,
                            b_feature_offset,
                            a_read.sequence,
                            b_read.sequence,
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
                    f"0-based offset {offset}: {a_nt}->{b_nt} but no features "
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
    for change in sorted(changes):
        details = changes[change]
        a_nt, b_nt = change

        change_type_count = defaultdict(int)
        for _, feature_info, _ in details:
            if feature_info:
                for _, change_type in feature_info:
                    change_type_count[change_type] += 1
            else:
                change_type_count["non-coding"] += 1

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
                for (feature_name, change_type), summary in zip(
                    feature_info, summaries
                ):
                    if use_rich:
                        rich.print(f"      [bold green]{feature_name} ({change_type})")
                    else:
                        print(f"      {feature_name} ({change_type})")
                    if verbose > 1 and summary:
                        print("      " + "\n      ".join(summary))
