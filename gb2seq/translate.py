from Bio.Seq import translate as bpTranslate
from itertools import groupby
from typing import Dict, List, Optional

from gb2seq import Gb2SeqError

from dark.aa import CODONS, STOP_CODONS
from dark.reads import AARead


class TranslationError(Gb2SeqError):
    "Error when using custom translation of sequences."


class NoSlipperySequenceError(TranslationError):
    "No slippery sequence could be found in a genome."


class NoStopCodonError(TranslationError):
    "No stop codon was found downstream from the slippery sequence."


class StopCodonTooDistantError(TranslationError):
    "The stop codon following the slippery sequence was too far away."


class TranslatedSequenceLengthError(TranslationError):
    "A sequence to be translated has an incorrect length."


class TranslatedReferenceAndGenomeLengthError(TranslationError):
    "The translated reference and genome have different lengths."


class TranslatedGapLengthError(TranslationError):
    "The sequence to translate has stretches of gaps that are not"
    "multiples of three."


codons = dict(
    [(codon, aa) for aa, cods in CODONS.items() for codon in cods]
    + [(codon, "*") for codon in STOP_CODONS]
    + [("---", "-")]
)

# The maximum difference (number of nucleotides) to allow between the
# offset of the start of the slippery sequence and the downstream stop
# codon.
_MAX_DISTANCE_TO_STOP = 20

SLIPPERY_SEQUENCE = "TTTAAAC"

_SLIPPERY_LEN = len(SLIPPERY_SEQUENCE)


def translate(
    seq: str,
    name: Optional[str] = None,
    sars2: bool = False,
    untranslatable: Optional[Dict[str, str]] = None,
) -> str:
    """
    Translate a sequence.

    @param seq: A C{str} nucelotide sequence.
    @param name: A C{str} feature name.
    @param untranslatable: A C{dict} with C{str} keys and values. If any of
        the keys appears in a codon, the corresponding value is added to the
        translation. This can be used e.g., to make occurrences of '?'
        translate into '-' or 'X'.
    @return: A translated C{str} amino acid sequence.
    """
    if name == "ORF1ab polyprotein":
        # See Fields Virology (figure 10.6a on page 421, 7th edition or
        # figure 28.7a on page 836, 6th edition) plus
        # https://www.ncbi.nlm.nih.gov/nuccore/NC_045512 for details of
        # what happens below. Note that the nucelotide sequence we are
        # passed is the one that's made from the alignment with the
        # reference ORF1ab nucleotide sequence (in sequence.py) and so is
        # just that ORF and does not include the leading ~265 nucleotides
        # of the 5' UTR. As a result, the offset used to begin the search
        # for the slippery sequence is 13000, which is chosen to be a bit
        # before 13462 - 265.  There are various occurrences of the
        # slippery sequence in the reference genome (and hence probably in
        # other CoV genomes), but only one in this region and with a stop
        # codon shortly (up to _MAX_DISTANCE_TO_STOP nucleotides) downstream.
        offset = seq.find(SLIPPERY_SEQUENCE, 13000)
        stop = seq.find("TAA", offset + _SLIPPERY_LEN)
        if offset == -1:
            raise NoSlipperySequenceError("No slippery sequence found.")
        if stop == -1:
            raise NoStopCodonError(
                f"Could not find a stop codon downstream from the start of "
                f"the slippery sequence at location {offset + 1}."
            )
        if stop - offset > _MAX_DISTANCE_TO_STOP:
            raise StopCodonTooDistantError(
                f"The stop codon was too far ({stop - offset} nucleotides) "
                f"downstream (max allowed distance is "
                f"{_MAX_DISTANCE_TO_STOP}) from the start of the slippery "
                f"sequence at location {offset + 1}."
            )

        seq = seq[: offset + _SLIPPERY_LEN] + seq[offset + _SLIPPERY_LEN - 1 :]

    # Pad with 'N' to avoid a 'BiopythonWarning: Partial codon' warning.
    remainder = len(seq) % 3
    seq += "N" * (3 - remainder if remainder else 0)

    if untranslatable:
        translation = []
        for index in range(0, len(seq), 3):
            codon = seq[index : index + 3]
            assert len(codon) == 3
            for char, replacement in untranslatable.items():
                if char in codon:
                    translation.append(replacement)
                    break
            else:
                translation.append(bpTranslate(codon))

        return "".join(translation)
    else:
        return bpTranslate(seq)


def translateSARS2Spike(seq: str) -> str:
    """
    Translate a SARS-CoV-2 Spike sequence, taking into account gaps introduced
    in the nucleotide alignment. This means that the amino acid sequences do
    not have to be re-aligned after translating, which avoids the introduction
    of gaps at the amino acid level that may be different from gaps at the
    nucleotide level.

    @param seq: A C{str} nucelotide sequence.
    @return: A translated C{str} amino acid sequence that retains the gaps.
    """
    sequence = ""
    current = 0
    seqLen = len(seq)

    if seqLen % 3:
        raise TranslatedSequenceLengthError(
            f"The length of a sequence to be translated must "
            f"be a multiple of 3 but is {seqLen!r}."
        )

    groups = groupby(seq)
    result = [(label, len(list(group))) for label, group in groups if label == "-"]

    for group, length in result:
        if length % 3 != 0:
            raise TranslatedGapLengthError(
                f"Length of stretch of gaps in {seq!r} is not divisible by 3.")

    # if any(length % 3 != 0 for g, length in result):
    # raise TranslatedGapLengthError("Length of stretch of gaps not divisible by 3.")

    while current + 3 <= seqLen:
        codon = seq[current : current + 3]
        codonGaps = codon.count("-")
        if codon in codons:
            # This is a codon that corresponds to a normal amino acid.
            sequence += codons[codon]
            current += 3
        elif codon not in codons and codonGaps == 0:
            # This is a codon that contains ambiguities.
            sequence += "X"
            current += 3
        elif codonGaps > 0:
            # Count how many gaps there are after the current codon.
            c = 3
            while seq[current + c] == "-":
                c += 1
            subsequentGaps = c - 3

            # Find the next nucleotide after the gap.
            index = 3 + codonGaps
            nextNt = seq[
                current + subsequentGaps + 3 : current + subsequentGaps + index
            ]

            newCodon = codon.strip("-") + nextNt

            sequence += codons.get(newCodon, "X")

            # Add the correct number of gaps after the amino acid.
            totalGaps = subsequentGaps + codonGaps
            gapAA = totalGaps // 3
            sequence += gapAA * "-"

            # Set current so it starts at the right place after the gap.
            current += subsequentGaps + index

    return sequence


_NO_COVERAGE = "no coverage "


def _summarizeNoCoverage(changes: List[str]) -> List[str]:
    """
    Summarize genome regions that have no coverage.

    @param changes: A C{list} of C{str} changes, as returned by
        getSubstitutionsString.
    @return: A C{list} of C{str} summarized changes.
    """
    result = []
    current = []
    for change in changes:
        if change.startswith(_NO_COVERAGE):
            current.append(change[len(_NO_COVERAGE) :])
        else:
            if current:
                result.append(_NO_COVERAGE + ", ".join(current))
                current = []
            result.append(change)

    if current:
        result.append(_NO_COVERAGE + ", ".join(current))

    return result


def getSubstitutionsString(
    referenceAa: AARead, genomeAa: AARead, summarizeNoCoverage: bool = False
) -> str:
    """
    Get a string with the substitutions.

    @param referenceAa: A C(dark.AARead) reference sequence.
    @param genomeAa: A C(dark.AARead) aligned sequence.
    @param summarizeNoCoverage: Instead of individually reporting each section
        of no coverage locations, summarize them. So instead of reporting
        "no coverage 33; no coverage 40-42", produce "no coverage 33, 40-42".
    @return: A C{str} summary of the substitutions from the reference
        to the genome.
    """
    changes = []
    previousXPosition = firstXposition = None
    site = refInsertCount = 0
    refSeq = referenceAa.sequence
    refSeqLen = len(refSeq)
    genSeq = genomeAa.sequence
    genSeqLen = len(genSeq)

    # Make sure zip doesn't silently truncate in case the sequences are not
    # the same length.
    if refSeqLen != genSeqLen:
        raise TranslatedReferenceAndGenomeLengthError(
            f"Reference and genome lengths unequal " f"({refSeqLen} != {genSeqLen})."
        )

    for site, (a, b) in enumerate(zip(refSeq, genSeq), start=1):
        if a == "-":
            refInsertCount += 1
        else:
            site -= refInsertCount
        if a != b:
            if b == "X":
                if previousXPosition == site - 1:
                    # There is already a string of Xs.
                    previousXPosition = site
                else:
                    # This is a new string of Xs.
                    changes.append(f"{_NO_COVERAGE}{site}")
                    previousXPosition = firstXposition = site
            else:
                if previousXPosition == site - 1:
                    assert changes[-1].startswith(_NO_COVERAGE)
                    if int(changes[-1][len(_NO_COVERAGE) :]) != site - 1:
                        changes[-1] += f"-{site - 1}"
                    # This is the first non-X after a string of Xs.
                    changes.append(f"{a}{site}{b}")
                    previousXPosition = site - 2
                else:
                    changes.append(f"{a}{site}{b}")

        elif previousXPosition == site - 1 and firstXposition != site - 1:
            assert changes[-1].startswith(_NO_COVERAGE)
            if int(changes[-1][len(_NO_COVERAGE) :]) != site - 1:
                changes[-1] += f"-{site - 1}"

    if previousXPosition is not None and previousXPosition == site:
        assert changes[-1].startswith(_NO_COVERAGE)
        if int(changes[-1][len(_NO_COVERAGE) :]) != site:
            changes[-1] += f"-{site}"

    if summarizeNoCoverage:
        changes = _summarizeNoCoverage(changes)

    return "; ".join(changes)


KNOWN_INSERTIONS = (
    ("QTKGIALSPR", 650, 700, 679, 683),
    ("NLVRTDRDLP", 200, 230, 214, 217),
    ("LVRAAGYLPQ", 200, 230, 214, 217),
    ("HRSSHQNLTP", 240, 255, 247, 250),  # Lineage B.3, -249H -250Q -251N
    ("GAAAYYYVG", 260, 270, 264, 265),  # Lineage B.1.1.257 -265Y
    ("SALLSDLQGTIT", 870, 900, 879, 882),  # -880D -881L -882Q
    ("TEKSKAENIIR", 93, 120, 98, 101),  # -99K -100A -101E
    ("MFVFFFVL", 0, 10, 3, 4),  # -4F L6F
    ("TWLLGSMHIS", 55, 80, 64, 65),  # -65L F66L H67G A68S I69M
    # ('LHRSSHQNLTP', 239, 254, 248, 251),     # -249H -250Q -251N
    ("NLVRAKKNDLP", 210, 240, 214, 218),  # -215A -216K -217K -218N
    ("NLVRVPWIDLP", 210, 240, 214, 218),  # -215V -216P -217W -218I
    ("VSQPFFFMD", 167, 182, 174, 175),  # -175F
    ("NLVRKLGPDLPQG", 207, 224, 214, 218),  # -215K -216L -217G -218P
    ("LHAPPATV", 515, 525, 520, 521),  # -521P
    ("LLHAAPPATV", 515, 540, 520, 522),  # -521A -522P
    ("LLPAAPPATV", 515, 540, 520, 522),  # -521A -522P
    ("LHAHPATV", 515, 540, 520, 521),  # -521H
    ("LLHAPPPATV", 515, 540, 520, 522),  # -521P -522P
    ("PGDSSSS", 248, 258, 253, 254),  # -254S
    ("AWNRKRISKRIS", 350, 390, 355, 359),  # -356K -357R -358I -359S
    ("LIGAAEH", 645, 660, 652, 653),  # -653A
    ("ALHRDSWGSY", 240, 260, 246, 250),  # -247D -248S -249W -250G
    ("CASYQTQTQT", 667, 680, 674, 676),  # -675Q -676T
    ("GSCCKFKFD", 1250, 1270, 1254, 1256),  # -1255K -1256F
    ("MLVFFF", 0, 10, 3, 4),  # -1M M2L
    ("SQCATLRINLT", 10, 30, 16, 20),  # -17T -18L -19R -20I
    ("FNDGVCVY", 84, 100, 89, 91),  # -90V -91C
    ("NPVLPLPFN", 78, 89, 84, 86),  # -85P -86L
    ("VLPFNVND", 79, 91, 86, 88),  # -87N -88V
    ("VLPFNDDDG", 80, 92, 87, 89),  # -88D -89D
    ("FDNPVLPLPF", 75, 90, 83, 85),  # -84L -85P
    ("LHRSSSLTYLT", 240, 255, 247, 251),  # -248S -249S -250L -251T
    ("KVCEFQFQ", 125, 137, 132, 134),  # -133F -134Q
    ("NLVRAQERDLP", 207, 225, 213, 217),  # -214R -215A -216Q -217E
    ("YHKNNNK", 140, 152, 147, 148),  # -148N
    ("IYKTPPP", 784, 797, 791, 792),  # -792P
    ("NLVRKRIDLP", 207, 220, 214, 217),  # -215K -216R -217I
    ("NLVRHSGDLP", 207, 220, 214, 217),  # -215H -216S -217G
    ("DEDDRCIDSE", 1250, 1270, 1259, 1263),  # -1260D -1261R -1262C -1263I
    ("FRVLDSS", 153, 166, 160, 161),  # -161D
    ("NLVRDLADLP", 206, 223, 214, 217),  # -215D -216L -217A
    ("LNDLCFTFTN", 380, 400, 391, 393),  # -392F -393T
    ("NLVRKFHDLP", 207, 222, 214, 217),  # -215K -216F -217H
    ("HRSSRWESVHLT", 240, 260, 247, 253),  # -248S -249R -250W -251E -252S
    # -253V
    ("AIHVVS", 64, 75, 69, 70),  # -70V
    ("NLVRANRNDLP", 207, 250, 214, 218),  # -215A -216N -217R -218N
    ("SGWTAASAG", 250, 265, 259, 262),  # -260A -261A -262S
    ("LVRDRSNLP", 207, 220, 215, 218),  # -216R -217S -218N
    ("FQTLHLL", 230, 250, 240, 242),  # -241L -242H
    ("LLACTPAT", 510, 530, 519, 520),  # -520C
    ("HRSFKTYLT", 240, 255, 247, 250),  # -248F -249K -250T
    ("THTNKAVRSPR", 670, 690, 679, 683),  # -680K -681A -682V -683R
    ("AYTGWLNSF", 20, 35, 29, 32),  # -30G -31W -32L
    ("SVITLTPG", 590, 605, 598, 600),  # -599T -600L
    ("LGVTSNH", 135, 150, 145, 146),  # Y144T Y145S -146N
    ("LGVTFNH", 135, 150, 145, 146),  # Y144T Y145F -146N
    ("LGVYSNH", 135, 150, 145, 146),  # Y145S -146N
    ("FLVVTSNH", 135, 150, 145, 146),  # G142V Y144T Y145S -146N
    ("FLGVISNH", 135, 150, 145, 146),  # Y144I Y145S -146N
    ("FLGVTSXH", 135, 150, 145, 146),  # Y144T Y145S -146X
    ("FLGVTXXX", 135, 150, 143, 144),  # -144T XXX
    # ('FLGVTSNHK', 130, 160, 143, 144),       # -144T Y145S Y146N
    ("YNSASSF", 365, 380, 372, 373),  # -373S
    ("LVRAAGAAGYLP", 205, 230, 215, 221),  # -216A -217G -218A -219A -220G
    # -221Y
    ("LVRASGYLP", 205, 230, 215, 218),  # -216S -217G -218Y
    ("HRSEIEAYL", 240, 255, 247, 251),  # -248E -249I -250E -251A
    ("SQPFFFM", 165, 180, 174, 175),  # -175F L177F
    ("LVRDNFGLPQG", 205, 222, 215, 218),  # -216N -217F -218G
    ("NLVRQASDL", 205, 222, 214, 217),  # -215Q -216A -217S
    ("QTQTFSRTNS", 650, 690, 677, 681),  # -678T -679F -680S -681R
    ("VGYLLQ", 260, 280, 269, 270),  # -270L
    ("CVADYYS", 355, 375, 364, 365),  # -365Y
    ("LVRERGLPQ", 207, 230, 215, 217),  # -216R -217G
    ("NRKRRIS", 320, 370, 356, 357),  # -357R
    ("GIVNNNT", 1120, 1140, 1133, 1134),  # -1134N
    ("PFNEDGV", 75, 100, 87, 88),  # -88E
    ("LEGAWLKQ", 170, 190, 182, 184),  # -183W -184L
    ("LVRAAGYLP", 200, 220, 215, 218),  # D215A -216A -217G -218Y
    ("KSNIIIR", 90, 110, 99, 100),  # -100I
    ("HADSQL", 610, 640, 627, 628),  # -628S
    ("HAIKHV", 60, 80, 68, 69),  # -69K
    ("FNCLHFP", 480, 500, 489, 490),  # Y489L -490H
    ("WTAGAAG", 250, 270, 259, 262),  # -260A -261G -262A
    ("SSSGWTGW", 250, 270, 256, 259),  # -257G -258W -259T
    ("ALHHGDRS", 240, 260, 245, 248),  # -246H -247G -248D
    ("HAPSAT", 510, 530, 521, 522),  # -522S
    ("FGTAVTLD", 100, 120, 108, 110),  # -109A -110V
    ("FRVINTTCYS", 150, 170, 159, 164),  # -160I -161N -162T -163T -164C
    ("LGVTYNN", 130, 160, 145, 146),  # Y144T -146N
    ("LGVTYNH", 130, 160, 145, 146),  # Y144T -146N
    ("MFVFFFF", 0, 10, 2, 4),  # -1M -2F M3V V5F
    ("NLVRQIDDLP", 210, 220, 214, 217),  # -215Q -216I -217D
    ("YNYLYFLRLF", 445, 465, 453, 455),  # -454F -455L
    ("FSNSVTW", 55, 70, 61, 62),  # -62S
    ("LVRKGEDLP", 205, 230, 214, 217),  # -215K -216G -217E
    ("SQSSIIA", 680, 700, 691, 692),  # -692S
    ("LVRDKLRSLPQG", 205, 230, 215, 219),  # -216K -217L -218R -219S
    ("MLTMLVFL", 0, 10, 1, 4),  # -1M -2L -3T F5L
    ("NLVRAPRDLPQ", 200, 230, 214, 217),  # -215A -216P -217R
    ("DEMINFTISAQY", 860, 880, 870, 875),  # -871N -872F -873T -874I -875S
    ("SQCVVREVRNLT", 0, 30, 16, 21),  # -17V -18R -19E -20V -21R
    ("SQPFLGECLMD", 165, 185, 175, 179),  # -176L -177G -178E -179C
    ("LVRKAFKQDLP", 205, 250, 214, 219),  # -215K -216A -217F -218K -219Q
    ("FLGVTXYH", 130, 185, 143, 144),  # -144T, X
    ("FLGVTYYH", 130, 185, 143, 144),  # -144T
    ("DPLYYPET", 290, 310, 297, 299),  # -298Y -299P
    ("LFRHYKYFKSN", 450, 490, 457, 462),  # -458H -459Y -460K -461Y -462F
    ("NLVRAVGYLP", 210, 225, 215, 218),  # -216V -217G -218Y
    ("GDSSSGSS", 250, 270, 254, 257),  # -255S -256S -257G
    ("GNFRQSKNL", 180, 200, 186, 189),  # -187R -188Q -189S
    ("CASYSLSYQT", 665, 690, 672, 676),  # -673S -674Y -675S -676L
    ("FPLPMCYESY", 480, 505, 492, 496),  # -493P -494M -495C -496Y
    ("NLVRDGRLLPQ", 205, 525, 215, 218),  # -216G -217R -218L
    ("FLGVTSY", 135, 155, 144, 145),  # Y144T -145S
    ("NLVRIDRDLP", 205, 230, 213, 216),  # -214R -215I -216D
    ("NLVRAAG", 205, 230, 215, 217),  # D215A -216A -217G
    ("LVRTDRDLP", 205, 230, 213, 216),  # -214R -215T -216D
    ("SSANNNC", 150, 180, 163, 164),  # -164N
    ("LVRDTRSLPQ", 200, 230, 215, 218),  # -216T -217R -218S
    ("VQPTLESIV", 310, 340, 323, 324),  # -324L
    ("MFVFFFFFFV", 0, 20, 3, 7),  # -4F -5F -6F -7F L9F
    ("TAGDGSDKSAAY", 250, 270, 261, 266),  # A262D -263G -264S -265D -266K
    # -267S
    ("LVRDQMRLPQG", 205, 225, 215, 218),  # -216Q -217M -218R
    ("HRSYPVGLT", 230, 260, 248, 251),  # -249P -250V -251G
    ("FLGVTSDH", 135, 150, 145, 146),  # Y144T Y145S -146D
    ("FLGVXSNH", 135, 150, 145, 146),  # no coverage 144 Y145S -146D
    ("NCVQPDY", 350, 380, 362, 364),  # -363Q A364P
    ("FLGVTXNH", 130, 185, 143, 144),  # -144T, X
    ("NLVRAAVYLP", 205, 230, 215, 218),  # -216A -217V -218Y
    ("VIHVISG", 50, 90, 70, 71),  # A67V -71I
    ("QTNGDSTSP", 670, 690, 679, 683),  # -680G -681D -682S -683T
    ("NLVREAGDLP", 205, 230, 214, 217),  # -215E -216A -217G
    ("NLVRWRRDLP", 205, 230, 214, 217),  # -215W -216R 217R
    ("NLVRNFGDLP", 205, 230, 214, 217),  # -215N -216F 217G
    ("NLVRWRRDLPQ", 205, 230, 213, 216),  # -214R -215W 216R
    ("VREPEDLP", 205, 230, 214, 217),  # -215E -216P -217E
    ("SNVRPGFQ", 630, 650, 642, 645),  # -643R -644P -645G
    ("QPFLMDLEGKQGNFKNLM", 170, 200, 175, 188),  # -176L -177M -178D -179L
    # -180E -181G -182K -183Q
    # -184G -185N -186F -187K
    # -188N
    ("RSFRRGRQGDL", 244, 260, 248, 255),  # Y248F -249R -250R -251G -252R
    # -253Q -254G -255D
    ("FASTDTE", 90, 110, 94, 96),  # -95T -96D
    ("IFGTAVXX", 100, 115, 110, 112),  # T109A L110V -111X -112X
    ("IHAISGT", 60, 80, 70, 71),  # V70A -71I
    ("FHARLHVS", 60, 80, 69, 70),  # I68R -69L
    ("GTNGTXXXXXXX", 60, 80, 76, 77),  # -77X
    ("HRSSGENLTP", 240, 260, 248, 251),  # Y248S -249G -250E -251N
)


def checkSpikeInsertions(accession: str, referenceAa: AARead, genomeAa: AARead) -> str:
    """
    Check and parse out known insertions in the Spike protein.

    @param accession: The C{str} accession number of the sequence.
    @param referenceAa: A C(dark.AARead) reference sequence.
    @param genomeAa: A C(dark.AARead) aligned genome sequence.
    @return: The C{str} corrected translated sequence.
    """
    seq = genomeAa.sequence
    seqLen = len(seq)
    if seqLen == 1274:
        # There are no insertions.
        return seq
    elif seqLen > 1274:
        for (target, findStart, findStop, sliceStart, sliceStop) in KNOWN_INSERTIONS:
            if seq.find(target, findStart, findStop) > -1:
                return seq[:sliceStart] + seq[sliceStop:]
        if seqLen < 1290:
            substitutions = getSubstitutionsString(referenceAa, genomeAa)
            raise TranslatedSequenceLengthError(
                f"Sequence with accession {accession} is too long "
                f"(length: {seqLen}) and does not have a known insertion. "
                f"Substitutions: {substitutions}"
            )
        else:
            raise TranslatedSequenceLengthError(
                f"Sequence with accession {accession} is too long "
                f"(length: {seqLen}) and does not have a known insertion."
            )
    else:
        raise TranslatedSequenceLengthError(
            f"Sequence with accession {accession} is too short " f"(length: {seqLen}."
        )
