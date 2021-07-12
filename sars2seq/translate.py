from Bio.Seq import Seq
from itertools import groupby

from dark.aa import CODONS, STOP_CODONS


class TranslationError(Exception):
    'Error when using custom translation of sequences.'


class NoSlipperySequenceError(TranslationError):
    'No slippery sequence could be found in a genome.'


class NoStopCodonError(TranslationError):
    'No stop codon was found downstream from the slippery sequence.'


class StopCodonTooDistantError(TranslationError):
    'The stop codon following the slippery sequence was too far away.'


class TranslatedSequenceLengthError(TranslationError):
    'The translated reference and genome have different lengths.'


class TranslatedReferenceAndGenomeLengthError(TranslationError):
    'The translated reference and genome have different lengths.'


class TranslatedGapLengthError(TranslationError):
    'The sequence to translate has stretches of gaps that are not '
    'multiples of three.'


codons = dict([(codon, aa) for aa, cods in CODONS.items() for codon in cods] +
              [(codon, '*') for codon in STOP_CODONS] + [('---', '-')])


# The maximum difference (number of nucleotides) to allow between the
# offset of the start of the slippery sequence and the downstream stop
# codon.
_MAX_DISTANCE_TO_STOP = 20

SLIPPERY_SEQUENCE = 'TTTAAAC'

_SLIPPERY_LEN = len(SLIPPERY_SEQUENCE)


def translate(seq, name=None):
    """
    Translate a sequence.

    @param seq: A C{str} nucelotide sequence.
    @param name: A C{str} feature name.
    @return: A translated C{str} amino acid sequence.
    """
    if name == 'ORF1ab polyprotein':
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
        stop = seq.find('TAA', offset + _SLIPPERY_LEN)
        if offset == -1:
            raise NoSlipperySequenceError('No slippery sequence found.')
        if stop == -1:
            raise NoStopCodonError(
                f'Could not find a stop codon downstream from the start of '
                f'the slippery sequence at location {offset + 1}.')
        if stop - offset > _MAX_DISTANCE_TO_STOP:
            raise StopCodonTooDistantError(
                f'The stop codon was too far ({stop - offset} nucleotides) '
                f'downstream (max allowed distance is '
                f'{_MAX_DISTANCE_TO_STOP}) from the start of the slippery '
                f'sequence at location {offset + 1}.')

        seq = seq[:offset + _SLIPPERY_LEN] + seq[offset + _SLIPPERY_LEN - 1:]

    # Pad with 'N' to avoid a 'BiopythonWarning: Partial codon' warning.
    remainder = len(seq) % 3
    seq += 'N' * (3 - remainder if remainder else 0)

    return Seq(seq).translate()


def translateSpike(seq):
    """
    Translate a Spike sequence, taking into account gaps introduced in the
    nucleotide alignment. This means that the amino acid sequences do not
    have to be re-aligned after translating, which avoids the introduction of
    gaps at the amino acid level that may be different from gaps at the
    nucleotide level.

    @param seq: A C{str} nucelotide sequence.
    @return: A translated C{str} amino acid sequence that retains the gaps.
    """
    sequence = ''
    current = 0
    seqLen = len(seq)

    if not seqLen % 3 == 0:
        raise TranslatedSequenceLengthError(
            f'The length of a sequence to be translated must '
            f'be a multiple of 3 but is {seqLen!r}.')

    groups = groupby(seq)
    result = [(label, len(list(group))) for label, group in groups if
              label == '-']

    if any(length % 3 != 0 for g, length in result):
        raise TranslatedGapLengthError(
            'Length of stretch of gaps not divisible by 3.')

    while current + 3 <= seqLen:
        codon = seq[current:current + 3]
        codonGaps = codon.count('-')
        if codon in codons:
            # This is a codon that corresponds to a normal amino acid.
            sequence += codons[codon]
            current += 3
        elif codon not in codons and codonGaps == 0:
            # This is a codon that contains ambiguities.
            sequence += 'X'
            current += 3
        elif codonGaps > 0:
            # Count how many gaps there are after the current codon.
            c = 3
            while seq[current + c] == '-':
                c += 1
            subsequentGaps = c - 3

            # Find the next nucleotide after the gap.
            index = 3 + codonGaps
            nextNt = seq[current + subsequentGaps + 3:
                         current + subsequentGaps + index]

            newCodon = codon.strip('-') + nextNt

            sequence += codons.get(newCodon, 'X')

            # Add the correct number of gaps after the amino acid.
            totalGaps = subsequentGaps + codonGaps
            gapAA = totalGaps // 3
            sequence += (gapAA * '-')

            # Set current so it starts at the right place after the gap.
            current += (subsequentGaps + index)

    return sequence


def checkSpikeInsertions(accession, seq):
    """
    Check and parse out known insertions in the Spike protein.

    @param accession: The C{str} accession number of the sequence.
    @param seq: A C{str} amino acid sequence.

    @return: The C{str} corrected translated sequence.
    """
    seqLen = len(seq)
    if seqLen == 1274:
        # There are no insertions
        return seq
    if seqLen > 1274:
        if seq.find('QTKGIALSPR', 650, 700) > -1:
            return seq[:679] + seq[683:]
        elif seq.find('NLVRTDRDLP', 200, 230) > -1:
            return seq[:214] + seq[217:]
        elif seq.find('LVRAAGYLPQ', 200, 230) > -1:
            return seq[:214] + seq[217:]
        elif seq.find('HRSSHQNLTP', 240, 255) > -1:  # Lineage B.3
            return seq[:247] + seq[250:]
        elif seq.find('GAAAYYYVG', 260, 270) > -1:   # Lineage B.1.1.257 -265Y
            return seq[:264] + seq[265:]
        elif seq.find('SALLSDLQGTIT', 870, 890) > -1:   # -880D -881L -882Q
            return seq[:879] + seq[882:]
        elif seq.find('TEKSKAENIIR', 93, 103) > -1:  # -99K -100A -101E
            return seq[:98] + seq[101:]
        elif seq.find('MFVFFFVL', 0, 10) > -1:   # -4F L6F
            return seq[:3] + seq[4:]
        elif seq.find('TWLLGSMHIS', 55, 80) > -1:  # -65L F66L H67G A68S I69M
            return seq[:64] + seq[65:]
        elif seq.find('LHRSSHQNLTP', 239, 254) > -1:  # -249H -250Q -251N
            return seq[:248] + seq[251:]
        elif seq.find('NLVRAKKNDLP', 210, 220) > -1:  # -215A -216K -217K -218N
            return seq[:214] + seq[218:]
        elif seq.find('VSQPFFFMD', 167, 182) > -1:  # -175F
            return seq[:174] + seq[175:]
        elif seq.find('NLVRKLGPDLPQG', 207, 221) > -1:  # -215K -216L -217G
            return seq[:214] + seq[218:]                # -218P
        elif seq.find('LHAPPATV', 515, 525) > -1:   # -521P
            return seq[:520] + seq[521:]
        elif seq.find('LLHAAPPATV', 515, 525) > -1:  # -521A -522P
            return seq[:520] + seq[522:]
        elif seq.find('LHAHPATV', 515, 525) > -1:   # -521H
            return seq[:520] + seq[521:]
        elif seq.find('LLHAPPPATV', 515, 525) > -1:  # -521P -522P
            return seq[:520] + seq[522:]
        elif seq.find('PGDSSSS', 248, 258) > -1:  # -254S
            return seq[:253] + seq[254:]
        elif seq.find('AWNRKRISKRIS', 350, 360) > -1:  # -356K -357R -358I
            return seq[:355] + seq[359:]               # -359S
        elif seq.find('LIGAAEH', 645, 660) > -1:  # -653A
            return seq[:652] + seq[653]
        elif seq.find('ALHRDSWGSY', 240, 250) > -1:  # -247D -248S -249W -250G
            return seq[:246] + seq[250:]
        elif seq.find('CASYQTQTQT', 667, 680) > -1:  # -675Q -676T
            return seq[:674] + seq[676:]
        elif seq.find('GSCCKFKFD', 1250, 1260) > -1:  # -1255K -1256F
            return seq[:1254] + seq[1256]
        elif seq.find('MLVFFF', 0, 10) > -1:  # -1M M2L
            return seq[:3] + seq[4:]
        elif seq.find('SQCATLRINLT', 10, 20) > -1:  # -17T -18L -19R -20I
            return seq[:16] + seq[20:]
        elif seq.find('FNDGVCVY', 84, 95) > -1:  # -90V -91C
            return seq[:89] + seq[91:]
        elif seq.find('NPVLPLPFN', 78, 89) > -1:  # -85P -86L
            return seq[:84] + seq[86:]
        elif seq.find('VLPFNVND', 79, 91) > -1:  # -87N -88V
            return seq[:86] + seq[88:]
        elif seq.find('VLPFNDDDG', 80, 92) > -1:  # -88D -89D
            return seq[:87] + seq[:89]
        elif seq.find('FDNPVLPLPF', 75, 90) > -1:  # -84L -85P
            return seq[:83] + seq[85:]
        elif seq.find('LHRSSSLTYLT', 240, 255) > -1:  # -248S -249S -250L -251T
            return seq[:247] + seq[251:]
        elif seq.find('KVCEFQFQ', 125, 137) > -1:  # -133F -134Q
            return seq[:132] + seq[134:]
        elif seq.find('NLVRAQERDLP', 207, 220) > -1:  # -214R -215A -216Q -217E
            return seq[:213] + seq[217:]
        elif seq.find('YHKNNNK', 140, 152) > -1:  # -148N
            return seq[:147] + seq[148:]
        elif seq.find('IYKTPPP', 784, 797) > -1:  # -792P
            return seq[:791] + seq[792:]
        elif seq.find('NLVRKRIDLP', 207, 220) > -1:  # -215K -216R -217I
            return seq[:214] + seq[217:]
        elif seq.find('DEDDRCIDSE', 1250, 1265) > -1:  # -1260D -1261R -1262C
            return seq[:1259] + seq[1263:]             # -1263I
        elif seq.find('FRVLDSS', 153, 166) > -1:  # -161D
            return seq[:160] + seq[161:]
        elif seq.find('NLVRDLADLP', 206, 223) > -1:  # -215D -216L -217A
            return seq[:214] + seq[217:]
        elif seq.find('LNDLCFTFTN', 380, 400) > -1:  # -392F -393T
            return seq[:391] + seq[393:]
        elif seq.find('NLVRKFHDLP', 207, 222) > -1:  # -215K -216F -217H
            return seq[:214] + seq[217:]
        elif seq.find('HRSSRWESVHLT', 240, 260) > -1:  # -248S -249R -250W
            return seq[:247] + seq[253:]               # -251E -252S -253V
        elif seq.find('AIHVVS', 64, 75) > -1:  # -70V
            return seq[:69] + seq[70:]
        elif seq.find('NLVRANRNDLP', 207, 220) > -1:  # -215A -216N -217R -218N
            return seq[:214] + seq[218:]
        elif seq.find('SGWTAASAG', 250, 265) > -1:  # -260A -261A -262S
            return seq[:259] + seq[262:]
        elif seq.find('LVRDRSNLP', 207, 220) > -1:  # -216R -217S -218N
            return seq[:215] + seq[218:]
        elif seq.find('FQTLHLL', 230, 250) > -1:  # -241L -242H
            return seq[:240] + seq[242:]
        elif seq.find('LLACTPAT', 230, 250) > -1:  # -520C
            return seq[:519] + seq[520:]
        elif seq.find('HRSFKTYLT', 240, 255) > -1:  # -248F -249K -250T
            return seq[:247] + seq[250:]
        elif seq.find('THTNKAVRSPR', 670, 690) > -1:  # -680K -681A -682V -683R
            return seq[:679] + seq[683:]
        elif seq.find('AYTGWLNSF', 20, 35) > -1:  # -30G -31W -32L
            return seq[:29] + seq[32:]
        elif seq.find('SVITLTPG', 590, 605) > -1:  # -599T -600L
            return seq[:598] + seq[600:]
        elif seq.find('FLGVTSNH', 135, 150) > -1:  # Y144T Y145S -146N
            return seq[:145] + seq[146:]
        elif seq.find('YNSASSF', 365, 380) > -1:  # -373S
            return seq[:372] + seq[373:]
        elif seq.find('LVRAAGAAGYLP', 205, 230) > -1:  # -216A -217G -218A
            return seq[:215] + seq[221:]               # -219A -220G -221Y
        elif seq.find('LVRASGYLP', 205, 230) > -1:  # -216S -217G -218Y
            return seq[:215] + seq[218:]
        elif seq.find('HRSEIEAYL', 240, 255) > -1:  # -248E -249I -250E -251A
            return seq[:247] + seq[251:]
        elif seq.find('SQPFFFM', 165, 180) > -1:  # -175F L177F
            return seq[:174] + seq[175:]
        elif seq.find('LVRDNFGLPQG', 205, 222) > -1:  # -216N -217F -218G
            return seq[:215] + seq[218:]
        elif seq.find('NLVRQASDL', 205, 222) > -1:  # -215Q -216A -217S
            return seq[:214] + seq[217:]
        elif seq.find('QTQTFSRTNS', 205, 222) > -1:  # -678T -679F -680S -681R
            return seq[:677] + seq[681:]
        elif seq.find('VGYLLQ', 260, 280) > -1:  # -270L
            return seq[:269] + seq[270:]
        elif seq.find('CVADYYS', 355, 375) > -1:  # -365Y
            return seq[:364] + seq[365:]
        elif seq.find('LVRERGLPQ', 207, 230) > -1:  # -216R -217G
            return seq[:215] + seq[217:]
        elif seq.find('NRKRRIS', 207, 230) > -1:  # -357R
            return seq[:356] + seq[357:]
        elif seq.find('GIVNNNT', 1120, 1140) > -1:  # -1134N
            return seq[:1133] + seq[1134:]
        elif seq.find('PFNEDGV', 75, 100) > -1:  # -88E
            return seq[:87] + seq[88:]
        elif seq.find('LEGAWLKQ', 170, 190) > -1:  # -183W -184L
            return seq[:182] + seq[184:]
        elif seq.find('LVRAAGYLP', 200, 220) > -1:  # D215A -216A -217G -218Y
            return seq[:215] + seq[218:]
        elif seq.find('KSNIIIR', 90, 110) > -1:  # -100I
            return seq[:99] + seq[100:]
        elif seq.find('HADSQL', 610, 640) > -1:  # -628S
            return seq[:627] + seq[628:]
        elif seq.find('HAIKHV', 60, 80) > -1:  # -69K
            return seq[:68] + seq[69:]
        elif seq.find('FNCLHFP', 480, 500) > -1:  # Y489L -490H
            return seq[:489] + seq[490:]
        elif seq.find('WTAGAAG', 250, 270) > -1:  # -260A -261G -262A
            return seq[:259] + seq[262:]
        elif seq.find('SSSGWTGW', 250, 270) > -1:  # -257G -258W -259T
            return seq[:256] + seq[260:]
        elif seq.find('ALHHGDRS', 240, 260) > -1:  # -246H -247G -248D
            return seq[:245] + seq[248:]
        elif seq.find('HAPSAT', 510, 530) > -1:  # -522S
            return seq[:521] + seq[522:]
        elif seq.find('FGTAVTLD', 100, 120) > -1:  # -109A -110V
            return seq[:108] + seq[110:]
        elif seq.find('FRVINTTCYS', 150, 170) > -1:  # -160I -161N -162T
            return seq[:159] + seq[164:]             # -163T -164C
        elif seq.find('LGVTYNN', 130, 160) > -1:  # Y144T -146N
            return seq[:145] + seq[146:]
        elif seq.find('MFVFFFF', 0, 10) > -1:  # -1M -2F M3V V5F
            return seq[:2] + seq[4:]
        elif seq.find('NLVRQIDDLP', 210, 220) > -1:  # -215Q -216I -217D
            return seq[:214] + seq[217:]
        elif seq.find('YNYLYFLRLF', 445, 465) > -1:  # -454F -455L
            return seq[:453] + seq[455:]
        elif seq.find('FSNSVTW', 55, 70) > -1:  # -62S
            return seq[:61] + seq[62:]
        elif seq.find('LVRKGEDLP', 205, 230) > -1:  # -215K -216G -217E
            return seq[:214] + seq[217:]
        elif seq.find('SQSSIIA', 680, 700) > -1:  # -692S
            return seq[:691] + seq[692:]
        elif seq.find('LVRDKLRSLPQG', 680, 700) > -1:  # -216K -217L -218R
            return seq[:215] + seq[219:]               # -219S
        elif seq.find('MLTMLVFL', 0, 10) > -1:  # -1M -2L -3T F5L
            return seq[:1] + seq[4:]
        elif seq.find('NLVRAPRDLPQ', 200, 230) > -1:  # -215A -216P -217R
            return seq[:214] + seq[217:]
        elif seq.find('DEMINFTISAQY', 860, 880) > -1:  # -871N -872F -873T
            return seq[:870] + seq[875:]               # -874I -875S
        elif seq.find('SQCVVREVRNLT', 0, 30) > -1:  # -17V -18R -19E -20V -21R
            return seq[:16] + seq[21:]
        elif seq.find('SQPFLGECLMD', 165, 185) > -1:  # -176L -177G -178E -179C
            return seq[:175] + seq[179:]
        elif seq.find('LVRKAFKQDLP', 205, 220) > -1:  # -215K -216A -217F
            return seq[:214] + seq[219:]              # -218K -219Q
        elif seq.find('FLGVTXYH', 165, 185) > -1:  # -144T, X
            return seq[:143] + seq[144:]
        elif seq.find('FLGVTYYH', 165, 185) > -1:  # -144T
            return seq[:143] + seq[144:]
        elif seq.find('DPLYYPET', 290, 310) > -1:  # -298Y -299P
            return seq[:297] + seq[299:]
        elif seq.find('LVRHYKYFKSN', 450, 470) > -1:  # -458H -459Y -460K -461Y
            return seq[:457] + seq[462:]              # -462F
        elif seq.find('NLVRAVGYLP', 210, 225) > -1:  # -216V -217G -218Y
            return seq[:215] + seq[218:]
        elif seq.find('GDSSSGSS', 250, 270) > -1:  # -255S -256S -257G
            return seq[:254] + seq[257:]
        elif seq.find('GNFRQSKNL', 180, 200) > -1:  # -187R -188Q -189S
            return seq[:186] + seq[189:]
        elif seq.find('CASYSLSYQT', 665, 690) > -1:  # -673S -674Y -675S -676L
            return seq[:672] + seq[676:]
        elif seq.find('FPLPMCYESY', 480, 505) > -1:  # -493P -494M -495C -496Y
            return seq[:492] + seq[496:]
        elif seq.find('NLVRDGRLLPQ', 205, 525) > -1:  # -216G -217R -218L
            return seq[:215] + seq[218:]
        elif seq.find('FLGVTSY', 135, 155) > -1:  # Y144T -145S
            return seq[:144] + seq[145:]
        elif seq.find('INLVRIDFDLP', 205, 225) > -1:  # -214R -215I -216D
            return seq[:213] + seq[216:]
        elif seq.find('SSANNNC', 150, 180) > -1:  # -164N
            return seq[:163] + seq[164:]
        elif seq.find('LVRDTRSLPQ', 200, 230) > -1:  # -216T -217R -218S
            return seq[:215] + seq[218:]
        elif seq.find('VQPTLESIV', 310, 340) > -1:  # -324L
            return seq[:323] + seq[324:]
        elif seq.find('MFVFFFFFFV', 0, 20) > -1:  # -4F -5F -6F -7F L9F
            return seq[:3] + seq[7:]
        elif seq.find('TAGDGSDKSAAY', 250, 270) > -1:  # A262D -263G -264S
            return seq[:261] + seq[266:]               # -265D -266K -267S
        elif seq.find('LVRDQMRLPQG', 205, 225) > -1:  # -216Q -217M -218R
            return seq[:215] + seq[218:]
        elif seq.find('HRSYPVGLT', 230, 260) > -1:  # -249P -250V -251G
            return seq[:248] + seq[251:]
        elif seq.find('FLGVTSDH', 135, 150) > -1:  # Y144T Y145S -146D
            return seq[:145] + seq[146:]
        elif seq.find('NCVQPDY', 350, 380) > -1:  # -363Q A364P
            return seq[:362] + seq[364:]
        elif seq.find('FLGVTXNH', 165, 185) > -1:  # -144T, X
            return seq[:143] + seq[144:]
        elif seq.find('NLVRAAVYLP', 205, 230) > -1:  # -216A -217V -218Y
            return seq[:215] + seq[218:]
        elif seq.find('VIHVISG', 205, 230) > -1:  # A67V -71I
            return seq[:70] + seq[71:]
        elif seq.find('QTNGDSTSP', 670, 690) > -1:  # -680G -681D -682S -683T
            return seq[:679] + seq[683:]
        elif seq.find('NLVREAGDLP', 205, 230) > -1:  # -215E -216A -217G
            return seq[:214] + seq[217:]
        elif seq.find('FLGVTSNHK', 130, 160) > -1:  # -144T Y145S Y146N
            return seq[:143] + seq[144:]
        elif seq.find('FLGVTSNHK', 130, 160) > -1:  # -144T Y145S Y146N
            return seq[:143] + seq[144:]
        else:
            raise TranslatedSequenceLengthError(
                f'Sequence with accession {accession} is too long '
                f'(length: {seqLen}) and does not have a known insertion.',
            )
    raise TranslatedSequenceLengthError(
        f'Sequence with accession {accession} is too short '
        f'(length: {seqLen}.'
    )
