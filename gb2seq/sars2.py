from gb2seq import Gb2SeqError


class NoSlipperySequenceError(Gb2SeqError):
    "No slippery sequence could be found in a genome."


class NoStopCodonError(Gb2SeqError):
    "No stop codon was found downstream from the slippery sequence."


class StopCodonTooDistantError(Gb2SeqError):
    "The stop codon following the slippery sequence was too far away."


# The maximum difference (number of nucleotides) to allow between the
# offset of the start of the slippery sequence and the downstream stop
# codon.
_MAX_DISTANCE_TO_STOP = 20

SLIPPERY_SEQUENCE = "TTTAAAC"

_SLIPPERY_LEN = len(SLIPPERY_SEQUENCE)


def getORF1abSequence(seq):
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
            f"the slippery sequence at site {offset + 1}."
        )
    if stop - offset > _MAX_DISTANCE_TO_STOP:
        raise StopCodonTooDistantError(
            f"The stop codon was too far ({stop - offset} nucleotides) "
            f"downstream (max allowed distance is "
            f"{_MAX_DISTANCE_TO_STOP}) from the start of the slippery "
            f"sequence at site {offset + 1}."
        )

    return seq[: offset + _SLIPPERY_LEN] + seq[offset + _SLIPPERY_LEN - 1 :]


# Provide convenient aliases for SARS-CoV-2 feature names. The alias is the
# key, the canonical name (as found in the GenBank file) is the value.
#
# Alphanumeric feature aliases must have lower case keys. If not they will not
# be detected (and the test suite will fail).


SARS_COV_2_ALIASES = {
    "2": "2'-O-ribose methyltransferase",
    "3clpro": "3C-like proteinase",
    "3utr": "3'UTR",
    "5utr": "5'UTR",
    "e": "envelope protein",
    "endornase": "endoRNAse",
    "envelope": "envelope protein",
    "exon": "3'-to-5' exonuclease",
    "exonuclease": "3'-to-5' exonuclease",
    "fcs": "furin cleavage site",
    "leader": "leader protein",
    "m": "membrane glycoprotein",
    "membrane": "membrane glycoprotein",
    "mpro": "3C-like proteinase",
    "n": "nucleocapsid phosphoprotein",
    "nsp1": "leader protein",
    "nsp5": "3C-like proteinase",
    "nsp12": "RNA-dependent RNA polymerase",
    "nsp13": "helicase",
    "nsp14": "3'-to-5' exonuclease",
    "nsp15": "endoRNAse",
    "nsp16": "2'-O-ribose methyltransferase",
    "orf4": "envelope protein",
    "orf5": "membrane glycoprotein",
    "orf1a": "ORF1a polyprotein",
    "orf1ab": "ORF1ab polyprotein",
    "orf3a": "ORF3a protein",
    "orf6": "ORF6 protein",
    "orf7a": "ORF7a protein",
    "orf7b": "ORF7b",
    "orf8": "ORF8 protein",
    "orf9": "nucleocapsid phosphoprotein",
    "orf10": "ORF10 protein",
    "rdrp": "RNA-dependent RNA polymerase",
    "s": "surface glycoprotein",
    "sl1": "stem loop 1",
    "sl2": "stem loop 2",
    "sl3": "stem loop 3",
    "sl4": "stem loop 4",
    "sl5": "stem loop 5",
    "s": "surface glycoprotein",
    "spike": "surface glycoprotein",
}

# Name of translated features, with (case sensitive!) names matching those in
# the GenBank file ../data/NC_045512.2.gb  These are in genome offset order.

SARS_COV_2_TRANSLATED = {
    "ORF1a polyprotein",
    "ORF1ab polyprotein",
    "leader protein",  # nsp1
    "nsp2",
    "nsp3",
    "nsp4",
    "3C-like proteinase",  # nsp5
    "nsp6",
    "nsp7",
    "nsp8",
    "nsp9",
    "nsp10",
    "nsp11",
    "RNA-dependent RNA polymerase",  # nsp12
    "helicase",  # nsp13
    "3'-to-5' exonuclease",  # nsp14
    "endoRNAse",  # nsp15
    "nsp16",
    "surface glycoprotein",
    "furin cleavage site",
    "ORF3a protein",
    "envelope protein",  # ORF4
    "membrane glycoprotein",  # ORF5
    "ORF6 protein",
    "ORF7a protein",
    "ORF7b",
    "ORF8 protein",
    "nucleocapsid phosphoprotein",  # ORF9
    "ORF10 protein",
}
