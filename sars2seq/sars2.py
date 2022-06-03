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
# the GenBank file ../data/NC_045512.2.gb

SARS_COV_2_TRANSLATED = {
    "3'-to-5' exonuclease",  # nsp14
    "3C-like proteinase",  # nsp5
    "endoRNAse",  # nsp15
    "envelope protein",  # ORF4
    "helicase",  # nsp13
    "leader protein",  # nsp1
    "membrane glycoprotein",  # ORF5
    "nsp2",
    "nsp3",
    "nsp4",
    "nsp6",
    "nsp7",
    "nsp8",
    "nsp9",
    "nsp10",
    "nsp11",
    "nucleocapsid phosphoprotein",  # ORF9
    "ORF1a polyprotein",
    "ORF1ab polyprotein",
    "ORF3a protein",
    "ORF6 protein",
    "ORF7a protein",
    "ORF7b",
    "ORF8 protein",
    "ORF10 protein",
    "RNA-dependent RNA polymerase",  # nsp12
    "surface glycoprotein",
}
