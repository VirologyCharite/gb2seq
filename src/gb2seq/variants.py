VARIANTS = {
    "spikeDeletion": {
        "description": "Spike amino acid deletions at locations 69 and 70",
        "changes": {
            "spike": {
                "aa": "69- 70-",
            },
        },
    },
    "N501Y": {
        "description": "Spike N501Y substitution",
        "changes": {
            "spike": {
                "aa": "N501Y",
            },
        },
    },
    "B117-typing": {
        "description": (
            "UK variant, based on just two tests (as in typing-PCR assays)."
        ),
        "changes": {
            "spike": {
                "aa": "69- 70- N501Y",
            },
        },
    },
    "VOC_20201201_UK": {
        "description": "UK variant of concern (VOC) 202012/01",
        "comment": (
            "From Table 1 of https://www.gov.uk/government/"
            "publications/investigation-of-novel-sars-cov-2"
            "-variant-variant-of-concern-20201201"
        ),
        "changes": {
            "orf1ab": {
                "aa": "1001I 1708D 2230T 3675- 3676- 3677-",
            },
            "spike": {
                "aa": "69- 70- 144- 501Y 570D 681H 716I 982A 1118H",
            },
            "orf8": {
                "aa": "Q27* 52I 73C",
            },
            "n": {
                "aa": "3L 235F",
            },
        },
    },
    "501Y.V2": {
        "description": "South African variant",
        "changes": {
            "e": {
                "aa": "P71L",
            },
            "orf1a": {
                "aa": "K1655N",
            },
            "spike": {
                "aa": "D80A D215G K417N E484K N501Y",
            },
            "n": {
                "aa": "T205I",
            },
        },
    },
    "victorSpike": {
        "description": "Spike amino acid changes of interest to Victor",
        "changes": {
            "spike": {
                "aa": (
                    "H69- V70- K417T K417N N439K Y453F E484K N501Y D614G "
                    "P681H V1176F"
                ),
            },
        },
    },
}
