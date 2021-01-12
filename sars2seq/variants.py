VARIANTS = {
    'spikeDeletion': {
        'description': 'Spike amino acid deletions at locations 69 and 70',
        'changes': {
            'spike': {
                'aa': '69- 70-',
            },
        },
    },

    'N501Y': {
        'description': 'Spike N501Y substitution',
        'changes': {
            'spike': {
                'aa': 'N501Y',
            },
        },
    },

    'VOC_20201201_UK': {
        'description': 'UK variant of concern (VOC) 202012/01',
        'comment': ('From Table 1 of https://www.gov.uk/government/'
                    'publications/investigation-of-novel-sars-cov-2'
                    '-variant-variant-of-concern-20201201'),
        # del:11288:9; del:21765:6; del:21991:3;
        'changes': {
            # orf1ab:T1001I; orf1ab:A1708D; orf1ab:I2230T;
            'orf1ab': {
                'aa': '1001I 1708D 2230T 3675- 3676- 3677-',
            },
            # S:N501Y; S:A570D; S:P681H; S:T716I; S:S982A; S:D1118H;
            'spike': {
                'aa': '69- 70- 144- 501Y 570D 681H 716I 982A 1118H',
            },
            # Orf8:Q27*; Orf8:R52I; Orf8:Y73C;
            'orf8': {
                'aa': '27* 52I 73C',
            },
            # N:D3L; N:S235F;
            'n': {
                'aa': '3L 235F',
            },
        },
    },

    '501Y.V2': {
        'description': 'South African variant',
        'changes': {
            'e': {
                'aa': 'P71L',
            },
            'orf1a': {
                'aa': 'K1655N',
            },
            'spike': {
                'aa': 'D80A D215G K417N E484K N501Y',
            },
            'n': {
                'aa': 'T205I',
            },
        },
    },
}


def isVOC_20201201Lineage(genome):
    """
    Does a genome fall into the VOC_20201201_UK lineage.
    """
    testCount, errorCount, _ = genome.checkVariant('VOC_20201201_UK')
    assert testCount == 20
