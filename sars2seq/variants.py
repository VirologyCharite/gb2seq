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
        'changes': {
            'orf1ab': {
                'aa': '1001I 1708D 2230T 3675- 3676- 3677-',
            },
            'spike': {
                'aa': '69- 70- 144- 501Y 570D 681H 716I 982A 1118H',
            },
            'orf8': {
                'aa': '27* 52I 73C',
            },
            'n': {
                'aa': '3L 235F',
            },
        },
    },
}
