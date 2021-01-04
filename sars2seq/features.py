from os import environ
from os.path import dirname, exists, join

from Bio import Entrez, SeqIO

import sars2seq

# Set ENTREZ_EMAIL in your environment to have your requests to NCBI Entez
# be accompanied by your address. If you don't do this you'll see warning
# messages and be limited to a lower rate of querying.
Entrez.email = environ.get('ENTREZ_EMAIL')

_DATA_DIR = join(dirname(dirname(sars2seq.__file__)), 'data')


# Feature aliases should have a lower case key.
ALIASES = {
    '2': "2'-O-ribose methyltransferase",
    '3clpro': '3C-like proteinase',
    '3utr': "3'UTR",
    '5utr': "5'UTR",
    'e': 'envelope protein',
    'endornase': 'endoRNAse',
    'envelope': 'envelope protein',
    'exon': "3'-to-5' exonuclease",
    'exonuclease': "3'-to-5' exonuclease",
    'leader': 'leader protein',
    'm': 'membrane glycoprotein',
    'membrane': 'membrane glycoprotein',
    'mpro': '3C-like proteinase',
    'n': 'nucleocapsid phosphoprotein',
    'orf10': 'ORF10 protein',
    'orf1a': 'ORF1a polyprotein',
    'orf1ab': 'ORF1ab polyprotein',
    'orf3a': 'ORF3a protein',
    'orf6': 'ORF6 protein',
    'orf7a': 'ORF7a protein',
    'orf7b': 'ORF7b',
    'orf8': 'ORF8 protein',
    'rdrp': 'RNA-dependent RNA polymerase',
    's': 'surface glycoprotein',
    'sl1': 'stem loop 1',
    'sl2': 'stem loop 2',
    'sl3': 'stem loop 3',
    'sl4': 'stem loop 4',
    'sl5': 'stem loop 5',
    'spike': 'surface glycoprotein',
    'surface glycoprotein': 'surface glycoprotein',
}


class Features:
    """
    Manage sequence features from the information in C{gbFile}.

    @param gbFile: The C{str} Genbank file containing the features.
    """
    REF_GB = join(_DATA_DIR, 'NC_045512.2.gb')

    def __init__(self, gbFile=None):
        gbFile = gbFile or self.REF_GB

        if exists(gbFile):
            with open(gbFile) as fp:
                record = SeqIO.read(fp, 'genbank')
        else:
            import sys
            print(f'Fetching Genbank record for {gbFile!r}.',
                  file=sys.stderr)
            client = Entrez.efetch(db='nucleotide', rettype='gb',
                                   retmode='text', id=gbFile)
            record = SeqIO.read(client, 'gb')
            client.close()

        self.features = self._extractFeatures(record)
        self.id = record.id

    def _extractFeatures(self, record):
        """
        Extract features.

        @param record: A BioPython C{SeqRecord} sequence record.
        """
        result = {}

        for feature in record.features:
            type_ = feature.type

            if type_ == "3'UTR" or type_ == "5'UTR":
                start = int(feature.location.start)
                end = int(feature.location.end)
                key = start, end
                assert key not in result
                result[key] = {
                    'sequence': str(record.seq)[start:end],
                    'name': type_,
                }

            elif type_ == 'stem_loop':
                start = int(feature.location.start)
                end = int(feature.location.end)
                key = start, end
                assert key not in result
                result[key] = {
                    'sequence': str(record.seq)[start:end],
                    'function': feature.qualifiers['function'][0],
                    'name': 'stem loop',
                }

            elif type_ == 'CDS':
                start = int(feature.location.start)
                end = int(feature.location.end)
                key = start, end
                assert key not in result
                product = feature.qualifiers['product'][0]
                # print('CDS:', product, 'location', feature.location)
                result[key] = {
                    'product': product,
                    'translation': feature.qualifiers['translation'][0],
                    'sequence': str(record.seq)[start:end],
                    'name': product,
                }

                for optional in ('note',):
                    try:
                        result[key][optional] = feature.qualifiers[optional][0]
                    except KeyError:
                        pass

            elif type_ == 'mat_peptide':
                start = int(feature.location.start)
                end = int(feature.location.end)
                key = start, end
                product = feature.qualifiers['product'][0]
                # print('mat:', product, 'location', feature.location)
                sequence = str(record.seq)[start:end]
                info = {
                    'product': product,
                    'sequence': sequence,
                    'name': product,
                }

                for optional in ('note',):
                    try:
                        info[key][optional] = feature.qualifiers[optional][0]
                    except KeyError:
                        pass

                if key in result:
                    assert result[key] == info
                else:
                    result[key] = info

            elif type_ not in {'source', 'gene'}:
                raise ValueError(f'Unknown feature type {type_!r}.')

        return result

    def featuresDict(self):
        """
        Return features as a dictionary.

        @return: A C{dict} keyed by C{str} product name.
        """
        result = {}

        def _key(name):
            return len(name), name

        namesSeen = set()

        for ((start, stop), values) in self.features.items():
            name = values['name']
            if name in namesSeen:
                try:
                    result[f'{name} 1'] = result[name]
                except KeyError:
                    pass
                else:
                    del result[name]
                for i in range(2, 1000000):
                    key = f'{name} {i}'
                    if key not in result:
                        break
                else:
                    raise ValueError(f'Name {name} found a million times!?')
            else:
                namesSeen.add(name)
                key = name

            result[key] = {
                'name': name,
                'start': start,
                'stop': stop,
            }
            result[key].update(values)

        return result

    def canonicalName(self, name):
        """
        Get the canonical name for a feature.

        @param name: A C{str} feature name to look up.
        @raise KeyError: If the name is unknown.
        @return: A C{str} canonical name.
        """
        featuresDict = self.featuresDict()
        if name in featuresDict:
            return name

        nameLower = name.lower()
        for featureName in featuresDict:
            if nameLower == featureName.lower():
                return featuresDict[featureName]

        alias = ALIASES.get(nameLower)
        if alias:
            assert alias in featuresDict
            return alias

        raise KeyError(name)

    def getFeature(self, name):
        """
        Find a feature by name.

        @param name: A C{str} feature name to look up.
        @return: A C{dict} for the feature.
        """
        featuresDict = self.featuresDict()
        try:
            return featuresDict[name]
        except KeyError:
            nameLower = name.lower()
            for featureName in featuresDict:
                if nameLower == featureName.lower():
                    name = featureName
                    break
            else:
                alt = ALIASES.get(nameLower)
                if alt is None:
                    raise KeyError(name)
                else:
                    name = alt

        return featuresDict[name]
