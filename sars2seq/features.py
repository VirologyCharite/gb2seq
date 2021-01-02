from os import environ
from os.path import dirname, exists, join

from Bio import Entrez, SeqIO
# from Bio.Seq import Seq

import sars2seq

Entrez.email = environ.get('ENTREZ_EMAIL')

_DATA_DIR = join(dirname(dirname(sars2seq.__file__)), 'data')


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
                'start': start,
                'stop': stop,
            }
            result[key].update(values)

        return result
