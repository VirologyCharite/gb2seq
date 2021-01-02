from os.path import dirname, join
import gffutils
import urllib

from Bio.Seq import Seq

import sars2seq

_DATA_DIR = join(dirname(dirname(sars2seq.__file__)), 'data')


class Annotations:
    """
    Manage sequence annotations from the information in C{gffFile}.

    @param fastaFile: The C{str} file containing a SARS-CoV-2 genome.
    @param id_: The C{str} id of the annotated sequence. If C{None}, the
        first sequence in C{fastaFile} will be used.
    @param gffFile: The C{str} GFF file containing the annotations for the
        sequence in C{fastaFile}.
    """
    REF_ID = 'NC_045512.2'
    REF_FASTA = join(_DATA_DIR, 'NC_045512.2.fasta')
    REF_GFF = join(_DATA_DIR, 'GCF_009858895.2_ASM985889v3_genomic.gff')

    def __init__(self, fastaFile=None, gffFile=None):
        self.fastaFile = fastaFile or self.REF_FASTA
        self.products = self.readProducts(gffFile or self.REF_GFF)

    def readProducts(self, gffFile):
        """
        Read annotation products from C{gffFile}.

        @param gffFile: A C{str} GFF file.
        """
        db = gffutils.create_db(
            gffFile, dbfn=':memory:', keep_order=False,
            merge_strategy='merge', sort_attribute_values=False)

        result = {}

        for attr in ('CDS',):
            for ft in db.features_of_type(attr):
                print(f'{attr}:', ft)
                for child in db.children(ft):
                    key = (child.start - 1, child.stop)
                    print(f'  {attr} CHILD {key}:\n    {ft}')
                    try:
                        info = result[key]
                    except KeyError:
                        sequence = child.sequence(self.fastaFile)
                        info = result[key] = {
                            'names': set(),
                            'sequence': sequence,
                            'translation': str(Seq(sequence).translate()),
                        }
                    for subattr in 'gene', 'product', 'note':
                        try:
                            value = child.attributes[subattr][0]
                        except KeyError:
                            pass
                        else:
                            print(f'    {attr} attr {subattr}', value)
                            info['names'].add(value)

                for grandchild in db.children(child):
                    print('    GRANDCHILD:', grandchild)

        for attr in 'three_prime_UTR', 'five_prime_UTR', 'stem_loop':
            for ft in db.features_of_type(attr):
                print(ft.attributes)
                key = (ft.start - 1, ft.stop)
                assert key not in result
                result[key] = {
                    'names': {attr},
                    'sequence': sequence,
                }
                try:
                    result[key]['function'] = urllib.parse.unquote(
                        ft.attributes['function'][0])
                except KeyError:
                    pass

        print('ALL FEATURES:')
        for ft in db.all_features():
            print('--->', ft)

        return result

    def productDict(self):
        """
        Turn the products into a dictionary.

        @return: A C{dict} keyed by C{str} product name.
        """
        result = {}

        def _key(name):
            return len(name), name

        namesSeen = set()

        for ((start, stop), values) in self.products.items():
            names = sorted(values['names'], key=_key)
            name = names[0]
            if name in namesSeen:
                try:
                    result[f'{name} 1'] = result[name]
                    del result[name]
                except KeyError:
                    pass
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
            result[key]['names'] = names

        return result
