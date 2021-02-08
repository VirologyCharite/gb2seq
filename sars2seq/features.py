import sys
from os import environ
from os.path import dirname, exists, join
import itertools
from collections import UserDict

from Bio import Entrez, SeqIO

from dark.aa import STOP_CODONS
from dark.reads import DNARead

import sars2seq

# Set ENTREZ_EMAIL in your environment to have your requests to NCBI Entez
# be accompanied by your address. If you don't do this you'll see warning
# messages and be limited to a lower rate of querying.
Entrez.email = environ.get('ENTREZ_EMAIL')

_DATA_DIR = join(dirname(dirname(sars2seq.__file__)), 'data')


class ReferenceWithGapError(Exception):
    'A GenBank reference sequence had a gap.'


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
    'nsp1': 'leader protein',
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


class Features(UserDict):
    """
    Manage sequence features from the information in C{gbFile}.

    @param spec: Either:
        * A C{str} name of a Genbank file containing the features.
        * A C{str} Genbank accession id.
        * A C{dict} of pre-prepared features, in which case C{reference}
              must not be C{None}. Passing a C{dict} is provided for testing.
        * C{None}, in which case the default reference, NC_045512.2.gb, is
              loaded.
    @param reference: A C{dark.reads.DNARead} instance if C{spec} is a C{dict},
        or C{None}.
    @raise ReferenceWithGapError: If the reference or one of its
        features has a gap in its nucleotide sequence.
    """
    REF_GB = join(_DATA_DIR, 'NC_045512.2.gb')

    def __init__(self, spec=None, reference=None):
        super().__init__()
        spec = self.REF_GB if spec is None else spec

        if isinstance(spec, str):
            if exists(spec):
                assert reference is None
                with open(spec) as fp:
                    record = SeqIO.read(fp, 'genbank')
            else:
                assert reference is None
                print(f'Fetching Genbank record for {spec!r}.',
                      file=sys.stderr)
                client = Entrez.efetch(db='nucleotide', rettype='gb',
                                       retmode='text', id=spec)
                record = SeqIO.read(client, 'gb')
                client.close()
            self._initializeFromGenBankRecord(record)
        elif isinstance(spec, dict):
            self.update(spec)
            self.reference = reference
        else:
            raise ValueError(f'Unrecognized specification {spec!r}')

    def _initializeFromGenBankRecord(self, record):
        """
        Initialize from a GenBank record.

        @param record: A BioPython C{SeqRecord} sequence record.
        """
        self.reference = DNARead(record.id, record.seq)

        for feature in record.features:
            type_ = feature.type
            start = int(feature.location.start)
            stop = int(feature.location.end)
            value = {}

            if type_ == "3'UTR" or type_ == "5'UTR":
                name = type_

            elif type_ == 'stem_loop':
                for n in itertools.count(1):
                    name = f'stem loop {n}'
                    if name not in self:
                        break
                value['function'] = feature.qualifiers['function'][0]

            elif type_ in {'CDS', 'mat_peptide'}:
                name = feature.qualifiers['product'][0]
                value['product'] = name

            elif type_ in {'source', 'gene'}:
                continue

            else:
                raise ValueError(f'Unknown feature type {type_!r}.')

            value.update({
                'name': name,
                'sequence': str(record.seq)[start:stop],
                'start': start,
                'stop': stop,
            })

            for optional in 'translation', 'note':
                try:
                    value[optional] = feature.qualifiers[optional][0]
                except KeyError:
                    pass

            # If there is a translation, add an amino acid '*' stop
            # indicator if there is not one already and the sequence ends
            # with a stop codon.
            try:
                translation = value['translation']
            except KeyError:
                pass
            else:
                if (not translation.endswith('*') and
                        value['sequence'][-3:].upper() in STOP_CODONS):
                    value['translation'] += '*'

            if name in self:
                assert self[name] == value
            else:
                self[name] = value

        self._checkForGaps()

    def canonicalName(self, name):
        """
        Get the canonical name for a feature.

        @param name: A C{str} feature name to look up.
        @raise KeyError: If the name is unknown.
        @return: A C{str} canonical name.
        """
        if name in self:
            return name

        nameLower = name.lower()
        for featureName in self:
            if nameLower == featureName.lower():
                return nameLower

        alias = ALIASES.get(nameLower)
        if alias:
            assert alias in self
            return alias

        raise KeyError(name)

    def __getitem__(self, name):
        """
        Find a feature by name.

        @param name: A C{str} feature name to look up.
        @return: A C{dict} for the feature.
        """
        try:
            return self.data[name]
        except KeyError:
            nameLower = name.lower()
            for featureName in self.data:
                if nameLower == featureName.lower():
                    name = featureName
                    break
            else:
                alt = ALIASES.get(nameLower)
                if alt is None:
                    raise KeyError(name)
                else:
                    name = alt

        return self.data[name]

    def _checkForGaps(self):
        """
        Check there are no gaps in the reference or and feature sequence.

        @raise ReferenceWithGapError: If the reference sequence or the
            sequence of any of its features has a gap.
        """
        referenceId = self.reference.id
        if self.reference.sequence.find('-') > -1:
            raise ReferenceWithGapError(
                f'Reference sequence {referenceId!r} has a gap!')

        for featureInfo in self.values():
            if featureInfo['sequence'].find('-') > -1:
                raise ReferenceWithGapError(
                    f'Feature {featureInfo["name"]!r} sequence in '
                    f'{referenceId!r} has a gap!')
