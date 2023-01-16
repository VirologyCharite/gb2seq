# A Python library and command-line scripts for working with GenBank genome sequences and annotations

`gb2seq` provides utility scripts for examining GenBank genome sequences
and their annotations, and a Python library API for working with both.

The library code knows how to extract genes at the nucleotide or amino acid
level, how to identify or check for expected genome changes, how to
translate between aa and nt offsets, and about all genome features.

## Conventions

### SARS-CoV-2 genomes

This code was written in late 2020 for processing SARS-CoV-2 genomes. It
has since been generalized to deal with any GenBank file. The command-line
scripts can be given a `--sars2` to have them use the Wuhan reference, in
which case there is no need to provide a GenBank file. Similarly, the
Python classes can take a `sars2` argument with the same effect when it is
set to `True`. If you're not working with SARS-CoV-2 genomes, just pass the
GenBank file for your genome to the command-line scripts with `--reference`
or to the Python functions. If you are working with SARS-CoV-2, and you
want to use the Wuhan reference, you can use the shortcut just described.

### Locations vs offsets

When giving locations within a genome or a genome feature (e.g., a
particular protein), a "site" refers to a 1-based location, as would
normally be used by a regular person, whereas an "offset" refers to a
0-based location, as would be used by a Python programmer.  The utility
scripts, run from the command line, expect the former, whereas library code
called by a programmer expects the latter. This terminology is used
throughout the code as well.

### Argument and return value ordering

In functions that deal with both a reference and a genome sequence, the
reference is always given or returned first.

## Alignment

`gb2seq` makes an alignment between each genome you give it and a
reference sequence.

The default alignment algorithm is
[MAFFT](https://mafft.cbrc.jp/alignment/software/), which can be slow but
is reliable. You can also run scripts with `--aligner edlib` (or pass
`aligner='edlib'` to library functions) to use [the Python
wrapper](https://pypi.org/project/edlib/) for the extremely fast
[edlib](https://github.com/Martinsos/edlib) library.

# Utility scripts

The scripts described below all accept a `--help` option. Below is
some representative usage.

## describe-feature.py

The simplest script is `describe-feature.py`, which can be used to get
information about features in a reference genome specified by passing the
name of a GenBank file with `--gbFile` (or just use `--sars2` to if you
want the Wuhan reference
[NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512)).

If you run `describe-feature.py` with no arguments, it will print summary
details of all features, with their (0-based) nucleotide offsets. E.g.,

```sh
$ describe-feature.py --sars2
Features for NC_045512.2:
2'-O-ribose methyltransferase:
  start: 20658
  stop: 21552
  length: 894
  product: 2'-O-ribose methyltransferase
  sequence    (len   894 nt): TCTAGTCAAGCGTGGCAACCGGGTGTTGCTATGCCTAATCTTTACAAAA...
3'-to-5' exonuclease:
  start: 18039
  stop: 19620
  length: 1581
  product: 3'-to-5' exonuclease
  sequence    (len  1581 nt): GCTGAAAATGTAACAGGACTCTTTAAAGATTGTAGTAAGGTAATCACTG...
3'UTR:
  start: 29674
  stop: 29903
  length: 229
  sequence    (len   229 nt): CAATCTTTAATCAGTGTGTAACATTAGGGAGGACTTGAAAGAGCCACCA...
# [Many additional output lines omitted here.]
```

You can also pass it a feature name (or multiple feature names):

```sh
$ describe-feature.py --sars2 --name spike
surface glycoprotein:
  start: 21562
  stop: 25384
  length: 3822
  product: surface glycoprotein
  sequence    (len  3822 nt): ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTA...
  translation (len  1274 aa): MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLH...
```

And you can pass in a genome to get information on the feature in both the
reference and the passed genome (which will be aligned to the
reference). The offset and length of the feature in the genome may of
course differ from the reference:

```sh
$ describe-feature.py --name s --sars2 --genome CSpecVir9290.fasta
Reference:
  surface glycoprotein:
    start: 21562
    stop: 25384
    length (nt): 3822
    product: surface glycoprotein
    note: structural protein; spike protein
    feature is translated left-to-right.
    sequence: ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGC...
    length (aa): 1274
    translation: MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFD...
  Genome BetaCoV/Rendsburg-Eckernfoerde/ChVir9290/2020:
    start: 21550
    stop: 25366
    length (nt): 3822
    sequence: ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGC...
    translation: MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAISGTNGTKRFDNP...
```

Or ask to see all known feature names (using `--names`). If you pass
`--sars2`, each name is printed followed by a colon and a (possibly empty)
list of aliases:

```sh
$ describe-feature.py --sars2 --names
2'-O-ribose methyltransferase: 2
3'-to-5' exonuclease: exon, exonuclease, nsp14
3'UTR: 3utr
3C-like proteinase: 3clpro, mpro, nsp5
5'UTR: 5utr
endoRNAse: endornase, nsp15
envelope protein: e, envelope, orf4
helicase: nsp13
leader protein: leader, nsp1
membrane glycoprotein: m, membrane, orf5
nsp2:
nsp3:
nsp4:
nsp6:
nsp7:
nsp8:
nsp9:
nsp10:
nsp11:
nucleocapsid phosphoprotein: n, orf9
ORF1ab polyprotein: orf1ab
ORF1a polyprotein: orf1a
ORF3a protein: orf3a
ORF6 protein: orf6
ORF7a protein: orf7a
ORF7b: orf7b
ORF8 protein: orf8
ORF10 protein: orf10
RNA-dependent RNA polymerase: nsp12, rdrp
stem loop 1: sl1
stem loop 2: sl2
stem loop 3: sl3
stem loop 4: sl4
stem loop 5: sl5
surface glycoprotein: s, spike
```

## describe-genome.py

`describe-genome.py` has many uses. It can extract multiple features from
multiple given genomes, as amino acids or nucleotides (or both). It will
print to standard output by default, but if you use the `--outDir` option
to provide a directory, individual output files with (hopefully)
self-explanatory names will be created in that directory. The directory
will be created for you if it doesn't exist.

A small example is perhaps best. Here we pull out the spike nucleotide and
amino acid sequence from B.1.1.7 and also ask for a summary of the amino
acid differences:

```sh
$ describe-genome.py --genome data/EPI_ISL_601443.fasta --outDir /tmp/out \
    --feature spike --printNtSequence --printAaSequence --printAaMatch
Examined 1 genome.

$ ls -l /tmp/out
total 24
-rw-rw-r--  1 terry  wheel   748 Apr 11 15:25 EPI_ISL_601443-spike-aa-match.txt
-rw-rw-r--  1 terry  wheel  1344 Apr 11 15:25 EPI_ISL_601443-spike-aa-sequence.fasta
-rw-rw-r--  1 terry  wheel  3886 Apr 11 15:25 EPI_ISL_601443-spike-nt-sequence.fasta
```

If the file given by `--genome` had contained more than one SARS-CoV-2
genome, the spike would have been extracted for all of them. Similarly, if
`--feature` had been repeated, multiple features would have been extracted
and compared. If no `--feature` is given, you'll get them all.

```sh
$ cat /tmp/out/EPI_ISL_601443-spike-aa-match.txt
Feature: spike amino acid match
  Matches: 1264/1274 (99.22%)
  Mismatches: 10/1274 (0.78%)
    Not involving gaps (i.e., conflicts): 7/1274 (0.55%)
    Involving a gap in one sequence: 3/1274 (0.24%)
    Involving a gap in both sequences: 0
    Id: NC_045512.2 (surface glycoprotein)
      Length: 1274
      Gaps: 0
    Id: EPI_ISL_601443 hCoV-19/England/MILK-9E05B3/2020 (surface glycoprotein)
      Length: 1274
      Gaps: 3/1274 (0.24%)
      Gap locations (1-based): 69, 70, 144
    Differences: site, aa1, aa2, ref nt codon start
        69 H - 21767
        70 V - 21770
       144 Y - 21992
       501 N Y 23063
       570 A D 23270
       614 D G 23402
       681 P H 23603
       716 T I 23708
       982 S A 24506
      1118 D H 24914
```

You can also ask for the changes in a variant to be checked.
      
```sh
$ describe-genome.py --sars2 --genome data/EPI_ISL_601443.fasta --checkVariant VOC_20201201_UK
EPI_ISL_601443 hCoV-19/England/MILK-9E05B3/2020
Variant summary:
  UK variant of concern (VOC) 202012/01:
  20 checks, 20 passed.
    orf1ab aa: PASS: 1001I, 1708D, 2230T, 3675-, 3676-, 3677-
    spike aa: PASS: 1118H, 144-, 501Y, 570D, 681H, 69-, 70-, 716I, 982A
    orf8 aa: PASS: 52I, 73C, Q27*
    n aa: PASS: 235F, 3L
```

There are some known variants, and you can also provide your own via a JSON
file and the `--variantFile` option. E.g.,:

```
VARIANTS = {
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
                'aa': 'Q27* 52I 73C',
            },
            'n': {
                'aa': '3L 235F',
            }
        }
    }
}
```

## describe-site.py

Will print information about a given location (a "site") in the genome,
showing you what's in the reference and in the genome you (optionally)
pass.

In the simplest case, just give a 1-based site and you'll see what's in the
reference (with 0-based Python offsets):

```sh
$ describe-site.py --sars2 --site 26000
{
    "alignmentOffset": 25999,
    "featureName": "ORF3a protein",
    "featureNames": [
        "ORF3a protein"
    ],
    "reference": {
        "aa": "L",
        "aaOffset": 202,
        "codon": "TTA",
        "frame": 1,
        "id": "NC_045512.2",
        "ntOffset": 25999
    }
}
```

At the moment, the default output is a JSON object, though this may someday
change to a more verbose/readable form. Use `--json` to make sure you
always get JSON.

You can also specify the site relative to a feature:

```sh
$ describe-site.py --sars2 --site 1501 --feature spike --relative
{
    "alignmentOffset": 23062,
    "featureName": "surface glycoprotein",
    "featureNames": [
        "surface glycoprotein"
    ],
    "reference": {
        "aa": "N",
        "aaOffset": 500,
        "codon": "AAT",
        "frame": 0,
        "id": "NC_045512.2",
        "ntOffset": 23062
    }
}
```

Or pass an amino acid site number (via `--aa`):

```sh
$ describe-site.py --sars2 --site 501 --feature spike --relative --aa
{
    "alignmentOffset": 23062,
    "featureName": "surface glycoprotein",
    "featureNames": [
        "surface glycoprotein"
    ],
    "reference": {
        "aa": "N",
        "aaOffset": 500,
        "codon": "AAT",
        "frame": 0,
        "id": "NC_045512.2",
        "ntOffset": 23062
    }
}
```

Of course it's more fun if you also provide a genome to compare the reference to.
Here's the `N501Y` change in B.1.1.7 (Alpha):

```sh
$ describe-site.py --sars2 --site 501 --relative --genome EPI_ISL_601443.fasta --feature spike --aa
{
    "alignmentOffset": 23062,
    "featureName": "surface glycoprotein",
    "featureNames": [
        "surface glycoprotein"
    ],
    "genome": {
        "aa": "Y",
        "aaOffset": 497,
        "codon": "TAT",
        "frame": 0,
        "id": "EPI_ISL_601443 hCoV-19/England/MILK-9E05B3/2020",
        "ntOffset": 22990
    },
    "reference": {
        "aa": "N",
        "aaOffset": 500,
        "codon": "AAT",
        "frame": 0,
        "id": "NC_045512.2",
        "ntOffset": 23062
    }
}
```

Other options include `--genomeAaOnly` to just print the amino acid at a
location in the genome, `--includeFeature` to also receive information
about the feature at the site, and `--minReferenceCoverage` to exclude
low-coverage genomes or features from the results.

## annotate-genome.py

Produces JSON output with annotation information for a genome. Here's an
example of the start of output for a Monkeypox genome (note that the
nucleotide sequences have been truncated to reduce output width):

```sh
$ annotate-genome.py --aligner edlib --reference ON676708.1.gb --genome ChVir28389.fasta
{
    "features": {
        "A15.5L": {
            "genome": {
                "sequence": "ATGATAAGTAATTACGAGCCGTTGCTGCTGTTAGTTATAAC...",
                "start": 121823,
                "stop": 121985,
                "translation": "MISNYEPLLLLVITCCVLLFNFTISSKTKIDIIFAVQTIVFIWFIFHFVYSAI*"
            },
            "reference": {
                "forward": false,
                "name": "A15.5L",
                "note": "Non-essential IMV membrane protein (Cop-A14.5L)",
                "product": "A15.5L",
                "sequence": "TTAAATCGCCGAATAAACAAAGTGGAATATAAACCATATAA...",
                "start": 121759,
                "stop": 121921,
                "translation": "MISNYEPLLLLVITCCVLLFNFTISSKTKIDIIFAVQTIVFIWFIFHFVYSAI*"
            }
        },
        "A32.5L": {
            "genome": {
                "sequence": "ATGTTAAAAATGTCAGCTGCCGACTTTTTGGAACGTTTGAT...",
                "start": 139696,
                "stop": 139825,
                "translation": "MLKMSAADFLERLIKAGIYIYVLRTKYVITALLVKNYPIKDE*"
            },
            "reference": {
                "forward": false,
                "name": "A32.5L",
                "note": "Viral membrane assembly proteins (VMAP) (Cop-A 30.5L)",
                "product": "A32.5L",
                "sequence": "TTATTCGTCTTTTATGGGATAGTTTTTAACTAGTAAAGCTGTAA...",
                "start": 139635,
                "stop": 139764,
                "translation": "MLKMSAADFLERLIKAGIYIYVLRTKYVITALLVKNYPIKDE*"
            }
        },
}
```

# Python API

There are two main Python classes provided by `gb2seq`: `Features` and
`Gb2Alignment`.

## Features

The `Features` class provides methods for accessing information about
genome features obtained from
[a GenBank flat file](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html).

You've likely run across these records before. E.g., here's the SARS-CoV-2
Wuhan reference,
[NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). To download
the GenBank file for that reference, click on "Send to" and select
"Complete Record", "File", and "GenBank" format. Then click "Create File".

You can pass the path to a GenBank file to the `Features` class. You can
also just pass an accession number, and the file will be downloaded for
you. If you pass `referenceSpecification` as `None` and `sars2` as `True`,
the Wuhan reference (version
[NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)) will be
used.

You can use a `Features` instance like a dictionary:

```py
from pprint import pprint as pp
from gb2seq.features import Features

>>> f = Features("NC_045512.2.gb")
>>> pp(f['e'])
{'name': 'envelope protein',
 'note': 'ORF4; structural protein; E protein',
 'product': 'envelope protein',
 'sequence': 'ATGTACTCATTCGTTTCGGAAGAGACAGGTACGTTAATAGTTAATAGCGTACTTCTTTTTCTTGCTTTCGTGGTATT...',
 'start': 26244,
 'stop': 26472,
 'translation': 'MYSFVSEETGTLIVNSVLLFLAFVVFLLVTLAILTALRLCAYCCNIVNVSLVKPSFYVYSRVKNLNSSRVPDLLV*'}

# You can use abbreviated names, and get the canonical names:
>>> f.canonicalName('s')
'surface glycoprotein'

# Get the list of aliases for a name:
>>> f.aliases('s')
{'spike', 'surface glycoprotein', 's'}

# Given an offset relative to a feature, get the offset in the genome:
>>> f.referenceOffset('spike', 1503)
23065

# Same thing, but with an amino acid offset.
>>> f.referenceOffset('spike', 501, aa=True)
23065

# What features are present at an offset?
>>> f.getFeatureNames(13450)
{'ORF1ab polyprotein', 'nsp11', 'ORF1a polyprotein', 'RNA-dependent RNA polymerase'}

# What features are present at an offset, including those that are not translated?
>>> f.getFeatureNames(21000, includeUntranslated=True)
{'ORF1ab polyprotein', "2'-O-ribose methyltransferase"}
```

## Gb2Alignment

The `Gb2Alignment` class can be used to extract and compare features from
the reference and the given genome sequence.

You pass it a `Read` instance from the
[dark-matter](https://github.com/acorg/dark-matter) module (which is
installed for you when you install `gb2seq`).

```py
from gb2seq.alignment import Gb2Alignment
from dark.reads import Read

alignment = Gb2Alignment(Read('id', 'AGCT...'))
```

These can also be read from a FASTA file:

```py
from gb2seq.alignment import Gb2Alignment
from dark.fasta import FastaReads

for read in FastaReads('sequences.fasta'):
    alignment = Gb2Alignment(read)
```

Once you have a `Gb2Alignment` instance, you can ask it for the aligned
sequences or features.

Below I'll use the GISAID `EPI_ISL_601443` (B.1.1.7, or Alpha, variant)
sequence, which you can find in
[data/EPI_ISL_601443.fasta](data/EPI_ISL_601443.fasta) in this repo:

```py
>>> from pathlib import Path
>>> from pprint import pprint as pp
>>> from gb2seq.alignment import Gb2Alignment
>>> from dark.fasta import FastaReads

>>> alpha = list(FastaReads(Path('data/EPI_ISL_601443.fasta')))[0]
>>> len(alpha)
29764
>>> alpha.id
'EPI_ISL_601443 hCoV-19/England/MILK-9E05B3/2020'
>>> alpha.sequence[:50]
'AGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCT'

>>> alignment = Gb2Alignment(alpha)
```

You'll find the aligned reference and genome in `alignment.referenceAligned`
and `alignment.genomeAligned`, both of which are  `Read` instances:

```py
>>> len(alignment.referenceAligned)
29903
>>> len(alignment.genomeAligned)
29903

# Get the nucleotide sequence for the spike protein for the reference and genome.
>>> referenceSpikeNt, genomeSpikeNt = alignment.ntSequences('spike')
>>> len(referenceSpikeNt)
3822

# Get the amino acid sequence for the spike protein for the reference and genome.
>>> referenceSpikeAa, genomeSpikeAa = alignment.aaSequences('spike')
>>> len(referenceSpikeAa)
1274

# There were three deletions in the alpha spike, so it only covers 3813 of the
# 3822 bases in the reference spike.
>>> alignment.coverage('s')
(3813, 3822)

# Get information about what's at an offset (0-based). This is what the describe-site.py
# utility does (though it takes a 1-based site).
#
# Here is the Alpha N501Y change:
>>> pp(alignment.offsetInfo(500, relativeToFeature=True, aa=True, featureName='s'))
{'featureName': 'surface glycoprotein',
 'featureNames': {'surface glycoprotein'},
 'alignmentOffset': 23062,
 'reference': {'aa': 'N',
  'codon': 'AAT',
  'frame': 0,
  'id': 'NC_045512.2',
  'aaOffset': 500,
  'ntOffset': 23062},
 'genome': {'aa': 'Y',
  'codon': 'TAT',
  'frame': 0,
  'id': 'EPI_ISL_601443 hCoV-19/England/MILK-9E05B3/2020',
  'aaOffset': 497,
  'ntOffset': 22990}}
```

You can check a feature for an expected change:

``` py
>>> alignment.checkFeature('spike', 'N501Y', aa=True)
(1, 0, {'N501Y': (True, 'N', True, 'Y')})
```

The return value gives the number of tests done (1), the number that failed
(0), and a `dict` with information about each test. The `dict` values are
4-tuples, indicating whether what was in the reference was as expected
(`True`), the value in the reference (`N`), whether what was in the genome was
as expected (`True`), the value in the genome (`Y`).

You can test multiple things:

``` py
# Note that we use 501Y in the following, not N501Y, since we might just
# want to check that something is in the genome without knowing or caring
# about what's in the reference.
>>> pp(alignment.checkFeature('spike', '501Y 69- 70-', aa=True))
(3,
 0,
 {'501Y': (True, 'N', True, 'Y'),
  '69-': (True, 'H', True, '-'),
  '70-': (True, 'V', True, '-')})
```


There is also a convenience `Checker` class that can check whether logical
combinations of amino acid and nucleotide changes are satisfied for a genome.
Continuing from the above:

``` py
>>> from gb2seq.checker import AAChecker, NTChecker

# Make a Boolean checker function to test whether a genome has the N501Y
# and A570D spike changes seen in Alpha.
>>> checker = AAChecker('spike', 'N501Y') & AAChecker('spike', 'A570D')
>>> checker(alignment)
True

# Check for some nucleotide changes in the nucleocapsid and some amino
# acid changes in the spike.
checker = (NTChecker('N', 'G7C A8T T9A G608A G609A G610C C704T') &
           AAChecker('S', 'N501Y H69- V70- Y144-'))

# You can also use `|` to check an OR condition.
>>> checker = AAChecker('spike', 'E484K') | AAChecker('spike', 'E484Q')
>>> checker(alignment)
False
```

Note that because the `Checker` class takes a string argument that likely
originates from a human-readable source, e.g., "N501Y", the changes
specified for the checker are in terms of 1-based sites.

Similiar to the `--checkVariant` argument to `describe-genome.py` (above),
you can also check a pre-defined variant:

```py
>>> pp(alignment.checkVariant('VOC_20201201_UK'))
(20,
 0,
 {'n': {'aa': {'235F': (True, 'S', True, 'F'), '3L': (True, 'D', True, 'L')}},
  'orf1ab': {'aa': {'1001I': (True, 'T', True, 'I'),
                    '1708D': (True, 'A', True, 'D'),
                    '2230T': (True, 'I', True, 'T'),
                    '3675-': (True, 'S', True, '-'),
                    '3676-': (True, 'G', True, '-'),
                    '3677-': (True, 'F', True, '-')}},
  'orf8': {'aa': {'52I': (True, 'R', True, 'I'),
                  '73C': (True, 'Y', True, 'C'),
                  'Q27*': (True, 'Q', True, '*')}},
  'spike': {'aa': {'1118H': (True, 'D', True, 'H'),
                   '144-': (True, 'Y', True, '-'),
                   '501Y': (True, 'N', True, 'Y'),
                   '570D': (True, 'A', True, 'D'),
                   '681H': (True, 'P', True, 'H'),
                   '69-': (True, 'H', True, '-'),
                   '70-': (True, 'V', True, '-'),
                   '716I': (True, 'T', True, 'I'),
                   '982A': (True, 'S', True, 'A')}}})
```

with the results summarizing the number of tests done, the number that
failed, and then the features checked, with `aa` and `nt` dictionaries for
the amino acid and nucleotide checks.

You can pass your own variant dictionary specifying what you want checked.

## To learn more

See the `Features` and `Gb2Alignment` classes in
[gb2seq/feature.py](gb2seq/feature.py) and
[gb2seq/alignment.py](gb2seq/alignment.py). Also, the tests (e.g., in
[test/test_feature.py](test/test_feature.py),
[test/test_alignment.py](test/test_alignment.py),
[test/test_checker.py](test/test_checker.py)
[test/test_variants.py](test/test_variants.py)) show you example uses of
these classes and their methods.  You can also look to see how the three
utility scripts described above (which can all be found in the [bin
directory](bin)) call the library functions and use the results.

# Developing

Run the tests via 

```sh
$ pytest
```

## Notes

* The tests will run `MAFFT`, so you need that installed and in
    your shell's path.
* The tests do not all pass if you set `DEFAULT_ALIGNER` to be "edlib" in
    `gb2seq/alignment.py`. This is because `MAFFT` and `edlib` produce
    slightly different results and some tests of the SARS-CoV-2 spike
    translation tests fail as a result. I would like to post-process the
    `edlib` output to match that of `MAFFT` (this is in cases where there
    are multiple equally-good alignments).
