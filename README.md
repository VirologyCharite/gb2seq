## A library and utility scripts for working with SARS-CoV-2 sequences

sars2seq provides three general utility scripts for examining SARS-CoV-2
genome sequences and a Python library API to allow you to roll your own.

The library code knows how to extract genes at the nucleotide or amino acid
level, how to identify or check for expected substitutions, how to
translate between aa and nt offsets, and about all genome features.

### Alignment

`sars2seq` makes an alignment between each genome you give it and a
reference sequence
([NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) by
default). The default alignment algorithm is
[MAFFT](https://mafft.cbrc.jp/alignment/software/), which can be slow but
is reliable. You can also run scripts with `--aligner edlib` (or pass
`aligner='edlib' to library functions) to use [the Python
wrapper](https://pypi.org/project/edlib/) for the extremely fast
[edlib](https://github.com/Martinsos/edlib) library.

## Utility scripts

The three scripts described below all accept a `--help` option. Below is
some representative usage.

### describe-feature.py

The simplest script is `describe-feature.py`, which can be used to get
information about features in a SARS-CoV-2 reference genome (you can
specify a different reference by passing the name of a GenBank file with
`--gbFile`).

If you run `describe-feature.py` with no arguments, it will print summary
details of all features, with their (0-based) nucleotide offsets:

```sh
$ describe-feature.py
Features for NC_045512.2:
2'-O-ribose methyltransferase:
  start: 20658
  stop: 21552
  length: 894
  product: 2'-O-ribose methyltransferase
  sequence    (len   894 nt): TCTAGTCAAGCGTGGCAACCGGGTGTTGCTATGCCTAATCTTTACAAAATGCAAAGAATGCTATTAGAAAAGTGTGACCT...
3'-to-5' exonuclease:
  start: 18039
  stop: 19620
  length: 1581
  product: 3'-to-5' exonuclease
  sequence    (len  1581 nt): GCTGAAAATGTAACAGGACTCTTTAAAGATTGTAGTAAGGTAATCACTGGGTTACATCCTACACAGGCACCTACACACCT...
3'UTR:
  start: 29674
  stop: 29903
  length: 229
  sequence    (len   229 nt): CAATCTTTAATCAGTGTGTAACATTAGGGAGGACTTGAAAGAGCCACCACATTTTCACCGAGGCCACGCGGAGTACGATC...
# [Many additional output lines omitted here.]
```

You can also pass it a feature name:

``` sh
$ describe-feature.py --name spike
surface glycoprotein:
  start: 21562
  stop: 25384
  length: 3822
  product: surface glycoprotein
  sequence    (len  3822 nt): ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGC...
  translation (len  1274 aa): MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFD...
```

Or ask to see all known feature names. Each is printed followed by a colon
and a (possibly empty) list of aliases:

``` sh
$ describe-feature.py --names
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

### describe-genome.py

`describe-genome.py` has many uses. It can extract multiple features from
multiple give genomes, as amino acids or nucleotides (or both). It will
print to standard error by default, but if you use the `--outDir` option to
provide a directory, individual output files with (hopefully)
self-explanatory names will be created in that directory. The directory
will be created for you if it doesn't exist.

A small example is perhaps best. Here we pull out the spike nucleotide and
amino acid sequence from B.1.1.7 and also ask for a summary of the amino
acid differences:

``` sh
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

``` sh
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
      
``` sh
$ describe-genome.py --genome data/EPI_ISL_601443.fasta --checkVariant VOC_20201201_UK
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

``` json
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

### describe-site.py

Will print information about a given location (a "site") in the genome,
showing you what's in the reference and in the genome you (optionally)
pass.

In the simplest case, just give a 1-based site and you'll see what's in the
reference:

``` sh
$ describe-site.py --site 26000
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

At the moment, the output is a JSON object, with 0-based genome offsets
(suitable for working in Python).

You can also specify the site relative to a feature:

```sh
$ describe-site.py --site 1501 --feature spike --relative
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

Or pass an amino acid site number (via `-aa`):

```sh
$ describe-site.py --site 501 --feature spike --relative --aa
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

Of course it's more fun if you also provide a genome to compare the reference to:

```sh
$ describe-site.py --site 501 --relative --genome EPI_ISL_601443.fasta --feature spike --aa
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
low-coverage genomes from the results.

## Python API

To be written!

For now, see the `SARS2Genome` in [sars2seq/genome.py](sars2seq/genome.py)
and the tests (e.g., in [test/test_genome.py], [test/test_checker.py]
[test/test_variants.py]).  You can also look to see how the three utility
scripts above call the library functions and use the results.

## Developing

Run the tests via 

``` sh
$ make pytest`
```
