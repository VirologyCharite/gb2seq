## 0.3.5 November 27, 2025

Fixed error in file opening mode (now needs to be binary due to changes in
dark-matter FASTA reader) in `bin/describe-genome.py`.

## 0.3.4 August 19, 2025

Added 2'-O-ribose methyltransferase to the list of translated SARS-CoV-2
features.

## 0.3.3 August 17, 2025

Switched to `uv` project organization with `pyproject.toml`.  Many small
(mainly cosmetic) changes to `bin/compare-genomes.py`.

## 0.3.2 August 15, 2025

Not released.

## 0.3.1 August 14, 2025

Added `bin/compare-genomes.py` script. This allows you to find differences
between two unannotated sequences via aligning both with an annotated
reference.

## 0.3.0 August 13, 2025

Improved the `Gb2Alignment.offsetInfo` method to be more careful about the nt
offset it returns in the genome (it was possible for it to return an index
that was longer than the genome). The 'aa', 'ntOffset' values may now come
back as `None`. See the long comment in that method for details on what it
returns.

## 0.2.27 August 12, 2025

Marked the SARS-CoV-2 furin cleavage site as a translated feature and an "fcs"
abbreviation for it.

## 0.2.26 August 12, 2025

Added CDS section for furin cleavage site with product name and translation
to SARS-2 genbank file. This doesn't seem to break anything (try running
`bin/describe-feature.py --sars2 --alsoInclude misc_feature`) but may be
violating some principle.

## 0.2.25 June 27, 2025

Added furin cleavage site miscellaneous feature to SARS-2 genbank file. Added
`--alsoInclude` to feature options.

## 0.2.24 January 25, 2025

Improve annotating an unknown genome so untranslated features aren't treated
as if they were translated. This code needs work!

## 0.2.23 September 22, 2024

Change Features class to no longer inherit from UserDict. Fix unbound local bug in Features.

## 0.2.22 August 27, 2024

Improved dealing with features that have multiple genome ranges to fix
https://github.com/VirologyCharite/gb2seq/issues/11

## 0.2.21 August 10, 2023
Added `--translated` option to `bin/describe-genome.py`.
This will extract all the features that are translated.

## 0.2.20 August 7, 2023

'function' may not be present on stem_loop features. Added
`--slashReplacement` and `--spaceReplacement` to `bin/describe-genome.py`
for more control over creating output filenames.

## 0.2.19 August 7, 2023

More robust treatment of feature names. After processing everything we want
to treat specially, just look for things that have a "product" and ignore
those that don't.

## 0.2.18 August 7, 2023

Added `ncRNA` to known feature types. Probably I should just ignore all
unknown feature types but at this point I still want to see them.

## 0.2.17 August 7, 2023

Change to use `dark.aaVars` to be compatible with dark-matter 4.0.68

## 0.2.4 January 16, 2023

Added `--genome` to `describe-feature.py`.

## 0.2.3 November 27, 2022

Added `nsp16` as a sars-2 abbreviation for `2'-O-ribose
methyltransferase`. Made `TranslationError` a subclass of `Gb2SeqError`.

## 0.2.2 August 4, 2022

Remove some SARS-2 defaults.

## 0.2.1 August 4, 2022

More small changes for the `gb2seq` name change.

## 0.2.0 August 4, 2022

Changed to `gb2seq`.

## 0.1.2 May 1, 2022

Switch to use [black](https://black.readthedocs.io/en/stable/index.html) formatting.

## 0.1.1 April 19, 2022

Renamed the `SARS2Genome` to `SARS2Alignment`. The old name is still
supported, so existing code will continue to work. A deprecation warning is
issued.

## 0.1.0 April 19, 2022

Two small backward-incompatible changes were introduced.

1. The `nt` argument to some checker functions was renamed to `aa` so you
   now need to pass the opposite Boolean value. This is to be consistent
   with other functions or scripts that already used the `aa` name.
1. The `featuresAt` method of the `Features` class was renamed to
   `getFeatureNames`.
