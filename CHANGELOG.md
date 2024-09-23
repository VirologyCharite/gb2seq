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
