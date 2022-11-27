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
