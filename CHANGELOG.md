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
