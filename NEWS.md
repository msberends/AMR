## 0.2.0
#### New
- Full support for Windows, Linux and macOS
- Function `guess_bactid` to determine the ID of a microorganism based on genus/species or known abbreviations like MRSA
- Functions `clipboard_import` and `clipboard_export` as helper functions to quickly copy and paste from/to software like Excel and SPSS
- Function `MDRO` to determine Multi Drug Resistant Organisms (MDRO) with support for country-specific guidelines. Suggest your own via https://github.com/msberends/AMR/issues/new. Functions `BRMO` and `MRGN` are wrappers for Dutch and German guidelines, respectively
- Function `freq` to create frequency tables, with additional info in a header
- New algorithm to determine weighted isolates, can now be `"points"` or `"keyantibiotics"`, see `?first_isolate`
- New print format for tibbles and data.tables

#### Changed
- Support for older R versions, only R 3.1.0 and later is needed
- Renamed dataset `ablist` to `antibiotics`
- Renamed dataset `bactlist` to `microorganisms`
- Added more microorganisms to `bactlist`
- Added analysis examples on help page of dataset `septic_patients`
- Added support for character vector in `join` functions
- Added warnings when a join results in more rows after than before the join
- Altered `%like%` to make it case insensitive
- For parameters of functions `first_isolate` and `EUCAST_rules` the column names are now case-insensitive
- Functions `as.rsi` and `as.mic` now add the package name and version as attribute

#### Other
- Expanded README.md with more examples
- Added ORC IDs of authors to DESCRIPTION file
- Added unit testing with the `testthat` package
- Added build tests for Linux and macOS using Travis CI (https://travis-ci.org/msberends/AMR)
- Added Line coverage checking using CodeCov (https://codecov.io/gh/msberends/AMR/tree/master/R)

## 0.1.1
- `EUCAST_rules` applies for amoxicillin even if ampicillin is missing
- Edited column names to comply with GLIMS, the laboratory information system
- Added more valid MIC values
- Renamed 'Daily Defined Dose' to 'Defined Daily Dose'
- Added barplots for `rsi` and `mic` classes

## 0.1.0
- First submission to CRAN.
