# 0.2.0.90xx (development version)
#### New
* Function `top_freq` function to get the top/below *n* items of frequency tables
* Vignette about frequency tables
* Possibility to globally set the default for the amount of items to print in frequency tables (`freq` function), with `options(max.print.freq = n)`
* Functions `clipboard_import` and `clipboard_export` as helper functions to quickly copy and paste from/to software like Excel and SPSS
* Function `g.test` to perform the Χ<sup>2</sup> distributed [*G*-test](https://en.wikipedia.org/wiki/G-test)
* Function `p.symbol` to transform p value to their related symbol: `0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`
* Function `vector2ratio` to transform a vector of values to a preset ratio. For example:
```r
vector2ratio(c(772, 1611, 737), ratio = "1:2:1")
# [1]  780 1560  780
```

#### Changed
* Frequency tables (function `freq`) now supports quasiquotation: `freq(mydata, mycolumn)`, or `mydata %>% freq(mycolumn)`
* Frequency tables are now actual `data.frame`s with altered console printing to make it look like a frequency table. Because of this, the parameter `toConsole` is not longer needed.
* Small translational improvements to the `septic_patients` dataset
* Combined MIC/RSI values will now be coerced by the `rsi` and `mic` functions:
  * `as.rsi("<=0.002; S")` will return `S`
  * `as.mic("<=0.002; S")` will return `<=0.002`
* Now possible to coerce MIC values with a space between operator and value, i.e. `as.mic("<= 0.002")` now works
* Added `"groups"` option for `atc_property(..., property)`. It will return a vector of the ATC hierarchy as defined by the [WHO](https://www.whocc.no/atc/structure_and_principles/). The new function `atc_groups` is a convenient wrapper around this.
* Build-in host check for `atc_property` as it requires the host set by `url` to be responsive
* Improved `first_isolate` algorithm to exclude isolates where bacteria ID or genus is unavailable
* Fix for warning *hybrid evaluation forced for row_number* ([`924b62`](https://github.com/tidyverse/dplyr/commit/924b62)) from the `dplyr` package v0.7.5 and above
* Support for 1 or 2 columns as input for `guess_bactid`
* Fix for printing tibbles where characters would be accidentally transformed to factors

# 0.2.0 (latest stable version)
#### New
* Full support for Windows, Linux and macOS
* Full support for old R versions, only R-3.0.0 (April 2013) or later is needed (needed packages may have other dependencies)
* Function `n_rsi` to count cases where antibiotic test results were available, to be used in conjunction with `dplyr::summarise`, see ?rsi
* Function `guess_bactid` to **determine the ID** of a microorganism based on genus/species or known abbreviations like MRSA
* Function `guess_atc` to **determine the ATC** of an antibiotic based on name, trade name, or known abbreviations
* Function `freq` to create **frequency tables**, with additional info in a header
* Function `MDRO` to **determine Multi Drug Resistant Organisms (MDRO)** with support for country-specific guidelines.
  * Suggest your own via [https://github.com/msberends/AMR/issues/new](https://github.com/msberends/AMR/issues/new?title=New%20guideline%20for%20MDRO&body=%3C--%20Please%20add%20your%20country%20code,%20guideline%20name,%20version%20and%20source%20below%20and%20remove%20this%20line--%3E)
  * [Exceptional resistances defined by EUCAST](http://www.eucast.org/expert_rules_and_intrinsic_resistance) are also supported instead of countries alone
  * Functions `BRMO` and `MRGN` are wrappers for Dutch and German guidelines, respectively
* New algorithm to determine weighted isolates, can now be `"points"` or `"keyantibiotics"`, see `?first_isolate`
* New print format for `tibble`s and `data.table`s

#### Changed
* Fixed `rsi` class for vectors that contain only invalid antimicrobial interpretations
* Renamed dataset `ablist` to `antibiotics`
* Renamed dataset `bactlist` to `microorganisms`
* Added common abbreviations and trade names to the `antibiotics` dataset
* Added more microorganisms to the `microorganisms` dataset
* Added analysis examples on help page of dataset `septic_patients`
* Added support for character vector in `join` functions
* Added warnings when a join results in more rows after than before the join
* Altered `%like%` to make it case insensitive
* For parameters of functions `first_isolate` and `EUCAST_rules` column names are now case-insensitive
* Functions `as.rsi` and `as.mic` now add the package name and version as attributes

#### Other
* Expanded `README.md` with more examples
* Added [ORCID](https://orcid.org) of authors to DESCRIPTION file
* Added unit testing with the `testthat` package
* Added build tests for Linux and macOS using Travis CI (https://travis-ci.org/msberends/AMR)
* Added line coverage checking using CodeCov (https://codecov.io/gh/msberends/AMR/tree/master/R)

# 0.1.1
* `EUCAST_rules` applies for amoxicillin even if ampicillin is missing
* Edited column names to comply with GLIMS, the laboratory information system
* Added more valid MIC values
* Renamed 'Daily Defined Dose' to 'Defined Daily Dose'
* Added barplots for `rsi` and `mic` classes

# 0.1.0
* First submission to CRAN.
