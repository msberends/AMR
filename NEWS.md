# 0.3.0.90xx (latest development version)

#### New
* Functions `count_R`, `count_IR`, `count_I`, `count_SI` and `count_S` to selectively count resistant or susceptible isolates
  * Extra function `count_df` (which works like `portion_df`) to get all counts of S, I and R of a data set with antibiotic columns, with support for grouped variables
* Function `is.rsi.eligible` to check for columns that have valid antimicrobial results, but do not have the `rsi` class yet. Transform the columns of your raw data with: `data %>% mutate_if(is.rsi.eligible, as.rsi)`
* Functions `as.mo` and `is.mo` as replacements for `as.bactid` and `is.bactid` (since the `microoganisms` data set not only contains bacteria). These last two functions are deprecated and will be removed in a future release.
* Renamed all previous references to `bactid` to `mo`, like:
  * Column names inputs of `EUCAST_rules`, `first_isolate` and `key_antibiotics`
  * Column names of datasets `microorganisms` and `septic_patients`
  * All old syntaxes will still work with this version, but will throw warnings
* Functions `as.atc` and `is.atc` to transform/look up antibiotic ATC codes as defined by the WHO. The existing function `guess_atc` is now an alias of `as.atc`.
* Aliases for existing function `mo_property`: `mo_family`, `mo_genus`, `mo_species`, `mo_subspecies`, `mo_fullname`, `mo_shortname`, `mo_aerobic`, `mo_type` and `mo_gramstain`. The last two functions have a `language` parameter, with support for Spanish, German and Dutch:
  ```r
  mo_gramstain("E. coli")
  # [1] "Negative rods"
  mo_gramstain("E. coli", language = "de") # "de" = Deutsch / German
  # [1] "Negative Staebchen"
  mo_gramstain("E. coli", language = "es") # "es" = Español / Spanish
  # [1] "Bacilos negativos"
  ```
* Function `ab_property` and its aliases: `ab_official`, `ab_tradenames`, `ab_certe`, `ab_umcg`, `ab_official_nl` and `ab_trivial_nl`
* Introduction to AMR as a vignette

#### Changed
* Added 231 microorganisms to the `microorganisms` data set and removed the few viruses it contained, now *n* = 2,669 (2,230 bacteria, 285 fungi/yeasts, 153 parasites, 1 other)
* Added three antimicrobial agents to the `antibiotics` data set: Terbinafine (D01BA02), Rifaximin (A07AA11) and Isoconazole (D01AC05)
* Added 163 trade names to the `antibiotics` data set, it now contains 298 different trade names in total, e.g.:
  ```r
  ab_official("Bactroban")
  # [1] "Mupirocin"
  ab_official(c("Bactroban", "Amoxil", "Zithromax", "Floxapen"))
  # [1] "Mupirocin" "Amoxicillin" "Azithromycin" "Flucloxacillin"
  ab_atc(c("Bactroban", "Amoxil", "Zithromax", "Floxapen"))
  # [1] "R01AX06" "J01CA04" "J01FA10" "J01CF05"
  ```
* Function `ratio` is now deprecated and will be removed in a future release, as it is not really the scope of this package
* Fix for `as.mic` for values ending in zeroes after a real number
* Tremendous speed improvement for `as.bactid` (now `as.mo`)
* Added parameters `minimum` and `as_percent` to `portion_df`
* Support for quasiquotation in the functions series `count_*` and `portions_*`, and `n_rsi`. This allows to check for more than 2 vectors or columns.
  ```r
  septic_patients %>% select(amox, cipr) %>% count_IR()
  # which is the same as:
  septic_patients %>% count_IR(amox, cipr)
  
  septic_patients %>% portion_S(amcl)
  septic_patients %>% portion_S(amcl, gent)
  septic_patients %>% portion_S(amcl, gent, pita)
  ```
* Edited `ggplot_rsi` and `geom_rsi` so they can cope with `count_df`. The new `fun` parameter has value `portion_df` at default, but can be set to `count_df`.
* Fix for `ggplot_rsi` when the `ggplot2` package was not loaded
* Added possibility to set any parameter to `geom_rsi` (and `ggplot_rsi`) so you can set your own preferences
* Fix for joins, where predefined suffices would not be honoured
* Support for types list and matrix for `freq`
  ```r
  my_matrix = with(septic_patients, matrix(c(age, sex), ncol = 2))
  freq(my_matrix)
  ```
  * Subsetting also possible for lists:
    ```r
    my_list = list(age = septic_patients$age, sex = septic_patients$sex)
    my_list %>% freq(age)
    my_list %>% freq(sex)
    ```

#### Other
* More unit tests to ensure better integrity of functions

# 0.3.0 (latest stable version)
**Published on CRAN: 2018-08-14**

#### New
* **BREAKING**: `rsi_df` was removed in favour of new functions `portion_R`, `portion_IR`, `portion_I`, `portion_SI` and `portion_S` to selectively calculate resistance or susceptibility. These functions are 20 to 30 times faster than the old `rsi` function. The old function still works, but is deprecated.
  * New function `portion_df` to get all portions of S, I and R of a data set with antibiotic columns, with support for grouped variables
* **BREAKING**: the methodology for determining first weighted isolates was changed. The antibiotics that are compared between isolates (call *key antibiotics*) to include more first isolates (afterwards called first *weighted* isolates) are now as follows:
  * Universal: amoxicillin, amoxicillin/clavlanic acid, cefuroxime, piperacillin/tazobactam, ciprofloxacin,  trimethoprim/sulfamethoxazole
  * Gram-positive: vancomycin, teicoplanin, tetracycline, erythromycin, oxacillin, rifampicin
  * Gram-negative: gentamicin, tobramycin, colistin, cefotaxime, ceftazidime, meropenem
* Support for `ggplot2`
  * New functions `geom_rsi`, `facet_rsi`, `scale_y_percent`, `scale_rsi_colours` and `theme_rsi`
  * New wrapper function `ggplot_rsi` to apply all above functions on a data set:
    * `septic_patients %>% select(tobr, gent) %>% ggplot_rsi` will show portions of S, I and R immediately in a pretty plot
    * Support for grouped variables, see `?ggplot_rsi`
* Determining bacterial ID:
  * New functions `as.bactid` and `is.bactid` to transform/ look up microbial ID's.
  * The existing function `guess_bactid` is now an alias of `as.bactid`
  * New Becker classification for *Staphylococcus* to categorise them into Coagulase Negative *Staphylococci* (CoNS) and Coagulase Positve *Staphylococci* (CoPS)
  * New Lancefield classification for *Streptococcus* to categorise them into Lancefield groups
* For convience, new descriptive statistical functions `kurtosis` and `skewness` that are lacking in base R - they are generic functions and have support for vectors, data.frames and matrices
* Function `g.test` to perform the Χ<sup>2</sup> distributed [*G*-test](https://en.wikipedia.org/wiki/G-test), which use is the same as `chisq.test`
* ~~Function `ratio` to transform a vector of values to a preset ratio~~
  * ~~For example: `ratio(c(10, 500, 10), ratio = "1:2:1")` would return `130, 260, 130`~~
* Support for Addins menu in RStudio to quickly insert `%in%` or `%like%` (and give them keyboard shortcuts), or to view the datasets that come with this package
* Function `p.symbol` to transform p values to their related symbols: `0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`
* Functions `clipboard_import` and `clipboard_export` as helper functions to quickly copy and paste from/to software like Excel and SPSS. These functions use the `clipr` package, but are a little altered to also support headless Linux servers (so you can use it in RStudio Server)
* New for frequency tables (function `freq`):
  * A vignette to explain its usage
  * Support for `rsi` (antimicrobial resistance) to use as input
  * Support for `table` to use as input: `freq(table(x, y))`
  * Support for existing functions `hist` and `plot` to use a frequency table as input: `hist(freq(df$age))`
  * Support for `as.vector`, `as.data.frame`, `as_tibble` and `format`
  * Support for quasiquotation: `freq(mydata, mycolumn)` is the same as `mydata %>% freq(mycolumn)`
  * Function `top_freq` function to return the top/below *n* items as vector
  * Header of frequency tables now also show Mean Absolute Deviaton (MAD) and Interquartile Range (IQR)
  * Possibility to globally set the default for the amount of items to print, with `options(max.print.freq = n)` where *n* is your preset value

#### Changed
* Improvements for forecasting with `resistance_predict` and added more examples
* More antibiotics added as parameters for EUCAST rules
* Updated version of the `septic_patients` data set to better reflect the reality
* Pretty printing for tibbles removed as it is not really the scope of this package
* Printing of `mic` and `rsi` classes now returns all values - use `freq` to check distributions
* Improved speed of key antibiotics comparison for determining first isolates
* Column names for the `key_antibiotics` function are now generic: 6 for broadspectrum ABs, 6 for Gram-positive specific and 6 for Gram-negative specific ABs
* Speed improvement for the `abname` function
* `%like%` now supports multiple patterns
* Frequency tables are now actual `data.frame`s with altered console printing to make it look like a frequency table. Because of this, the parameter `toConsole` is not longer needed.
* Fix for `freq` where the class of an item would be lost
* Small translational improvements to the `septic_patients` dataset and the column `bactid` now has the new class `"bactid"`
* Small improvements to the `microorganisms` dataset (especially for *Salmonella*) and the column `bactid` now has the new class `"bactid"`
* Combined MIC/RSI values will now be coerced by the `rsi` and `mic` functions:
  * `as.rsi("<=0.002; S")` will return `S`
  * `as.mic("<=0.002; S")` will return `<=0.002`
* Now possible to coerce MIC values with a space between operator and value, i.e. `as.mic("<= 0.002")` now works
* Classes `rsi` and `mic` do not add the attribute `package.version` anymore
* Added `"groups"` option for `atc_property(..., property)`. It will return a vector of the ATC hierarchy as defined by the [WHO](https://www.whocc.no/atc/structure_and_principles/). The new function `atc_groups` is a convenient wrapper around this.
* Build-in host check for `atc_property` as it requires the host set by `url` to be responsive
* Improved `first_isolate` algorithm to exclude isolates where bacteria ID or genus is unavailable
* Fix for warning *hybrid evaluation forced for row_number* ([`924b62`](https://github.com/tidyverse/dplyr/commit/924b62)) from the `dplyr` package v0.7.5 and above
* Support for empty values and for 1 or 2 columns as input for `guess_bactid` (now called `as.bactid`)
  * So `yourdata %>% select(genus, species) %>% as.bactid()` now also works
* Other small fixes

#### Other
* Added integration tests (check if everything works as expected) for all releases of R 3.1 and higher
  * Linux and macOS: https://travis-ci.org/msberends/AMR
  * Windows: https://ci.appveyor.com/project/msberends/amr
* Added thesis advisors to DESCRIPTION file

# 0.2.0
**Published on CRAN: 2018-05-03**

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
**Published on CRAN: 2018-03-14**

* `EUCAST_rules` applies for amoxicillin even if ampicillin is missing
* Edited column names to comply with GLIMS, the laboratory information system
* Added more valid MIC values
* Renamed 'Daily Defined Dose' to 'Defined Daily Dose'
* Added barplots for `rsi` and `mic` classes

# 0.1.0
**Published on CRAN: 2018-02-22**

* First submission to CRAN.
