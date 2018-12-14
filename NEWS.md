# 0.5.0.90xx (latest development version)

#### New
* Function `mo_failures` to review values that could not be coerced to a valid MO code, using `as.mo`. This latter function will now only show a maximum of 25 uncoerced values.
* Function `mo_renamed` to get a list of all returned values from `as.mo` that have had taxonomic renaming

#### Changed
* Improvements for `as.mo`:
  * Finds better results when input is in other languages
  * Better handling for subspecies
  * Better handling for *Salmonellae*
  * There will be looked for uncertain results at default - these results will be returned with an informative warning
  * Manual now contains more info about the algorithms
  * Progress bar will be shown when it takes more than 3 seconds to get results
* Function `first_isolate`:
  * Will now use a column named like "patid" for the patient ID (parameter `col_patientid`), when this parameter was left blank
  * Will now use a column named like "key(...)ab" or "key(...)antibiotics" for the key antibiotics (parameter `col_keyantibiotics`), when this parameter was left blank
* A note to the manual pages of the `portion` functions, that low counts can infuence the outcome and that the `portion` functions may camouflage this, since they only return the portion (albeit being dependent on the `minimum` parameter)
* Function `mo_taxonomy` now contains the kingdom too
* Function `first_isolate` will now use a column named like "patid" for the patient ID, when this parameter was left blank
* Reduce false positives for `is.rsi.eligible`
* Summaries of class `mo` will now return the top 3 and the unique count, e.g. using `summary(mo)`
* Small text updates to summaries of class `rsi` and `mic`
* Frequency tables (`freq` function):
  * Added header info for class `mo` to show unique count of families, genera and species
  * Now honours the `decimal.mark` setting, which just like `format` defaults to `getOption("OutDec")`
  * The new `big.mark` parameter will at default be `","` when `decimal.mark = "."` and `"."` otherwise
  * Fix for header text where all observations are `NA`


# 0.5.0 (latest stable release)
**Published on CRAN: 2018-11-30**

#### New
* Repository moved to GitLab: https://gitlab.com/msberends/AMR
* Function `count_all` to get all available isolates (that like all `portion_*` and `count_*` functions also supports `summarise` and `group_by`), the old `n_rsi` is now an alias of `count_all`
* Function `get_locale` to determine language for language-dependent output for some `mo_*` functions. This is now the default value for their `language` parameter, by which the system language will be used at default.
* Data sets `microorganismsDT`, `microorganisms.prevDT`, `microorganisms.unprevDT` and `microorganisms.oldDT` to improve the speed of `as.mo`. They are for reference only, since they are primarily for internal use of `as.mo`.
* Function `read.4D` to read from the 4D database of the MMB department of the UMCG
* Functions `mo_authors` and `mo_year` to get specific values about the scientific reference of a taxonomic entry

#### Changed
* Functions `MDRO`, `BRMO`, `MRGN` and `EUCAST_exceptional_phenotypes` were renamed to `mdro`, `brmo`, `mrgn` and `eucast_exceptional_phenotypes`
* `EUCAST_rules` was renamed to `eucast_rules`, the old function still exists as a deprecated function
* Big changes to the `eucast_rules` function:
  * Now also applies rules from the EUCAST 'Breakpoint tables for bacteria', version 8.1, 2018, http://www.eucast.org/clinical_breakpoints/ (see Source of the function)
  * New parameter `rules` to specify which rules should be applied (expert rules, breakpoints, others or all)
  * New parameter `verbose` which can be set to `TRUE` to get very specific messages about which columns and rows were affected
  * Better error handling when rules cannot be applied (i.e. new values could not be inserted)
  * The number of affected values will now only be measured once per row/column combination
  * Data set `septic_patients` now reflects these changes
  * Added parameter `pipe` for piperacillin (J01CA12), also to the `mdro` function
  * Small fixes to EUCAST clinical breakpoint rules
* Added column `kingdom` to the microorganisms data set, and function `mo_kingdom` to look up values
* Tremendous speed improvement for `as.mo` (and subsequently all `mo_*` functions), as empty values wil be ignored *a priori*
* Fewer than 3 characters as input for `as.mo` will return NA
* Function `as.mo` (and all `mo_*` wrappers) now supports genus abbreviations with "species" attached
  ```r
  as.mo("E. species")        # B_ESCHR
  mo_fullname("E. spp.")     # "Escherichia species"
  as.mo("S. spp")            # B_STPHY
  mo_fullname("S. species")  # "Staphylococcus species"
  ```
* Added parameter `combine_IR` (TRUE/FALSE) to functions `portion_df` and `count_df`, to indicate that all values of I and R must be merged into one, so the output only consists of S vs. IR (susceptible vs. non-susceptible)
* Fix for `portion_*(..., as_percent = TRUE)` when minimal number of isolates would not be met
* Added parameter `also_single_tested` for `portion_*` and `count_*` functions to also include cases where not all antibiotics were tested but at least one of the tested antibiotics includes the target antimicribial interpretation, see `?portion`
* Using `portion_*` functions now throws a warning when total available isolate is below parameter `minimum`
* Functions `as.mo`, `as.rsi`, `as.mic`, `as.atc` and `freq` will not set package name as attribute anymore
* Frequency tables - `freq()`:
  * Support for grouping variables, test with:
    ```r
    septic_patients %>% 
      group_by(hospital_id) %>% 
      freq(gender)
    ```
  * Support for (un)selecting columns:
    ```r
    septic_patients %>% 
      freq(hospital_id) %>% 
      select(-count, -cum_count) # only get item, percent, cum_percent
    ```
  * Check for `hms::is.hms`
  * Now prints in markdown at default in non-interactive sessions
  * No longer adds the factor level column and sorts factors on count again
  * Support for class `difftime`
  * New parameter `na`, to choose which character to print for empty values
  * New parameter `header` to turn the header info off (default when `markdown = TRUE`)
  * New parameter `title` to manually setbthe title of the frequency table
* `first_isolate` now tries to find columns to use as input when parameters are left blank
* Improvements for MDRO algorithm (function `mdro`)
* Data set `septic_patients` is now a `data.frame`, not a tibble anymore
* Removed diacritics from all authors (columns `microorganisms$ref` and `microorganisms.old$ref`) to comply with CRAN policy to only allow ASCII characters
* Fix for `mo_property` not working properly
* Fix for `eucast_rules` where some Streptococci would become ceftazidime R in EUCAST rule 4.5
* Support for named vectors of class `mo`, useful for `top_freq()`
* `ggplot_rsi` and `scale_y_percent` have `breaks` parameter
* AI improvements for `as.mo`:
  * `"CRS"` -> *Stenotrophomonas maltophilia*
  * `"CRSM"` -> *Stenotrophomonas maltophilia*
  * `"MSSA"` -> *Staphylococcus aureus*
  * `"MSSE"` -> *Staphylococcus epidermidis*
* Fix for `join` functions
* Speed improvement for `is.rsi.eligible`, now 15-20 times faster
* In `g.test`, when `sum(x)` is below 1000 or any of the expected values is below 5, Fisher's Exact Test will be suggested
* `ab_name` will try to fall back on `as.atc` when no results are found
* Removed the addin to view data sets
* Percentages will now will rounded more logically (e.g. in `freq` function)

#### Other
* New dependency on package `crayon`, to support formatted text in the console
* Dependency `tidyr` is now mandatory (went to `Import` field) since `portion_df` and `count_df` rely on it
* Updated vignettes to comply with README


# 0.4.0
**Published on CRAN: 2018-10-01**

#### New
* The data set `microorganisms` now contains **all microbial taxonomic data from ITIS** (kingdoms Bacteria, Fungi and Protozoa), the Integrated Taxonomy Information System, available via https://itis.gov. The data set now contains more than 18,000 microorganisms with all known bacteria, fungi and protozoa according ITIS with genus, species, subspecies, family, order, class, phylum and subkingdom. The new data set `microorganisms.old` contains all previously known taxonomic names from those kingdoms.
* New functions based on the existing function `mo_property`:
  * Taxonomic names: `mo_phylum`, `mo_class`, `mo_order`, `mo_family`, `mo_genus`, `mo_species`, `mo_subspecies`
  * Semantic names: `mo_fullname`, `mo_shortname`
  * Microbial properties: `mo_type`, `mo_gramstain`
  * Author and year: `mo_ref`
  
  They also come with support for German, Dutch, French, Italian, Spanish and Portuguese:
  ```r
  mo_gramstain("E. coli")
  # [1] "Gram negative"
  mo_gramstain("E. coli", language = "de") # German
  # [1] "Gramnegativ"
  mo_gramstain("E. coli", language = "es") # Spanish
  # [1] "Gram negativo"
  mo_fullname("S. group A", language = "pt") # Portuguese
  # [1] "Streptococcus grupo A"
  ```
  
  Furthermore, former taxonomic names will give a note about the current taxonomic name:
  ```r
  mo_gramstain("Esc blattae")
  # Note: 'Escherichia blattae' (Burgess et al., 1973) was renamed 'Shimwellia blattae' (Priest and Barker, 2010)
  # [1] "Gram negative"
  ```
* Functions `count_R`, `count_IR`, `count_I`, `count_SI` and `count_S` to selectively count resistant or susceptible isolates
  * Extra function `count_df` (which works like `portion_df`) to get all counts of S, I and R of a data set with antibiotic columns, with support for grouped variables
* Function `is.rsi.eligible` to check for columns that have valid antimicrobial results, but do not have the `rsi` class yet. Transform the columns of your raw data with: `data %>% mutate_if(is.rsi.eligible, as.rsi)`
* Functions `as.mo` and `is.mo` as replacements for `as.bactid` and `is.bactid` (since the `microoganisms` data set not only contains bacteria). These last two functions are deprecated and will be removed in a future release. The `as.mo` function determines microbial IDs using Artificial Intelligence (AI):
  ```r
  as.mo("E. coli")
  # [1] B_ESCHR_COL
  as.mo("MRSA")
  # [1] B_STPHY_AUR
  as.mo("S group A")
  # [1] B_STRPTC_GRA
  ```
  And with great speed too - on a quite regular Linux server from 2007 it takes us less than 0.02 seconds to transform 25,000 items:
  ```r
  thousands_of_E_colis <- rep("E. coli", 25000)
  microbenchmark::microbenchmark(as.mo(thousands_of_E_colis), unit = "s")
  # Unit: seconds
  #         min       median         max  neval
  #  0.01817717  0.01843957  0.03878077    100
  ```
* Added parameter `reference_df` for `as.mo`, so users can supply their own microbial IDs, name or codes as a reference table
* Renamed all previous references to `bactid` to `mo`, like:
  * Column names inputs of `EUCAST_rules`, `first_isolate` and `key_antibiotics`
  * Column names of datasets `microorganisms` and `septic_patients`
  * All old syntaxes will still work with this version, but will throw warnings
* Function `labels_rsi_count` to print datalabels on a RSI `ggplot2` model
* Functions `as.atc` and `is.atc` to transform/look up antibiotic ATC codes as defined by the WHO. The existing function `guess_atc` is now an alias of `as.atc`.

* Function `ab_property` and its aliases: `ab_name`, `ab_tradenames`, `ab_certe`, `ab_umcg` and `ab_trivial_nl`
* Introduction to AMR as a vignette
* Removed clipboard functions as it violated the CRAN policy
* Renamed `septic_patients$sex` to `septic_patients$gender`

#### Changed
* Added three antimicrobial agents to the `antibiotics` data set: Terbinafine (D01BA02), Rifaximin (A07AA11) and Isoconazole (D01AC05)
* Added 163 trade names to the `antibiotics` data set, it now contains 298 different trade names in total, e.g.:
  ```r
  ab_official("Bactroban")
  # [1] "Mupirocin"
  ab_name(c("Bactroban", "Amoxil", "Zithromax", "Floxapen"))
  # [1] "Mupirocin" "Amoxicillin" "Azithromycin" "Flucloxacillin"
  ab_atc(c("Bactroban", "Amoxil", "Zithromax", "Floxapen"))
  # [1] "R01AX06" "J01CA04" "J01FA10" "J01CF05"
  ```
* For `first_isolate`, rows will be ignored when there's no species available
* Function `ratio` is now deprecated and will be removed in a future release, as it is not really the scope of this package
* Fix for `as.mic` for values ending in zeroes after a real number
* Small fix where *B. fragilis* would not be found in the `microorganisms.umcg` data set
* Added `prevalence` column to the `microorganisms` data set
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
* Added datalabels function `labels_rsi_count` to `ggplot_rsi`
* Added possibility to set any parameter to `geom_rsi` (and `ggplot_rsi`) so you can set your own preferences
* Fix for joins, where predefined suffices would not be honoured
* Added parameter `quote` to the `freq` function
* Added generic function `diff` for frequency tables
* Added longest en shortest character length in the frequency table (`freq`) header of class `character`
* Support for types (classes) list and matrix for `freq`
  ```r
  my_matrix = with(septic_patients, matrix(c(age, gender), ncol = 2))
  freq(my_matrix)
  ```
  For lists, subsetting is possible:
  ```r
  my_list = list(age = septic_patients$age, gender = septic_patients$gender)
  my_list %>% freq(age)
  my_list %>% freq(gender)
  ```

#### Other
* More unit tests to ensure better integrity of functions

# 0.3.0
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
