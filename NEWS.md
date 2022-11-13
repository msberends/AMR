# AMR 1.8.2.9049

This version will eventually become v2.0! We're happy to reach a new major milestone soon!

### Breaking
* Removed all species of the taxonomic kingdom Chromista from the package. This was done for multiple reasons:
  * CRAN allows packages to be around 5 MB maximum, some packages are exempted but this package is not one of them
  * Chromista are not relevant when it comes to antimicrobial resistance, thus lacking the primary scope of this package
  * Chromista are almost never clinically relevant, thus lacking the secondary scope of this package
* The `microorganisms` no longer relies on the Catalogue of Life, but now primarily on the List of Prokaryotic names with Standing in Nomenclature (LPSN) and is supplemented with the Global Biodiversity Information Facility (GBIF). The structure of this data set has changed to include separate LPSN and GBIF identifiers. Almost all previous MO codes were retained. It contains over 1,000 taxonomic names from 2022 already.
* **The `microorganisms.old` data set was removed**, and all previously accepted names are now included in the `microorganisms` data set. A new column `status` contains `"accepted"` for currently accepted names and `"synonym"` for taxonomic synonyms; currently invalid names. All previously accepted names now have a microorganisms ID and - if available - an LPSN, GBIF and SNOMED CT identifier.
* The MO matching score algorithm (`mo_matching_score()`) now counts deletions and substitutions as 2 instead of 1, which impacts the outcome of `as.mo()` and any `mo_*()` function
* **Argument `combine_IR` has been removed** from this package (affecting functions `count_df()`, `proportion_df()`, and `rsi_df()` and some plotting functions), since it was replaced with `combine_SI` three years ago
* Interpretation **guidelines older than 10 years were removed**, the oldest now included guidelines of EUCAST and CLSI are from 2013
* Using `units` in `ab_ddd(..., units = "...")` had been deprecated and is now not supported anymore. Use `ab_ddd_units()` instead.

### New
* **EUCAST 2022 and CLSI 2022 guidelines** have been added for `as.rsi()`. EUCAST 2022 is now the new default guideline for all MIC and disks diffusion interpretations.
* Support for the following languages: Chinese, Greek, Japanese, Polish, Turkish and Ukrainian. We are very grateful for the valuable input by our colleagues from other countries. The `AMR` package is now available in 16 languages. The automatic language determination will give a note at start-up on systems in supported languages.
* **All new algorithm for `as.mo()`** (and thus all `mo_*()` functions) while still following our original set-up as described in our recently submitted JSS paper (DOI [10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)).
  * A new argument `keep_synonyms` allows to *not* correct for updated taxonomy, in favour of the now deleted argument `allow_uncertain`
  * It has increased tremendously in speed and returns generally more consequent results
  * Sequential coercion is now extremely fast as results are stored to the package environment, although coercion of unknown values must be run once per session. Previous results can be reset/removed with the new `mo_reset_session()` function.
  * Support for microorganism codes of the ASIan Antimicrobial Resistance Surveillance Network (ASIARS-Net)
* **Extensive support for antiviral agents!** For the first time, the `AMR` package has extensive support for antiviral drugs and to work with their names, codes and other data in any way.
  * The `antivirals` data set has been extended with 18 new drugs (also from the [new J05AJ ATC group](https://www.whocc.no/atc_ddd_index/?code=J05AJ&showdescription=no)) and now also contains antiviral identifiers and LOINC codes
  * A new data type `av` (*antivirals*) has been added, which is functionally similar to `ab` for antibiotics
  * Functions `as.av()`, `av_name()`, `av_atc()`, `av_synonyms()`, `av_from_text()` have all been added as siblings to their `ab_*()` equivalents
* **Other new functions!**
  * Function `rsi_confidence_interval()` to add confidence intervals in AMR calculation. This is also included in `rsi_df()` and `proportion_df()`
  * Function `mean_amr_distance()` to calculate the mean AMR distance. The mean AMR distance is a normalised numeric value to compare AMR test results and can help to identify similar isolates, without comparing antibiograms by hand.
  * Function `rsi_interpretation_history()` to view the history of previous runs of `as.rsi()`. This returns a 'logbook' with the selected guideline, reference table and specific interpretation of each row in a data set on which `as.rsi()` was run.
  * Function `mo_current()` to get the currently valid taxonomic name of a microorganism
  * Function `add_custom_antimicrobials()` to add custom antimicrobial codes and names to the `AMR` package
* New and updated entries for the `antibiotics` data set
  * The following **20 antibiotics have been added** (also includes the [new J01RA ATC group](https://www.whocc.no/atc_ddd_index/?code=J01RA&showdescription=no)): azithromycin/fluconazole/secnidazole (AFC), cefepime/amikacin (CFA), cefixime/ornidazole (CEO), ceftriaxone/beta-lactamase inhibitor (CEB), ciprofloxacin/metronidazole (CIM), ciprofloxacin/ornidazole (CIO), ciprofloxacin/tinidazole (CIT), furazidin (FUR), isoniazid/sulfamethoxazole/trimethoprim/pyridoxine (IST), lascufloxacin (LSC), levofloxacin/ornidazole (LEO), nemonoxacin (NEM), norfloxacin/metronidazole (NME), norfloxacin/tinidazole (NTI), ofloxacin/ornidazole (OOR), oteseconazole (OTE), rifampicin/ethambutol/isoniazid (REI), sarecycline (SRC), tetracycline/oleandomycin (TOL), and thioacetazone (TAT)
  * Added some missing ATC codes
  * Updated DDDs and PubChem Compound IDs
  * Updated some antibiotic name spelling, now used by WHOCC (such as cephalexin -> cefalexin, and phenethicillin -> pheneticillin)
  * Antibiotic code "CEI" for ceftolozane/tazobactam has been replaced with "CZT" to comply with EARS-Net and WHONET 2022. The old code will still work in all cases when using `as.ab()` or any of the `ab_*()` functions.
  * Support for antimicrobial interpretation of anaerobic bacteria, by adding a 'placeholder' code `B_ANAER` to the `microorganisms` data set and add the breakpoints of anaerobics to the `rsi_interpretation` data set, which is used by `as.rsi()` when interpreting MIC and disk diffusion values
* Support for `data.frame`-enhancing R packages, more specifically: `data.table::data.table`, `janitor::tabyl`, `tibble::tibble`, and `tsibble::tsibble`. AMR package functions that have a data set as output (such as `rsi_df()` and `bug_drug_combinations()`), will now return the same data type as the input.
* All data sets in this package are **now exported as `tibble`**, instead of base R `data.frame`s. Older R versions are still supported.
* Our data sets are now also continually exported to **Apache Feather and Apache Parquet formats**. You can find more info [in this article on our website](https://msberends.github.io/AMR/articles/datasets.html).
* Support for using antibiotic selectors in scoped `dplyr` verbs (with or without `vars()`), such as in: `... %>% summarise_at(aminoglycosides(), resistance)`, see `resistance()`

### Changed
* Fix for using `as.rsi()` on certain EUCAST breakpoints for MIC values
* Fix for using `as.rsi()` on `NA` values (e.g. `as.rsi(as.disk(NA), ...)`)
* Fix for using `as.rsi()` on bug-drug combinations with multiple breakpoints for different body sites
* Removed `as.integer()` for MIC values, since MIC are not integer values and running `table()` on MIC values  consequently failed for not being able to retrieve the level position (as that's how normally `as.integer()` on `factor`s work)
* `droplevels()` on MIC will now return a common `factor` at default and will lose the `mic` class. Use `droplevels(..., as.mic = TRUE)` to keep the `mic` class.
* Small fix for using `ab_from_text()`
* Fixes for reading in text files using `set_mo_source()`, which now also allows the source file to contain valid taxonomic names instead of only valid microorganism ID of this package
* Using any `random_*()` function (such as `random_mic()`) is now possible by directly calling the package without loading it first: `AMR::random_mic(10)`
* Added *Toxoplasma gondii* (`P_TXPL_GOND`) to the `microorganisms` data set, together with its genus, family, and order
* Changed value in column `prevalence` of the `microorganisms` data set from 3 to 2 for these genera: *Acholeplasma*, *Alistipes*, *Alloprevotella*, *Bergeyella*, *Borrelia*, *Brachyspira*, *Butyricimonas*, *Cetobacterium*, *Chlamydia*, *Chlamydophila*, *Deinococcus*, *Dysgonomonas*, *Elizabethkingia*, *Empedobacter*, *Haloarcula*, *Halobacterium*, *Halococcus*, *Myroides*, *Odoribacter*, *Ornithobacterium*, *Parabacteroides*, *Pedobacter*, *Phocaeicola*, *Porphyromonas*, *Riemerella*, *Sphingobacterium*, *Streptobacillus*, *Tenacibaculum*, *Terrimonas*, *Victivallis*, *Wautersiella*, *Weeksella*
* Extended support for the `vctrs` package, used internally by the tidyverse. This allows to change values of class `mic`, `disk`, `rsi`, `mo` and `ab` in tibbles, and to use antibiotic selectors for selecting/filtering, e.g. `df[carbapenems() == "R", ]`
* Fix for using `info = FALSE` in `mdro()`
* For all interpretation guidelines using `as.rsi()` on amoxicillin, the rules for ampicillin will be used if amoxicillin rules are not available
* Fix for using `ab_atc()` on non-existing ATC codes
* Black and white message texts are now reversed in colour if using an RStudio dark theme
* `mo_snomed()` now returns class `character`, not `numeric` anymore (to make long SNOMED codes readable)
* Fix for using `as.ab()` on `NA` values
* Updated support for all WHONET 2022 microorganism codes 
* Antimicrobial interpretation 'SDD' (susceptible dose-dependent, coined by CLSI) will be interpreted as 'I' to comply with EUCAST's 'I' in `as.rsi()`
* Fix for `mo_shortname()` in case of higher taxonomic ranks (order, class, phylum)

### Other
* New website to make use of the new Bootstrap 5 and pkgdown 2.0. The website now contains results for all examples and will be automatically regenerated with every change to our repository, using GitHub Actions
* Added Peter Dutey-Magni, Dmytro Mykhailenko, Anton Mymrikov, and Jonas Salm as contributors, to thank them for their valuable input
* All R and Rmd files in this project are now styled using the `styler` package
* Set scalar conditional expressions (`&&` and `||`) where possible to comply with the upcoming R 4.3
* An enormous lot of code cleaning, fixing some small bugs on the way


# AMR 1.8.2

This is a small intermediate update to include the reference to our publication in the Journal of Statistical Software, DOI 10.18637/jss.v104.i03.

A major update will be released by the end of 2022 or early 2023 to include the most recent EUCAST and CLSI guidelines, updated microbial taxonomy, and support for 16 languages.


# AMR 1.8.1

### Changed
* Fix for using `as.rsi()` on values containing capped values (such as `>=`), sometimes leading to `NA`
* Support for antibiotic interpretations of the MIPS laboratory system: `"U"` for S ('susceptible urine'), `"D"` for I ('susceptible dose-dependent')
* Improved algorithm of `as.mo()`, especially for ignoring non-taxonomic text, such as:
  ```r
  mo_name("methicillin-resistant S. aureus (MRSA)")
  #> [1] "Staphylococcus aureus"
  ```
* More informative warning messages
* Added 192 as valid MIC
* Updated MIC printing in tibbles
* Increased speed for loading the package

### Other
* Fix for unit testing on R 3.3
* Fix for size of some image elements, as requested by CRAN


# AMR 1.8.0

### Breaking changes
* Removed `p_symbol()` and all `filter_*()` functions (except for `filter_first_isolate()`), which were all deprecated in a previous package version
* Removed the `key_antibiotics()` and `key_antibiotics_equal()` functions, which were deprecated and superseded by `key_antimicrobials()` and `antimicrobials_equal()`
* Removed all previously implemented `ggplot2::ggplot()` generics for classes `<mic>`, `<disk>`, `<rsi>` and `<resistance_predict>` as they did not follow the `ggplot2` logic. They were replaced with `ggplot2::autoplot()` generics.
* Renamed function `get_locale()` to `get_AMR_locale()` to prevent conflicts with other packages

### New
* Support for the CLSI 2021 guideline for interpreting MIC/disk diffusion values, which are incorporated in the `rsi_translation` data set. This data set now more strictly follows the WHONET software as well.
* Support for EUCAST Intrinsic Resistance and Unusual Phenotypes v3.3 (October 2021). This is now the default EUCAST guideline in the package (all older guidelines are still available) for `eucast_rules()`, `mo_is_intrinsic_resistant()` and `mdro()`. The `intrinsic_resistant` data set was also updated accordingly.
* Support for all antimicrobial drug (group) names and colloquial microorganism names in Danish, Dutch, English, French, German, Italian, Portuguese, Russian, Spanish and Swedish
* Function `set_ab_names()` to rename data set columns that resemble antimicrobial drugs. This allows for quickly renaming columns to official names, ATC codes, etc. Its second argument can be a tidyverse way of selecting:
  ```r
  example_isolates %>% set_ab_names(where(is.rsi))
  example_isolates %>% set_ab_names(AMC:GEN, property = "atc")
  ```
* Function `mo_lpsn()` to retrieve the [LPSN](https://lpsn.dsmz.de) record ID
* Function `ab_ddd_units()` to get units of DDDs (daily defined doses), deprecating the use of `ab_ddd(..., units = TRUE)` to be more consistent in data types of function output

### Changed
* Updated the bacterial taxonomy to 5 October 2021 (according to [LPSN](https://lpsn.dsmz.de)), including all 11 new staphylococcal species named since 1 January last year
* The `antibiotics` data set now contains **all ATC codes** that are available through the [WHOCC website](https://www.whocc.no), regardless of drugs being present in more than one ATC group. This means that:
  * Some drugs now contain multiple ATC codes (e.g., metronidazole contains 5)
  * `antibiotics$atc` is now a `list` containing `character` vectors, and this `atc` column was moved to the 5th position of the `antibiotics` data set
  * `ab_atc()` does not always return a character vector of length 1, and returns a `list` if the input is larger than length 1
  * `ab_info()` has a slightly different output
  * Some DDDs (daily defined doses) were added or updated according to newly included ATC codes
* Antibiotic selectors
  * They now also work in R-3.0 and R-3.1, supporting every version of R since 2013 like the rest of the package
  * Added more selectors for antibiotic classes: `aminopenicillins()`, `antifungals()`, `antimycobacterials()`, `lincosamides()`, `lipoglycopeptides()`, `polymyxins()`, `quinolones()`, `streptogramins()`, `trimethoprims()` and `ureidopenicillins()`
  * Added specific selectors for certain types for treatment: `administrable_per_os()` and `administrable_iv()`, which are based on available Defined Daily Doses (DDDs), as defined by the WHOCC. These are ideal for e.g. analysing pathogens in primary care where IV treatment is not an option. They can be combined with other AB selectors, e.g. to select penicillins that are only administrable per os (i.e., orally):
    ```r
    example_isolates[, penicillins() & administrable_per_os()]          # base R
    example_isolates %>% select(penicillins() & administrable_per_os()) # dplyr
    ```
  * Added the selector `ab_selector()`, which accepts a filter to be used internally on the `antibiotics` data set, yielding great flexibility on drug properties, such as selecting antibiotic columns with an oral DDD of at least 1 gram:
    ```r
    example_isolates[, ab_selector(oral_ddd > 1 & oral_units == "g")]          # base R
    example_isolates %>% select(ab_selector(oral_ddd > 1 & oral_units == "g")) # dplyr
    ```
  * Added the selector `not_intrinsic_resistant()`, which only keeps antibiotic columns that are not intrinsic resistant for all microorganisms in a data set, based on the latest EUCAST guideline on intrinsic resistance. For example, if a data set contains only microorganism codes or names of *E. coli* and *K. pneumoniae* and contains a column "vancomycin", this column will be removed (or rather, unselected) using this function.
  * Added argument `only_treatable`, which defaults to `TRUE` and will exclude drugs that are only for laboratory tests and not for treating patients (such as imipenem/EDTA and gentamicin-high)
  * Fix for using selectors multiple times in one call (e.g., using them in `dplyr::filter()` and immediately after in `dplyr::select()`)
  * Fix for using having multiple columns that are coerced to the same antibiotic agent
  * Fixed for using `all()` or `any()` on antibiotic selectors in an R Markdown file
* Added the following antimicrobial agents that are now covered by the WHO: aztreonam/nacubactam (ANC), cefepime/nacubactam (FNC), exebacase (EXE), ozenoxacin (OZN), zoliflodacin (ZFD), manogepix (MGX), ibrexafungerp (IBX), and rezafungin (RZF). None of these agents have an ATC code yet.
* Fixed the Gram stain (`mo_gramstain()`) determination of the taxonomic class Negativicutes within the phylum of Firmicutes - they were considered Gram-positives because of their phylum but are actually Gram-negative. This impacts 137 taxonomic species, genera and families, such as *Negativicoccus* and *Veillonella*.
* Dramatic speed improvement for `first_isolate()`
* Fix to prevent introducing `NA`s for old MO codes when running `as.mo()` on them
* Added more informative error messages when any of the `proportion_*()` and `count_*()` functions fail
* When printing a tibble with any old MO code, a warning will be thrown that old codes should be updated using `as.mo()`
* Improved automatic column selector when `col_*` arguments are left blank, e.g. in `first_isolate()`
* The right input types for `random_mic()`, `random_disk()` and `random_rsi()` are now enforced
* `as.rsi()` has an improved algorithm and can now also correct for textual input (such as "Susceptible", "Resistant") in all supported languages
* `as.mic()` has an improved algorithm
* When warnings are thrown because of too few isolates in any `count_*()`, `proportion_*()` function (or `resistant()` or `susceptible()`), the `dplyr` group will be shown, if available
* Fix for legends created with `scale_rsi_colours()` when using `ggplot2` v3.3.4 or higher (this is ggplot2 bug 4511, soon to be fixed)
* Fix for minor translation errors
* Fix for the MIC interpretation of *Morganellaceae* (such as *Morganella* and *Proteus*) when using the EUCAST 2021 guideline
* Improved algorithm of `as.mo()`
* Improved algorithm for generating random MICs with `random_mic()`
* Improved plot legends for MICs and disk diffusion values
* Improved speed of `as.ab()` and all `ab_*()` functions
* Added `fortify()` extensions for plotting methods
* `NA` values of the classes `<mic>`, `<disk>` and `<rsi>` are now exported objects of this package, e.g. `NA_mic_` is an `NA` of class `mic` (just like the base R `NA_character_` is an `NA` of class `character`)
* The `proportion_df()`, `count_df()` and `rsi_df()` functions now return with the additional S3 class 'rsi_df' so they can be extended by other packages
* The `mdro()` function now returns `NA` for all rows that have no test results
* The `species_id` column in the `microorganisms` data set now only contains LPSN record numbers. For this reason, this column is now numeric instead of a character, and `mo_url()` has been updated to reflect this change.
* Fixed a small bug in the functions `get_episode()` and `is_new_episode()`
* `get_episode()` and `is_new_episode()` can now cope with `NA`s

### Other
* This package is now being maintained by two epidemiologists and a data scientist from two different non-profit healthcare organisations.


# AMR older versions

For the changelog of older versions, please see [our archive](https://github.com/msberends/AMR/blob/v1.8.0/NEWS.md).
