# AMR 2.0.0.9024

## New
* Clinical breakpoints and intrinsic resistance of EUCAST 2023 and CLSI 2023 have been added for `as.sir()`. EUCAST 2023 (v13.0) is now the new default guideline for all MIC and disks diffusion interpretations, and for `eucast_rules()` to apply EUCAST Expert Rules.
* The EUCAST dosage guideline of v13.0 has been added to the `dosage` data set
* ECOFF: the `clinical_breakpoints` data set now contains the new column `ecoff`, in which the epidemiological cut-off (ECOFF) are available. These ECOFFs can be used for MIC/disk interpretation using `as.sir(..., ecoff = TRUE)`, which is an important new addition for veterinary microbiology.
* Added oxygen tolerance from BacDive to over 25,000 bacteria in the `microorganisms` data set
  * Added `mo_oxygen_tolerance()` to retrieve the values
  * Added `mo_is_anaerobic()` to determine which genera/species are obligate anaerobic bacteria
* Added LPSN and GBIF identifiers, and oxygen tolerance to `mo_info()`
* Added SAS Transport files (file extension `.xpt`) to [our download page](https://msberends.github.io/AMR/articles/datasets.html) to use in SAS software
* Added microbial codes for Gram-negative/positive anaerobic bacteria

## Changed
* `mo_rank()` now returns `NA` for 'unknown' microorganisms (`B_ANAER`, `B_ANAER-NEG`, `B_ANAER-POS`, `B_GRAMN`, `B_GRAMP`, `F_FUNGUS`, `F_YEAST`, and `UNKNOWN`)
* Fixed formatting for `sir_interpretation_history()`
* Fixed some WHONET codes for microorganisms and consequently a couple of entries in `clinical_breakpoints`
* Fixed a bug for `as.mo()` that led to coercion of `NA` values when using custom microorganism codes
* Fixed usage of `icu_exclude` in `first_isolates()`
* Improved `as.mo()` algorithm:
  * Now allows searching on only species names
  * Fix for using the `keep_synonyms` argument when using MO codes as input
  * Fix for using the `minimum_matching_score` argument
* Updated the code table in `microorganisms.codes`


# AMR 2.0.0

This is a new major release of the AMR package, with great new additions but also some breaking changes for current users. These are all listed below.

**[TL;DR](https://en.wikipedia.org/wiki/TL;DR)**

* All functions and arguments with 'rsi' were replaced with 'sir', such as the interpretation of MIC values (now `as.sir()` instead of `as.rsi()`) - all old functions still work for now
* Many new interesting functions, such as `antibiogram()` (for generating traditional/combined/syndromic/WISCA antibiograms), `sir_confidence_interval()` and `mean_amr_distance()`, and `add_custom_microorganisms()` to add custom microorganisms to this package
* Clinical breakpoints added for EUCAST 2022 and CLSI 2022
* Microbiological taxonomy (`microorganisms` data set) updated to 2022 and now based on LPSN and GBIF
* Much increased algorithms to translate user input to valid taxonomy, e.g. by using [recent scientific work](https://doi.org/10.1099/mic.0.001269) about per-species human pathogenicity
* 20 new antibiotics added and updated all DDDs and ATC codes
* Extended support for antiviral agents (`antivirals` data set), with many new functions
* Now available in 20 languages
* Many small bug fixes

## New

### SIR vs. RSI

For this milestone version, we replaced all mentions of RSI with SIR, to comply with what is actually being commonly used in the field of clinical microbiology when it comes to this tri-form regarding AMR.

While existing functions such as `as.rsi()`, `rsi_df()` and `ggplot_rsi()` still work, their replacements `as.sir()`,  `sir_df()`, `ggplot_sir()` are now the current functions for AMR data analysis. A warning will be thrown once a session to remind users about this. The data set `rsi_translation` is now called `clinical_breakpoints` to better reflect its content.

The 'RSI functions' will be removed in a future version, but not before late 2023 / early 2024.

### New antibiogram function

With the new `antibiogram()` function, users can now generate traditional, combined, syndromic, and even weighted-incidence syndromic combination antibiograms (WISCA). With this, we follow the logic in the previously described work of Klinker *et al.* (2021, DOI [10.1177/20499361211011373](https://doi.org/10.1177/20499361211011373)) and Barbieri *et al.* (2021, DOI [10.1186/s13756-021-00939-2](https://doi.org/10.1186/s13756-021-00939-2)).

The help page for `antibiogram()` extensively elaborates on use cases, and `antibiogram()` also supports printing in R Markdown and Quarto, with support for 20 languages.

Furthermore, different plotting methods were implemented to allow for graphical visualisations as well.

### Interpretation of MIC and disk diffusion values

The clinical breakpoints and intrinsic resistance of EUCAST 2022 and CLSI 2022 have been added for `as.sir()`. EUCAST 2022 (v12.0) is now the new default guideline for all MIC and disks diffusion interpretations, and for `eucast_rules()` to apply EUCAST Expert Rules. The default guideline (EUCAST) can now be changed with the new `AMR_guideline` option, such as: `options(AMR_guideline = "CLSI 2020")`.

With the new arguments `include_PKPD` (default: `TRUE`) and `include_screening` (default: `FALSE`), users can now specify whether breakpoints for screening and from the PK/PD table should be included when interpreting MICs and disks diffusion values. These options can be set globally, which can be read in [our new manual](https://msberends.github.io/AMR/reference/AMR-options.html).
 
Interpretation guidelines older than 10 years were removed, the oldest now included guidelines of EUCAST and CLSI are from 2013.

### Supported languages

We added support for the following ten languages: Chinese (simplified), Czech, Finnish, Greek, Japanese, Norwegian (bokmål), Polish, Romanian, Turkish and Ukrainian. All antibiotic names are now available in these languages, and the AMR package will automatically determine a supported language based on the user's system language.

We are very grateful for the valuable input by our colleagues from other countries. The `AMR` package is now available in 20 languages in total, and according to download stats used in almost all countries in the world!

### Outbreak management

For analysis in outbreak management, we updated the `get_episode()` and `is_new_episode()` functions: they now contain an argument `case_free_days`. This argument can be used to quantify the duration of case-free days (the inter-epidemic interval), after which a new episode will start. 

This is common requirement in outbreak management, e.g. when determining the number of norovirus outbreaks in a hospital. The case-free period could then be 14 or 28 days, so that new norovirus cases after that time will be considered a different (or new) episode.

### Microbiological taxonomy

The `microorganisms` data set no longer relies on the Catalogue of Life, but on the List of Prokaryotic names with Standing in Nomenclature (LPSN) and is supplemented with the 'backbone taxonomy' from the Global Biodiversity Information Facility (GBIF). The structure of this data set has changed to include separate LPSN and GBIF identifiers. Almost all previous MO codes were retained. It contains over 1,400 taxonomic names from 2022.

We previously relied on our own experience to categorise species into pathogenic groups, but we were very happy to encounter the very recent work of Bartlett *et al.* (2022, DOI [10.1099/mic.0.001269](https://doi.org/10.1099/mic.0.001269)) who extensively studied medical-scientific literature to categorise all bacterial species into groups. See `mo_matching_score()` on how their work was incorporated into the `prevalence` column of the `microorganisms` data set. Using their results, the `as.mo()` and all `mo_*()` functions are now much better capable of converting user input to valid taxonomic records.

The new function `add_custom_microorganisms()` allows users to add custom microorganisms to the `AMR` package.

We also made the following changes regarding the included taxonomy or microorganisms functions:

* Updated full microbiological taxonomy according to the latest daily LPSN data set (December 2022) and latest yearly GBIF taxonomy backbone (November 2022)
* Added function `mo_current()` to get the currently valid taxonomic name of a microorganism
* Support for all 1,516 city-like serovars of *Salmonella*, such as *Salmonella* Goldcoast. Formally, these are serovars belonging to the *S. enterica* species, but they are reported with only the name of the genus and the city. For this reason, the serovars are in the `subspecies` column of the `microorganisms` data set and "enterica" is in the `species` column, but the full name does not contain the species name (*enterica*).
* All new algorithm for `as.mo()` (and thus all `mo_*()` functions) while still following our original set-up as described in our recently published JSS paper (DOI [10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)).
  * A new argument `keep_synonyms` allows to *not* correct for updated taxonomy, in favour of the now deleted argument `allow_uncertain`
  * It has increased tremendously in speed and returns generally more consequent results
  * Sequential coercion is now extremely fast as results are stored to the package environment, although coercion of unknown values must be run once per session. Previous results can be reset/removed with the new `mo_reset_session()` function.
  * Support for microorganism codes of the ASIan Antimicrobial Resistance Surveillance Network (ASIARS-Net)
  * The MO matching score algorithm (`mo_matching_score()`) now counts deletions and substitutions as 2 instead of 1, which impacts the outcome of `as.mo()` and any `mo_*()` function
* **Removed all species of the taxonomic kingdom Chromista** from the package. This was done for multiple reasons:
  * CRAN allows packages to be around 5 MB maximum, some packages are exempted but this package is not one of them
  * Chromista are not relevant when it comes to antimicrobial resistance, thus lacking the primary scope of this package
  * Chromista are almost never clinically relevant, thus lacking the secondary scope of this package
* The `microorganisms.old` data set was removed, and all previously accepted names are now included in the `microorganisms` data set. A new column `status` contains `"accepted"` for currently accepted names and `"synonym"` for taxonomic synonyms; currently invalid names. All previously accepted names now have a microorganisms ID and - if available - an LPSN, GBIF and SNOMED CT identifier.

### Antibiotic agents and selectors

The new function `add_custom_antimicrobials()` allows users to add custom antimicrobial codes and names to the `AMR` package.

The `antibiotics` data set was greatly updated:

  * The following 20 antibiotics have been added (also includes the [new J01RA ATC group](https://www.whocc.no/atc_ddd_index/?code=J01RA&showdescription=no)): azithromycin/fluconazole/secnidazole (AFC), cefepime/amikacin (CFA), cefixime/ornidazole (CEO), ceftriaxone/beta-lactamase inhibitor (CEB), ciprofloxacin/metronidazole (CIM), ciprofloxacin/ornidazole (CIO), ciprofloxacin/tinidazole (CIT), furazidin (FUR), isoniazid/sulfamethoxazole/trimethoprim/pyridoxine (IST), lascufloxacin (LSC), levofloxacin/ornidazole (LEO), nemonoxacin (NEM), norfloxacin/metronidazole (NME), norfloxacin/tinidazole (NTI), ofloxacin/ornidazole (OOR), oteseconazole (OTE), rifampicin/ethambutol/isoniazid (REI), sarecycline (SRC), tetracycline/oleandomycin (TOL), and thioacetazone (TAT)
  * Added some missing ATC codes
  * Updated DDDs and PubChem Compound IDs
  * Updated some antibiotic name spelling, now used by WHOCC (such as cephalexin -> cefalexin, and phenethicillin -> pheneticillin)
  * Antibiotic code "CEI" for ceftolozane/tazobactam has been replaced with "CZT" to comply with EARS-Net and WHONET 2022. The old code will still work in all cases when using `as.ab()` or any of the `ab_*()` functions.
  * Support for antimicrobial interpretation of anaerobic bacteria, by adding a 'placeholder' code `B_ANAER` to the `microorganisms` data set and adding the breakpoints of anaerobics to the `clinical_breakpoints` data set, which is used by `as.sir()` for interpretion of MIC and disk diffusion values

Also, we added support for using antibiotic selectors in scoped `dplyr` verbs (with or without using `vars()`), such as in: `... %>% summarise_at(aminoglycosides(), resistance)`, please see `resistance()` for examples.

### Antiviral agents

We now added extensive support for antiviral agents! For the first time, the `AMR` package has extensive support for antiviral drugs and to work with their names, codes and other data in any way.

* The `antivirals` data set has been extended with 18 new drugs (also from the [new J05AJ ATC group](https://www.whocc.no/atc_ddd_index/?code=J05AJ&showdescription=no)) and now also contains antiviral identifiers and LOINC codes
* A new data type `av` (*antivirals*) has been added, which is functionally similar to `ab` for antibiotics
* Functions `as.av()`, `av_name()`, `av_atc()`, `av_synonyms()`, `av_from_text()` have all been added as siblings to their `ab_*()` equivalents

### Other new functions

* Function `sir_confidence_interval()` to add confidence intervals in AMR calculation. This is now also included in `sir_df()` and `proportion_df()`.
* Function `mean_amr_distance()` to calculate the mean AMR distance. The mean AMR distance is a normalised numeric value to compare AMR test results and can help to identify similar isolates, without comparing antibiograms by hand.
* Function `sir_interpretation_history()` to view the history of previous runs of `as.sir()` (previously `as.rsi()`). This returns a 'logbook' with the selected guideline, reference table and specific interpretation of each row in a data set on which `as.sir()` was run.


## Changes

* `get_episode()` (and its wrapper `is_new_episode()`):
  * Fix for working with `NA` values
  * Fix for unsorted dates of length 2
  * Now returns class `integer` instead of `numeric` since they are always whole numbers
* Argument `combine_IR` has been removed from this package (affecting functions `count_df()`, `proportion_df()`, and `sir_df()` and some plotting functions), since it was replaced with `combine_SI` three years ago
* Using `units` in `ab_ddd(..., units = "...")` had been deprecated for some time and is now not supported anymore. Use `ab_ddd_units()` instead.
* Support for `data.frame`-enhancing R packages, more specifically: `data.table::data.table`, `janitor::tabyl`, `tibble::tibble`, and `tsibble::tsibble`. AMR package functions that have a data set as output (such as `sir_df()` and `bug_drug_combinations()`), will now return the same data type as the input.
* All data sets in this package are now a `tibble`, instead of base R `data.frame`s. Older R versions are still supported, even if they do not support `tibble`s.
* Our data sets are now also continually exported to **Apache Feather and Apache Parquet formats**. You can find more info [in this article on our website](https://msberends.github.io/AMR/articles/datasets.html).
* For `as.sir()`:
  * Fixed certain EUCAST breakpoints for MIC values
  * Allow `NA` values (e.g. `as.sir(as.disk(NA), ...)`)
  * Fix for bug-drug combinations with multiple breakpoints for different body sites
  * Interpretation from MIC and disk zones is now more informative about availability of breakpoints and more robust
* Removed the `as.integer()` method for MIC values, since MIC are not integer values and running `table()` on MIC values consequently failed for not being able to retrieve the level position (as that's how normally `as.integer()` on `factor`s work)
* Fixed determination of Gram stains (`mo_gramstain()`), since the taxonomic phyla Actinobacteria, Chloroflexi, Firmicutes, and Tenericutes have been renamed to respectively Actinomycetota, Chloroflexota, Bacillota, and Mycoplasmatota in 2021
* `droplevels()` on MIC will now return a common `factor` at default and will lose the `mic` class. Use `droplevels(..., as.mic = TRUE)` to keep the `mic` class.
* Small fix for using `ab_from_text()`
* Fixes for reading in text files using `set_mo_source()`, which now also allows the source file to contain valid taxonomic names instead of only valid microorganism ID of this package
* Fixed a bug for `mdro()` when using similar column names with the Magiorakos guideline
* Using any `random_*()` function (such as `random_mic()`) is now possible by directly calling the package without loading it first: `AMR::random_mic(10)`
* Extended support for the `vctrs` package, used internally by the tidyverse. This allows to change values of class `mic`, `disk`, `sir`, `mo` and `ab` in tibbles, and to use antibiotic selectors for selecting/filtering, e.g. `df[carbapenems() == "R", ]`
* Fix for using `info = FALSE` in `mdro()`
* For all interpretation guidelines using `as.sir()` on amoxicillin, the rules for ampicillin will be used if amoxicillin rules are not available
* Fix for using `ab_atc()` on non-existing ATC codes
* Black and white message texts are now reversed in colour if using an RStudio dark theme
* `mo_snomed()` now returns class `character`, not `numeric` anymore (to make long SNOMED codes readable)
* Fix for using `as.ab()` on `NA` values
* Updated support for all WHONET 2022 microorganism codes 
* Antimicrobial interpretation 'SDD' (susceptible dose-dependent, coined by CLSI) will be interpreted as 'I' to comply with EUCAST's 'I' in `as.sir()`
* Fix for `mo_shortname()` in case of higher taxonomic ranks (order, class, phylum)
* Cleaning columns with `as.sir()`, `as.mic()`, or `as.disk()` will now show the column name in the warning for invalid results
* Fix for using `g.test()` with zeroes in a 2x2 table
* `mo_synonyns()` now contains the scientific reference as names

## Other

* Added Peter Dutey-Magni, Dmytro Mykhailenko, Anton Mymrikov, Andrew Norgan, Jonas Salm, and Anita Williams as contributors, to thank them for their valuable input
* New website to make use of the new Bootstrap 5 and pkgdown 2.0. The website now contains results for all examples and will be automatically regenerated with every change to our repository, using GitHub Actions
* All R and Rmd files in this project are now styled using the `styler` package
* Set scalar conditional expressions (`&&` and `||`) where possible to comply with the upcoming R 4.3
* An enormous lot of code cleaning, fixing some small bugs along the way

----

This changelog only contains changes from AMR v2.0 and later. For prior versions, please see [our archive](https://github.com/msberends/AMR/blob/v1.8.2/NEWS.md).
