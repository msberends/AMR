# AMR 3.0.1.9067

Planned as v3.1.0, end of June 2026.

### Breaking Changes
* The former *kingdoms* Bacteria and Archaea are now each divided into four kingdoms with new top-level *domains* 'Bacteria' and 'Archaea' (Göker and Oren, 2024, DOI: 10.1099/ijsem.0.006242). Following this, a new `domain` column in the `microorganisms` data set was added, and more importantly, `mo_kingdom()` now returns the formal kingdom (e.g. `"Pseudomonadati"` instead of `"Bacteria"`). Use `mo_domain()` for the old behaviour. For non-prokaryotic kingdoms (Fungi, Protozoa, etc.), `kingdom` and `domain` are identical.
* Faster parallel computing via the `future` package for `as.sir()` and `wisca()`: a non-sequential plan (e.g. `future::plan(future::multisession)`) must be active before using `parallel = TRUE`.

### New
* EUCAST 2026 and CLSI 2026 breakpoints: over 5,700 new breakpoints added to the `clinical_breakpoints` data set; EUCAST 2026 is now the default for all MIC and disk diffusion interpretations
* Wildtype/Non-wildtype (WT/NWT) output when using ECOFF-based interpretation, by setting `breakpoint_type = "ECOFF"` in `as.sir()`; WT/NWT results are fully supported in all resistance/susceptibility functions and plots (#254)
* *tidymodels* integration for using SIR, MIC and disk data in modelling pipelines: `step_mic_log2()`, `step_sir_numeric()`, and new column selectors `all_sir()`, `all_mic()`, `all_disk()`
* New `esbl_isolates` data set for practising AMR modelling
* New antimicrobial selectors: `ionophores()`, `peptides()`, `phosphonics()`, `spiropyrimidinetriones()`
* New antimicrobials: cefepime/taniborbactam (`FTA`), ceftibuten/avibactam (`CTA`), clorobiocin (`CLB`), kasugamycin (`KAS`), ostreogrycin (`OST`), taniborbactam (`TAN`), thiostrepton (`THS`), xeruborbactam (`XER`), zorbamycin (`ZOR`)
* New `interpretive_rules()`, a unified function for EUCAST and CLSI interpretive rules; `eucast_rules()` is now a wrapper around it (#235, #259)
* New `morphology` column in the `microorganisms` data set and corresponding `mo_morphology()` function, returning the cell shape of bacteria. Data sourced from BacDive; values prefixed with "likely" are extrapolated from genus-level consensus. New `add_morphology` argument was added to `mo_gramstain()` to return combined results such as `"Gram-negative rods"`.
* New `amr_course()` to download and unpack course or webinar materials from GitHub in one call
* Typed missing value constants `NA_ab_` and `NA_mo_`, for use in pipelines that need missing values of a specific class
* New `wisca_plot()` to assess the susceptibility and incidence distributions from the Monte Carlo simulations

### Fixed
* `as.sir()`
  * On data frames: already-converted SIR columns no longer dropped on re-run (#278)
  * Metadata columns (e.g. `patient`, `ward`) no longer misidentified as antibiotic columns
  * `info = FALSE` now suppresses all messages, including for columns without breakpoints
  * Assumption of disk zones are now preferred over MIC values when input is only whole numbers (#291)
* `as.mic()`: values in scientific notation (e.g. `1e-3`) now handled correctly
* `as.ab()`: codes containing "PH" or "TH" (e.g. `ETH`, `PHE`) no longer return `NA` when mixed with unrecognised input (#245)
* Combined MIC/SIR input values (e.g. `"<= 0.002; S"` or `"S; 0.002"`) now parsed correctly (#252)
* `as.mo()`: 
  * Input of the form `"X complex"` now falls back to `"X"` when the complex is not a distinct taxon in the database, preventing `NA` results for valid clinical descriptions such as `"Proteus vulgaris complex"` (#287)
  * Abbreviated-genus input (e.g. `"S. apiospermum"`) now correctly ranks candidates whose species epithet exactly matches the input above more-prevalent organisms whose species does not match; fixes `"S. apiospermum"` resolving to *Staphylococcus* instead of *Scedosporium apiospermum* (#288)
* `get_author_year()` in the microorganism reproduction script now strips `emend.` and everything after it, so `ref` reflects the combination authority rather than the emendation author (e.g. *Rhodococcus equi* now returns "Goodfellow et al., 1977" instead of "Nouioui et al., 2018")
* BRMO classification now includes bacterial complexes (#275)
* Translation fixes for Italian CoNS/CoPS names (#256), Dutch antimicrobials, and `sir_df()` foreign-language output (#272)
* Fixed some EUCAST Expert Rules, mostly on *S. pneumoniae*

### Updated
* Taxonomic update for all microorganisms, now updated to June 2026
* `mo_kingdom()` now returns the formal taxonomic kingdom; a one-time note per session explains the change when querying bacterial or archaeal records.
* `mo_taxonomy()` and `mo_info()` gained `domain` for the list output
* `antibiogram()` and `wisca()` now also support parallel computing via the argument `parallel = TRUE` (#281)
* `custom_eucast_rules()` renamed to `custom_interpretive_rules()`; old name deprecated but still works (#268)
* `mdro()` can now infer resistance from a drug+inhibitor combination when the base drug column is absent (e.g. piperacillin inferred from piperacillin/tazobactam); controlled via new `infer_from_combinations` argument (default `TRUE`) (#209)
* `wisca()` now more strictly follows Bielicki et al. (2016) by using $\text{Beta}(1, 9999)$ for intrinsically resistant pairs, forcing near-zero susceptibility regardless of observed data (based on EUCAST Expected Resistant Phenotypes)
* `susceptibility()` / `resistance()`: new `guideline` argument (default EUCAST) to ensure the 'I' category is interpreted correctly per guideline
* Capped MIC handling in `as.sir()` reworked into four clearly defined options: `"none"`, `"conservative"` (new default), `"standard"`, `"lenient"` (#243)
* `as.mic()` / `rescale_mic()`: new `round_to_next_log2` argument to round values up to the nearest log2 dilution level (#255)
* `antimicrobials$group` is now a `list`, so that drugs belonging to multiple groups are fully represented; use `ab_group(all_groups = TRUE)` to retrieve all groups for a drug (#246)
* Improved console messages with clickable links throughout, powered by `cli` if it is installed (#191, #265)
* `as.disk()`: input validation is now more strict, rejecting values that are not recognisable as a numeric disk zone diameter


# AMR 3.0.1

This is a bugfix release following the release of v3.0.0 in June 2025.

### Changed
* Fixed bugs introduced by `ggplot2` v4.0.0 (#236)
  * MIC scale functions (such as `scale_y_mic()`) will now be applied automatically when plotting values of class `mic`
  * SIR scale functions (such as `scale_x_sir()`) will now be applied automatically when plotting values of class `sir`
* Fixed a bug in `antibiogram()` for when no antimicrobials are set
* Fixed a bug in `antibiogram()` to allow column names containing the `+` character (#222)
* Fixed a bug in `as.ab()` for antimicrobial codes with a number in it if they are preceded by a space
* Fixed a bug in `eucast_rules()` for using specific custom rules
* Fixed a bug in `as.sir()` to allow any tidyselect language (#220)
* Fixed a bug in `as.sir()` to pick right breakpoint when `uti = FALSE` (#216)
* Fixed a bug in `ggplot_sir()` when using `combine_SI = FALSE` (#213)
* Fixed a bug in `mdro()` to make sure all genes specified in arguments are acknowledged
* Fixed a bug the `antimicrobials` data set to remove statins (#229)
* Fixed a bug the `microorganisms` data set for MycoBank IDs and synonyms (#233)
* Fixed ATC J01CR05 to map to piperacillin/tazobactam rather than piperacillin/sulbactam (#230)
* Fixed skimmers (`skimr` package) of class `ab`, `sir`, and `disk` (#234)
* Fixed all plotting to contain a separate colour for SDD (susceptible dose-dependent) (#223)
* Fixed some specific Dutch translations for antimicrobials
* Added a warning to `as.ab()` if input resembles antiviral codes or names (#232)
* Added all reasons in verbose output of `mdro()` (#227)
* Added `names` to `age_groups()` so that custom names can be given (#215)
* Added note to `as.sir()` to make it explicit when higher-level taxonomic breakpoints are used (#218)
* Added antibiotic codes from the Comprehensive Antibiotic Resistance Database (CARD) to the `antimicrobials` data set (#225)
* Updated Fosfomycin to be of antibiotic class Phosphonics (#225)
* Updated `random_mic()` and `random_disk()` to set skewedness of the distribution and allow multiple microorganisms


# AMR 3.0.0

This package now supports not only tools for AMR data analysis in clinical settings, but also for veterinary and environmental microbiology. This was made possible through a collaboration with the [University of Prince Edward Island's Atlantic Veterinary College](https://www.upei.ca/avc), Canada. To celebrate this great improvement of the package, we also updated the package logo to reflect this change.

### Breaking
* Data set `antibiotics` has been renamed to `antimicrobials` as the data set contains more than just antibiotics. Using `antibiotics` will still work, but now returns a warning.
* Removed all functions and references that used the deprecated `rsi` class, which were all replaced with their `sir` equivalents over two years ago.
* Functions `resistance_predict()` and `sir_predict()` are now deprecated and will be removed in a future version. Use the `tidymodels` framework instead, for which we [wrote a basic introduction](https://amr-for-r.org/articles/AMR_with_tidymodels.html).

### New
* **One Health implementation**
  * Function `as.sir()` now has extensive support for veterinary breakpoints from CLSI. Use `breakpoint_type = "animal"` and set the `host` argument to a variable that contains animal species names.
  * The `clinical_breakpoints` data set contains all these breakpoints, and can be downloaded on our [download page](https://amr-for-r.org/articles/datasets.html).
  * The (new) `antimicrobials` data set contains all veterinary antimicrobials, such as pradofloxacin and enrofloxacin. All WHOCC codes for veterinary use have been added as well.
  * `ab_atc()` now supports ATC codes of veterinary antimicrobials (that all start with "Q")
  * `ab_url()` now supports retrieving the WHOCC url of their ATCvet pages
* **Support for WISCA antibiograms**  
  * The `antibiogram()` function now supports creating true Weighted-Incidence Syndromic Combination Antibiograms (WISCA), a powerful Bayesian method for estimating regimen coverage probabilities using pathogen incidence and antimicrobial susceptibility data. WISCA offers improved precision for syndrome-specific treatment, even in data sets with sparse data. A dedicated `wisca()` function is also available for easy usage.
* **More global coverage of languages**
  * Added full support for 8 new languages: Arabic, Bengali, Hindi, Indonesian, Korean, Swahili, Urdu, and Vietnamese. The `AMR` package is now available in 28 languages.
* **Major update to fungal taxonomy and tools for mycologists**
  * MycoBank has now been integrated as the primary taxonomic source for fungi. The `microorganisms` data set has been enriched with new columns (`mycobank`, `mycobank_parent`, and `mycobank_renamed_to`) that provide detailed information for fungal species.
  * A remarkable addition of over 20,000 new fungal records
  * New function `mo_mycobank()` to retrieve the MycoBank record number, analogous to existing functions such as `mo_lpsn()` and `mo_gbif()`.
  * The `as.mo()` function and all `mo_*()` functions now include an `only_fungi` argument, allowing users to restrict results solely to fungal species. This ensures fungi are prioritised over bacteria during microorganism identification. This can also be set globally with the new `AMR_only_fungi` option.
  * Also updated other kingdoms, welcoming a total of 2,149 new records from 2023 and 927 from 2024.
* **Updated clinical breakpoints**
  * Breakpoint of 2024 and 2025 of both CLSI and EUCAST are now supported, by adding all of their over 10,000 new clinical breakpoints to the `clinical_breakpoints` data set for usage in `as.sir()`. EUCAST 2025 is now the new default guideline for all MIC and disk diffusion interpretations.
  * Added all Expected Resistant Phenotypes from EUCAST (v1.2). The default `rules` for `eucast_rules()` are now: `c("breakpoints", "expected_phenotypes")`.
  * Updated the `intrinsic_resistant` data set, which is now based on EUCAST Expected Resistant Phenotypes v1.2
  * `as.sir()` now brings additional factor levels: "NI" for non-interpretable and "SDD" for susceptible dose-dependent. Currently, the `clinical_breakpoints` data set contains 24 breakpoints that can return the value "SDD" instead of "I".
  * EUCAST interpretive rules (using `eucast_rules()`) are now available for EUCAST 12 (2022), 13 (2023), 14 (2024), and 15 (2025).
  * EUCAST dosage tables (`dosage` data set) are now available for EUCAST 13 (2023), 14 (2024), and 15 (2025).
* **New advanced ggplot2 extensions for MIC and SIR plotting and transforming**
  * New function group `scale_*_mic()`, namely: `scale_x_mic()`, `scale_y_mic()`, `scale_colour_mic()` and `scale_fill_mic()`. They allow easy plotting of MIC values. They allow for manual range definition and plotting missing intermediate log2 levels.
  * New function group `scale_*_sir()`, namely: `scale_x_sir()`, `scale_colour_sir()` and `scale_fill_sir()`. They allow to plot the `sir` class, and translates into the system language at default. They also set colourblind-safe colours to the plots.
  * New function `rescale_mic()`, which allows users to rescale MIC values to a manually set range. This is the powerhouse behind the `scale_*_mic()` functions, but it can be used independently to, for instance, compare equality in MIC distributions by rescaling them to the same range first.
* **Support for Python**
  * While using R for the heavy lifting, [our 'AMR' Python Package](https://pypi.org/project/AMR/) was developed to run the AMR R package natively in Python. The Python package will always have the same version number as the R package, as it is built automatically with every code change.
* **Support for `tidymodels`**
  * All antimicrobial selectors (such as `aminoglycosides()` and `betalactams()`) are now supported in `tidymodels` packages such as `recipe` and `parsnip`. See for more info [our tutorial](https://amr-for-r.org/articles/AMR_with_tidymodels.html) on using these AMR functions for predictive modelling.
* **Other**
  * New function `top_n_microorganisms()` to filter a data set to the top *n* of any taxonomic property, e.g., filter to the top 3 species, filter to any species in the top 5 genera, or filter to the top 3 species in each of the top 5 genera
  * New function `mo_group_members()` to retrieve the member microorganisms of a microorganism group. For example, `mo_group_members("Strep group C")` returns a vector of all microorganisms that belong to that group.
  * New functions `mic_p50()` and `mic_p90()` to retrieve the 50th and 90th percentile of MIC values.

### Changed
* SIR interpretation
  * Support for parallel computing to greatly improve speed using the `parallel` package (part of base R). Use `as.sir(your_data, parallel = TRUE)` to run SIR interpretation using multiple cores.
  * It is now possible to use column names for arguments `guideline`, `ab`, `mo`, and `uti`: `as.sir(..., ab = "column1", mo = "column2", uti = "column3")`. This greatly improves the flexibility for users.
  * Users can now set their own criteria (using regular expressions) as to what should be considered S, I, R, SDD, and NI.
  * To get quantitative values, `as.double()` on a `sir` object will return 1 for S, 2 for SDD/I, and 3 for R (NI will become `NA`). Other functions using `sir` classes (e.g., `summary()`) are updated to reflect the change to contain NI and SDD.
  * Following CLSI interpretation rules, values outside the log2-dilution range will be rounded upwards to the nearest log2-level before interpretation. Only if using a CLSI guideline.
  * Combined MIC values (e.g., from CLSI) are now supported
  * The argument `conserve_capped_values` in `as.sir()` has been replaced with `capped_mic_handling`, which allows greater flexibility in handling capped MIC values (`<`, `<=`, `>`, `>=`). The four available options (`"standard"`, `"strict"`, `"relaxed"`, `"inverse"`) provide full control over whether these values should be interpreted conservatively or ignored. Using `conserve_capped_values` is now deprecated and returns a warning.
  * Added argument `info` to silence all console messages
* `antibiogram()` function
  * Argument `antibiotics` has been renamed to `antimicrobials`. Using `antibiotics` will still work, but now returns a warning.
  * Added argument `formatting_type` to set any of the 22 options for the formatting of all 'cells'. This defaults to `18` for non-WISCA and `14` for WISCA, changing the output of antibiograms to cells with more info.
  * For this reason, `add_total_n` is now deprecated and `FALSE` at default since the denominators are added to the cells dependent on the `formatting_type` setting
  * The `ab_transform` argument now defaults to `"name"`, displaying antibiotic column names instead of codes
* Antimicrobial selectors (previously: *antibiotic selectors*)
  * 'Antibiotic selectors' are now called 'antimicrobial selectors' since their scope is broader than just antibiotics. All documentation have been updated, and `ab_class()` and `ab_selector()` have been replaced with `amr_class()` and `amr_selector()`. The old functions are now deprecated and will be removed in a future version.
  * Added selectors `isoxazolylpenicillins()`, `monobactams()`, `nitrofurans()`, `phenicols()`, `rifamycins()`, and `sulfonamides()`
  * When using antimicrobial selectors that exclude non-treatable drugs (such as gentamicin-high when using `aminoglycosides()`), the function now always returns a warning that these can be included using `only_treatable = FALSE`
  * Added a new argument `return_all` to all selectors, which defaults to `TRUE` to include any match. With `FALSE`, the old behaviour, only the first hit for each unique antimicrobial is returned.
  * All selectors can now be run as a separate command to retrieve a vector of all possible antimicrobials that the selector can select
  * The selectors `lincosamides()` and `macrolides()` do not overlap anymore - each antibiotic is now classified as either of these and not both
  * Fixed selector `fluoroquinolones()`, which now really only selects second-generation quinolones and up (first-generation quinolones do not contain a fluorine group)
* `antimicrobials` data set
  * Added agents used for screening, with an ID all ending with `-S`: benzylpenicillin screening test (`PEN-S`), beta-lactamase screening test (`BLA-S`), cefotaxime screening test (`CTX-S`), clindamycin inducible screening test (`CLI-S`), nalidixic acid screening test (`NAL-S`), norfloxacin screening test (`NOR-S`), oxacillin screening test (`OXA-S`), pefloxacin screening test (`PEF-S`), and tetracycline screening test (`TCY-S`). The ID of cefoxitin screening was renamed from `FOX1` to `FOX-S`, while the old code remains to work.
  * For this reason, the antimicrobial selectors `cephalosporins()`, `cephalosporins_3rd()`, `lincosamides()`, `isoxazolylpenicillins()`, `quinolones()`, `fluoroquinolones()`, and `tetracyclines()` now contain the argument `only_treatable = TRUE` (similar to other antimicrobial selectors that contain non-treatable drugs)
  * Added amorolfine (`AMO`, D01AE16), an antimycotic, which is now also part of the `antifungals()` selector
  * Added cefepime/enmetazobactam (`FPE`), a 4th gen cephalosporin
  * Added tigemonam (`TNM`), a monobactam
  * Added bleomycin (`BLM`), a glycopeptide
  * Added efflux (`EFF`), to allow mapping to AMRFinderPlus
  * Updated all ATC codes, trade names, and DDDs
* MICs
  * Added as valid levels: 4096, 6 powers of 0.0625, and 5 powers of 192 (192, 384, 576, 768, 960)
  * Fixed a bug in `as.mic()` that failed translation of scientifically formatted numbers
  * Added new argument `keep_operators` to `as.mic()`. This can be `"all"` (default), `"none"`, or `"edges"`. This argument is also available in the new `rescale_mic()` and `scale_*_mic()` functions.
  * Comparisons of MIC values are now more strict. For example, `>32` is higher than (and never equal to) `32`. Thus, `as.mic(">32") == as.mic(32)` now returns `FALSE`, and `as.mic(">32") > as.mic(32)` now returns `TRUE`.
  * Sorting of MIC values (using `sort()`) was fixed in the same manner; `<0.001` now gets sorted before `0.001`, and `>0.001` gets sorted after `0.001`.
  * Intermediate log2 levels used for MIC plotting are now more common values instead of following a strict dilution range
  * `is.mic()` now returns a vector of `TRUE`/`FALSE` if the input is a `data.frame`, just like `as.sir()`
* `eucast_rules()` now has an argument `overwrite` (default: `FALSE`) to indicate whether non-`NA` values should be overwritten
* Disks of 0 to 5 mm are now allowed, the newly allowed range for disk diffusion (`as.disk()`) is now between 0 and 50 mm
* Updated `italicise_taxonomy()` to support HTML output
* `custom_eucast_rules()` now supports multiple antimicrobials and antimicrobial groups to be affected by a single rule
* `mo_info()` now contains an extra element `rank` and `group_members` (with the contents of the new `mo_group_members()` function)
* Updated all ATC codes from WHOCC
* Updated all antimicrobial DDDs from WHOCC
* Fix for using a manual value for `mo_transform` in `antibiogram()`
* Fixed a bug for when `antibiogram()` returns an empty data set
* Argument `only_sir_columns` now defaults to `TRUE` if any column of a data set contains a class 'sir' (functions `eucast_rules()`, `key_antimicrobials()`, `mdro()`, etc.)
* Added Sensititre codes for animals, antimicrobials and microorganisms
* Fix for mapping 'high level' antimicrobials in `as.ab()` (amphotericin B-high, gentamicin-high, kanamycin-high, streptomycin-high, tobramycin-high)
* Improved overall algorithm of `as.ab()` for better performance and accuracy, including the new function `as_reset_session()` to remove earlier coercions.
* Improved overall algorithm of `as.mo()` for better performance and accuracy, specifically:
  * More weight is given to genus and species combinations in cases where the subspecies is miswritten, so that the result will be the correct genus and species
  * Genera from the World Health Organization's (WHO) Priority Pathogen List now have the highest prevalence
* Fixed a bug for `sir_confidence_interval()` when there are no isolates available
* Updated the prevalence calculation to include genera from the World Health Organization's (WHO) Priority Pathogen List
* Improved algorithm of `first_isolate()` when using the phenotype-based method, to prioritise records with the highest availability of SIR values
* `scale_y_percent()` can now cope with ranges outside the 0-100% range
* MDRO determination (using `mdro()`)
  * The Verbose Mode (`verbose = TRUE`) now includes the guideline name
  * Implemented the new Dutch national MDRO guideline (SRI-richtlijn BRMO, Nov 2024)
  * Added arguments `esbl`, `carbapenemase`, `mecA`, `mecC`, `vanA`, `vanB` to denote column names or logical values indicating presence of these genes (or production of their proteins)
  * Added upport for antimicrobial selectors to use as as a custom rule (`custom_mdro_guideline()`)
* Added console colours support of `sir` class for Positron

### Other
* New website domain: <https://amr-for-r.org>! The old domain will remain to work.
* Added Dr. Larisse Bolton and Aislinn Cook as contributors for their fantastic implementation of WISCA in a mathematically solid way
* Added Matthew Saab, Dr. Jordan Stull, and Prof. Javier Sanchez as contributors for their tremendous input on veterinary breakpoints and interpretations
* Added Prof. Kathryn Holt, Dr. Jane Hawkey, and Dr. Natacha Couto as contributors for their many suggestions, ideas and bugfixes
* Greatly improved `vctrs` integration, a Tidyverse package working in the background for many Tidyverse functions. For users, this means that functions such as `dplyr`'s `bind_rows()`, `rowwise()` and `c_across()` are now supported for e.g. columns of class `mic`. Despite this, this `AMR` package is still zero-dependent on any other package, including `dplyr` and `vctrs`.
* Greatly updated and expanded documentation
* Stopped support for SAS (`.xpt`) files, since their file structure and extremely inefficient and requires more disk space than GitHub allows in a single commit.

## Older Versions

This changelog only contains changes from AMR v3.0 (June 2025) and later.

* For prior v2 versions, please see [our v2 archive](https://github.com/msberends/AMR/blob/v2.1.1/NEWS.md).
* For prior v1 versions, please see [our v1 archive](https://github.com/msberends/AMR/blob/v1.8.2/NEWS.md).
