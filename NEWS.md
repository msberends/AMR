# AMR 2.1.1.9083

*(this beta version will eventually become v3.0. We're happy to reach a new major milestone soon, which will be all about the new One Health support! Install this beta using [the instructions here](https://msberends.github.io/AMR/#latest-development-version).)*

#### A New Milestone: AMR v3.0 with One Health Support (= Human + Veterinary + Environmental)
This package now supports not only tools for AMR data analysis in clinical settings, but also for veterinary and environmental microbiology. This was made possible through a collaboration with the [University of Prince Edward Island's Atlantic Veterinary College](https://www.upei.ca/avc), Canada. To celebrate this great improvement of the package, we also updated the package logo to reflect this change.

## Breaking
* Removed all functions and references that used the deprecated `rsi` class, which were all replaced with their `sir` equivalents over a year ago

## New
* **One Health implementation**
  * Function `as.sir()` now has extensive support for veterinary breakpoints from CLSI. Use `breakpoint_type = "animal"` and set the `host` argument to a variable that contains animal species names.
  * The CLSI VET09 guideline has been implemented to address cases where veterinary breakpoints are missing (only applies when `guideline` is set to CLSI)
  * The `clinical_breakpoints` data set contains all these breakpoints, and can be downloaded on our [download page](https://msberends.github.io/AMR/articles/datasets.html).
  * The `antibiotics` data set contains all veterinary antibiotics, such as pradofloxacin and enrofloxacin. All WHOCC codes for veterinary use have been added as well.
  * `ab_atc()` now supports ATC codes of veterinary antibiotics (that all start with "Q")
  * `ab_url()` now supports retrieving the WHOCC url of their ATCvet pages
* **Major update to fungal taxonomy and tools for mycologists**
  * MycoBank has now been integrated as the primary taxonomic source for fungi. The `microorganisms` data set has been enriched with new columns (`mycobank`, `mycobank_parent`, and `mycobank_renamed_to`) that provide detailed information for fungal species.
  * A remarkable addition of over 20,000 new fungal records
  * New function `mo_mycobank()` to retrieve the MycoBank record number, analogous to existing functions such as `mo_lpsn()` and `mo_gbif()`.
  * The `as.mo()` function and all `mo_*()` functions now include an `only_fungi` argument, allowing users to restrict results solely to fungal species. This ensures fungi are prioritised over bacteria during microorganism identification. This can also be set globally with the new `AMR_only_fungi` option.
  * Also updated other kingdoms, welcoming a total of 2,149 new records from 2023 and 927 from 2024.
* **Updated clinical breakpoints**
  * EUCAST 2024 and CLSI 2024 are now supported, by adding all of their over 4,000 new clinical breakpoints to the `clinical_breakpoints` data set for usage in `as.sir()`. EUCAST 2024 is now the new default guideline for all MIC and disk diffusion interpretations.
  * `as.sir()` now brings additional factor levels: "NI" for non-interpretable and "SDD" for susceptible dose-dependent. Currently, the `clinical_breakpoints` data set contains 24 breakpoints that can return the value "SDD" instead of "I".
* **New forms for MIC plotting and transforming**
  * New function group `scale_*_mic()`, namely: `scale_x_mic()`, `scale_y_mic()`, `scale_colour_mic()` and `scale_fill_mic()`. They are advanced ggplot2 extensions to allow easy plotting of MIC values. They allow for manual range definition and plotting missing intermediate log2 levels.
  * New function `rescale_mic()`, which allows users to rescale MIC values to a manually set range. This is the powerhouse behind the `scale_*_mic()` functions, but it can be used independently to, for instance, compare equality in MIC distributions by rescaling them to the same range first.
* **Other**
  * New function `mo_group_members()` to retrieve the member microorganisms of a microorganism group. For example, `mo_group_members("Strep group C")` returns a vector of all microorganisms that belong to that group.

## Changed
* SIR interpretation
  * It is now possible to use column names for argument `ab`, `mo`, and `uti`: `as.sir(..., ab = "column1", mo = "column2", uti = "column3")`. This greatly improves the flexibility for users.
  * Users can now set their own criteria (using regular expressions) as to what should be considered S, I, R, SDD, and NI.
  * To get quantitative values, `as.double()` on a `sir` object will return 1 for S, 2 for SDD/I, and 3 for R (NI will become `NA`). Other functions using `sir` classes (e.g., `summary()`) are updated to reflect the change to contain NI and SDD.
* `antibiogram()` function
  * New argument `formatting_type` to set any of the 12 options for the formatting of all 'cells'. This defaults to `10`, changing the output of antibiograms to cells with `5% (15/300)` instead of the previous standard of just `5`.
  * For this reason, `add_total_n` is now `FALSE` at default since the denominators are added to the cells
  * The `ab_transform` argument now defaults to `"name"`, displaying antibiotic column names instead of codes
* `antibiotics` data set
  * Added "clindamycin inducible screening" as `CLI1`. Since clindamycin is a lincosamide, the antibiotic selector `lincosamides()` now contains the argument `only_treatable = TRUE` (similar to other antibiotic selectors that contain non-treatable drugs)
  * Added Amorolfine (`AMO`, D01AE16), which is now also part of the `antifungals()` selector
* Antibiotic selectors
  * Added selectors `nitrofurans()` and `rifamycins()`
  * When using antibiotic selectors such as `aminoglycosides()` that exclude non-treatable drugs like gentamicin-high, the function now always returns a warning that these can be included using `only_treatable = FALSE`
* MICs
  * Added as valid levels: 4096, 6 powers of 0.0625, and 5 powers of 192 (192, 384, 576, 768, 960)
  * Added new argument `keep_operators` to `as.mic()`. This can be `"all"` (default), `"none"`, or `"edges"`. This argument is also available in the new `rescale_mic()` and `scale_*_mic()` functions.
  * Comparisons of MIC values are now more strict. For example, `>32` is higher than (and never equal to) `32`. Thus, `as.mic(">32") == as.mic(32)` now returns `FALSE`, and `as.mic(">32") > as.mic(32)` now returns `TRUE`.
  * Sorting of MIC values (using `sort()`) was fixed in the same manner; `<0.001` now gets sorted before `0.001`, and `>0.001` gets sorted after `0.001`.
  * Intermediate log2 levels used for MIC plotting are now more common values instead of following a strict dilution range
* Updated `italicise_taxonomy()` to support HTML output
* `custom_eucast_rules()` now supports multiple antibiotics and antibiotic groups to be affected by a single rule
* `mo_info()` now contains an extra element `rank` and `group_members` (with the contents of the new `mo_group_members()` function)
* Greatly improved `vctrs` integration, a Tidyverse package working in the background for many Tidyverse functions. For users, this means that functions such as `dplyr`'s `bind_rows()`, `rowwise()` and `c_across()` are now supported for e.g. columns of class `mic`. Despite this, this `AMR` package is still zero-dependent on any other package, including `dplyr` and `vctrs`.
* Updated all ATC codes from WHOCC
* Updated all antibiotic DDDs from WHOCC
* Fix for using a manual value for `mo_transform` in `antibiogram()`
* Fix for mapping 'high level' antibiotics in `as.ab()` (amphotericin B-high, gentamicin-high, kanamycin-high, streptomycin-high, tobramycin-high)
* Improved overall algorithm of `as.ab()` for better performance and accuracy
* Improved overall algorithm of `as.mo()` for better performance and accuracy. Specifically:
  * More weight is given to genus and species combinations in cases where the subspecies is miswritten, so that the result will be the correct genus and species
  * Genera from the World Health Organization's (WHO) Priority Pathogen List now have the highest prevalence
* Fixed a bug for when `antibiogram()` returns an empty data set
* Updated the prevalence calculation to include genera from the World Health Organization's (WHO) Priority Pathogen List
* Improved algorithm of `first_isolate()` when using the phenotype-based method, to prioritise records with the highest availability of SIR values

## Other
* Greatly updated and expanded documentation
* Added Larisse Bolton, Jordan Stull, Matthew Saab, and Javier Sanchez as contributors, to thank them for their valuable input
* Stopped support for SAS (`.xpt`) files, since their file structure and extremely inefficient and requires more disk space than GitHub allows in a single commit.

## Older Versions

This changelog only contains changes from AMR v3.0 (October 2024) and later.

* For prior v2 versions, please see [our v2 archive](https://github.com/msberends/AMR/blob/v2.1.1/NEWS.md).
* For prior v1 versions, please see [our v1 archive](https://github.com/msberends/AMR/blob/v1.8.2/NEWS.md).
