# Changelog

## AMR 3.0.1.9008

#### New

- Integration with the **tidymodels** framework to allow seamless use of
  SIR, MIC and disk data in modelling pipelines via `recipes`
  - [`step_mic_log2()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    to transform `<mic>` columns with log2, and
    [`step_sir_numeric()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    to convert `<sir>` columns to numeric
  - New `tidyselect` helpers:
    - [`all_sir()`](https://amr-for-r.org/reference/amr-tidymodels.md),
      [`all_sir_predictors()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    - [`all_mic()`](https://amr-for-r.org/reference/amr-tidymodels.md),
      [`all_mic_predictors()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    - [`all_disk()`](https://amr-for-r.org/reference/amr-tidymodels.md),
      [`all_disk_predictors()`](https://amr-for-r.org/reference/amr-tidymodels.md)
- Data set `esbl_isolates` to practise with AMR modelling
- AMR selectors
  [`phosphonics()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  and
  [`spiropyrimidinetriones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
- [`ab_group()`](https://amr-for-r.org/reference/ab_property.md) gained
  an argument `all_groups` to return all groups the antimicrobial drug
  is in ([\#246](https://github.com/msberends/AMR/issues/246))

#### Changed

- Fixed a bug in
  [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) for
  when no antimicrobials are set
- Added taniborbactam (`TAN`) and cefepime/taniborbactam (`FTA`) to the
  `antimicrobials` data set
- Fixed a bug in [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)
  where for numeric input the arguments `S`, `i`, and `R` would not be
  considered ([\#244](https://github.com/msberends/AMR/issues/244))
- Added explaining message to
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md) when
  interpreting numeric values (e.g., 1 for S, 2 for I, 3 for R)
  ([\#244](https://github.com/msberends/AMR/issues/244))
- Updated handling of capped MIC values (`<`, `<=`, `>`, `>=`) in
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md) in the
  argument `capped_mic_handling`:
  ([\#243](https://github.com/msberends/AMR/issues/243))
  - Introduced four clearly defined options: `"none"`, `"conservative"`
    (default), `"standard"`, and `"lenient"`
  - Interpretation of capped MIC values now consistently returns `"NI"`
    (non-interpretable) when the true MIC could be at either side of a
    breakpoint, depending on the selected handling mode
  - This results in more reliable behaviour compared to previous
    versions for capped MIC values
  - Removed the `"inverse"` option, which has now become redundant
- Fixed some foreign translations of antimicrobial drugs
- `antimicrobials$group` is now a `list` instead of a `character`, to
  contain any group the drug is in
  ([\#246](https://github.com/msberends/AMR/issues/246))
- [`ab_group()`](https://amr-for-r.org/reference/ab_property.md) now
  returns values consist with the AMR selectors
  ([\#246](https://github.com/msberends/AMR/issues/246))

## AMR 3.0.1

CRAN release: 2025-09-20

This is a bugfix release following the release of v3.0.0 in June 2025.

#### Changed

- Fixed bugs introduced by `ggplot2` v4.0.0
  ([\#236](https://github.com/msberends/AMR/issues/236))
  - MIC scale functions (such as
    [`scale_y_mic()`](https://amr-for-r.org/reference/plot.md)) will now
    be applied automatically when plotting values of class `mic`
  - SIR scale functions (such as
    [`scale_x_sir()`](https://amr-for-r.org/reference/plot.md)) will now
    be applied automatically when plotting values of class `sir`
- Fixed a bug in
  [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) for
  when no antimicrobials are set
- Fixed a bug in
  [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) to
  allow column names containing the `+` character
  ([\#222](https://github.com/msberends/AMR/issues/222))
- Fixed a bug in [`as.ab()`](https://amr-for-r.org/reference/as.ab.md)
  for antimicrobial codes with a number in it if they are preceded by a
  space
- Fixed a bug in
  [`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md)
  for using specific custom rules
- Fixed a bug in [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)
  to allow any tidyselect language
  ([\#220](https://github.com/msberends/AMR/issues/220))
- Fixed a bug in [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)
  to pick right breakpoint when `uti = FALSE`
  ([\#216](https://github.com/msberends/AMR/issues/216))
- Fixed a bug in
  [`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md) when
  using `combine_SI = FALSE`
  ([\#213](https://github.com/msberends/AMR/issues/213))
- Fixed a bug in [`mdro()`](https://amr-for-r.org/reference/mdro.md) to
  make sure all genes specified in arguments are acknowledged
- Fixed a bug the `antimicrobials` data set to remove statins
  ([\#229](https://github.com/msberends/AMR/issues/229))
- Fixed a bug the `microorganisms` data set for MycoBank IDs and
  synonyms ([\#233](https://github.com/msberends/AMR/issues/233))
- Fixed ATC J01CR05 to map to piperacillin/tazobactam rather than
  piperacillin/sulbactam
  ([\#230](https://github.com/msberends/AMR/issues/230))
- Fixed skimmers (`skimr` package) of class `ab`, `sir`, and `disk`
  ([\#234](https://github.com/msberends/AMR/issues/234))
- Fixed all plotting to contain a separate colour for SDD (susceptible
  dose-dependent) ([\#223](https://github.com/msberends/AMR/issues/223))
- Fixed some specific Dutch translations for antimicrobials
- Added a warning to
  [`as.ab()`](https://amr-for-r.org/reference/as.ab.md) if input
  resembles antiviral codes or names
  ([\#232](https://github.com/msberends/AMR/issues/232))
- Added all reasons in verbose output of
  [`mdro()`](https://amr-for-r.org/reference/mdro.md)
  ([\#227](https://github.com/msberends/AMR/issues/227))
- Added `names` to
  [`age_groups()`](https://amr-for-r.org/reference/age_groups.md) so
  that custom names can be given
  ([\#215](https://github.com/msberends/AMR/issues/215))
- Added note to [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)
  to make it explicit when higher-level taxonomic breakpoints are used
  ([\#218](https://github.com/msberends/AMR/issues/218))
- Added antibiotic codes from the Comprehensive Antibiotic Resistance
  Database (CARD) to the `antimicrobials` data set
  ([\#225](https://github.com/msberends/AMR/issues/225))
- Updated Fosfomycin to be of antibiotic class Phosphonics
  ([\#225](https://github.com/msberends/AMR/issues/225))
- Updated [`random_mic()`](https://amr-for-r.org/reference/random.md)
  and [`random_disk()`](https://amr-for-r.org/reference/random.md) to
  set skewedness of the distribution and allow multiple microorganisms

## AMR 3.0.0

CRAN release: 2025-06-02

This package now supports not only tools for AMR data analysis in
clinical settings, but also for veterinary and environmental
microbiology. This was made possible through a collaboration with the
[University of Prince Edward Island’s Atlantic Veterinary
College](https://www.upei.ca/avc), Canada. To celebrate this great
improvement of the package, we also updated the package logo to reflect
this change.

#### Breaking

- Data set `antibiotics` has been renamed to `antimicrobials` as the
  data set contains more than just antibiotics. Using `antibiotics` will
  still work, but now returns a warning.
- Removed all functions and references that used the deprecated `rsi`
  class, which were all replaced with their `sir` equivalents over two
  years ago.
- Functions
  [`resistance_predict()`](https://amr-for-r.org/reference/resistance_predict.md)
  and
  [`sir_predict()`](https://amr-for-r.org/reference/resistance_predict.md)
  are now deprecated and will be removed in a future version. Use the
  `tidymodels` framework instead, for which we [wrote a basic
  introduction](https://amr-for-r.org/articles/AMR_with_tidymodels.html).

#### New

- **One Health implementation**
  - Function [`as.sir()`](https://amr-for-r.org/reference/as.sir.md) now
    has extensive support for veterinary breakpoints from CLSI. Use
    `breakpoint_type = "animal"` and set the `host` argument to a
    variable that contains animal species names.
  - The `clinical_breakpoints` data set contains all these breakpoints,
    and can be downloaded on our [download
    page](https://amr-for-r.org/articles/datasets.html).
  - The (new) `antimicrobials` data set contains all veterinary
    antimicrobials, such as pradofloxacin and enrofloxacin. All WHOCC
    codes for veterinary use have been added as well.
  - [`ab_atc()`](https://amr-for-r.org/reference/ab_property.md) now
    supports ATC codes of veterinary antimicrobials (that all start with
    “Q”)
  - [`ab_url()`](https://amr-for-r.org/reference/ab_property.md) now
    supports retrieving the WHOCC url of their ATCvet pages
- **Support for WISCA antibiograms**
  - The
    [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md)
    function now supports creating true Weighted-Incidence Syndromic
    Combination Antibiograms (WISCA), a powerful Bayesian method for
    estimating regimen coverage probabilities using pathogen incidence
    and antimicrobial susceptibility data. WISCA offers improved
    precision for syndrome-specific treatment, even in data sets with
    sparse data. A dedicated
    [`wisca()`](https://amr-for-r.org/reference/antibiogram.md) function
    is also available for easy usage.
- **More global coverage of languages**
  - Added full support for 8 new languages: Arabic, Bengali, Hindi,
    Indonesian, Korean, Swahili, Urdu, and Vietnamese. The `AMR` package
    is now available in 28 languages.
- **Major update to fungal taxonomy and tools for mycologists**
  - MycoBank has now been integrated as the primary taxonomic source for
    fungi. The `microorganisms` data set has been enriched with new
    columns (`mycobank`, `mycobank_parent`, and `mycobank_renamed_to`)
    that provide detailed information for fungal species.
  - A remarkable addition of over 20,000 new fungal records
  - New function
    [`mo_mycobank()`](https://amr-for-r.org/reference/mo_property.md) to
    retrieve the MycoBank record number, analogous to existing functions
    such as
    [`mo_lpsn()`](https://amr-for-r.org/reference/mo_property.md) and
    [`mo_gbif()`](https://amr-for-r.org/reference/mo_property.md).
  - The [`as.mo()`](https://amr-for-r.org/reference/as.mo.md) function
    and all `mo_*()` functions now include an `only_fungi` argument,
    allowing users to restrict results solely to fungal species. This
    ensures fungi are prioritised over bacteria during microorganism
    identification. This can also be set globally with the new
    `AMR_only_fungi` option.
  - Also updated other kingdoms, welcoming a total of 2,149 new records
    from 2023 and 927 from 2024.
- **Updated clinical breakpoints**
  - Breakpoint of 2024 and 2025 of both CLSI and EUCAST are now
    supported, by adding all of their over 10,000 new clinical
    breakpoints to the `clinical_breakpoints` data set for usage in
    [`as.sir()`](https://amr-for-r.org/reference/as.sir.md). EUCAST 2025
    is now the new default guideline for all MIC and disk diffusion
    interpretations.
  - Added all Expected Resistant Phenotypes from EUCAST (v1.2). The
    default `rules` for
    [`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md)
    are now: `c("breakpoints", "expected_phenotypes")`.
  - Updated the `intrinsic_resistant` data set, which is now based on
    EUCAST Expected Resistant Phenotypes v1.2
  - [`as.sir()`](https://amr-for-r.org/reference/as.sir.md) now brings
    additional factor levels: “NI” for non-interpretable and “SDD” for
    susceptible dose-dependent. Currently, the `clinical_breakpoints`
    data set contains 24 breakpoints that can return the value “SDD”
    instead of “I”.
  - EUCAST interpretive rules (using
    [`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md))
    are now available for EUCAST 12 (2022), 13 (2023), 14 (2024), and 15
    (2025).
  - EUCAST dosage tables (`dosage` data set) are now available for
    EUCAST 13 (2023), 14 (2024), and 15 (2025).
- **New advanced ggplot2 extensions for MIC and SIR plotting and
  transforming**
  - New function group `scale_*_mic()`, namely:
    [`scale_x_mic()`](https://amr-for-r.org/reference/plot.md),
    [`scale_y_mic()`](https://amr-for-r.org/reference/plot.md),
    [`scale_colour_mic()`](https://amr-for-r.org/reference/plot.md) and
    [`scale_fill_mic()`](https://amr-for-r.org/reference/plot.md). They
    allow easy plotting of MIC values. They allow for manual range
    definition and plotting missing intermediate log2 levels.
  - New function group `scale_*_sir()`, namely:
    [`scale_x_sir()`](https://amr-for-r.org/reference/plot.md),
    [`scale_colour_sir()`](https://amr-for-r.org/reference/plot.md) and
    [`scale_fill_sir()`](https://amr-for-r.org/reference/plot.md). They
    allow to plot the `sir` class, and translates into the system
    language at default. They also set colourblind-safe colours to the
    plots.
  - New function
    [`rescale_mic()`](https://amr-for-r.org/reference/as.mic.md), which
    allows users to rescale MIC values to a manually set range. This is
    the powerhouse behind the `scale_*_mic()` functions, but it can be
    used independently to, for instance, compare equality in MIC
    distributions by rescaling them to the same range first.
- **Support for Python**
  - While using R for the heavy lifting, [our ‘AMR’ Python
    Package](https://pypi.org/project/AMR/) was developed to run the AMR
    R package natively in Python. The Python package will always have
    the same version number as the R package, as it is built
    automatically with every code change.
- **Support for `tidymodels`**
  - All antimicrobial selectors (such as
    [`aminoglycosides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    and
    [`betalactams()`](https://amr-for-r.org/reference/antimicrobial_selectors.md))
    are now supported in `tidymodels` packages such as `recipe` and
    `parsnip`. See for more info [our
    tutorial](https://amr-for-r.org/articles/AMR_with_tidymodels.html)
    on using these AMR functions for predictive modelling.
- **Other**
  - New function
    [`top_n_microorganisms()`](https://amr-for-r.org/reference/top_n_microorganisms.md)
    to filter a data set to the top *n* of any taxonomic property, e.g.,
    filter to the top 3 species, filter to any species in the top 5
    genera, or filter to the top 3 species in each of the top 5 genera
  - New function
    [`mo_group_members()`](https://amr-for-r.org/reference/mo_property.md)
    to retrieve the member microorganisms of a microorganism group. For
    example, `mo_group_members("Strep group C")` returns a vector of all
    microorganisms that belong to that group.
  - New functions
    [`mic_p50()`](https://amr-for-r.org/reference/as.mic.md) and
    [`mic_p90()`](https://amr-for-r.org/reference/as.mic.md) to retrieve
    the 50th and 90th percentile of MIC values.

#### Changed

- SIR interpretation
  - Support for parallel computing to greatly improve speed using the
    `parallel` package (part of base R). Use
    `as.sir(your_data, parallel = TRUE)` to run SIR interpretation using
    multiple cores.
  - It is now possible to use column names for arguments `guideline`,
    `ab`, `mo`, and `uti`:
    `as.sir(..., ab = "column1", mo = "column2", uti = "column3")`. This
    greatly improves the flexibility for users.
  - Users can now set their own criteria (using regular expressions) as
    to what should be considered S, I, R, SDD, and NI.
  - To get quantitative values,
    [`as.double()`](https://rdrr.io/r/base/double.html) on a `sir`
    object will return 1 for S, 2 for SDD/I, and 3 for R (NI will become
    `NA`). Other functions using `sir` classes (e.g.,
    [`summary()`](https://rdrr.io/r/base/summary.html)) are updated to
    reflect the change to contain NI and SDD.
  - Following CLSI interpretation rules, values outside the
    log2-dilution range will be rounded upwards to the nearest
    log2-level before interpretation. Only if using a CLSI guideline.
  - Combined MIC values (e.g., from CLSI) are now supported
  - The argument `conserve_capped_values` in
    [`as.sir()`](https://amr-for-r.org/reference/as.sir.md) has been
    replaced with `capped_mic_handling`, which allows greater
    flexibility in handling capped MIC values (`<`, `<=`, `>`, `>=`).
    The four available options (`"standard"`, `"strict"`, `"relaxed"`,
    `"inverse"`) provide full control over whether these values should
    be interpreted conservatively or ignored. Using
    `conserve_capped_values` is now deprecated and returns a warning.
  - Added argument `info` to silence all console messages
- [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md)
  function
  - Argument `antibiotics` has been renamed to `antimicrobials`. Using
    `antibiotics` will still work, but now returns a warning.
  - Added argument `formatting_type` to set any of the 22 options for
    the formatting of all ‘cells’. This defaults to `18` for non-WISCA
    and `14` for WISCA, changing the output of antibiograms to cells
    with more info.
  - For this reason, `add_total_n` is now deprecated and `FALSE` at
    default since the denominators are added to the cells dependent on
    the `formatting_type` setting
  - The `ab_transform` argument now defaults to `"name"`, displaying
    antibiotic column names instead of codes
- Antimicrobial selectors (previously: *antibiotic selectors*)
  - ‘Antibiotic selectors’ are now called ‘antimicrobial selectors’
    since their scope is broader than just antibiotics. All
    documentation have been updated, and
    [`ab_class()`](https://amr-for-r.org/reference/AMR-deprecated.md)
    and
    [`ab_selector()`](https://amr-for-r.org/reference/AMR-deprecated.md)
    have been replaced with
    [`amr_class()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    and
    [`amr_selector()`](https://amr-for-r.org/reference/antimicrobial_selectors.md).
    The old functions are now deprecated and will be removed in a future
    version.
  - Added selectors
    [`isoxazolylpenicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    [`monobactams()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    [`nitrofurans()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    [`phenicols()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    [`rifamycins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    and
    [`sulfonamides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  - When using antimicrobial selectors that exclude non-treatable drugs
    (such as gentamicin-high when using
    [`aminoglycosides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)),
    the function now always returns a warning that these can be included
    using `only_treatable = FALSE`
  - Added a new argument `return_all` to all selectors, which defaults
    to `TRUE` to include any match. With `FALSE`, the old behaviour,
    only the first hit for each unique antimicrobial is returned.
  - All selectors can now be run as a separate command to retrieve a
    vector of all possible antimicrobials that the selector can select
  - The selectors
    [`lincosamides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    and
    [`macrolides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    do not overlap anymore - each antibiotic is now classified as either
    of these and not both
  - Fixed selector
    [`fluoroquinolones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    which now really only selects second-generation quinolones and up
    (first-generation quinolones do not contain a fluorine group)
- `antimicrobials` data set
  - Added agents used for screening, with an ID all ending with `-S`:
    benzylpenicillin screening test (`PEN-S`), beta-lactamase screening
    test (`BLA-S`), cefotaxime screening test (`CTX-S`), clindamycin
    inducible screening test (`CLI-S`), nalidixic acid screening test
    (`NAL-S`), norfloxacin screening test (`NOR-S`), oxacillin screening
    test (`OXA-S`), pefloxacin screening test (`PEF-S`), and
    tetracycline screening test (`TCY-S`). The ID of cefoxitin screening
    was renamed from `FOX1` to `FOX-S`, while the old code remains to
    work.
  - For this reason, the antimicrobial selectors
    [`cephalosporins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    [`cephalosporins_3rd()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    [`lincosamides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    [`isoxazolylpenicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    [`quinolones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    [`fluoroquinolones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
    and
    [`tetracyclines()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    now contain the argument `only_treatable = TRUE` (similar to other
    antimicrobial selectors that contain non-treatable drugs)
  - Added amorolfine (`AMO`, D01AE16), an antimycotic, which is now also
    part of the
    [`antifungals()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    selector
  - Added cefepime/enmetazobactam (`FPE`), a 4th gen cephalosporin
  - Added tigemonam (`TNM`), a monobactam
  - Added bleomycin (`BLM`), a glycopeptide
  - Added efflux (`EFF`), to allow mapping to AMRFinderPlus
  - Updated all ATC codes, trade names, and DDDs
- MICs
  - Added as valid levels: 4096, 6 powers of 0.0625, and 5 powers of 192
    (192, 384, 576, 768, 960)
  - Fixed a bug in
    [`as.mic()`](https://amr-for-r.org/reference/as.mic.md) that failed
    translation of scientifically formatted numbers
  - Added new argument `keep_operators` to
    [`as.mic()`](https://amr-for-r.org/reference/as.mic.md). This can be
    `"all"` (default), `"none"`, or `"edges"`. This argument is also
    available in the new
    [`rescale_mic()`](https://amr-for-r.org/reference/as.mic.md) and
    `scale_*_mic()` functions.
  - Comparisons of MIC values are now more strict. For example, `>32` is
    higher than (and never equal to) `32`. Thus,
    `as.mic(">32") == as.mic(32)` now returns `FALSE`, and
    `as.mic(">32") > as.mic(32)` now returns `TRUE`.
  - Sorting of MIC values (using
    [`sort()`](https://rdrr.io/r/base/sort.html)) was fixed in the same
    manner; `<0.001` now gets sorted before `0.001`, and `>0.001` gets
    sorted after `0.001`.
  - Intermediate log2 levels used for MIC plotting are now more common
    values instead of following a strict dilution range
  - [`is.mic()`](https://amr-for-r.org/reference/as.mic.md) now returns
    a vector of `TRUE`/`FALSE` if the input is a `data.frame`, just like
    [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)
- [`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md)
  now has an argument `overwrite` (default: `FALSE`) to indicate whether
  non-`NA` values should be overwritten
- Disks of 0 to 5 mm are now allowed, the newly allowed range for disk
  diffusion ([`as.disk()`](https://amr-for-r.org/reference/as.disk.md))
  is now between 0 and 50 mm
- Updated
  [`italicise_taxonomy()`](https://amr-for-r.org/reference/italicise_taxonomy.md)
  to support HTML output
- [`custom_eucast_rules()`](https://amr-for-r.org/reference/custom_eucast_rules.md)
  now supports multiple antimicrobials and antimicrobial groups to be
  affected by a single rule
- [`mo_info()`](https://amr-for-r.org/reference/mo_property.md) now
  contains an extra element `rank` and `group_members` (with the
  contents of the new
  [`mo_group_members()`](https://amr-for-r.org/reference/mo_property.md)
  function)
- Updated all ATC codes from WHOCC
- Updated all antimicrobial DDDs from WHOCC
- Fix for using a manual value for `mo_transform` in
  [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md)
- Fixed a bug for when
  [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md)
  returns an empty data set
- Argument `only_sir_columns` now defaults to `TRUE` if any column of a
  data set contains a class ‘sir’ (functions
  [`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md),
  [`key_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md),
  [`mdro()`](https://amr-for-r.org/reference/mdro.md), etc.)
- Added Sensititre codes for animals, antimicrobials and microorganisms
- Fix for mapping ‘high level’ antimicrobials in
  [`as.ab()`](https://amr-for-r.org/reference/as.ab.md) (amphotericin
  B-high, gentamicin-high, kanamycin-high, streptomycin-high,
  tobramycin-high)
- Improved overall algorithm of
  [`as.ab()`](https://amr-for-r.org/reference/as.ab.md) for better
  performance and accuracy, including the new function
  `as_reset_session()` to remove earlier coercions.
- Improved overall algorithm of
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md) for better
  performance and accuracy, specifically:
  - More weight is given to genus and species combinations in cases
    where the subspecies is miswritten, so that the result will be the
    correct genus and species
  - Genera from the World Health Organization’s (WHO) Priority Pathogen
    List now have the highest prevalence
- Fixed a bug for
  [`sir_confidence_interval()`](https://amr-for-r.org/reference/proportion.md)
  when there are no isolates available
- Updated the prevalence calculation to include genera from the World
  Health Organization’s (WHO) Priority Pathogen List
- Improved algorithm of
  [`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md)
  when using the phenotype-based method, to prioritise records with the
  highest availability of SIR values
- [`scale_y_percent()`](https://amr-for-r.org/reference/plot.md) can now
  cope with ranges outside the 0-100% range
- MDRO determination (using
  [`mdro()`](https://amr-for-r.org/reference/mdro.md))
  - The Verbose Mode (`verbose = TRUE`) now includes the guideline name
  - Implemented the new Dutch national MDRO guideline (SRI-richtlijn
    BRMO, Nov 2024)
  - Added arguments `esbl`, `carbapenemase`, `mecA`, `mecC`, `vanA`,
    `vanB` to denote column names or logical values indicating presence
    of these genes (or production of their proteins)
  - Added upport for antimicrobial selectors to use as as a custom rule
    ([`custom_mdro_guideline()`](https://amr-for-r.org/reference/custom_mdro_guideline.md))
- Added console colours support of `sir` class for Positron

#### Other

- New website domain: <https://amr-for-r.org>! The old domain will
  remain to work.
- Added Dr. Larisse Bolton and Aislinn Cook as contributors for their
  fantastic implementation of WISCA in a mathematically solid way
- Added Matthew Saab, Dr. Jordan Stull, and Prof. Javier Sanchez as
  contributors for their tremendous input on veterinary breakpoints and
  interpretations
- Added Prof. Kathryn Holt, Dr. Jane Hawkey, and Dr. Natacha Couto as
  contributors for their many suggestions, ideas and bugfixes
- Greatly improved `vctrs` integration, a Tidyverse package working in
  the background for many Tidyverse functions. For users, this means
  that functions such as `dplyr`’s
  [`bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html),
  [`rowwise()`](https://dplyr.tidyverse.org/reference/rowwise.html) and
  [`c_across()`](https://dplyr.tidyverse.org/reference/c_across.html)
  are now supported for e.g. columns of class `mic`. Despite this, this
  `AMR` package is still zero-dependent on any other package, including
  `dplyr` and `vctrs`.
- Greatly updated and expanded documentation
- Stopped support for SAS (`.xpt`) files, since their file structure and
  extremely inefficient and requires more disk space than GitHub allows
  in a single commit.

### Older Versions

This changelog only contains changes from AMR v3.0 (June 2025) and
later.

- For prior v2 versions, please see [our v2
  archive](https://github.com/msberends/AMR/blob/v2.1.1/NEWS.md).
- For prior v1 versions, please see [our v1
  archive](https://github.com/msberends/AMR/blob/v1.8.2/NEWS.md).
