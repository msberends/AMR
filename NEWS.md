# AMR 3.1.0

This substantial release makes most core functions considerably faster, adds EUCAST 2026 and CLSI 2026 breakpoints, adopts the new prokaryotic kingdoms of Göker and Oren (2024), introduces WT/NWT interpretation for ECOFFs, unifies the EUCAST and CLSI interpretive rules, and supports `tidymodels` pipelines. As it contains several breaking changes, please read the *Breaking* section below before updating.

### Breaking
* The former *kingdoms* Bacteria and Archaea are now each divided into four kingdoms, with new top-level *domains* Bacteria and Archaea (Göker and Oren, 2024, DOI: 10.1099/ijsem.0.006242). Following this, a new `domain` column was added to the `microorganisms` data set, and `mo_kingdom()` now returns the formal kingdom (e.g., `"Pseudomonadati"` instead of `"Bacteria"`). Use `mo_domain()` for the old behaviour. For non-prokaryotic kingdoms (Fungi, Protozoa, etc.), `kingdom` and `domain` are identical.

  | `mo_kingdom()` < 3.1.0 | `mo_kingdom()` now  | `mo_domain()` (unchanged) |
  |------------------------|---------------------|---------------------------|
  | Bacteria               | Bacillati           | Bacteria                  |
  |                        | Fusobacteriati      | Bacteria                  |
  |                        | Pseudomonadati      | Bacteria                  |
  |                        | Thermotogati        | Bacteria                  |
  |                        |                     |                           |
  | Archaea                | Methanobacteriati   | Archaea                   |
  |                        | Nanobdellati        | Archaea                   |
  |                        | Promethearchaeati   | Archaea                   |
  |                        | Thermoproteati      | Archaea                   |
  |                        |                     |                           |
  | Fungi                  | Fungi               | Fungi                     |
  |                        |                     |                           |
  | Protozoa               | Protozoa            | Protozoa                  |

  Previously, `mo_domain()` was an alias of `mo_kingdom()`; it is now a distinct function. Its output is therefore unchanged, while `mo_kingdom()` now returns the new formal kingdom. A one-time note per session explains this change when bacterial or archaeal records are queried.
* Parallel computing in `as.sir()`, `antibiogram()`, and `wisca()` now uses the `future` and `future.apply` packages, which is faster. A non-sequential plan (e.g., `future::plan(future::multisession)`) must be active before using `parallel = TRUE` (#281).
* The options of `capped_mic_handling` in `as.sir()` have been redefined as `"none"`, `"conservative"`, `"standard"` (default), and `"lenient"`. The option `"inverse"` is no longer available (#243).
* `as.sir()` now returns WT/NWT (Wildtype/Non-wildtype) instead of S/R when using ECOFF-based interpretation with `breakpoint_type = "ECOFF"`. WT/NWT results are fully supported in all resistance and susceptibility functions and in all plots (#254).
* Aggregated analyses now follow the guideline for the 'I' category, set with the new `guideline` argument or the `AMR_guideline` option: EUCAST considers 'I' susceptible and CLSI considers 'I' resistant, while SDD and WT are always considered susceptible, and NWT and NS always resistant (see section *Interpretation of SIR* in `?as.sir`). This applies to `count_df()`, `proportion_df()`, `sir_df()`, `ggplot_sir()`, `geom_sir()`, and `labels_sir_count()`, which with CLSI now return `"S"` vs. `"IR"` instead of `"SI"` vs. `"R"`, and to `antibiogram()`, `wisca()`, and `format()` on the output of `bug_drug_combinations()`. Previously, these functions always considered 'I' susceptible, as with EUCAST. Use `guideline = "EUCAST"` for the old behaviour.
* The `combine_SI` argument of these functions is deprecated in favour of `guideline`. While set, it is followed literally and takes precedence over `guideline`. Use `guideline = NULL` instead of `combine_SI = FALSE` to keep S, I, and R separate.

### New
* **Updated clinical breakpoints**
  * EUCAST 2026 and CLSI 2026 are now supported, by adding over 5,700 new breakpoints to the `clinical_breakpoints` data set. EUCAST 2026 is now the default guideline for all MIC and disk diffusion interpretations.
  * Human EUCAST breakpoints from 2019 (v9.0) onwards are now retrieved directly from the EUCAST Excel files, instead of relying on the WHONET repository. As a result:
    * `clinical_breakpoints` now contains blocking rows (`NA` breakpoints, with the reason in `note`) where EUCAST gives no breakpoint for an organism, so that `as.sir()` no longer falls back to the breakpoint of a broader taxon. This includes infections other than UTI where EUCAST gives breakpoints in brackets or no breakpoint, but does give a breakpoint for UTI, e.g., aminoglycosides for *Acinetobacter* spp.
    * Restrictions in EUCAST notes are applied, e.g., zone diameter breakpoints for *E. coli* only
    * EUCAST PK-PD (non-species related) breakpoints of 2019-2023 and *H. parainfluenzae* MIC breakpoints are now included
* **Unified interpretive rules**
  * New function `interpretive_rules()` applies the interpretive rules of both EUCAST and CLSI. `eucast_rules()` is now a wrapper around this function (#235, #259).
  * Species named in a rule now also match all other names of the same accepted taxon (e.g., 'Clostridium ramosum' for *Thomasclavelia ramosa*). This also updated the `intrinsic_resistant` data set.
* **Support for `tidymodels` pipelines**
  * New recipe steps `step_mic_log2()` and `step_sir_numeric()` to use MIC and SIR data in modelling pipelines
  * New column selectors `all_sir()`, `all_mic()`, and `all_disk()`
  * New `esbl_isolates` data set for practising AMR modelling
* **Morphology of bacteria**
  * New `morphology` column in the `microorganisms` data set, with the corresponding function `mo_morphology()` to return the cell shape of bacteria. Data are sourced from BacDive; values prefixed with "likely" are extrapolated from the consensus at genus level.
  * New argument `add_morphology` in `mo_gramstain()` to return combined results, such as `"Gram-negative rods"`
* **Antimicrobials**
  * New function `sir_combine()` to combine SIR results of multiple antimicrobials into one variable, for combination therapy (`method = "most_susceptible"`) or for resistance within a group (`method = "most_resistant"`) (#321)
  * New function `add_custom_antimicrobial_synonyms()` to add local names, such as trade names, to existing antimicrobials (#307)
  * New antimicrobial selectors `ionophores()`, `peptides()`, `phosphonics()`, and `spiropyrimidinetriones()`
  * New antimicrobials: cefepime/taniborbactam (`FTA`), ceftibuten/avibactam (`CTA`), clorobiocin (`CLB`), kasugamycin (`KAS`), ostreogrycin (`OST`), sulbactam/durlobactam (`SUD1`, #308), taniborbactam (`TAN`), thiostrepton (`THS`), xeruborbactam (`XER`), and zorbamycin (`ZOR`)
* **Other**
  * New function `wisca_plot()` to assess the distributions of susceptibility and incidence from the Monte Carlo simulations of WISCA
  * New function `amr_course()` to download and unpack course or webinar materials from GitHub in one call
  * New typed missing value constants `NA_ab_` and `NA_mo_`, for use in pipelines that need missing values of a specific class
  * A registry of all MO codes since v2.0.0 now prevents existing codes from being reused or mixed up. Codes of earlier releases are translated to the current code of the same taxon.
  * New column `source` in the `microorganisms.codes` data set, stating the source of a code and the name it refers to, e.g., `"WHONET (Enterobacter aerogenes)"`

### Changed
* SIR interpretation (`as.sir()`)
  * New argument `enforce_method` (`"auto"`, `"mic"`, or `"disk"`) to force the interpretation method when S3 class information is lost, e.g., when called from Python (#291)
  * When input consists of whole numbers only, disk zones are now assumed rather than MIC values (#291)
  * Capped MICs (e.g., `<=0.5`) now return `I` instead of `S` or `NI` if there is no S category, i.e., in case of an off-scale S breakpoint such as EUCAST's 0.001 mg/L (#295)
  * Capped MICs are no longer forced to `S` or `R` if there is no breakpoint for the drug-bug combination
  * Fixed the documented default of `capped_mic_handling`, which is `"standard"`
  * On data frames, already-converted SIR columns are no longer dropped on a re-run (#278)
  * Metadata columns (e.g., `patient`, `ward`) are no longer misidentified as antimicrobial columns
  * `info = FALSE` now suppresses all messages, including those for columns without breakpoints
  * Combined MIC/SIR input values (e.g., `"<= 0.002; S"` or `"S; 0.002"`) are now parsed correctly (#252)
  * Setting a year-qualified guideline, such as `options(AMR_guideline = "EUCAST 2012")`, no longer causes errors or silently wrong results in `interpretive_rules()`, `resistance()`, `susceptibility()`, `count_resistant()`, `count_susceptible()`, and the SIR plotting and printing functions (#298)
  * With an unknown `uti` (`NA`, the default), non-UTI breakpoints are now consistently preferred over UTI breakpoints. Previously, this preference was only applied together with its note, which is shown once. UTI breakpoints are only used if there are no non-UTI breakpoints, which is now noted in `sir_interpretation_history()` for every interpretation. A blocking row of a broader taxon does not override a UTI breakpoint of the organism itself, e.g., Enterobacterales and nitroxoline for *E. coli*.
* `clinical_breakpoints` data set
  * Fixed EUCAST breakpoint rules for *Aerococcus* (v12 onwards), which were never applied, and for *Pasteurella*, which were partly filed under *Neisseria meningitidis* (#290)
  * Fixed streptomycin breakpoints imported from WHONET, which were coded as streptoduocin (`STR`) instead of streptomycin (`STR1`)
  * Fixed CLSI breakpoints for sulbactam/durlobactam, which were coded as sulfadimethoxine (`SUD`), and added the missing EUCAST 2011 breakpoints for roxithromycin (#308)
  * Fixed some EUCAST Expert Rules, mostly on *S. pneumoniae*
  * Zone diameters for which the guideline has no S category (coded by WHONET as an S breakpoint of 80 or 99 mm) are now blocking rows for EUCAST, e.g. cefuroxime oral in *H. influenzae* in EUCAST 2011, and keep only the R breakpoint for CLSI, e.g. fluconazole in *C. glabrata*. Previously, these S breakpoints were stored as 50 mm.
  * Breakpoints and ECOFFs for *Salmonella* serovars Typhi, Typhimurium, and Paratyphi, and ECOFFs for the *Bacteroides fragilis* complex, are now assigned to these taxa instead of to *S. enterica* subsp. *enterica* and *B. fragilis*, and ECOFFs of amphotericin B for *Neocosmospora solani* were added
* Microorganism coercion (`as.mo()`) and properties (`mo_*()`)
  * Considerably faster through one consistent session cache design, which is cleared when lookup tables change, e.g., after `add_custom_microorganisms()` (#313)
  * `mo_*()` functions now retrieve properties of valid `mo` input directly, and only for unique values
  * Subsetting, combining, and binding rows of `mo` vectors (e.g., in `rbind()`) is now much faster, as codes are checked against a hashed set instead of the whole `microorganisms` data set
  * The mapping of outdated taxonomic names to their current names is now remembered per code
  * Previously coerced input is now recognised again, e.g., for an unknown species of a known genus, instead of being matched again on every call (#313)
  * `mo_uncertainties()` now also lists uncertainties of previously coerced input
  * Translated names in the set language are now recognised, such as the output of `mo_name(x, language = "fr")`, and cached per language
  * Input of the form `"X complex"` now falls back to `"X"` when the complex is not a distinct taxon, preventing `NA` for valid clinical descriptions such as `"Proteus vulgaris complex"`, while a complex named after several of its species returns the complex of one of them, e.g., `"Acinetobacter calcoaceticus-baumannii complex"` returns *A. baumannii* complex (#287)
  * Abbreviated genus input now ranks candidates with an exactly matching species epithet above more prevalent organisms, so that e.g. `"S. apiospermum"` returns *Scedosporium apiospermum* instead of a *Staphylococcus* species (#288)
  * Abbreviated genus input for species with subspecies (e.g., `"P. ovale"`) now returns the species instead of a more prevalent organism, while explicit subspecies (e.g., `"P. ovale curtisi"`) are preserved (#288)
  * *Salmonella* serovars written with the species, such as `"Salmonella enterica serovar Typhi"`, no longer return *Salmonella bongori* or another subspecies
  * `as.mo()` and `mo_current()` now find the current name of synonyms with multiple sources (e.g., *Chaetomium abuense*)
  * `mo_taxonomy()` and `mo_info()` now contain an extra element `domain`
* Own microorganism codes (`set_mo_source()` and `get_mo_source()`)
  * Fixed `get_mo_source()` removing the `mo_source` file at the default location in interactive sessions, when called with a `destination` that does not exist
  * Fixed `get_mo_source()` returning the `mo_source` already loaded instead of the one at the requested `destination`
  * Fixed `get_mo_source()`, and therefore `as.mo()` and all `mo_*()` functions, throwing an error in interactive sessions after the original file of the `mo_source` was moved or deleted. The `mo_source` file is now removed, as documented.
  * An updated original file now updates the `mo_source` file at its own `destination` instead of at the default location
* `microorganisms` data set
  * Taxonomic update of all microorganisms, now based on COL 2026-04-18 XR, LPSN and BacDive of October 2026, and MycoBank of January 2026
  * Names listed twice by a source now take the correct record (e.g., *Eggerthella lenta*, *Gordonia amarae*, and *Anisakis simplex* were synonyms without a current name)
  * Synonyms without a current name are now linked via LPSN or removed, with their codes retired. The *M. tuberculosis* complex is always kept.
  * Higher taxonomy and parent records are now consistent, and all data sets refer only to current names (except `intrinsic_resistant` and `microorganisms.codes`)
  * *Erysipelatoclostridium* and its species are now synonyms of *Thomasclavelia*; *E. merdavium* is retired
  * Author references now strip `emend.` and everything after it, so that `ref` reflects the combination authority rather than the emendation author (e.g., *Rhodococcus equi* now returns "Goodfellow et al., 1977" instead of "Nouioui et al., 2018")
* `microorganisms.codes` data set
  * Fixed WHONET codes that referred to an unrelated microorganism, e.g., `VAG` ("Vaginal flora") to *Vagococcus*, `PAR` ("Parasite") to *Paracoccus*, the *Leptospira interrogans* serovars to *Legionella*, and *Salmonella* O-groups to *Salmonella* Group A. Codes for organisms that are not in the `microorganisms` data set, such as viruses, now refer to `UNKNOWN`, and names that only partly exist refer to the genus.
  * WHONET codes for *Salmonella* serovars now refer to the serovar instead of *S. enterica* subsp. *enterica*, e.g., `SAT` to *Salmonella* Typhi
  * Codes now refer to the name they were named after, also if that name is outdated, so that the code and the microorganism code remain relatable, e.g., `EAE` refers to *Enterobacter aerogenes* (`B_ENTRBC_AERG`). `as.mo()` still returns the current name, here *Klebsiella aerogenes*.
  * Added 271 WHONET codes
* Antimicrobial coercion (`as.ab()`) and `antimicrobials` data set
  * Considerably faster through the same session cache design as `as.mo()` (#313)
  * Input that is not a known name or code is now coerced up to 70 times faster, with identical results, as LOINC codes, synonyms, and abbreviations are looked up in one step instead of per antimicrobial, and translated names are built once per language
  * Fixed a bug where results of `fast_mode` (used internally) were remembered and reused by later regular calls, which could return `NA` for valid input (#313)
  * Codes containing "PH" or "TH" (e.g., `ETH`, `PHE`) no longer return `NA` when mixed with unrecognised input (#245)
  * Input that is not valid UTF-8 (e.g., text exported in a local encoding such as CP949) now returns `NA` with a warning instead of an error
  * `ab_reset_session()` no longer clears the uncertainties of `as.mo()` (#313)
  * `antimicrobials$group` is now a `list`, so that drugs belonging to multiple groups are fully represented. Use `ab_group(all_groups = TRUE)` to retrieve all groups of a drug (#246).
  * Trade names of combinations (Tazocin, Timentin, Recarbrio, Exblifep) were listed under a single component, and names of screening libraries and suppliers (e.g., Prestwick, Spectrum) as synonyms. These have been corrected (#315).
* Antibiograms (`antibiogram()` and `wisca()`)
  * `wisca()` now follows Bielicki *et al.* (2016) more strictly by using $\text{Beta}(1, 9999)$ for intrinsically resistant pairs, forcing near-zero susceptibility regardless of observed data (based on EUCAST Expected Resistant Phenotypes)
  * Fixed a bug where both functions failed without the `antimicrobials` argument if `tidyselect` was not available, e.g., for binary installations without `dplyr`
  * Isolates of unknown species are now included in grouped antibiograms and in `wisca()`, as they already were in traditional antibiograms. They were previously excluded silently, with the warning "invalid microorganism code, NA generated". This can change WISCA coverage estimates.
  * Grouped traditional antibiograms now show the (translated) pathogen names, as ungrouped antibiograms do, instead of microorganism codes
  * `wisca()` and `wisca_plot()` are considerably faster, especially with syndromic groups, by determining intrinsic resistance and pathogen names only once per regimen and pathogen
* Resistance and susceptibility calculation
  * New argument `guideline` (default: the `AMR_guideline` option, otherwise EUCAST) in `susceptibility()`, `resistance()`, `count_df()`, `proportion_df()`, `sir_df()`, `ggplot_sir()`, `geom_sir()`, `labels_sir_count()`, `antibiogram()`, `wisca()`, and `bug_drug_combinations()` to interpret the 'I' category correctly per guideline. Functions that return results per category also accept `guideline = NULL` to keep S, I, and R separate, while functions that return a single proportion or count give an informative error.
  * When EUCAST is assumed because neither the `guideline` argument nor the `AMR_guideline` option is set, a message explains how to change this, once per session and not in R Markdown or Quarto documents
  * `count_*()` and `proportion_*()` now handle WT, NWT, NS, and NI results correctly (#321)
* MICs and disk diffusion
  * New argument `round_to_next_log2` in `as.mic()` and `rescale_mic()` to round values up to the nearest log2 dilution level (#255)
  * `as.mic()` now handles values in scientific notation (e.g., `1e-3`) correctly
  * Input validation of `as.disk()` is now more strict, rejecting values that are not recognisable as a numeric disk zone diameter. Values higher than 50 mm now return `NA` with a warning, instead of being set to 50 mm.
* Plotting (`plot()`, `barplot()`, and `autoplot()` for MIC, disk diffusion and SIR values, and the `scale_*_mic()` and `scale_*_sir()` functions)
  * The colours and legend labels of interpreted MIC and disk values are now based on the interpretation itself instead of matching colours, so that the legend of `plot()` no longer shows SDD when `colours_SIR` gives SDD the colour of S (e.g., with one or three colours), and NI and NS are labelled too
  * Fixed `autoplot()` for disk diffusion values, which left SDD bars without colour and legend entry
  * Fixed `barplot()` for MIC values, which ignored the `language` and `expand` arguments, and which now also has the `include_PKPD` and `breakpoint_type` arguments
  * Fixed `barplot()` for SIR values, which gave NI, WT, and NWT the wrong colours, and always showed empty bars for WT, NWT, and NS. `barplot()` and `autoplot()` for SIR values now always show S and R, and all other values (such as SDD, I, WT, NWT, and NS) only when they occur. Data with only WT and NWT shows these two.
  * `autoplot()` for disk diffusion values now validates `include_PKPD` and `breakpoint_type`, and `plot()` and `autoplot()` for SIR values now validate `language`
  * The y axis of `plot()` and `barplot()` now has a readable number of tick marks for large counts, instead of one per count
  * A named `colours_SIR` is now matched by name in all plotting functions, and S, I and R get their default colour if not set. SDD gets the colour of S if not set. Previously, the colours of S/I and I/R were taken by position, and a named vector of three colours gave missing colours.
  * `keep_operators = TRUE` and `FALSE` in the `scale_*_mic()` functions now behave like `"all"` and `"none"`, as in `as.mic()`
  * Fixed the `scale_*_sir()` functions failing when `values` was set. A named `values` (e.g., `values = c(S = "green")`) now replaces the colours of the SIR values it names, and `scale_x_sir()` warns that it ignores `values`.
  * Fixed `scale_x_sir()` ignoring other arguments, such as `name`
  * The `scale_*_sir()` functions now call 'I' "Intermediate" by default if the `AMR_guideline` option is set to EUCAST before 2019
  * `facet_sir()` now accepts unquoted column names, e.g., `facet_sir(antibiotic)`
* MDRO determination (`mdro()`)
  * Resistance can now be inferred from a drug/inhibitor combination when the column of the base drug is absent (e.g., piperacillin inferred from piperacillin/tazobactam). This is controlled with the new argument `infer_from_combinations` (default: `TRUE`) (#209).
  * BRMO classification now includes bacterial complexes (#275)
  * Fixed a bug in the CMI 2012 guideline where pandrug-resistant (PDR) isolates were classified as extensively drug-resistant (XDR), and another isolate could be marked as PDR instead, if the isolate was not among the first rows of the data. This can change the PDR and XDR results of `guideline = "CMI 2012"` (the default).
  * Considerably faster, as the classes per isolate are now determined per column instead of per cell
* `top_n_microorganisms()`
  * New argument `property_for_each` for sub-grouping within the top *n* groups. Rank ordering is enforced, so that only lower taxonomic ranks are allowed.
  * `property = NULL` is now accepted
  * The inner filter now tracks the original row indices, to prevent contamination across groups
* Translations
  * Improved French translations, kindly provided by Mathilde Garé (#322)
  * Improved translations of combination names and antimicrobial groups in 20 languages
  * Translating is now around 50 times faster, with identical results, as the translations are prepared once per language and applied to all values at once. This also makes `ab_group()` faster in other languages.
  * Translations are no longer translated again (e.g., "Ticarcillinee" in Dutch and "groupee" in French), and parts of combination names are no longer left untranslated
  * Fixed Italian names of CoNS and CoPS (#256), Dutch names of antimicrobials, and output of `sir_df()` in other languages (#272)
  * Added translations of "Susceptible dose-dependent", "Non-susceptible", "Wildtype", and "Non-wildtype", which were previously shown half translated in plots
  * Corrected the translations of the interpretation categories (S, I, R, and NI) in 18 languages, which were partly incorrect (e.g., Ukrainian words in Russian, "temporary" for Intermediate in Finnish, and "intermediate product" in Czech), and made them consistent with the newly added categories. Where available, the terms of national antimicrobial susceptibility testing committees and laboratories are used, e.g., "Sensibel bei erhöhter Exposition" (German), "Sensible, cuando se incrementa la exposición" (Spanish), and "Känslig vid ökad exponering" (Swedish).
* Interpretive rules (`interpretive_rules()`, `eucast_rules()`, and `clsi_rules()`)
  * Considerably faster, around 11 to 18 times, with identical results
  * The number of added and changed values now only counts values that were actually written, so that values left untouched due to `overwrite = FALSE` or `add_if_missing = FALSE` are no longer reported as changes or additions, also not in the output of `verbose = TRUE` (#259)
  * The changes of custom rules are now also printed when no other rules were applied
  * The `rule` column of the output of `verbose = TRUE` no longer contains formatting markup (e.g., `{.field AMC}`) for the inhibitor-inheritance rules
* `custom_eucast_rules()` has been renamed to `custom_interpretive_rules()`. The old name is now deprecated but still works (#268).
*  The package option `AMR_eucastrules` has been replaced by `AMR_interpretive_rules`, which applies to `interpretive_rules()`, `eucast_rules()`, and `clsi_rules()`. Its default is now `c("breakpoints", "expected_phenotypes")`. A set `AMR_eucastrules` option is now ignored, with a warning in `eucast_rules()`.
* Logical comparison (`==`, `!=`) of the `sir` class with a `factor` is now possible on R >= 4.3
* New argument `guideline` in `random_mic()` and `random_disk()`, defaulting to the `AMR_guideline` option (otherwise EUCAST), so that generated values are based on the latest breakpoints of one guideline instead of those of all guidelines and years combined
* Fixed `format()` on the output of `bug_drug_combinations()`, which counted I and SDD as resistant with the default `combine_SI = TRUE`, contrary to its documentation. It now reads the guideline used in `bug_drug_combinations()`, which can be overridden with its own `guideline` argument, and also counts NS as resistant.
* Fixed a bug in `g.test()` that returned the degrees of freedom as `argument` instead of `parameter`, so that they were not printed
* Fixed a bug where R could crash after unloading the `data.table` package, and where suggested packages available on the machine that built the package were used even if they were not available to the user

### Other
* Improved console messages with clickable links throughout the package, powered by `cli` if it is installed (#191, #265)
* The section *Interpretation of SIR* (in `?as.sir`) now explains for all functions with a `guideline` argument for aggregated analyses how S, SDD, I, R, WT, NWT, NS, and NI are categorised per guideline. The handling of 'I' per MDRO guideline in `mdro()` will be revised in a future version (#323).
* `labels_sir_count()` is now documented together with `ggplot_sir()`
* Corrected the documentation of `colours_SIR` and `eucast_I` in the plotting functions, and an example of `scale_colour_sir()` in which the title was not applied
* The documentation of `get_AMR_locale()` now describes the actual order in which the language is determined: the `AMR_locale` option, then the `LANG` system variable, then `LANGUAGE` only if `LANG` is not set, and then the system locale
* Revamped console output of `interpretive_rules()`, with section headers, S/I/R values styled as in the printing of `sir` columns in tibbles, and the summary of edits sorted by SIR level, using `cli` if it is installed (#259)
* The printing of `custom_interpretive_rules()` and `custom_mdro_guideline()` now styles SIR values as in tibbles, i.e., the colour of `I` depends on the guideline. `WT`, `NWT`, and `NS` were previously shown as `SDD`, `I`, and `R`, and values other than `S` and `R` as `I`.
* The 'AMR for Python' vignette now contains sections on installation channels (stable from CRAN, or development from GitHub via `AMR.beta`) and on using `enforce_method` in `as_sir()` from Python


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
  * For `as.sir()`, it is now possible to use column names for arguments `guideline`, `ab`, `mo`, and `uti`
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
