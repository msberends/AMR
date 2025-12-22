# Package index

## Introduction to the package

Please find the introduction to (and some general information about) our
package here.

- [`AMR-package`](https://amr-for-r.org/reference/AMR.md)
  [`AMR`](https://amr-for-r.org/reference/AMR.md) :

  The `AMR` Package

## Preparing data: microorganisms

These functions are meant to get taxonomically valid properties of
microorganisms from any input, but also properties derived from
taxonomy, such as the Gram stain
([`mo_gramstain()`](https://amr-for-r.org/reference/mo_property.md)) ,
or [`mo_is_yeast()`](https://amr-for-r.org/reference/mo_property.md).
Use [`mo_source()`](https://amr-for-r.org/reference/mo_source.md) to
teach this package how to translate your own codes to valid
microorganisms, and use
[`add_custom_microorganisms()`](https://amr-for-r.org/reference/add_custom_microorganisms.md)
to add your own custom microorganisms to this package.

- [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)
  [`is.mo()`](https://amr-for-r.org/reference/as.mo.md)
  [`mo_uncertainties()`](https://amr-for-r.org/reference/as.mo.md)
  [`mo_renamed()`](https://amr-for-r.org/reference/as.mo.md)
  [`mo_failures()`](https://amr-for-r.org/reference/as.mo.md)
  [`mo_reset_session()`](https://amr-for-r.org/reference/as.mo.md)
  [`mo_cleaning_regex()`](https://amr-for-r.org/reference/as.mo.md) :
  Transform Arbitrary Input to Valid Microbial Taxonomy
- [`mo_name()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_fullname()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_shortname()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_subspecies()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_species()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_genus()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_family()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_order()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_class()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_phylum()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_kingdom()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_domain()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_type()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_status()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_pathogenicity()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_gramstain()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_is_gram_negative()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_is_gram_positive()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_is_yeast()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_is_intrinsic_resistant()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_oxygen_tolerance()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_is_anaerobic()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_snomed()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_ref()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_authors()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_year()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_lpsn()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_mycobank()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_gbif()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_rank()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_taxonomy()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_synonyms()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_current()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_group_members()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_info()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_url()`](https://amr-for-r.org/reference/mo_property.md)
  [`mo_property()`](https://amr-for-r.org/reference/mo_property.md) :
  Get Properties of a Microorganism
- [`add_custom_microorganisms()`](https://amr-for-r.org/reference/add_custom_microorganisms.md)
  [`clear_custom_microorganisms()`](https://amr-for-r.org/reference/add_custom_microorganisms.md)
  : Add Custom Microorganisms
- [`set_mo_source()`](https://amr-for-r.org/reference/mo_source.md)
  [`get_mo_source()`](https://amr-for-r.org/reference/mo_source.md) :
  User-Defined Reference Data Set for Microorganisms

## Preparing data: antimicrobials

Use these functions to get valid properties of antimicrobials from any
input or to clean your input. You can even retrieve drug names and doses
from clinical text records, using
[`ab_from_text()`](https://amr-for-r.org/reference/ab_from_text.md).

- [`as.ab()`](https://amr-for-r.org/reference/as.ab.md)
  [`is.ab()`](https://amr-for-r.org/reference/as.ab.md)
  [`ab_reset_session()`](https://amr-for-r.org/reference/as.ab.md) :
  Transform Input to an Antibiotic ID
- [`ab_name()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_cid()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_synonyms()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_tradenames()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_group()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_atc()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_atc_group1()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_atc_group2()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_loinc()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_ddd()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_ddd_units()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_info()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_url()`](https://amr-for-r.org/reference/ab_property.md)
  [`ab_property()`](https://amr-for-r.org/reference/ab_property.md)
  [`set_ab_names()`](https://amr-for-r.org/reference/ab_property.md) :
  Get Properties of an Antibiotic
- [`ab_from_text()`](https://amr-for-r.org/reference/ab_from_text.md) :
  Retrieve Antimicrobial Drug Names and Doses from Clinical Text
- [`atc_online_property()`](https://amr-for-r.org/reference/atc_online.md)
  [`atc_online_groups()`](https://amr-for-r.org/reference/atc_online.md)
  [`atc_online_ddd()`](https://amr-for-r.org/reference/atc_online.md)
  [`atc_online_ddd_units()`](https://amr-for-r.org/reference/atc_online.md)
  : Get ATC Properties from WHOCC Website
- [`add_custom_antimicrobials()`](https://amr-for-r.org/reference/add_custom_antimicrobials.md)
  [`clear_custom_antimicrobials()`](https://amr-for-r.org/reference/add_custom_antimicrobials.md)
  : Add Custom Antimicrobials

## Preparing data: antimicrobial results

With [`as.mic()`](https://amr-for-r.org/reference/as.mic.md) and
[`as.disk()`](https://amr-for-r.org/reference/as.disk.md) you can
transform your raw input to valid MIC or disk diffusion values. Use
[`as.sir()`](https://amr-for-r.org/reference/as.sir.md) for cleaning raw
data to let it only contain “R”, “I” and “S”, or to interpret MIC or
disk diffusion values as SIR based on the lastest EUCAST and CLSI
guidelines. Afterwards, you can extend antibiotic interpretations by
applying [EUCAST
rules](https://www.eucast.org/expert_rules_and_intrinsic_resistance/)
with
[`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md).

- [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)
  [`NA_sir_`](https://amr-for-r.org/reference/as.sir.md)
  [`is.sir()`](https://amr-for-r.org/reference/as.sir.md)
  [`is_sir_eligible()`](https://amr-for-r.org/reference/as.sir.md)
  [`sir_interpretation_history()`](https://amr-for-r.org/reference/as.sir.md)
  : Interpret MIC and Disk Diffusion as SIR, or Clean Existing SIR Data
- [`as.mic()`](https://amr-for-r.org/reference/as.mic.md)
  [`is.mic()`](https://amr-for-r.org/reference/as.mic.md)
  [`NA_mic_`](https://amr-for-r.org/reference/as.mic.md)
  [`rescale_mic()`](https://amr-for-r.org/reference/as.mic.md)
  [`mic_p50()`](https://amr-for-r.org/reference/as.mic.md)
  [`mic_p90()`](https://amr-for-r.org/reference/as.mic.md)
  [`droplevels(`*`<mic>`*`)`](https://amr-for-r.org/reference/as.mic.md)
  : Transform Input to Minimum Inhibitory Concentrations (MIC)
- [`as.disk()`](https://amr-for-r.org/reference/as.disk.md)
  [`NA_disk_`](https://amr-for-r.org/reference/as.disk.md)
  [`is.disk()`](https://amr-for-r.org/reference/as.disk.md) : Transform
  Input to Disk Diffusion Diameters
- [`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md)
  [`eucast_dosage()`](https://amr-for-r.org/reference/eucast_rules.md) :
  Apply EUCAST Rules
- [`custom_eucast_rules()`](https://amr-for-r.org/reference/custom_eucast_rules.md)
  : Define Custom EUCAST Rules

## Analysing data

Use these function for the analysis part. You can use
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md) or
[`resistance()`](https://amr-for-r.org/reference/proportion.md) on any
antibiotic column. With
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md), you
can generate a traditional, combined, syndromic, or weighted-incidence
syndromic combination antibiogram (WISCA). This function also comes with
support for R Markdown and Quarto. Be sure to first select the isolates
that are appropiate for analysis, by using
[`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md) or
[`is_new_episode()`](https://amr-for-r.org/reference/get_episode.md).
You can also filter your data on certain resistance in certain
antibiotic classes
([`carbapenems()`](https://amr-for-r.org/reference/antimicrobial_selectors.md),
[`aminoglycosides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)),
or determine multi-drug resistant microorganisms (MDRO,
[`mdro()`](https://amr-for-r.org/reference/mdro.md)).

- [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md)
  [`wisca()`](https://amr-for-r.org/reference/antibiogram.md)
  [`retrieve_wisca_parameters()`](https://amr-for-r.org/reference/antibiogram.md)
  [`plot(`*`<antibiogram>`*`)`](https://amr-for-r.org/reference/antibiogram.md)
  [`autoplot(`*`<antibiogram>`*`)`](https://amr-for-r.org/reference/antibiogram.md)
  [`knit_print(`*`<antibiogram>`*`)`](https://amr-for-r.org/reference/antibiogram.md)
  : Generate Traditional, Combination, Syndromic, or WISCA Antibiograms

- [`resistance()`](https://amr-for-r.org/reference/proportion.md)
  [`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
  [`sir_confidence_interval()`](https://amr-for-r.org/reference/proportion.md)
  [`proportion_R()`](https://amr-for-r.org/reference/proportion.md)
  [`proportion_IR()`](https://amr-for-r.org/reference/proportion.md)
  [`proportion_I()`](https://amr-for-r.org/reference/proportion.md)
  [`proportion_SI()`](https://amr-for-r.org/reference/proportion.md)
  [`proportion_S()`](https://amr-for-r.org/reference/proportion.md)
  [`proportion_df()`](https://amr-for-r.org/reference/proportion.md)
  [`sir_df()`](https://amr-for-r.org/reference/proportion.md) :
  Calculate Antimicrobial Resistance

- [`count_resistant()`](https://amr-for-r.org/reference/count.md)
  [`count_susceptible()`](https://amr-for-r.org/reference/count.md)
  [`count_S()`](https://amr-for-r.org/reference/count.md)
  [`count_SI()`](https://amr-for-r.org/reference/count.md)
  [`count_I()`](https://amr-for-r.org/reference/count.md)
  [`count_IR()`](https://amr-for-r.org/reference/count.md)
  [`count_R()`](https://amr-for-r.org/reference/count.md)
  [`count_all()`](https://amr-for-r.org/reference/count.md)
  [`n_sir()`](https://amr-for-r.org/reference/count.md)
  [`count_df()`](https://amr-for-r.org/reference/count.md) : Count
  Available Isolates

- [`get_episode()`](https://amr-for-r.org/reference/get_episode.md)
  [`is_new_episode()`](https://amr-for-r.org/reference/get_episode.md) :
  Determine Clinical or Epidemic Episodes

- [`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md)
  [`filter_first_isolate()`](https://amr-for-r.org/reference/first_isolate.md)
  : Determine First Isolates

- [`key_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md)
  [`all_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md)
  [`antimicrobials_equal()`](https://amr-for-r.org/reference/key_antimicrobials.md)
  : (Key) Antimicrobials for First Weighted Isolates

- [`mdro()`](https://amr-for-r.org/reference/mdro.md)
  [`brmo()`](https://amr-for-r.org/reference/mdro.md)
  [`mrgn()`](https://amr-for-r.org/reference/mdro.md)
  [`mdr_tb()`](https://amr-for-r.org/reference/mdro.md)
  [`mdr_cmi2012()`](https://amr-for-r.org/reference/mdro.md)
  [`eucast_exceptional_phenotypes()`](https://amr-for-r.org/reference/mdro.md)
  : Determine Multidrug-Resistant Organisms (MDRO)

- [`custom_mdro_guideline()`](https://amr-for-r.org/reference/custom_mdro_guideline.md)
  [`c(`*`<custom_mdro_guideline>`*`)`](https://amr-for-r.org/reference/custom_mdro_guideline.md)
  : Define Custom MDRO Guideline

- [`bug_drug_combinations()`](https://amr-for-r.org/reference/bug_drug_combinations.md)
  [`format(`*`<bug_drug_combinations>`*`)`](https://amr-for-r.org/reference/bug_drug_combinations.md)
  : Determine Bug-Drug Combinations

- [`aminoglycosides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`aminopenicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`antifungals()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`antimycobacterials()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`betalactams()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`betalactams_with_inhibitor()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`carbapenems()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`cephalosporins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`cephalosporins_1st()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`cephalosporins_2nd()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`cephalosporins_3rd()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`cephalosporins_4th()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`cephalosporins_5th()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`fluoroquinolones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`glycopeptides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`isoxazolylpenicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`lincosamides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`lipoglycopeptides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`macrolides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`monobactams()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`nitrofurans()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`oxazolidinones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`penicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`phenicols()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`polymyxins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`quinolones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`rifamycins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`streptogramins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`sulfonamides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`tetracyclines()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`trimethoprims()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`ureidopenicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`amr_class()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`amr_selector()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`administrable_per_os()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`administrable_iv()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  [`not_intrinsic_resistant()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  : Antimicrobial Selectors

- [`top_n_microorganisms()`](https://amr-for-r.org/reference/top_n_microorganisms.md)
  :

  Filter Top *n* Microorganisms

- [`mean_amr_distance()`](https://amr-for-r.org/reference/mean_amr_distance.md)
  [`amr_distance_from_row()`](https://amr-for-r.org/reference/mean_amr_distance.md)
  : Calculate the Mean AMR Distance

- [`resistance_predict()`](https://amr-for-r.org/reference/resistance_predict.md)
  [`sir_predict()`](https://amr-for-r.org/reference/resistance_predict.md)
  [`plot(`*`<resistance_predict>`*`)`](https://amr-for-r.org/reference/resistance_predict.md)
  [`ggplot_sir_predict()`](https://amr-for-r.org/reference/resistance_predict.md)
  [`autoplot(`*`<resistance_predict>`*`)`](https://amr-for-r.org/reference/resistance_predict.md)
  : Predict Antimicrobial Resistance

- [`guess_ab_col()`](https://amr-for-r.org/reference/guess_ab_col.md) :
  Guess Antibiotic Column

## Plotting data

Use these functions for the plotting part. The `scale_*_mic()` functions
extend the ggplot2 package to allow plotting of MIC values, even within
a manually set range. If using
[`plot()`](https://amr-for-r.org/reference/plot.md) (base R) or
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
(ggplot2) on MIC values or disk diffusion values, the user can set the
interpretation guideline to give the bars the right SIR colours. The
[`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md) function
is a short wrapper for users not much accustomed to ggplot2 yet. The
[`ggplot_pca()`](https://amr-for-r.org/reference/ggplot_pca.md) function
is a specific function to plot so-called biplots for PCA (principal
component analysis).

- [`scale_x_mic()`](https://amr-for-r.org/reference/plot.md)
  [`scale_y_mic()`](https://amr-for-r.org/reference/plot.md)
  [`scale_colour_mic()`](https://amr-for-r.org/reference/plot.md)
  [`scale_fill_mic()`](https://amr-for-r.org/reference/plot.md)
  [`scale_x_sir()`](https://amr-for-r.org/reference/plot.md)
  [`scale_colour_sir()`](https://amr-for-r.org/reference/plot.md)
  [`scale_fill_sir()`](https://amr-for-r.org/reference/plot.md)
  [`plot(`*`<mic>`*`)`](https://amr-for-r.org/reference/plot.md)
  [`autoplot(`*`<mic>`*`)`](https://amr-for-r.org/reference/plot.md)
  [`plot(`*`<disk>`*`)`](https://amr-for-r.org/reference/plot.md)
  [`autoplot(`*`<disk>`*`)`](https://amr-for-r.org/reference/plot.md)
  [`plot(`*`<sir>`*`)`](https://amr-for-r.org/reference/plot.md)
  [`autoplot(`*`<sir>`*`)`](https://amr-for-r.org/reference/plot.md)
  [`facet_sir()`](https://amr-for-r.org/reference/plot.md)
  [`scale_y_percent()`](https://amr-for-r.org/reference/plot.md)
  [`scale_sir_colours()`](https://amr-for-r.org/reference/plot.md)
  [`theme_sir()`](https://amr-for-r.org/reference/plot.md)
  [`labels_sir_count()`](https://amr-for-r.org/reference/plot.md) :
  Plotting Helpers for AMR Data Analysis

- [`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md)
  [`geom_sir()`](https://amr-for-r.org/reference/ggplot_sir.md) :

  AMR Plots with `ggplot2`

- [`ggplot_pca()`](https://amr-for-r.org/reference/ggplot_pca.md) :

  PCA Biplot with `ggplot2`

## AMR-specific options

The AMR package is customisable, by providing settings that can be set
per user or per team. For example, the default interpretation guideline
can be changed from EUCAST to CLSI, or a supported language can be set
for the whole team (system-language independent) for antibiotic names in
a foreign language.

- [`AMR-options`](https://amr-for-r.org/reference/AMR-options.md) :
  Options for the AMR package

## Other: antiviral drugs

This package also provides extensive support for antiviral agents, even
though it is not the primary scope of this package. Working with data
containing information about antiviral drugs was never easier. Use these
functions to get valid properties of antiviral drugs from any input or
to clean your input. You can even retrieve drug names and doses from
clinical text records, using
[`av_from_text()`](https://amr-for-r.org/reference/av_from_text.md).

- [`as.av()`](https://amr-for-r.org/reference/as.av.md)
  [`is.av()`](https://amr-for-r.org/reference/as.av.md) : Transform
  Input to an Antiviral Drug ID
- [`av_name()`](https://amr-for-r.org/reference/av_property.md)
  [`av_cid()`](https://amr-for-r.org/reference/av_property.md)
  [`av_synonyms()`](https://amr-for-r.org/reference/av_property.md)
  [`av_tradenames()`](https://amr-for-r.org/reference/av_property.md)
  [`av_group()`](https://amr-for-r.org/reference/av_property.md)
  [`av_atc()`](https://amr-for-r.org/reference/av_property.md)
  [`av_loinc()`](https://amr-for-r.org/reference/av_property.md)
  [`av_ddd()`](https://amr-for-r.org/reference/av_property.md)
  [`av_ddd_units()`](https://amr-for-r.org/reference/av_property.md)
  [`av_info()`](https://amr-for-r.org/reference/av_property.md)
  [`av_url()`](https://amr-for-r.org/reference/av_property.md)
  [`av_property()`](https://amr-for-r.org/reference/av_property.md) :
  Get Properties of an Antiviral Drug
- [`av_from_text()`](https://amr-for-r.org/reference/av_from_text.md) :
  Retrieve Antiviral Drug Names and Doses from Clinical Text

## Other: background information on included data

Some pages about our package and its external sources. Be sure to read
our [How To’s](https://amr-for-r.org/articles/index.md) for more
information about how to work with functions in this package.

- [`microorganisms`](https://amr-for-r.org/reference/microorganisms.md)
  : Data Set with 78 679 Taxonomic Records of Microorganisms
- [`antimicrobials`](https://amr-for-r.org/reference/antimicrobials.md)
  [`antibiotics`](https://amr-for-r.org/reference/antimicrobials.md)
  [`antivirals`](https://amr-for-r.org/reference/antimicrobials.md) :
  Data Sets with 618 Antimicrobial Drugs
- [`clinical_breakpoints`](https://amr-for-r.org/reference/clinical_breakpoints.md)
  : Data Set with Clinical Breakpoints for SIR Interpretation
- [`example_isolates`](https://amr-for-r.org/reference/example_isolates.md)
  : Data Set with 2 000 Example Isolates
- [`esbl_isolates`](https://amr-for-r.org/reference/esbl_isolates.md) :
  Data Set with 500 ESBL Isolates
- [`microorganisms.codes`](https://amr-for-r.org/reference/microorganisms.codes.md)
  : Data Set with 6 036 Common Microorganism Codes
- [`microorganisms.groups`](https://amr-for-r.org/reference/microorganisms.groups.md)
  : Data Set with 534 Microorganisms In Species Groups
- [`intrinsic_resistant`](https://amr-for-r.org/reference/intrinsic_resistant.md)
  : Data Set Denoting Bacterial Intrinsic Resistance
- [`dosage`](https://amr-for-r.org/reference/dosage.md) : Data Set with
  Treatment Dosages as Defined by EUCAST
- [`WHOCC`](https://amr-for-r.org/reference/WHOCC.md) : WHOCC: WHO
  Collaborating Centre for Drug Statistics Methodology
- [`example_isolates_unclean`](https://amr-for-r.org/reference/example_isolates_unclean.md)
  : Data Set with Unclean Data
- [`WHONET`](https://amr-for-r.org/reference/WHONET.md) : Data Set with
  500 Isolates - WHONET Example

## Other: miscellaneous functions

These functions are mostly for internal use, but some of them may also
be suitable for your analysis. Especially the ‘like’ function can be
useful: `if (x %like% y) {...}`.

- [`age_groups()`](https://amr-for-r.org/reference/age_groups.md) :
  Split Ages into Age Groups
- [`age()`](https://amr-for-r.org/reference/age.md) : Age in Years of
  Individuals
- [`export_ncbi_biosample()`](https://amr-for-r.org/reference/export_ncbi_biosample.md)
  : Export Data Set as NCBI BioSample Antibiogram
- [`availability()`](https://amr-for-r.org/reference/availability.md) :
  Check Availability of Columns
- [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md)
  [`set_AMR_locale()`](https://amr-for-r.org/reference/translate.md)
  [`reset_AMR_locale()`](https://amr-for-r.org/reference/translate.md)
  [`translate_AMR()`](https://amr-for-r.org/reference/translate.md) :
  Translate Strings from the AMR Package
- [`italicise_taxonomy()`](https://amr-for-r.org/reference/italicise_taxonomy.md)
  [`italicize_taxonomy()`](https://amr-for-r.org/reference/italicise_taxonomy.md)
  : Italicise Taxonomic Families, Genera, Species, Subspecies
- [`inner_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
  [`left_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
  [`right_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
  [`full_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
  [`semi_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
  [`anti_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
  : Join microorganisms to a Data Set
- [`like()`](https://amr-for-r.org/reference/like.md)
  [`` `%like%` ``](https://amr-for-r.org/reference/like.md)
  [`` `%unlike%` ``](https://amr-for-r.org/reference/like.md)
  [`` `%like_case%` ``](https://amr-for-r.org/reference/like.md)
  [`` `%unlike_case%` ``](https://amr-for-r.org/reference/like.md) :
  Vectorised Pattern Matching with Keyboard Shortcut
- [`mo_matching_score()`](https://amr-for-r.org/reference/mo_matching_score.md)
  : Calculate the Matching Score for Microorganisms
- [`pca()`](https://amr-for-r.org/reference/pca.md) : Principal
  Component Analysis (for AMR)
- [`random_mic()`](https://amr-for-r.org/reference/random.md)
  [`random_disk()`](https://amr-for-r.org/reference/random.md)
  [`random_sir()`](https://amr-for-r.org/reference/random.md) : Random
  MIC Values/Disk Zones/SIR Generation

## Other: statistical tests

Some statistical tests or methods are not part of base R and were added
to this package for convenience.

- [`g.test()`](https://amr-for-r.org/reference/g.test.md) :

  *G*-test for Count Data

- [`kurtosis()`](https://amr-for-r.org/reference/kurtosis.md) : Kurtosis
  of the Sample

- [`skewness()`](https://amr-for-r.org/reference/skewness.md) : Skewness
  of the Sample

## Other: deprecated functions/arguments/datasets

These objects are deprecated, meaning that they will still work but show
a warning that they will be removed in a future version.

- [`ab_class()`](https://amr-for-r.org/reference/AMR-deprecated.md)
  [`ab_selector()`](https://amr-for-r.org/reference/AMR-deprecated.md) :
  Deprecated Functions, Arguments, or Datasets
