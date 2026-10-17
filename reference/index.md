# Package index

## Getting started

An overview of the package, its options for users and teams, and its 28
supported languages.

[The `AMR` Package](https://amr-for-r.org/reference/AMR.md)
:   [`AMR-package`](https://amr-for-r.org/reference/AMR.md)
    [`AMR`](https://amr-for-r.org/reference/AMR.md)

[Options for the AMR Package](https://amr-for-r.org/reference/AMR-options.md)
:   [`AMR-options`](https://amr-for-r.org/reference/AMR-options.md)

[Translate Strings from the AMR Package](https://amr-for-r.org/reference/translate.md)
:   [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md)
    [`set_AMR_locale()`](https://amr-for-r.org/reference/translate.md)
    [`reset_AMR_locale()`](https://amr-for-r.org/reference/translate.md)
    [`translate_AMR()`](https://amr-for-r.org/reference/translate.md)

[Download and Unpack an AMR Course Repository](https://amr-for-r.org/reference/amr_course.md)
:   [`amr_course()`](https://amr-for-r.org/reference/amr_course.md)

## Cleaning data

Transform raw input of any kind into valid microorganisms with
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md), antimicrobials
with [`as.ab()`](https://amr-for-r.org/reference/as.ab.md), and MIC
values and disk diffusion zones with
[`as.mic()`](https://amr-for-r.org/reference/as.mic.md) and
[`as.disk()`](https://amr-for-r.org/reference/as.disk.md), and retrieve
their properties.

### Microorganisms

[Transform Arbitrary Input to Valid Microbial Taxonomy](https://amr-for-r.org/reference/as.mo.md)
:   [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)
    [`is.mo()`](https://amr-for-r.org/reference/as.mo.md)
    [`mo_uncertainties()`](https://amr-for-r.org/reference/as.mo.md)
    [`mo_renamed()`](https://amr-for-r.org/reference/as.mo.md)
    [`mo_failures()`](https://amr-for-r.org/reference/as.mo.md)
    [`mo_reset_session()`](https://amr-for-r.org/reference/as.mo.md)
    [`mo_cleaning_regex()`](https://amr-for-r.org/reference/as.mo.md)
    [`NA_mo_`](https://amr-for-r.org/reference/as.mo.md)

[Get Properties of a Microorganism](https://amr-for-r.org/reference/mo_property.md)
:   [`mo_name()`](https://amr-for-r.org/reference/mo_property.md)
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
    [`mo_morphology()`](https://amr-for-r.org/reference/mo_property.md)
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
    [`mo_property()`](https://amr-for-r.org/reference/mo_property.md)

[User-Defined Reference Data Set for Microorganisms](https://amr-for-r.org/reference/mo_source.md)
:   [`set_mo_source()`](https://amr-for-r.org/reference/mo_source.md)
    [`get_mo_source()`](https://amr-for-r.org/reference/mo_source.md)

[Add Custom Microorganisms](https://amr-for-r.org/reference/add_custom_microorganisms.md)
:   [`add_custom_microorganisms()`](https://amr-for-r.org/reference/add_custom_microorganisms.md)
    [`clear_custom_microorganisms()`](https://amr-for-r.org/reference/add_custom_microorganisms.md)

[Calculate the Matching Score for Microorganisms](https://amr-for-r.org/reference/mo_matching_score.md)
:   [`mo_matching_score()`](https://amr-for-r.org/reference/mo_matching_score.md)

### Antimicrobials

[Transform Input to an Antimicrobial Drug ID](https://amr-for-r.org/reference/as.ab.md)
:   [`as.ab()`](https://amr-for-r.org/reference/as.ab.md)
    [`is.ab()`](https://amr-for-r.org/reference/as.ab.md)
    [`ab_reset_session()`](https://amr-for-r.org/reference/as.ab.md)
    [`NA_ab_`](https://amr-for-r.org/reference/as.ab.md)

[Get Properties of an Antimicrobial Drug](https://amr-for-r.org/reference/ab_property.md)
:   [`ab_name()`](https://amr-for-r.org/reference/ab_property.md)
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
    [`set_ab_names()`](https://amr-for-r.org/reference/ab_property.md)

[Retrieve Antimicrobial Drug Names and Doses from Clinical Text](https://amr-for-r.org/reference/ab_from_text.md)
:   [`ab_from_text()`](https://amr-for-r.org/reference/ab_from_text.md)

[Add Custom Antimicrobials](https://amr-for-r.org/reference/add_custom_antimicrobials.md)
:   [`add_custom_antimicrobials()`](https://amr-for-r.org/reference/add_custom_antimicrobials.md)
    [`add_custom_antimicrobial_synonyms()`](https://amr-for-r.org/reference/add_custom_antimicrobials.md)
    [`clear_custom_antimicrobials()`](https://amr-for-r.org/reference/add_custom_antimicrobials.md)

[Get ATC Properties from WHOCC Website](https://amr-for-r.org/reference/atc_online.md)
:   [`atc_online_property()`](https://amr-for-r.org/reference/atc_online.md)
    [`atc_online_groups()`](https://amr-for-r.org/reference/atc_online.md)
    [`atc_online_ddd()`](https://amr-for-r.org/reference/atc_online.md)
    [`atc_online_ddd_units()`](https://amr-for-r.org/reference/atc_online.md)

### Antivirals

[Transform Input to an Antiviral Drug ID](https://amr-for-r.org/reference/as.av.md)
:   [`as.av()`](https://amr-for-r.org/reference/as.av.md)
    [`is.av()`](https://amr-for-r.org/reference/as.av.md)

[Get Properties of an Antiviral Drug](https://amr-for-r.org/reference/av_property.md)
:   [`av_name()`](https://amr-for-r.org/reference/av_property.md)
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
    [`av_property()`](https://amr-for-r.org/reference/av_property.md)

[Retrieve Antiviral Drug Names and Doses from Clinical Text](https://amr-for-r.org/reference/av_from_text.md)
:   [`av_from_text()`](https://amr-for-r.org/reference/av_from_text.md)

### Test results

[Transform Input to Minimum Inhibitory Concentrations (MIC)](https://amr-for-r.org/reference/as.mic.md)
:   [`as.mic()`](https://amr-for-r.org/reference/as.mic.md)
    [`is.mic()`](https://amr-for-r.org/reference/as.mic.md)
    [`NA_mic_`](https://amr-for-r.org/reference/as.mic.md)
    [`rescale_mic()`](https://amr-for-r.org/reference/as.mic.md)
    [`mic_p50()`](https://amr-for-r.org/reference/as.mic.md)
    [`mic_p90()`](https://amr-for-r.org/reference/as.mic.md)
    [`droplevels(`*`<mic>`*`)`](https://amr-for-r.org/reference/as.mic.md)

[Transform Input to Disk Diffusion Diameters](https://amr-for-r.org/reference/as.disk.md)
:   [`as.disk()`](https://amr-for-r.org/reference/as.disk.md)
    [`NA_disk_`](https://amr-for-r.org/reference/as.disk.md)
    [`is.disk()`](https://amr-for-r.org/reference/as.disk.md)

## Interpreting results as SIR

Interpret MIC values and disk diffusion zones as SIR with
[`as.sir()`](https://amr-for-r.org/reference/as.sir.md), according to
EUCAST or CLSI, for every guideline year from 2011 to 2026, and apply
interpretive rules.
[`as.sir()`](https://amr-for-r.org/reference/as.sir.md) also cleans
existing SIR data.

[Interpret MIC and Disk Diffusion as SIR, or Clean Existing SIR Data](https://amr-for-r.org/reference/as.sir.md)
:   [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)
    [`NA_sir_`](https://amr-for-r.org/reference/as.sir.md)
    [`is.sir()`](https://amr-for-r.org/reference/as.sir.md)
    [`is_sir_eligible()`](https://amr-for-r.org/reference/as.sir.md)
    [`sir_interpretation_history()`](https://amr-for-r.org/reference/as.sir.md)

[Apply Interpretive Rules](https://amr-for-r.org/reference/interpretive_rules.md)
:   [`interpretive_rules()`](https://amr-for-r.org/reference/interpretive_rules.md)
    [`eucast_rules()`](https://amr-for-r.org/reference/interpretive_rules.md)
    [`clsi_rules()`](https://amr-for-r.org/reference/interpretive_rules.md)
    [`eucast_dosage()`](https://amr-for-r.org/reference/interpretive_rules.md)

[Define Custom Interpretive Rules](https://amr-for-r.org/reference/custom_interpretive_rules.md)
:   [`custom_interpretive_rules()`](https://amr-for-r.org/reference/custom_interpretive_rules.md)

## Analysing data

Select first isolates with
[`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md),
generate traditional, combination, syndromic, and WISCA antibiograms
with [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md)
and [`wisca()`](https://amr-for-r.org/reference/antibiogram.md),
calculate resistance, and determine multidrug-resistant organisms.

### Selecting isolates and antimicrobials

[Determine First Isolates](https://amr-for-r.org/reference/first_isolate.md)
:   [`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md)
    [`filter_first_isolate()`](https://amr-for-r.org/reference/first_isolate.md)

[Determine Clinical or Epidemic Episodes](https://amr-for-r.org/reference/get_episode.md)
:   [`get_episode()`](https://amr-for-r.org/reference/get_episode.md)
    [`is_new_episode()`](https://amr-for-r.org/reference/get_episode.md)

[(Key) Antimicrobials for First Weighted Isolates](https://amr-for-r.org/reference/key_antimicrobials.md)
:   [`key_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md)
    [`all_antimicrobials()`](https://amr-for-r.org/reference/key_antimicrobials.md)
    [`antimicrobials_equal()`](https://amr-for-r.org/reference/key_antimicrobials.md)

[Filter Top *n* Microorganisms](https://amr-for-r.org/reference/top_n_microorganisms.md)
:   [`top_n_microorganisms()`](https://amr-for-r.org/reference/top_n_microorganisms.md)

[Antimicrobial Selectors](https://amr-for-r.org/reference/antimicrobial_selectors.md)
:   [`aminoglycosides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
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
    [`ionophores()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`isoxazolylpenicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`lincosamides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`lipoglycopeptides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`macrolides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`monobactams()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`nitrofurans()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`oxazolidinones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`penicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`peptides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`phenicols()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`phosphonics()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`polymyxins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`quinolones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`rifamycins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
    [`spiropyrimidinetriones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
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

[Guess Antimicrobial Column](https://amr-for-r.org/reference/guess_ab_col.md)
:   [`guess_ab_col()`](https://amr-for-r.org/reference/guess_ab_col.md)

[Check Availability of Columns](https://amr-for-r.org/reference/availability.md)
:   [`availability()`](https://amr-for-r.org/reference/availability.md)

### Generating antibiograms and WISCA

[Generate Antibiograms (WISCA, Traditional, Combination, or Syndromic)](https://amr-for-r.org/reference/antibiogram.md)
:   [`wisca()`](https://amr-for-r.org/reference/antibiogram.md)
    [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md)
    [`retrieve_wisca_parameters()`](https://amr-for-r.org/reference/antibiogram.md)
    [`plot(`*`<antibiogram>`*`)`](https://amr-for-r.org/reference/antibiogram.md)
    [`autoplot(`*`<antibiogram>`*`)`](https://amr-for-r.org/reference/antibiogram.md)
    [`wisca_plot()`](https://amr-for-r.org/reference/antibiogram.md)
    [`knit_print(`*`<antibiogram>`*`)`](https://amr-for-r.org/reference/antibiogram.md)

### Calculating resistance and susceptibility

[Calculate Antimicrobial Resistance](https://amr-for-r.org/reference/proportion.md)
:   [`resistance()`](https://amr-for-r.org/reference/proportion.md)
    [`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
    [`sir_confidence_interval()`](https://amr-for-r.org/reference/proportion.md)
    [`proportion_R()`](https://amr-for-r.org/reference/proportion.md)
    [`proportion_IR()`](https://amr-for-r.org/reference/proportion.md)
    [`proportion_I()`](https://amr-for-r.org/reference/proportion.md)
    [`proportion_SI()`](https://amr-for-r.org/reference/proportion.md)
    [`proportion_S()`](https://amr-for-r.org/reference/proportion.md)
    [`proportion_df()`](https://amr-for-r.org/reference/proportion.md)
    [`sir_df()`](https://amr-for-r.org/reference/proportion.md)

[Count Available Isolates](https://amr-for-r.org/reference/count.md)
:   [`count_resistant()`](https://amr-for-r.org/reference/count.md)
    [`count_susceptible()`](https://amr-for-r.org/reference/count.md)
    [`count_S()`](https://amr-for-r.org/reference/count.md)
    [`count_SI()`](https://amr-for-r.org/reference/count.md)
    [`count_I()`](https://amr-for-r.org/reference/count.md)
    [`count_IR()`](https://amr-for-r.org/reference/count.md)
    [`count_R()`](https://amr-for-r.org/reference/count.md)
    [`count_all()`](https://amr-for-r.org/reference/count.md)
    [`n_sir()`](https://amr-for-r.org/reference/count.md)
    [`count_df()`](https://amr-for-r.org/reference/count.md)

[Combine SIR Results of Multiple Antimicrobials](https://amr-for-r.org/reference/sir_combine.md)
:   [`sir_combine()`](https://amr-for-r.org/reference/sir_combine.md)

[Determine Bug-Drug Combinations](https://amr-for-r.org/reference/bug_drug_combinations.md)
:   [`bug_drug_combinations()`](https://amr-for-r.org/reference/bug_drug_combinations.md)
    [`format(`*`<bug_drug_combinations>`*`)`](https://amr-for-r.org/reference/bug_drug_combinations.md)

### Determining multidrug resistance

[Determine Multidrug-Resistant Organisms (MDRO)](https://amr-for-r.org/reference/mdro.md)
:   [`mdro()`](https://amr-for-r.org/reference/mdro.md)
    [`brmo()`](https://amr-for-r.org/reference/mdro.md)
    [`mrgn()`](https://amr-for-r.org/reference/mdro.md)
    [`mdr_tb()`](https://amr-for-r.org/reference/mdro.md)
    [`mdr_cmi2012()`](https://amr-for-r.org/reference/mdro.md)
    [`eucast_exceptional_phenotypes()`](https://amr-for-r.org/reference/mdro.md)

[Define Custom MDRO Guideline](https://amr-for-r.org/reference/custom_mdro_guideline.md)
:   [`custom_mdro_guideline()`](https://amr-for-r.org/reference/custom_mdro_guideline.md)
    [`c(`*`<custom_mdro_guideline>`*`)`](https://amr-for-r.org/reference/custom_mdro_guideline.md)

## Visualising data

Plot MIC values, disk zones, and SIR results in base R or `ggplot2`,
with scales and colours tailored to AMR data.

[Plotting Helpers for AMR Data Analysis](https://amr-for-r.org/reference/plot.md)
:   [`scale_x_mic()`](https://amr-for-r.org/reference/plot.md)
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

[AMR Plots with `ggplot2`](https://amr-for-r.org/reference/ggplot_sir.md)
:   [`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md)
    [`geom_sir()`](https://amr-for-r.org/reference/ggplot_sir.md)
    [`labels_sir_count()`](https://amr-for-r.org/reference/ggplot_sir.md)

[PCA Biplot with `ggplot2`](https://amr-for-r.org/reference/ggplot_pca.md)
:   [`ggplot_pca()`](https://amr-for-r.org/reference/ggplot_pca.md)

## Modelling

Use AMR data types in `tidymodels` pipelines, and perform principal
component analysis or distance calculations on resistance profiles.

[AMR Extensions for Tidymodels](https://amr-for-r.org/reference/amr-tidymodels.md)
:   [`all_sir()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    [`all_sir_predictors()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    [`all_mic()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    [`all_mic_predictors()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    [`all_disk()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    [`all_disk_predictors()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    [`step_mic_log2()`](https://amr-for-r.org/reference/amr-tidymodels.md)
    [`step_sir_numeric()`](https://amr-for-r.org/reference/amr-tidymodels.md)

[Principal Component Analysis (for AMR)](https://amr-for-r.org/reference/pca.md)
:   [`pca()`](https://amr-for-r.org/reference/pca.md)

[Calculate the Mean AMR Distance](https://amr-for-r.org/reference/mean_amr_distance.md)
:   [`mean_amr_distance()`](https://amr-for-r.org/reference/mean_amr_distance.md)
    [`amr_distance_from_row()`](https://amr-for-r.org/reference/mean_amr_distance.md)

## Reference data

The data sets underlying this package, all of which can also be
[downloaded in other
formats](https://amr-for-r.org/articles/datasets.md).

### Microorganisms

[Data Set with 103 222 Taxonomic Records of Microorganisms](https://amr-for-r.org/reference/microorganisms.md)
:   [`microorganisms`](https://amr-for-r.org/reference/microorganisms.md)

[Data Set with 6 293 Common Microorganism Codes](https://amr-for-r.org/reference/microorganisms.codes.md)
:   [`microorganisms.codes`](https://amr-for-r.org/reference/microorganisms.codes.md)

[Data Set with 466 Microorganisms In Species Groups](https://amr-for-r.org/reference/microorganisms.groups.md)
:   [`microorganisms.groups`](https://amr-for-r.org/reference/microorganisms.groups.md)

### Antimicrobials

[Data Sets with 626 Antimicrobial Drugs](https://amr-for-r.org/reference/antimicrobials.md)
:   [`antimicrobials`](https://amr-for-r.org/reference/antimicrobials.md)
    [`antibiotics`](https://amr-for-r.org/reference/antimicrobials.md)
    [`antivirals`](https://amr-for-r.org/reference/antimicrobials.md)

[WHOCC: WHO Collaborating Centre for Drug Statistics Methodology](https://amr-for-r.org/reference/WHOCC.md)
:   [`WHOCC`](https://amr-for-r.org/reference/WHOCC.md)

[Data Set with Treatment Dosages as Defined by EUCAST](https://amr-for-r.org/reference/dosage.md)
:   [`dosage`](https://amr-for-r.org/reference/dosage.md)

### Breakpoints and intrinsic resistance

[Data Set with Clinical Breakpoints for SIR Interpretation](https://amr-for-r.org/reference/clinical_breakpoints.md)
:   [`clinical_breakpoints`](https://amr-for-r.org/reference/clinical_breakpoints.md)

[Data Set Denoting Bacterial Intrinsic Resistance](https://amr-for-r.org/reference/intrinsic_resistant.md)
:   [`intrinsic_resistant`](https://amr-for-r.org/reference/intrinsic_resistant.md)

### Example data

[Data Set with 2 000 Example Isolates](https://amr-for-r.org/reference/example_isolates.md)
:   [`example_isolates`](https://amr-for-r.org/reference/example_isolates.md)

[Data Set with Unclean Data](https://amr-for-r.org/reference/example_isolates_unclean.md)
:   [`example_isolates_unclean`](https://amr-for-r.org/reference/example_isolates_unclean.md)

[Data Set with 500 ESBL Isolates](https://amr-for-r.org/reference/esbl_isolates.md)
:   [`esbl_isolates`](https://amr-for-r.org/reference/esbl_isolates.md)

[Data Set with 500 Isolates - WHONET Example](https://amr-for-r.org/reference/WHONET.md)
:   [`WHONET`](https://amr-for-r.org/reference/WHONET.md)

## Utilities

Helper functions for pattern matching with `%like%`, ages, joining and
italicising taxonomy, generating random test results, and data export.

[Vectorised Pattern Matching with Keyboard Shortcut](https://amr-for-r.org/reference/like.md)
:   [`like()`](https://amr-for-r.org/reference/like.md)
    [`` `%like%` ``](https://amr-for-r.org/reference/like.md)
    [`` `%unlike%` ``](https://amr-for-r.org/reference/like.md)
    [`` `%like_case%` ``](https://amr-for-r.org/reference/like.md)
    [`` `%unlike_case%` ``](https://amr-for-r.org/reference/like.md)

[Age in Years of Individuals](https://amr-for-r.org/reference/age.md)
:   [`age()`](https://amr-for-r.org/reference/age.md)

[Split Ages into Age Groups](https://amr-for-r.org/reference/age_groups.md)
:   [`age_groups()`](https://amr-for-r.org/reference/age_groups.md)

[Join microorganisms to a Data Set](https://amr-for-r.org/reference/join.md)
:   [`inner_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
    [`left_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
    [`right_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
    [`full_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
    [`semi_join_microorganisms()`](https://amr-for-r.org/reference/join.md)
    [`anti_join_microorganisms()`](https://amr-for-r.org/reference/join.md)

[Italicise Taxonomic Families, Genera, Species, and Subspecies](https://amr-for-r.org/reference/italicise_taxonomy.md)
:   [`italicise_taxonomy()`](https://amr-for-r.org/reference/italicise_taxonomy.md)
    [`italicize_taxonomy()`](https://amr-for-r.org/reference/italicise_taxonomy.md)

[Generate Random MIC Values, Disk Zones, or SIR Results](https://amr-for-r.org/reference/random.md)
:   [`random_mic()`](https://amr-for-r.org/reference/random.md)
    [`random_disk()`](https://amr-for-r.org/reference/random.md)
    [`random_sir()`](https://amr-for-r.org/reference/random.md)

[Export Data Set as NCBI BioSample Antibiogram](https://amr-for-r.org/reference/export_ncbi_biosample.md)
:   [`export_ncbi_biosample()`](https://amr-for-r.org/reference/export_ncbi_biosample.md)

### Statistical tests

[*G*-test for Count Data](https://amr-for-r.org/reference/g.test.md)
:   [`g.test()`](https://amr-for-r.org/reference/g.test.md)

[Kurtosis of the Sample](https://amr-for-r.org/reference/kurtosis.md)
:   [`kurtosis()`](https://amr-for-r.org/reference/kurtosis.md)

[Skewness of the Sample](https://amr-for-r.org/reference/skewness.md)
:   [`skewness()`](https://amr-for-r.org/reference/skewness.md)

## Deprecated functions

These still work, but show a warning and will be removed in a future
version.

[Deprecated Functions, Arguments, or Data Sets](https://amr-for-r.org/reference/AMR-deprecated.md)
:   [`custom_eucast_rules()`](https://amr-for-r.org/reference/AMR-deprecated.md)
    [`ab_class()`](https://amr-for-r.org/reference/AMR-deprecated.md)
    [`ab_selector()`](https://amr-for-r.org/reference/AMR-deprecated.md)

[Predict Antimicrobial Resistance](https://amr-for-r.org/reference/resistance_predict.md)
:   [`resistance_predict()`](https://amr-for-r.org/reference/resistance_predict.md)
    [`sir_predict()`](https://amr-for-r.org/reference/resistance_predict.md)
    [`plot(`*`<resistance_predict>`*`)`](https://amr-for-r.org/reference/resistance_predict.md)
    [`ggplot_sir_predict()`](https://amr-for-r.org/reference/resistance_predict.md)
    [`autoplot(`*`<resistance_predict>`*`)`](https://amr-for-r.org/reference/resistance_predict.md)
