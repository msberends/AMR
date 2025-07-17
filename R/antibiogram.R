# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

#' Generate Traditional, Combination, Syndromic, or WISCA Antibiograms
#'
#' @description
#' Create detailed antibiograms with options for traditional, combination, syndromic, and Bayesian WISCA methods.
#'
#' Adhering to previously described approaches (see *Source*) and especially the Bayesian WISCA model (Weighted-Incidence Syndromic Combination Antibiogram) by Bielicki *et al.*, these functions provide flexible output formats including plots and tables, ideal for integration with R Markdown and Quarto reports.
#' @param x A [data.frame] containing at least a column with microorganisms and columns with antimicrobial results (class 'sir', see [as.sir()]).
#' @param antimicrobials A vector specifying the antimicrobials containing SIR values to include in the antibiogram (see *Examples*). Will be evaluated using [guess_ab_col()]. This can be:
#'   - Any antimicrobial name or code that could match (see [guess_ab_col()]) to any column in `x`
#'   - Any [antimicrobial selector][antimicrobial_selectors], such as [aminoglycosides()] or [carbapenems()]
#'   - A combination of the above, using `c()`, e.g.:
#'     - `c(aminoglycosides(), "AMP", "AMC")`
#'     - `c(aminoglycosides(), carbapenems())`
#'   - Column indices using numbers
#'   - Combination therapy, indicated by using `"+"`, with or without [antimicrobial selectors][antimicrobial_selectors], e.g.:
#'     - `"cipro + genta"`
#'     - `"TZP+TOB"`
#'     - `c("TZP", "TZP+GEN", "TZP+TOB")`
#'     - `carbapenems() + "GEN"`
#'     - `carbapenems() + c("", "GEN")`
#'     - `carbapenems() + c("", aminoglycosides())`
#' @param mo_transform A character to transform microorganism input - must be `"name"`, `"shortname"` (default), `"gramstain"`, or one of the column names of the [microorganisms] data set: `r vector_or(colnames(microorganisms), sort = FALSE, quotes = TRUE)`. Can also be `NULL` to not transform the input or `NA` to consider all microorganisms 'unknown'.
#' @param ab_transform A character to transform antimicrobial input - must be one of the column names of the [antimicrobials] data set (defaults to `"name"`): `r vector_or(colnames(antimicrobials), sort = FALSE, quotes = TRUE)`. Can also be `NULL` to not transform the input.
#' @param syndromic_group A column name of `x`, or values calculated to split rows of `x`, e.g. by using [ifelse()] or [`case_when()`][dplyr::case_when()]. See *Examples*.
#' @param add_total_n *(deprecated in favour of `formatting_type`)* A [logical] to indicate whether `n_tested` available numbers per pathogen should be added to the table (default is `TRUE`). This will add the lowest and highest number of available isolates per antimicrobial (e.g, if for *E. coli* 200 isolates are available for ciprofloxacin and 150 for amoxicillin, the returned number will be "150-200"). This option is unavailable when `wisca = TRUE`; in that case, use [retrieve_wisca_parameters()] to get the parameters used for WISCA.
#' @param only_all_tested (for combination antibiograms): a [logical] to indicate that isolates must be tested for all antimicrobials, see *Details*.
#' @param digits Number of digits to use for rounding the antimicrobial coverage, defaults to 1 for WISCA and 0 otherwise.
#' @param formatting_type Numeric value (1–22 for WISCA, 1-12 for non-WISCA) indicating how the 'cells' of the antibiogram table should be formatted. See *Details* > *Formatting Type* for a list of options.
#' @param col_mo Column name of the names or codes of the microorganisms (see [as.mo()]) - the default is the first column of class [`mo`]. Values will be coerced using [as.mo()].
#' @param language Language to translate text, which defaults to the system language (see [get_AMR_locale()]).
#' @param minimum The minimum allowed number of available (tested) isolates. Any isolate count lower than `minimum` will return `NA` with a warning. The default number of `30` isolates is advised by the Clinical and Laboratory Standards Institute (CLSI) as best practice, see *Source*.
#' @param combine_SI A [logical] to indicate whether all susceptibility should be determined by results of either S, SDD, or I, instead of only S (default is `TRUE`).
#' @param sep A separating character for antimicrobial columns in combination antibiograms.
#' @param sort_columns A [logical] to indicate whether the antimicrobial columns must be sorted on name.
#' @param wisca A [logical] to indicate whether a Weighted-Incidence Syndromic Combination Antibiogram (WISCA) must be generated (default is `FALSE`). This will use a Bayesian decision model to estimate regimen coverage probabilities using [Monte Carlo simulations](https://en.wikipedia.org/wiki/Monte_Carlo_method). Set `simulations`, `conf_interval`, and `interval_side` to adjust.
#' @param simulations (for WISCA) a numerical value to set the number of Monte Carlo simulations.
#' @param conf_interval A numerical value to set confidence interval (default is `0.95`).
#' @param interval_side The side of the confidence interval, either `"two-tailed"` (default), `"left"` or `"right"`.
#' @param info A [logical] to indicate info should be printed - the default is `TRUE` only in interactive mode.
#' @param object An [antibiogram()] object.
#' @param ... When used in [R Markdown or Quarto][knitr::kable()]: arguments passed on to [knitr::kable()] (otherwise, has no use).
#' @details These functions return a table with values between 0 and 100 for *susceptibility*, not resistance.
#'
#' **Remember that you should filter your data to let it contain only first isolates!** This is needed to exclude duplicates and to reduce selection bias. Use [first_isolate()] to determine them with one of the four available algorithms: isolate-based, patient-based, episode-based, or phenotype-based.
#'
#' For estimating antimicrobial coverage, especially when creating a WISCA, the outcome might become more reliable by only including the top *n* species encountered in the data. You can filter on this top *n* using [top_n_microorganisms()]. For example, use `top_n_microorganisms(your_data, n = 10)` as a pre-processing step to only include the top 10 species in the data.
#'
#' The numeric values of an antibiogram are stored in a long format as the [attribute][attributes()] `long_numeric`. You can retrieve them using `attributes(x)$long_numeric`, where `x` is the outcome of [antibiogram()] or [wisca()]. This is ideal for e.g. advanced plotting.
#'
#' ### Formatting Type
#'
#' The formatting of the 'cells' of the table can be set with the argument `formatting_type`. In these examples, `5` indicates the antimicrobial coverage (`4-6` the confidence level), `15` the number of susceptible isolates, and `300` the number of tested (i.e., available) isolates:
#'
#' 1. 5
#' 2. 15
#' 3. 300
#' 4. 15/300
#' 5. 5 (300)
#' 6. 5% (300)
#' 7. 5 (N=300)
#' 8. 5% (N=300)
#' 9. 5 (15/300)
#' 10. 5% (15/300)
#' 11. 5 (N=15/300)
#' 12. 5% (N=15/300)
#' 13. 5 (4-6)
#' 14. 5% (4-6%) - **default for WISCA**
#' 15. 5 (4-6,300)
#' 16. 5% (4-6%,300)
#' 17. 5 (4-6,N=300)
#' 18. 5% (4-6%,N=300) - **default for non-WISCA**
#' 19. 5 (4-6,15/300)
#' 20. 5% (4-6%,15/300)
#' 21. 5 (4-6,N=15/300)
#' 22. 5% (4-6%,N=15/300)
#'
#' The default can be set globally with the package option [`AMR_antibiogram_formatting_type`][AMR-options], e.g. `options(AMR_antibiogram_formatting_type = 5)`. Do note that for WISCA, the total numbers of tested and susceptible isolates are less useful to report, since these are included in the Bayesian model and apparent from the susceptibility and its confidence level.
#'
#' Set `digits` (defaults to `0`) to alter the rounding of the susceptibility percentages.
#'
#' ### Antibiogram Types
#'
#' There are various antibiogram types, as summarised by Klinker *et al.* (2021, \doi{10.1177/20499361211011373}), and they are all supported by [antibiogram()].
#'
#' For clinical coverage estimations, **use WISCA whenever possible**, since it provides more precise coverage estimates by accounting for pathogen incidence and antimicrobial susceptibility, as has been shown by Bielicki *et al.* (2020, \doi{10.1001/jamanetworkopen.2019.21124}). See the section *Explaining WISCA* on this page. Do note that WISCA is pathogen-agnostic, meaning that the outcome is not stratied by pathogen, but rather by syndrome.
#'
#' 1. **Traditional Antibiogram**
#'
#'    Case example: Susceptibility of *Pseudomonas aeruginosa* to piperacillin/tazobactam (TZP)
#'
#'    Code example:
#'
#'    ```r
#'    antibiogram(your_data,
#'                antimicrobials = "TZP")
#'    ```
#'
#' 2. **Combination Antibiogram**
#'
#'    Case example: Additional susceptibility of *Pseudomonas aeruginosa* to TZP + tobramycin versus TZP alone
#'
#'    Code example:
#'
#'    ```r
#'    antibiogram(your_data,
#'                antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"))
#'    ```
#'
#' 3. **Syndromic Antibiogram**
#'
#'    Case example: Susceptibility of *Pseudomonas aeruginosa* to TZP among respiratory specimens (obtained among ICU patients only)
#'
#'    Code example:
#'
#'    ```r
#'    antibiogram(your_data,
#'                antimicrobials = penicillins(),
#'                syndromic_group = "ward")
#'    ```
#'
#' 4. **Weighted-Incidence Syndromic Combination Antibiogram (WISCA)**
#'
#'    WISCA can be applied to any antibiogram, see the section *Explaining WISCA* on this page for more information.
#'
#'    Code example:
#'
#'    ```r
#'    antibiogram(your_data,
#'                antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
#'                wisca = TRUE)
#'
#'    # this is equal to:
#'    wisca(your_data,
#'          antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"))
#'    ```
#'
#'    WISCA uses a sophisticated Bayesian decision model to combine both local and pooled antimicrobial resistance data. This approach not only evaluates local patterns but can also draw on multi-centre datasets to improve regimen accuracy, even in low-incidence infections like paediatric bloodstream infections (BSIs).
#'
#' ### Grouped tibbles
#'
#' For any type of antibiogram, grouped [tibbles][tibble::tibble] can also be used to calculate susceptibilities over various groups.
#'
#' Code example:
#'
#' ```r
#' library(dplyr)
#' your_data %>%
#'   group_by(has_sepsis, is_neonate, sex) %>%
#'   wisca(antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"))
#' ```
#'
#' ### Stepped Approach for Clinical Insight
#'
#' In clinical practice, antimicrobial coverage decisions evolve as more microbiological data becomes available. This theoretical stepped approach ensures empirical coverage can continuously assessed to improve patient outcomes:
#'
#' 1. **Initial Empirical Therapy (Admission / Pre-Culture Data)**
#'
#'    At admission, no pathogen information is available.
#'
#'    - Action: broad-spectrum coverage is based on local resistance patterns and syndromic antibiograms. Using the pathogen-agnostic yet incidence-weighted WISCA is preferred.
#'    - Code example:
#'
#'      ```r
#'      antibiogram(your_data,
#'                  antimicrobials = selected_regimens,
#'                  mo_transform = NA) # all pathogens set to `NA`
#'
#'      # preferred: use WISCA
#'      wisca(your_data,
#'            antimicrobials = selected_regimens)
#'      ```
#'
#' 2. **Refinement with Gram Stain Results**
#'
#'    When a blood culture becomes positive, the Gram stain provides an initial and crucial first stratification (Gram-positive vs. Gram-negative).
#'
#'    - Action: narrow coverage based on Gram stain-specific resistance patterns.
#'    - Code example:
#'
#'      ```r
#'      antibiogram(your_data,
#'                  antimicrobials = selected_regimens,
#'                  mo_transform = "gramstain") # all pathogens set to Gram-pos/Gram-neg
#'      ```
#'
#' 3. **Definitive Therapy Based on Species Identification**
#'
#'    After cultivation of the pathogen, full pathogen identification allows precise targeting of therapy.
#'
#'    - Action: adjust treatment to pathogen-specific antibiograms, minimizing resistance risks.
#'    - Code example:
#'
#'      ```r
#'      antibiogram(your_data,
#'                  antimicrobials = selected_regimens,
#'                  mo_transform = "shortname") # all pathogens set to 'G. species', e.g., E. coli
#'      ```
#'
#' By structuring antibiograms around this stepped approach, clinicians can make data-driven adjustments at each stage, ensuring optimal empirical and targeted therapy while reducing unnecessary broad-spectrum antimicrobial use.
#'
#' ### Inclusion in Combination Antibiograms
#'
#' Note that for combination antibiograms, it is important to realise that susceptibility can be calculated in two ways, which can be set with the `only_all_tested` argument (default is `FALSE`). See this example for two antimicrobials, Drug A and Drug B, about how [antibiogram()] works to calculate the %SI:
#'
#' ```
#' --------------------------------------------------------------------
#'                     only_all_tested = FALSE  only_all_tested = TRUE
#'                     -----------------------  -----------------------
#'  Drug A    Drug B   considered   considered  considered   considered
#'                     susceptible    tested    susceptible    tested
#' --------  --------  -----------  ----------  -----------  ----------
#'  S or I    S or I        X            X           X            X
#'    R       S or I        X            X           X            X
#'   <NA>     S or I        X            X           -            -
#'  S or I      R           X            X           X            X
#'    R         R           -            X           -            X
#'   <NA>       R           -            -           -            -
#'  S or I     <NA>         X            X           -            -
#'    R        <NA>         -            -           -            -
#'   <NA>      <NA>         -            -           -            -
#' --------------------------------------------------------------------
#' ```
#'
#' ### Plotting
#'
#' All types of antibiograms as listed above can be plotted (using [ggplot2::autoplot()] or base \R's [plot()] and [barplot()]). As mentioned above, the numeric values of an antibiogram are stored in a long format as the [attribute][attributes()] `long_numeric`. You can retrieve them using `attributes(x)$long_numeric`, where `x` is the outcome of [antibiogram()] or [wisca()].
#'
#' The outcome of [antibiogram()] can also be used directly in R Markdown / Quarto (i.e., `knitr`) for reports. In this case, [knitr::kable()] will be applied automatically and microorganism names will even be printed in italics at default (see argument `italicise`).
#'
#' You can also use functions from specific 'table reporting' packages to transform the output of [antibiogram()] to your needs, e.g. with `flextable::as_flextable()` or `gt::gt()`.
#'
#' @section Explaining WISCA:
#'
#' WISCA (Weighted-Incidence Syndromic Combination Antibiogram) estimates the probability of empirical coverage for combination regimens.
#'
#' It weights susceptibility by pathogen prevalence within a clinical syndrome and provides credible intervals around the expected coverage.
#'
#' For more background, interpretation, and examples, see [the WISCA vignette](https://amr-for-r.org/articles/WISCA.html).
#' @source
#' * Bielicki JA *et al.* (2016). **Selecting appropriate empirical antibiotic regimens for paediatric bloodstream infections: application of a Bayesian decision model to local and pooled antimicrobial resistance surveillance data** *Journal of Antimicrobial Chemotherapy* 71(3); \doi{10.1093/jac/dkv397}
#' * Bielicki JA *et al.* (2020). **Evaluation of the coverage of 3 antibiotic regimens for neonatal sepsis in the hospital setting across Asian countries** *JAMA Netw Open.* 3(2):e1921124; \doi{10.1001/jamanetworkopen.2019.21124}
#' * Klinker KP *et al.* (2021). **Antimicrobial stewardship and antibiograms: importance of moving beyond traditional antibiograms**. *Therapeutic Advances in Infectious Disease*, May 5;8:20499361211011373; \doi{10.1177/20499361211011373}
#' * Barbieri E *et al.* (2021). **Development of a Weighted-Incidence Syndromic Combination Antibiogram (WISCA) to guide the choice of the empiric antibiotic treatment for urinary tract infection in paediatric patients: a Bayesian approach** *Antimicrobial Resistance & Infection Control* May 1;10(1):74; \doi{10.1186/s13756-021-00939-2}
#' * **M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 5th Edition**, 2022, *Clinical and Laboratory Standards Institute (CLSI)*. <https://clsi.org/standards/products/microbiology/documents/m39/>.
#' @author Implementation: Dr. Larisse Bolton and Dr. Matthijs Berends
#' @rdname antibiogram
#' @name antibiogram
#' @export
#' @examples
#' # example_isolates is a data set available in the AMR package.
#' # run ?example_isolates for more info.
#' example_isolates
#'
#' \donttest{
#' # Traditional antibiogram ----------------------------------------------
#'
#' antibiogram(example_isolates,
#'   antimicrobials = c(aminoglycosides(), carbapenems())
#' )
#'
#' antibiogram(example_isolates,
#'   antimicrobials = aminoglycosides(),
#'   ab_transform = "atc",
#'   mo_transform = "gramstain"
#' )
#'
#' antibiogram(example_isolates,
#'   antimicrobials = carbapenems(),
#'   ab_transform = "name",
#'   mo_transform = "name"
#' )
#'
#'
#' # Combined antibiogram -------------------------------------------------
#'
#' # combined antimicrobials yield higher empiric coverage
#' antibiogram(example_isolates,
#'   antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
#'   mo_transform = "gramstain"
#' )
#'
#' # you can use any antimicrobial selector with `+` too:
#' antibiogram(example_isolates,
#'   antimicrobials = ureidopenicillins() + c("", "GEN", "tobra"),
#'   mo_transform = "gramstain"
#' )
#'
#' # names of antimicrobials do not need to resemble columns exactly:
#' antibiogram(example_isolates,
#'   antimicrobials = c("Cipro", "cipro + genta"),
#'   mo_transform = "gramstain",
#'   ab_transform = "name",
#'   sep = " & "
#' )
#'
#'
#' # Syndromic antibiogram ------------------------------------------------
#'
#' # the data set could contain a filter for e.g. respiratory specimens
#' antibiogram(example_isolates,
#'   antimicrobials = c(aminoglycosides(), carbapenems()),
#'   syndromic_group = "ward"
#' )
#'
#' # now define a data set with only E. coli
#' ex1 <- example_isolates[which(mo_genus() == "Escherichia"), ]
#'
#' # with a custom language, though this will be determined automatically
#' # (i.e., this table will be in Spanish on Spanish systems)
#' antibiogram(ex1,
#'   antimicrobials = aminoglycosides(),
#'   ab_transform = "name",
#'   syndromic_group = ifelse(ex1$ward == "ICU",
#'     "UCI", "No UCI"
#'   ),
#'   language = "es"
#' )
#'
#'
#' # WISCA antibiogram ----------------------------------------------------
#'
#' # WISCA are not stratified by species, but rather on syndromes
#' antibiogram(example_isolates,
#'   antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
#'   syndromic_group = "ward",
#'   wisca = TRUE
#' )
#'
#'
#' # Print the output for R Markdown / Quarto -----------------------------
#'
#' ureido <- antibiogram(example_isolates,
#'   antimicrobials = ureidopenicillins(),
#'   syndromic_group = "ward",
#'   wisca = TRUE
#' )
#'
#' # in an Rmd file, you would just need to return `ureido` in a chunk,
#' # but to be explicit here:
#' if (requireNamespace("knitr")) {
#'   cat(knitr::knit_print(ureido))
#' }
#'
#'
#' # Generate plots with ggplot2 or base R --------------------------------
#'
#' ab1 <- antibiogram(example_isolates,
#'   antimicrobials = c("AMC", "CIP", "TZP", "TZP+TOB"),
#'   mo_transform = "gramstain"
#' )
#' ab2 <- antibiogram(example_isolates,
#'   antimicrobials = c("AMC", "CIP", "TZP", "TZP+TOB"),
#'   mo_transform = "gramstain",
#'   syndromic_group = "ward"
#' )
#'
#' if (requireNamespace("ggplot2")) {
#'   ggplot2::autoplot(ab1)
#' }
#' if (requireNamespace("ggplot2")) {
#'   ggplot2::autoplot(ab2)
#' }
#'
#' plot(ab1)
#' plot(ab2)
#' }
antibiogram <- function(x,
                        antimicrobials = where(is.sir),
                        mo_transform = "shortname",
                        ab_transform = "name",
                        syndromic_group = NULL,
                        add_total_n = FALSE,
                        only_all_tested = FALSE,
                        digits = ifelse(wisca, 1, 0),
                        formatting_type = getOption("AMR_antibiogram_formatting_type", ifelse(wisca, 14, 18)),
                        col_mo = NULL,
                        language = get_AMR_locale(),
                        minimum = 30,
                        combine_SI = TRUE,
                        sep = " + ",
                        sort_columns = TRUE,
                        wisca = FALSE,
                        simulations = 1000,
                        conf_interval = 0.95,
                        interval_side = "two-tailed",
                        info = interactive(),
                        ...) {
  UseMethod("antibiogram")
}

#' @method antibiogram default
#' @export
antibiogram.default <- function(x,
                                antimicrobials = where(is.sir),
                                mo_transform = "shortname",
                                ab_transform = "name",
                                syndromic_group = NULL,
                                add_total_n = FALSE,
                                only_all_tested = FALSE,
                                digits = ifelse(wisca, 1, 0),
                                formatting_type = getOption("AMR_antibiogram_formatting_type", ifelse(wisca, 14, 18)),
                                col_mo = NULL,
                                language = get_AMR_locale(),
                                minimum = 30,
                                combine_SI = TRUE,
                                sep = " + ",
                                sort_columns = TRUE,
                                wisca = FALSE,
                                simulations = 1000,
                                conf_interval = 0.95,
                                interval_side = "two-tailed",
                                info = interactive(),
                                ...) {
  meet_criteria(x, allow_class = "data.frame")
  x <- ascertain_sir_classes(x, "x")
  meet_criteria(wisca, allow_class = "logical", has_length = 1)
  if (isTRUE(wisca)) {
    if (!is.null(mo_transform) && !missing(mo_transform)) {
      warning_("WISCA must be based on the species level as WISCA parameters are based on this. For that reason, `mo_transform` will be ignored.")
    }
    mo_transform <- function(x) suppressMessages(suppressWarnings(paste(mo_genus(x, keep_synonyms = TRUE, language = NULL), mo_species(x, keep_synonyms = TRUE, language = NULL))))
  }
  if ("antibiotics" %in% names(list(...))) {
    deprecation_warning("antibiotics", "antimicrobials", fn = "antibiogram", is_argument = TRUE)
    antimicrobials <- list(...)$antibiotics
  }
  meet_criteria(antimicrobials, allow_class = c("character", "numeric", "integer"), allow_NA = FALSE, allow_NULL = FALSE)
  if (!is.function(mo_transform)) {
    meet_criteria(mo_transform, allow_class = "character", has_length = 1, is_in = c("name", "shortname", "gramstain", colnames(AMR::microorganisms)), allow_NULL = TRUE, allow_NA = TRUE)
  }
  if (!is.function(ab_transform)) {
    meet_criteria(ab_transform, allow_class = "character", has_length = 1, is_in = colnames(AMR::antimicrobials), allow_NULL = TRUE)
  }
  meet_criteria(syndromic_group, allow_class = "character", allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(add_total_n, allow_class = "logical", has_length = 1)
  if (isTRUE(add_total_n)) {
    deprecation_warning("add_total_n", "formatting_type", fn = "antibiogram", is_argument = TRUE)
  }
  meet_criteria(only_all_tested, allow_class = "logical", has_length = 1)
  meet_criteria(digits, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE)
  meet_criteria(formatting_type, allow_class = c("numeric", "integer"), has_length = 1, is_in = c(1:22))
  meet_criteria(col_mo, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  language <- validate_language(language)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(sep, allow_class = "character", has_length = 1)
  meet_criteria(sort_columns, allow_class = "logical", has_length = 1)
  meet_criteria(simulations, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE, is_positive = TRUE)
  meet_criteria(conf_interval, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE, is_positive = TRUE)
  meet_criteria(interval_side, allow_class = "character", has_length = 1, is_in = c("two-tailed", "left", "right"))
  meet_criteria(info, allow_class = "logical", has_length = 1)

  # try to find columns based on type
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = info)
    stop_if(is.null(col_mo), "`col_mo` must be set")
  }
  # transform MOs
  x$`.mo` <- x[, col_mo, drop = TRUE]
  if (is.null(mo_transform)) {
    # leave as is, no transformation
  } else if (is.function(mo_transform)) {
    x$`.mo` <- mo_transform(x$`.mo`)
  } else if (is.na(mo_transform)) {
    x$`.mo` <- NA_character_
  } else if (mo_transform == "gramstain") {
    x$`.mo` <- mo_gramstain(x$`.mo`, language = language)
  } else if (mo_transform == "shortname") {
    x$`.mo` <- mo_shortname(x$`.mo`, language = language)
  } else if (mo_transform == "name") {
    x$`.mo` <- mo_name(x$`.mo`, language = language)
  } else {
    x$`.mo` <- mo_property(x$`.mo`, property = mo_transform, language = language)
  }
  x$`.mo`[is.na(x$`.mo`)] <- "(??)"

  # get syndromic groups
  if (!is.null(syndromic_group)) {
    if (length(syndromic_group) == 1 && syndromic_group %in% colnames(x)) {
      x$`.syndromic_group` <- x[, syndromic_group, drop = TRUE]
    } else if (!is.null(syndromic_group)) {
      x$`.syndromic_group` <- syndromic_group
    }
    x$`.syndromic_group`[is.na(x$`.syndromic_group`) | x$`.syndromic_group` == ""] <- paste0("(", translate_AMR("unknown", language = language), ")")
    has_syndromic_group <- TRUE
  } else {
    has_syndromic_group <- FALSE
  }

  # get antimicrobials
  ab_trycatch <- tryCatch(colnames(suppressWarnings(x[, antimicrobials, drop = FALSE])), error = function(e) NULL)
  if (is.null(ab_trycatch)) {
    stop_ifnot(is.character(suppressMessages(antimicrobials)), "`antimicrobials` must be an antimicrobial selector, or a character vector.")
    antimicrobials.bak <- antimicrobials
    # split antimicrobials on separator and make it a list
    antimicrobials <- strsplit(gsub(" ", "", antimicrobials), "+", fixed = TRUE)
    # get available antimicrobials in data set
    df_ab <- get_column_abx(x, verbose = FALSE, info = FALSE)
    # get antimicrobials from user
    user_ab <- suppressMessages(suppressWarnings(lapply(antimicrobials, as.ab, flag_multiple_results = FALSE, info = FALSE)))
    non_existing <- character(0)
    user_ab <- lapply(user_ab, function(x) {
      out <- unname(df_ab[match(x, names(df_ab))])
      non_existing <<- c(non_existing, x[is.na(out) & !is.na(x)])
      # remove non-existing columns
      out[!is.na(out)]
    })
    user_ab <- user_ab[unlist(lapply(user_ab, length)) > 0]

    if (length(non_existing) > 0) {
      warning_("The following antimicrobials were not available and ignored: ", vector_and(ab_name(non_existing, language = NULL, tolower = TRUE), quotes = FALSE))
    }

    # make list unique
    antimicrobials <- unique(user_ab)
    # go through list to set AMR in combinations
    for (i in seq_along(antimicrobials)) {
      abx <- antimicrobials[[i]]
      for (ab in abx) {
        # make sure they are SIR columns
        x[, ab] <- as.sir(x[, ab, drop = TRUE])
        # set NI as NA
        x[[ab]][x[[ab]] == "NI"] <- NA_sir_
      }
      new_colname <- paste0(trimws(abx), collapse = sep)
      if (length(abx) == 1) {
        next
      } else {
        # determine whether this new column should contain S, I, R, or NA
        if (isTRUE(combine_SI)) {
          S_values <- c("S", "SDD", "I")
        } else {
          S_values <- "S"
        }
        other_values <- setdiff(c("S", "SDD", "I", "R"), S_values)
        x_transposed <- as.list(as.data.frame(t(x[, abx, drop = FALSE]), stringsAsFactors = FALSE))
        if (isTRUE(only_all_tested)) {
          x[new_colname] <- as.sir(vapply(FUN.VALUE = character(1), x_transposed, function(x) ifelse(anyNA(x), NA_character_, ifelse(any(x %in% S_values), "S", "R")), USE.NAMES = FALSE))
        } else {
          x[new_colname] <- as.sir(vapply(
            FUN.VALUE = character(1), x_transposed, function(x) ifelse(any(x %in% S_values, na.rm = TRUE), "S", ifelse(anyNA(x), NA_character_, "R")),
            USE.NAMES = FALSE
          ))
        }
      }
      antimicrobials[[i]] <- new_colname
    }
    antimicrobials <- unlist(antimicrobials)
  } else {
    existing_ab_combined_cols <- ab_trycatch[ab_trycatch %like% "[+]" & ab_trycatch %in% colnames(x)]
    if (length(existing_ab_combined_cols) > 0 && !is.null(ab_transform)) {
      ab_transform <- NULL
      warning_(
        "Detected column name(s) containing the '+' character, which conflicts with the expected syntax in `antibiogram()`: the '+' is used to combine separate antimicrobial agent columns (e.g., \"AMP+GEN\").\n\n",
        "To avoid incorrectly guessing which antimicrobials this represents, `ab_transform` was automatically set to `NULL`.\n\n",
        "If this is unintended, please rename the column(s) to avoid using '+' in the name, or set `ab_transform = NULL` explicitly to suppress this message."
      )
    }
    antimicrobials <- ab_trycatch
  }

  if (isTRUE(has_syndromic_group)) {
    out <- x %pm>%
      pm_select(.syndromic_group, .mo, antimicrobials) %pm>%
      pm_group_by(.syndromic_group)
  } else {
    out <- x %pm>%
      pm_select(.mo, antimicrobials)
  }

  # get numbers of S, I, R (per group)
  out <- out %pm>%
    bug_drug_combinations(
      col_mo = ".mo",
      FUN = function(x) x,
      include_n_rows = TRUE
    )
  colnames(out)[colnames(out) == "total"] <- "n_tested"
  colnames(out)[colnames(out) == "total_rows"] <- "n_total"
  out$ab <- factor(out$ab, levels = antimicrobials, ordered = TRUE)
  out <- out[order(out$mo, out$ab), , drop = FALSE]

  counts <- out

  if (isTRUE(combine_SI)) {
    out$n_susceptible <- out$S + out$I + out$SDD
  } else {
    out$n_susceptible <- out$S
  }
  if (all(out$n_tested < minimum, na.rm = TRUE) && wisca == FALSE) {
    warning_("All combinations had less than `minimum = ", minimum, "` results, returning an empty antibiogram")
    return(as_original_data_class(data.frame(), class(x), extra_class = "antibiogram"))
  } else if (any(out$n_tested < minimum, na.rm = TRUE)) {
    mins <- sum(out$n_tested < minimum, na.rm = TRUE)
    if (wisca == FALSE) {
      out <- out %pm>%
        subset(n_tested >= minimum)
      if (isTRUE(info) && mins > 0) {
        message_("NOTE: ", mins, " combinations had less than `minimum = ", minimum, "` results and were ignored", add_fn = font_red)
      }
    }
  }
  if (NROW(out) == 0) {
    return(as_original_data_class(data.frame(), class(x), extra_class = "antibiogram"))
  }

  out$p_susceptible <- out$n_susceptible / out$n_tested

  # add confidence levels
  out$lower_ci <- NA_real_
  out$upper_ci <- NA_real_
  for (r in seq_len(NROW(out))) {
    if (!is.na(out$n_susceptible[r]) && !is.na(out$n_tested[r]) && out$n_tested[r] > 0) {
      ci <- stats::binom.test(out$n_susceptible[r], out$n_tested[r], conf.level = conf_interval)$conf.int
      out$lower_ci[r] <- ci[1]
      out$upper_ci[r] <- ci[2]
    }
  }

  # regroup for summarising
  if (isTRUE(has_syndromic_group)) {
    colnames(out)[1] <- "syndromic_group"
    out <- out %pm>%
      pm_group_by(syndromic_group, mo, ab)
  } else {
    out <- out %pm>%
      pm_group_by(mo, ab)
  }

  long_numeric <- out %pm>%
    pm_summarise(
      coverage = p_susceptible,
      lower_ci = lower_ci,
      upper_ci = upper_ci,
      n_total = n_total,
      n_tested = n_tested,
      n_susceptible = n_susceptible
    )

  wisca_parameters <- data.frame()

  # WISCA START
  if (wisca == TRUE) {
    if (isTRUE(has_syndromic_group)) {
      colnames(out)[1] <- "syndromic_group"
      out_wisca <- out %pm>%
        pm_group_by(syndromic_group, ab)
    } else {
      out_wisca <- out %pm>%
        pm_group_by(ab)
    }

    out_wisca <- out_wisca %pm>%
      pm_summarise(
        coverage = NA_real_,
        lower_ci = NA_real_,
        upper_ci = NA_real_,
        n_total = sum(n_total, na.rm = TRUE),
        n_tested = sum(n_tested, na.rm = TRUE),
        n_susceptible = sum(n_susceptible, na.rm = TRUE)
      )

    if (any(out_wisca$n_tested < minimum, na.rm = TRUE) && message_not_thrown_before("antibiogram", wisca)) {
      warning_("Number of tested isolates should exceed ", minimum, " for each regimen (and group). WISCA coverage estimates might be inaccurate.", call = FALSE)
    }

    if (isTRUE(has_syndromic_group)) {
      out$group <- paste(out$syndromic_group, out$ab)
      out_wisca$group <- paste(out_wisca$syndromic_group, out_wisca$ab)
    } else {
      out$group <- out$ab
      out_wisca$group <- out_wisca$ab
    }

    wisca_parameters <- out

    progress <- progress_ticker(
      n = length(unique(wisca_parameters$group)) * simulations,
      n_min = 25,
      print = info,
      title = paste("Calculating WISCA for", length(unique(wisca_parameters$group)), "regimens")
    )
    on.exit(close(progress))

    # run WISCA per group
    for (group in unique(wisca_parameters$group)) {
      params_current <- wisca_parameters[wisca_parameters$group == group, , drop = FALSE]
      if (sum(params_current$n_tested, na.rm = TRUE) == 0) {
        next
      }

      # prepare priors
      priors_current <- create_wisca_priors(params_current)

      # Monte Carlo simulations
      coverage_simulations <- vapply(
        FUN.VALUE = double(1),
        seq_len(simulations), function(i) {
          progress$tick()
          simulate_coverage(priors_current)
        }
      )

      # summarise results
      coverage_mean <- mean(coverage_simulations)

      if (interval_side == "two-tailed") {
        probs <- c((1 - conf_interval) / 2, 1 - (1 - conf_interval) / 2)
      } else if (interval_side == "left") {
        probs <- c(0, conf_interval)
      } else if (interval_side == "right") {
        probs <- c(1 - conf_interval, 1)
      }

      coverage_ci <- unname(stats::quantile(coverage_simulations, probs = probs))

      out_wisca$coverage[out_wisca$group == group] <- coverage_mean
      out_wisca$lower_ci[out_wisca$group == group] <- coverage_ci[1]
      out_wisca$upper_ci[out_wisca$group == group] <- coverage_ci[2]
    }

    close(progress)

    # final output preparation
    out <- out_wisca
    wisca_parameters <- wisca_parameters[, colnames(wisca_parameters)[!colnames(wisca_parameters) %in% c(levels(NA_sir_), "lower_ci", "upper_ci", "group")], drop = FALSE]

    if (isTRUE(has_syndromic_group)) {
      long_numeric <- out_wisca %pm>%
        pm_ungroup() %pm>%
        pm_select(
          syndromic_group = syndromic_group,
          ab = ab,
          coverage = coverage,
          lower_ci = lower_ci,
          upper_ci = upper_ci,
          n_total = n_total,
          n_tested = n_tested,
          n_susceptible = n_susceptible
        )
    } else {
      long_numeric <- out_wisca %pm>%
        pm_ungroup() %pm>%
        pm_select(
          ab = ab,
          coverage = coverage,
          lower_ci = lower_ci,
          upper_ci = upper_ci,
          n_total = n_total,
          n_tested = n_tested,
          n_susceptible = n_susceptible
        )
    }
  }

  if (isFALSE(wisca)) {
    out$coverage <- out$p_susceptible
  }

  # formatting type:
  # 1. 5
  # 2. 15
  # 3. 300
  # 4. 15/300
  # 5. 5 (300)
  # 6. 5% (300)
  # 7. 5 (N=300)
  # 8. 5% (N=300)
  # 9. 5 (15/300)
  # 10. 5% (15/300)
  # 11. 5 (N=15/300)
  # 12. 5% (N=15/300)
  # 13. 5 (4-6)
  # 14. 5% (4-6%)
  # 15. 5 (4-6,300)
  # 16. 5% (4-6%,300)
  # 17. 5 (4-6,N=300)
  # 18. 5% (4-6%,N=300)
  # 19. 5 (4-6,15/300)
  # 20. 5% (4-6%,15/300)
  # 21. 5 (4-6,N=15/300)
  # 22. 5% (4-6%,N=15/300)
  if (wisca == TRUE && !formatting_type %in% c(1, 2, 13, 14) && info == TRUE && message_not_thrown_before("antibiogram", wisca, formatting_type)) {
    message_("Using WISCA with a `formatting_type` that includes the denominator is not useful")
  }
  out$digits <- digits # since pm_sumarise() cannot work with an object outside the current frame
  if (formatting_type == 1) out <- out %pm>% pm_summarise(out_value = round(coverage * 100, digits = digits))
  if (formatting_type == 2) out <- out %pm>% pm_summarise(out_value = n_susceptible)
  if (formatting_type == 3) out <- out %pm>% pm_summarise(out_value = n_tested)
  if (formatting_type == 4) out <- out %pm>% pm_summarise(out_value = paste0(n_susceptible, "/", n_tested))
  if (formatting_type == 5) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), " (", n_tested, ")"))
  if (formatting_type == 6) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), "% (", n_tested, ")"))
  if (formatting_type == 7) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), " (N=", n_tested, ")"))
  if (formatting_type == 8) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), "% (N=", n_tested, ")"))
  if (formatting_type == 9) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), " (", n_susceptible, "/", n_tested, ")"))
  if (formatting_type == 10) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), "% (", n_susceptible, "/", n_tested, ")"))
  if (formatting_type == 11) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), " (N=", n_susceptible, "/", n_tested, ")"))
  if (formatting_type == 12) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), "% (N=", n_susceptible, "/", n_tested, ")"))
  if (formatting_type == 13) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), " (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), ")"))
  if (formatting_type == 14) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), "% (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), "%)"))
  if (formatting_type == 15) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), " (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), ",", n_tested, ")"))
  if (formatting_type == 16) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), "% (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), "%,", n_tested, ")"))
  if (formatting_type == 17) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), " (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), ",N=", n_tested, ")"))
  if (formatting_type == 18) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), "% (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), "%,N=", n_tested, ")"))
  if (formatting_type == 19) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), " (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), ",", n_susceptible, "/", n_tested, ")"))
  if (formatting_type == 20) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), "% (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), "%,", n_susceptible, "/", n_tested, ")"))
  if (formatting_type == 21) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), " (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), ",N=", n_susceptible, "/", n_tested, ")"))
  if (formatting_type == 22) out <- out %pm>% pm_summarise(out_value = paste0(round(coverage * 100, digits = digits), "% (", round(lower_ci * 100, digits = digits), "-", round(upper_ci * 100, digits = digits), "%,N=", n_susceptible, "/", n_tested, ")"))
  if (formatting_type >= 4) {
    out$out_value[out$out_value %like% "^NA"] <- NA_character_
  }

  # transform names of antimicrobials
  ab_naming_function <- function(x, t, l, s) {
    x <- strsplit(as.character(x), s, fixed = TRUE)
    out <- character(length = length(x))
    for (i in seq_along(x)) {
      a <- x[[i]]
      if (is.null(t)) {
        # leave as is
      } else if (is.function(t)) {
        a <- t(a)
      } else if (t == "atc") {
        a <- ab_atc(a, only_first = TRUE, language = l)
      } else {
        a <- ab_property(a, property = t, language = l)
      }
      if (length(a) > 1) {
        a <- paste0(trimws(a), collapse = sep)
      }
      out[i] <- a
    }
    out
  }
  out$ab <- ab_naming_function(out$ab, t = ab_transform, l = language, s = sep)
  long_numeric$ab <- ab_naming_function(long_numeric$ab, t = ab_transform, l = language, s = sep)

  # transform long to wide
  long_to_wide <- function(object) {
    if (wisca == TRUE) {
      # column `mo` has already been removed, but we create here a surrogate to make the stats::reshape() work since it needs an identifier
      object$mo <- 1 # seq_len(NROW(object))
    }
    object <- object %pm>%
      # an unclassed data.frame is required for stats::reshape()
      as.data.frame(stringsAsFactors = FALSE) %pm>%
      stats::reshape(direction = "wide", idvar = "mo", timevar = "ab", v.names = "out_value")
    colnames(object) <- gsub("^out_value?[.]", "", colnames(object))
    if (wisca == TRUE) {
      object <- object[, colnames(object)[colnames(object) != "mo"], drop = FALSE]
    }
    return(object)
  }

  # ungroup for long -> wide transformation
  attr(out, "pm_groups") <- NULL
  attr(out, "groups") <- NULL
  class(out) <- class(out)[!class(out) %in% c("grouped_df", "grouped_data")]

  if (isTRUE(sort_columns)) {
    sort_fn <- base::sort
  } else {
    sort_fn <- function(x) x
  }

  if (isTRUE(has_syndromic_group)) {
    grps <- unique(out$syndromic_group)
    for (i in seq_along(grps)) {
      grp <- grps[i]
      if (i == 1) {
        new_df <- long_to_wide(out[which(out$syndromic_group == grp), , drop = FALSE])
      } else {
        new_df <- rbind_AMR(
          new_df,
          long_to_wide(out[which(out$syndromic_group == grp), , drop = FALSE])
        )
      }
    }
    if (wisca == TRUE) {
      # sort rows
      new_df <- new_df %pm>% pm_arrange(syndromic_group)
      # sort columns
      new_df <- new_df[, c("syndromic_group", sort_fn(colnames(new_df)[colnames(new_df) != "syndromic_group"])), drop = FALSE]
      colnames(new_df)[1] <- translate_AMR("Syndromic Group", language = language)
    } else {
      # sort rows
      new_df <- new_df %pm>% pm_arrange(mo, syndromic_group)
      # sort columns
      new_df <- new_df[, c("syndromic_group", "mo", sort_fn(colnames(new_df)[!colnames(new_df) %in% c("syndromic_group", "mo")])), drop = FALSE]
      colnames(new_df)[1:2] <- translate_AMR(c("Syndromic Group", "Pathogen"), language = language)
    }
  } else {
    new_df <- long_to_wide(out)
    if (wisca == TRUE) {
      # sort columns
      new_df <- new_df[, c(sort_fn(colnames(new_df))), drop = FALSE]
    } else {
      # sort rows
      new_df <- new_df %pm>% pm_arrange(mo)
      # sort columns
      new_df <- new_df[, c("mo", sort_fn(colnames(new_df)[colnames(new_df) != "mo"])), drop = FALSE]
      colnames(new_df)[1] <- translate_AMR("Pathogen", language = language)
    }
  }

  # add n_tested N if indicated
  if (isTRUE(add_total_n) && isFALSE(wisca)) {
    if (isTRUE(has_syndromic_group)) {
      n_per_mo <- counts %pm>%
        pm_group_by(mo, .syndromic_group) %pm>%
        pm_summarise(paste0(min(n_tested, na.rm = TRUE), "-", max(n_tested, na.rm = TRUE)))
      colnames(n_per_mo) <- c("mo", "syn", "count")
      count_group <- n_per_mo$count[match(paste(new_df[[2]], new_df[[1]]), paste(n_per_mo$mo, n_per_mo$syn))]
      edit_col <- 2
    } else {
      n_per_mo <- counts %pm>%
        pm_group_by(mo) %pm>%
        pm_summarise(paste0(min(n_tested, na.rm = TRUE), "-", max(n_tested, na.rm = TRUE)))
      colnames(n_per_mo) <- c("mo", "count")
      count_group <- n_per_mo$count[match(new_df[[1]], n_per_mo$mo)]
      edit_col <- 1
    }
    if (NCOL(new_df) == edit_col + 1) {
      # only 1 antibiotic
      new_df[[edit_col]] <- paste0(new_df[[edit_col]], " (", unlist(lapply(strsplit(x = count_group, split = "-", fixed = TRUE), function(x) x[1])), ")")
      colnames(new_df)[edit_col] <- paste(colnames(new_df)[edit_col], "(N)")
    } else {
      # more than 1 antibiotic
      new_df[[edit_col]] <- paste0(new_df[[edit_col]], " (", count_group, ")")
      colnames(new_df)[edit_col] <- paste(colnames(new_df)[edit_col], "(N min-max)")
    }
  }

  out <- structure(as_original_data_class(new_df, class(x), extra_class = "antibiogram"),
    has_syndromic_group = has_syndromic_group,
    combine_SI = combine_SI,
    wisca = wisca,
    conf_interval = conf_interval,
    formatting_type = formatting_type,
    wisca_parameters = as_original_data_class(wisca_parameters, class(x)),
    long_numeric = as_original_data_class(long_numeric, class(x))
  )
  rownames(out) <- NULL
  out
}

#' @method antibiogram grouped_df
#' @export
antibiogram.grouped_df <- function(x,
                                   antimicrobials = where(is.sir),
                                   mo_transform = NULL,
                                   ab_transform = "name",
                                   syndromic_group = NULL,
                                   add_total_n = FALSE,
                                   only_all_tested = FALSE,
                                   digits = ifelse(wisca, 1, 0),
                                   formatting_type = getOption("AMR_antibiogram_formatting_type", ifelse(wisca, 14, 18)),
                                   col_mo = NULL,
                                   language = get_AMR_locale(),
                                   minimum = 30,
                                   combine_SI = TRUE,
                                   sep = " + ",
                                   sort_columns = TRUE,
                                   wisca = FALSE,
                                   simulations = 1000,
                                   conf_interval = 0.95,
                                   interval_side = "two-tailed",
                                   info = interactive(),
                                   ...) {
  stop_ifnot(is.null(mo_transform), "`mo_transform` must not be set if creating an antibiogram using a grouped tibble. The groups will become the variables over which the antimicrobials are calculated, which could include the pathogen information (though not necessary). Nonetheless, this makes `mo_transform` redundant.", call = FALSE)
  stop_ifnot(is.null(syndromic_group), "`syndromic_group` must not be set if creating an antibiogram using a grouped tibble. The groups will become the variables over which the antimicrobials are calculated, making `syndromic_groups` redundant.", call = FALSE)
  groups <- attributes(x)$groups
  n_groups <- NROW(groups)
  progress <- progress_ticker(
    n = n_groups,
    n_min = 5,
    print = info,
    title = paste("Calculating AMR for", n_groups, "groups")
  )
  on.exit(close(progress))

  out <- NULL
  wisca_parameters <- NULL
  long_numeric <- NULL

  for (i in seq_len(n_groups)) {
    progress$tick()
    rows <- unlist(groups[i, ]$.rows)
    if (length(rows) == 0) {
      next
    }
    new_out <- antibiogram(as.data.frame(x)[rows, , drop = FALSE],
      antimicrobials = antimicrobials,
      mo_transform = NULL,
      ab_transform = ab_transform,
      syndromic_group = NULL,
      add_total_n = add_total_n,
      only_all_tested = only_all_tested,
      digits = digits,
      formatting_type = formatting_type,
      col_mo = col_mo,
      language = language,
      minimum = minimum,
      combine_SI = combine_SI,
      sep = sep,
      sort_columns = sort_columns,
      wisca = wisca,
      simulations = simulations,
      conf_interval = conf_interval,
      interval_side = interval_side,
      info = FALSE,
      ...
    )
    new_wisca_parameters <- attributes(new_out)$wisca_parameters
    new_long_numeric <- attributes(new_out)$long_numeric

    if (NROW(new_out) == 0) {
      next
    }

    # remove first column 'Pathogen' (in whatever language), except WISCA since that never has Pathogen column
    if (isFALSE(wisca)) {
      new_out <- new_out[, -1, drop = FALSE]
      new_long_numeric <- new_long_numeric[, -1, drop = FALSE]
    }

    # add group names to data set
    for (col in rev(seq_len(NCOL(groups) - 1))) {
      col_name <- colnames(groups)[col]
      col_value <- groups[i, col, drop = TRUE]
      new_out[, col_name] <- col_value
      new_out <- new_out[, c(col_name, setdiff(names(new_out), col_name))] # set place to 1st col

      if (isTRUE(wisca)) {
        new_wisca_parameters[, col_name] <- col_value
        new_wisca_parameters <- new_wisca_parameters[, c(col_name, setdiff(names(new_wisca_parameters), col_name))] # set place to 1st col
      }

      new_long_numeric[, col_name] <- col_value
      new_long_numeric <- new_long_numeric[, c(col_name, setdiff(names(new_long_numeric), col_name))] # set place to 1st col
    }

    if (i == 1) {
      # the first go
      out <- new_out
      wisca_parameters <- new_wisca_parameters
      long_numeric <- new_long_numeric
    } else {
      out <- rbind_AMR(out, new_out)
      wisca_parameters <- rbind_AMR(wisca_parameters, new_wisca_parameters)
      long_numeric <- rbind_AMR(long_numeric, new_long_numeric)
    }
  }

  close(progress)

  out <- structure(as_original_data_class(out, class(x), extra_class = "antibiogram"),
    has_syndromic_group = FALSE,
    combine_SI = isTRUE(combine_SI),
    wisca = isTRUE(wisca),
    conf_interval = conf_interval,
    formatting_type = formatting_type,
    wisca_parameters = as_original_data_class(wisca_parameters, class(x)),
    long_numeric = as_original_data_class(long_numeric, class(x))
  )
  rownames(out) <- NULL
  out
}

#' @export
#' @rdname antibiogram
wisca <- function(x,
                  antimicrobials = where(is.sir),
                  ab_transform = "name",
                  syndromic_group = NULL,
                  only_all_tested = FALSE,
                  digits = 1,
                  formatting_type = getOption("AMR_antibiogram_formatting_type", 14),
                  col_mo = NULL,
                  language = get_AMR_locale(),
                  combine_SI = TRUE,
                  sep = " + ",
                  sort_columns = TRUE,
                  simulations = 1000,
                  conf_interval = 0.95,
                  interval_side = "two-tailed",
                  info = interactive(),
                  ...) {
  antibiogram(
    x = x,
    antimicrobials = antimicrobials,
    ab_transform = ab_transform,
    mo_transform = NULL,
    syndromic_group = syndromic_group,
    add_total_n = FALSE,
    only_all_tested = only_all_tested,
    digits = digits,
    formatting_type = formatting_type,
    col_mo = col_mo,
    language = language,
    combine_SI = combine_SI,
    sep = sep,
    sort_columns = sort_columns,
    wisca = TRUE,
    simulations = simulations,
    conf_interval = conf_interval,
    interval_side = interval_side,
    info = info,
    ...
  )
}

create_wisca_priors <- function(data) {
  pathogens <- unique(data$mo)
  n_pathogens <- length(pathogens)

  # Dirichlet prior (gamma parameters)
  gamma_prior <- rep(1, times = n_pathogens)
  multinomial_obs <- data$n_total
  gamma_posterior <- gamma_prior + multinomial_obs

  # beta priors
  beta_prior_alpha <- rep(1, times = n_pathogens)
  beta_prior_beta <- rep(1, times = n_pathogens)

  r <- data$n_susceptible
  n <- data$n_tested
  diff_nr <- n - r

  beta_posterior_1 <- beta_prior_alpha + r
  beta_posterior_2 <- beta_prior_beta + diff_nr

  list(
    gamma_posterior = gamma_posterior,
    beta_posterior_1 = beta_posterior_1,
    beta_posterior_2 = beta_posterior_2
  )
}

simulate_coverage <- function(params) {
  n_pathogens <- length(params$gamma_posterior)

  # random draws per pathogen
  random_incidence <- stats::runif(n = n_pathogens)
  random_susceptibility <- stats::runif(n = n_pathogens)

  simulated_incidence <- stats::qgamma(
    p = random_incidence,
    shape = params$gamma_posterior,
    scale = 1
  )

  # normalise incidence
  simulated_incidence <- simulated_incidence / sum(simulated_incidence, na.rm = TRUE)

  simulated_susceptibility <- stats::qbeta(
    p = random_susceptibility,
    shape1 = params$beta_posterior_1,
    shape2 = params$beta_posterior_2
  )

  # weighted coverage
  sum(simulated_incidence * simulated_susceptibility, na.rm = TRUE)
}

#' @export
#' @param wisca_model The outcome of [wisca()] or [`antibiogram(..., wisca = TRUE)`][antibiogram()].
#' @rdname antibiogram
retrieve_wisca_parameters <- function(wisca_model, ...) {
  stop_ifnot(isTRUE(attributes(wisca_model)$wisca), "This function only applies to WISCA models. Use `wisca()` or `antibiogram(..., wisca = TRUE)` to create a WISCA model.")
  attributes(wisca_model)$wisca_parameters
}

# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(pillar::tbl_sum, antibiogram)
tbl_sum.antibiogram <- function(x, ...) {
  dims <- paste(format(NROW(x), big.mark = ","), AMR_env$cross_icon, format(NCOL(x), big.mark = ","))
  names(dims) <- "An Antibiogram"
  if (isTRUE(attributes(x)$wisca)) {
    dims <- c(dims, Type = paste0("WISCA with ", attributes(x)$conf_interval * 100, "% CI"))
  } else if (isTRUE(attributes(x)$formatting_type >= 13)) {
    dims <- c(dims, Type = paste0("Non-WISCA with ", attributes(x)$conf_interval * 100, "% CI"))
  } else {
    dims <- c(dims, Type = paste0("Non-WISCA without CI"))
  }
  dims
}

# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(pillar::tbl_format_footer, antibiogram)
tbl_format_footer.antibiogram <- function(x, ...) {
  footer <- NextMethod()
  if (NROW(x) == 0) {
    return(footer)
  }
  c(footer, font_subtle(paste0(
    "# Use `ggplot2::autoplot()` or base R `plot()` to create a plot of this antibiogram,\n",
    "# or use it directly in R Markdown or ",
    font_url("https://quarto.org", "Quarto"), ", see ", word_wrap("?antibiogram")
  )))
}

#' @export
#' @rdname antibiogram
plot.antibiogram <- function(x, ...) {
  df <- attributes(x)$long_numeric
  if (!"mo" %in% colnames(df)) {
    df$mo <- ""
  }
  if ("syndromic_group" %in% colnames(df)) {
    # barplot in base R does not support facets - paste columns together
    df$mo <- paste(df$mo, "-", df$syndromic_group)
    df$syndromic_group <- NULL
    df <- df[order(df$mo), , drop = FALSE]
  }
  mo_levels <- unique(df$mo)
  mfrow_old <- graphics::par()$mfrow
  sqrt_levels <- sqrt(length(mo_levels))
  graphics::par(mfrow = c(ceiling(sqrt_levels), floor(sqrt_levels)))

  for (i in seq_along(mo_levels)) {
    mo <- mo_levels[i]
    df_sub <- df[as.character(df$mo) == mo, , drop = FALSE]

    bp <- barplot(
      height = df_sub$coverage * 100,
      xlab = NULL,
      ylab = ifelse(isTRUE(attributes(x)$combine_SI), "%SI", "%S"),
      names.arg = df_sub$ab,
      col = "#aaaaaa",
      beside = TRUE,
      main = mo,
      legend = NULL
    )

    if (isTRUE(attributes(x)$wisca)) {
      lower_ci <- df_sub$lower_ci * 100
      upper_ci <- df_sub$upper_ci * 100
      arrows(
        x0 = bp, y0 = lower_ci, # Start of error bar (lower bound)
        x1 = bp, y1 = upper_ci, # End of error bar (upper bound)
        angle = 90, code = 3, length = 0.05, col = "black"
      )
    }
  }

  graphics::par(mfrow = mfrow_old)
}

#' @export
#' @noRd
barplot.antibiogram <- function(height, ...) {
  plot(height, ...)
}

#' @method autoplot antibiogram
#' @rdname antibiogram
# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(ggplot2::autoplot, antibiogram)
autoplot.antibiogram <- function(object, ...) {
  df <- attributes(object)$long_numeric
  if (!"mo" %in% colnames(df)) {
    df$mo <- ""
  }
  groups <- colnames(df)[seq_len(which(colnames(df) %in% c("mo", "ab"))[1] - 1)]
  group_name <- paste(groups, collapse = "/")
  if (length(groups) > 1) {
    df$syndromic_group <- apply(df[groups], 1, function(x) {
      paste(stats::na.omit(x), collapse = "/")
    })
  } else if ("syndromic_group" %in% colnames(df)) {
    group_name <- colnames(object)[1]
  }
  out <- ggplot2::ggplot(df,
    mapping = ggplot2::aes(
      x = ab,
      y = coverage * 100,
      fill = if ("syndromic_group" %in% colnames(df)) {
        syndromic_group
      } else {
        NULL
      }
    )
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::geom_errorbar(
      mapping = ggplot2::aes(ymin = lower_ci * 100, ymax = upper_ci * 100),
      position = ggplot2::position_dodge2(preserve = "single", width = 1)
    ) +
    ggplot2::labs(
      y = ifelse(isTRUE(attributes(object)$combine_SI), "%SI", "%S"),
      x = NULL,
      fill = if ("syndromic_group" %in% colnames(df)) {
        group_name
      } else {
        NULL
      }
    )
  if (!all(as.character(df$mo) == "", na.rm = TRUE)) {
    out <- out +
      ggplot2::facet_wrap("mo")
  }
  out
}

#' @method knit_print antibiogram
#' @param italicise A [logical] to indicate whether the microorganism names in the [knitr][knitr::kable()] table should be made italic, using [italicise_taxonomy()].
#' @param na Character to use for showing `NA` values.
#' @rdname antibiogram
# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(knitr::knit_print, antibiogram)
knit_print.antibiogram <- function(x, italicise = TRUE, na = getOption("knitr.kable.NA", default = ""), ...) {
  stop_ifnot_installed("knitr")
  meet_criteria(italicise, allow_class = "logical", has_length = 1)
  meet_criteria(na, allow_class = "character", has_length = 1, allow_NA = TRUE)

  add_MO_lookup_to_AMR_env()

  for (i in which(vapply(FUN.VALUE = logical(1), x, is.character))) {
    # make all microorganism names italic, according to nomenclature
    x[[i]] <- italicise_taxonomy(x[[i]], type = "markdown")
  }

  old_option <- getOption("knitr.kable.NA")
  options(knitr.kable.NA = na)
  on.exit(options(knitr.kable.NA = old_option))

  out <- paste(c("", "", knitr::kable(x, ..., output = FALSE)), collapse = "\n")
  knitr::asis_output(out)
}
