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

#' Generate Antibiograms (WISCA, Traditional, Combination, or Syndromic)
#'
#' @description
#' Generate antibiograms from antimicrobial susceptibility data, with support for traditional, combination, syndromic, and WISCA (Weighted-Incidence Syndromic Combination Antibiogram) methods.
#'
#' **For empirical therapy guidance, WISCA is the recommended approach.** When initiating empirical treatment, the causative pathogen is unknown, and the clinically relevant question is: *"what is the probability that this regimen will cover whatever pathogen turns out to cause the infection?"* WISCA answers that question directly by weighting susceptibility by pathogen incidence within a syndrome and providing credible intervals via Bayesian Monte Carlo simulation. Traditional antibiograms remain appropriate for tracking resistance per species for surveillance purposes. See the section *Explaining WISCA* on this page and the [WISCA vignette](https://amr-for-r.org/articles/WISCA.html) for details.
#'
#' All antibiogram types adhere to previously described approaches (see *Source*), and the WISCA method implements the Bayesian decision model by Bielicki *et al.* (2016, \doi{10.1093/jac/dkv397}). Output formats include plots and tables, ideal for integration with R Markdown and Quarto reports.
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
#' @param mo_transform A character to transform microorganism input - must be `"name"`, `"shortname"` (default), `"gramstain"`, or one of the column names of the [microorganisms] data set: `r vector_or(colnames(microorganisms), sort = FALSE, documentation = TRUE)`. Can also be `NULL` to not transform the input or `NA` to consider all microorganisms 'unknown'.
#' @param ab_transform A character to transform antimicrobial input - must be one of the column names of the [antimicrobials] data set (defaults to `"name"`): `r vector_or(colnames(antimicrobials), sort = FALSE, documentation = TRUE)`. Can also be `NULL` to not transform the input.
#' @param syndromic_group A column name of `x`, or values calculated to split rows of `x`, e.g. by using [ifelse()] or [`case_when()`][dplyr::case_when()]. See *Examples*.
#' @param add_total_n *(deprecated in favour of `formatting_type`)* A [logical] to indicate whether `n_tested` available numbers per pathogen should be added to the table (default is `TRUE`). This will add the lowest and highest number of available isolates per antimicrobial (e.g., if for *E. coli* 200 isolates are available for ciprofloxacin and 150 for amoxicillin, the returned number will be "150-200"). This option is unavailable when `wisca = TRUE`; in that case, use [retrieve_wisca_parameters()] to get the parameters used for WISCA.
#' @param only_all_tested (for combination antibiograms): a [logical] to indicate that isolates must be tested for all antimicrobials, see *Details*.
#' @param digits Number of digits to use for rounding the antimicrobial coverage, defaults to 1 for WISCA and 0 otherwise.
#' @param formatting_type Numeric value (1-22 for WISCA, 1-12 for non-WISCA) indicating how the 'cells' of the antibiogram table should be formatted. See *Details* > *Formatting Type* for a list of options.
#' @param col_mo Column name of the names or codes of the microorganisms (see [as.mo()]) - the default is the first column of class [`mo`]. Values will be coerced using [as.mo()].
#' @param language Language to translate text, which defaults to the system language (see [get_AMR_locale()]).
#' @param minimum The minimum allowed number of available (tested) isolates. Any isolate count lower than `minimum` will return `NA` with a warning. The default number of `30` isolates is advised by the Clinical and Laboratory Standards Institute (CLSI) as best practice, see *Source*.
#' @param combine_SI A [logical] to indicate whether all susceptibility should be determined by results of either S, SDD, or I, instead of only S (default is `TRUE`).
#' @param sep A separating character for antimicrobial columns in combination antibiograms.
#' @param sort_columns A [logical] to indicate whether the antimicrobial columns must be sorted on name.
#' @param wisca A [logical] to indicate whether a Weighted-Incidence Syndromic Combination Antibiogram (WISCA) must be generated (default is `FALSE`). This will use a Bayesian decision model to estimate regimen coverage probabilities using [Monte Carlo simulations](https://en.wikipedia.org/wiki/Monte_Carlo_method). Per \doi{10.1093/jac/dkv397}, susceptibility priors are \eqn{\beta(0.5, 0.5)} (Jeffreys) and intrinsically resistant pairs (based on [intrinsic_resistant]) use \eqn{\beta(1, 9999)}.
#'
#' Set `simulations`, `conf_interval`, and `interval_side` to adjust.
#' @param simulations (for WISCA) a numerical value to set the number of Monte Carlo simulations.
#' @param conf_interval A numerical value to set confidence interval (default is `0.95`).
#' @param interval_side The side of the confidence interval, either `"two-tailed"` (default), `"left"` or `"right"`.
#' @param parallel A [logical] to indicate if parallel computing must be used, defaults to `FALSE`. Requires the [`future.apply`][future.apply::future_lapply()] package. For WISCA, Monte Carlo simulations are distributed across workers; for grouped antibiograms, each group is processed by a separate worker. **A non-sequential [future::plan()] must already be active before setting `parallel = TRUE`** -- for example, `future::plan(future::multisession)`. An error is thrown if `parallel = TRUE` is used without a plan set by the user.
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
#' ### When to Use WISCA vs. Traditional Antibiograms
#'
#' There are various antibiogram types, as summarised by Klinker *et al.* (2021, \doi{10.1177/20499361211011373}), and they are all supported by [antibiogram()]: traditional, combination, syndromic, and WISCA.
#'
#' **If your goal is to guide empirical therapy, use WISCA.** Traditional antibiograms fragment susceptibility information by species, but at the point of prescribing, the clinician does not know which species is causing the infection. WISCA shifts the unit of analysis from the isolate to the patient: it estimates the probability that a regimen will cover the infection, given the local distribution of causative pathogens. It evaluates combination regimens, weights by pathogen incidence, and provides credible intervals that honestly communicate uncertainty. Hebert *et al.* (2012) demonstrated this concretely for the first time: ciprofloxacin showed 84% susceptibility against *E. coli* in the traditional antibiogram, but WISCA coverage was only 62% for UTI and 37% for abdominal infections, because other species (including intrinsically resistant enterococci) contribute substantially to these syndromes. Note that WISCA is pathogen-agnostic: the outcome is not stratified by species, but by syndrome.
#'
#' **Traditional, combination, and syndromic antibiograms remain appropriate for AMR surveillance**, i.e., tracking resistance trends per species over time. They are the right tool when the question is *"how resistant is species X to drug Y in our setting?"* rather than *"what regimen best covers this syndrome?"*.
#'
#' All four types are demonstrated in the *Examples* section below.
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
#' WISCA (Weighted-Incidence Syndromic Combination Antibiogram) estimates the probability that an empirical antimicrobial regimen will provide adequate coverage for a given infection syndrome, before the causative pathogen has been identified.
#'
#' It does so by combining two quantities: the relative incidence of each pathogen within the syndrome (modelled as a Dirichlet distribution) and the susceptibility of each pathogen to the regimen (modelled as Beta distributions). These are combined via Monte Carlo simulation to produce a coverage estimate with a credible interval.
#'
#' **Prior distributions:** Pathogen incidence uses a non-informative \eqn{Dirichlet(1, 1, \ldots, 1)} prior. Susceptibility proportions use the Jeffreys prior, \eqn{\beta(0.5, 0.5)}, except for pathogen-drug combinations with known intrinsic resistance, which use a strongly informative \eqn{\beta(1, 9999)} prior that forces near-zero susceptibility regardless of observed data. Intrinsic resistance is determined using the [intrinsic_resistant] data set, which is based on `r format_eucast_version_nr(names(EUCAST_VERSION_EXPECTED_PHENOTYPES[1]))`.
#'
#' **Interpreting the output:** Overlapping credible intervals between regimens indicate no significant difference in coverage; if a narrower-spectrum regimen overlaps with a broader one, the narrower-spectrum option may be preferred on stewardship grounds. Non-overlapping intervals indicate a clinically meaningful difference. For small sample sizes, consider pooling data from multiple sites to improve precision, provided pathogen distributions are sufficiently similar (Bielicki *et al.*, 2016).
#'
#' For the full mathematical derivation and worked examples, see the [WISCA vignette](https://amr-for-r.org/articles/WISCA.html).
#' @references
#' * Hebert C *et al.* (2012). **Demonstration of the weighted-incidence syndromic combination antibiogram: an empiric prescribing decision aid.** *Infection Control & Hospital Epidemiology* 33(4):381-388; \doi{10.1086/664768}
#' * Bielicki JA *et al.* (2016). **Selecting appropriate empirical antibiotic regimens for paediatric bloodstream infections: application of a Bayesian decision model to local and pooled antimicrobial resistance surveillance data.** *Journal of Antimicrobial Chemotherapy* 71(3):794-802; \doi{10.1093/jac/dkv397}
#' * Cook A *et al.* (2022). **Improving empiric antibiotic prescribing in pediatric bloodstream infections: a potential application of weighted-incidence syndromic combination antibiograms (WISCA).** *Expert Review of Anti-infective Therapy* 20(3):445-456; \doi{10.1080/14787210.2021.1967145}
#' * Klinker KP *et al.* (2021). **Antimicrobial stewardship and antibiograms: importance of moving beyond traditional antibiograms.** *Therapeutic Advances in Infectious Disease*, May 5;8:20499361211011373; \doi{10.1177/20499361211011373}
#' * Barbieri E *et al.* (2021). **Development of a Weighted-Incidence Syndromic Combination Antibiogram (WISCA) to guide the choice of the empiric antibiotic treatment for urinary tract infection in paediatric patients: a Bayesian approach.** *Antimicrobial Resistance & Infection Control* May 1;10(1):74; \doi{10.1186/s13756-021-00939-2}
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
#' # WISCA antibiogram (recommended for empirical therapy) -----------------
#'
#' # basic WISCA: empirical coverage per regimen, weighted by pathogen
#' # incidence, with 95% credible intervals
#' wisca(example_isolates,
#'   antimicrobials = c("AMC", "AMC+CIP", "AMC+GEN")
#' )
#'
#' # equivalent using antibiogram():
#' antibiogram(example_isolates,
#'   antimicrobials = c("AMC", "AMC+CIP", "AMC+GEN"),
#'   wisca = TRUE
#' )
#'
#' # stratified by syndrome or clinical group
#' wisca(example_isolates,
#'   antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
#'   syndromic_group = "ward"
#' )
#'
#' # stratified using grouped tibbles (e.g. by age and gender)
#' if (requireNamespace("dplyr")) {
#'   library(dplyr)
#'   example_isolates %>%
#'     top_n_microorganisms(n = 10) %>%
#'     group_by(
#'       age_group = age_groups(age, c(25, 50, 75)),
#'       gender) %>%
#'     wisca(antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"))
#' }
#'
#'
#' # Traditional antibiogram (for AMR surveillance) ------------------------
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
#'
#' # Combination antibiogram (for AMR surveillance) ------------------------
#'
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
#' # Syndromic antibiogram (for AMR surveillance) --------------------------
#'
#' antibiogram(example_isolates,
#'   antimicrobials = c(aminoglycosides(), carbapenems()),
#'   syndromic_group = "ward"
#' )
#'
#' # with a custom language, though this will be determined automatically
#' # (i.e., this table will be in Spanish on Spanish systems)
#' ex1 <- example_isolates[which(mo_genus() == "Escherichia"), ]
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
#' # Print the output for R Markdown / Quarto -----------------------------
#'
#' ureido <- wisca(example_isolates,
#'   antimicrobials = ureidopenicillins(),
#'   syndromic_group = "ward"
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
#' ab2 <- wisca(example_isolates,
#'   antimicrobials = c("AMC", "CIP", "TZP", "TZP+TOB"),
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
                  parallel = FALSE,
                  ...) {
  antibiogram(
    x = x,
    antimicrobials = antimicrobials,
    ab_transform = ab_transform,
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
    parallel = parallel,
    ...
  )
}

#' @export
#' @rdname antibiogram
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
                        parallel = FALSE,
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
                                parallel = FALSE,
                                ...) {
  meet_criteria(x, allow_class = "data.frame")
  x <- ascertain_sir_classes(x, "x")
  meet_criteria(wisca, allow_class = "logical", has_length = 1)
  if (wisca) {
    if (!is.null(mo_transform) && !missing(mo_transform)) {
      stop_("{.arg mo_transform} cannot be used when creating a WISCA. WISCA already integrates pathogen incidence into the coverage estimate, so the output is inherently pathogen-agnostic. To stratify results, use {.arg syndromic_group} instead.")
    }
    mo_transform <- function(x) suppressMessages(suppressWarnings(as.mo(x, keep_synonyms = TRUE, language = NULL, info = FALSE)))
  }
  if ("antibiotics" %in% names(list(...))) {
    deprecation_warning("antibiotics", "antimicrobials", fn = "antibiogram", is_argument = TRUE)
    antimicrobials <- list(...)$antibiotics
  }
  meet_criteria(antimicrobials, allow_class = c("character", "numeric", "integer", "function"), allow_NA = FALSE, allow_NULL = FALSE)
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
  meet_criteria(parallel, allow_class = "logical", has_length = 1)

  # get syndromic groups
  if (!is.null(syndromic_group)) {
    if (length(syndromic_group) == 1 && syndromic_group %in% colnames(x)) {
      x$`.syndromic_group` <- x[, syndromic_group, drop = TRUE]
    } else if (length(syndromic_group) > 1 && all(syndromic_group %in% colnames(x))) {
      x$`.syndromic_group` <- do.call(paste, c(x[syndromic_group], list(sep = "||")))
      attr(x, "antibiogram_groups") <- syndromic_group
    } else if (!is.null(syndromic_group) && length(syndromic_group) == 1) {
      x$`.syndromic_group` <- syndromic_group
    } else if (NCOL(syndromic_group) > 1) {
      stop_("{.arg syndromic_group} should be a 1-dimensional computed value, or 1 or more column names of {.arg x}.")
    } else {
      x$`.syndromic_group` <- syndromic_group
    }
    if (any(x$`.syndromic_group` %in% c(NA, ""))) {
      x$`.syndromic_group`[x$`.syndromic_group` %in% c(NA, "")] <- paste0("(", translate_AMR("unknown", language = language), ")")
    }
    has_syndromic_group <- TRUE
  } else {
    has_syndromic_group <- FALSE
  }

  # parallel gate - identical pattern to as.sir()
  if (requireNamespace("future.apply", quietly = TRUE)) {
    if (!inherits(future::plan(), "sequential")) {
      if (isFALSE(parallel)) {
        message_("Assuming {.code parallel = TRUE} since parallel computing has been set up using the {.pkg future} package before. Set {.help [{.fun plan}](future::plan)} to sequential to prevent this.")
      }
      parallel <- TRUE
    }
    if (wisca && interactive() && message_not_thrown_before("antibiogram", "wisca_parallel") && inherits(future::plan(), "sequential") && isFALSE(parallel) && simulations > 100) {
      advised_multi <- ifelse(.Platform$OS.type == "windows" || in_rstudio(), "multisession", "multicore")
      sims <- simulations * length(antimicrobials)
      if (has_syndromic_group) {
        sims <- sims * length(unique(x$`.syndromic_group`))
      }
      message_("Are you sure you want to run in non-parallel (=sequential) mode?", as_note = FALSE)
      message_("WISCA can take a ", ifelse(sims > 10000, font_bold("very "), ""), "long time for the ", format(sims, decimal.mark = ".", big.mark = " "), " simulations you require, and you already have the {.pkg future} package installed.", as_note = FALSE)
      q <- utils::menu(c(
        "Yes, still run in sequential mode",
        format_inline_("No, run in parallel mode and set {.help [future::plan(", advised_multi, ")](future::plan)}, and reset after WISCA finishes"),
        format_inline_("No, run in parallel mode and set {.help [future::plan(", advised_multi, ")](future::plan)}, and do not reset afterwards"),
        "Cancel WISCA calculation"
      ), graphics = FALSE, title = "")
      if (q %in% c(4, 0)) {
        return(invisible(NULL))
      } else if (q %in% c(2, 3)) {
        parallel <- TRUE
        AMR_env$wisca_parallel_choice <- "parallel"
        obj <- get(advised_multi, envir = asNamespace("future"))
        future::plan(obj)
        if (q == 2) {
          AMR_env$wisca_parallel_choice <- "parallel_reset"
        }
      } else {
        AMR_env$wisca_parallel_choice <- "sequential"
      }
    } else if (wisca && !is.null(AMR_env$wisca_parallel_choice)) {
      if (AMR_env$wisca_parallel_choice %in% c("parallel", "parallel_reset")) {
        parallel <- TRUE
      }
    }
    if (identical(AMR_env$wisca_parallel_choice, "parallel_reset") && inherits(future::plan(), "uniprocess", which = FALSE) == FALSE) {
      on.exit(
        {
          message_("Resetting {.fn future::plan}...", as_note = FALSE)
          future::plan(future::sequential)
          AMR_env$wisca_parallel_choice <- NULL
          message_("Parallel setting was reset to `future::plan(future::sequential)`.", as_check = TRUE)
        },
        add = TRUE
      )
    }
  }
  if (isTRUE(parallel)) {
    stop_ifnot(
      requireNamespace("future.apply", quietly = TRUE),
      "Setting {.code parallel = TRUE} requires the {.pkg future.apply} package.\n",
      "Install it with {.code install.packages(\"future.apply\")}."
    )
    stop_if(inherits(future::plan(), "sequential"),
      "Setting {.code parallel = TRUE} requires a non-sequential {.help [{.fun future::plan}](future::plan)} to be active.\n",
      "For your system, you could first run: {.code library(future); ",
      ifelse(.Platform$OS.type == "windows" || in_rstudio(),
        "plan(multisession)",
        "plan(multicore)"
      ),
      "}",
      call = FALSE
    )
    n_workers <- future::nbrOfWorkers()
  } else {
    n_workers <- 1L
  }

  # try to find columns based on type
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = info)
    stop_if(is.null(col_mo), "{.arg col_mo} must be set")
  }
  # transform MOs
  x$`.mo` <- x[, col_mo, drop = TRUE]
  if (is.null(mo_transform)) {
    # leave as is, no transformation, but do add backup
    x$`.mo.bak` <- x$`.mo`
  } else {
    x$`.mo` <- as.mo(x$`.mo`, keep_synonyms = TRUE, info = FALSE)
    x$`.mo.bak` <- x$`.mo`
    if (is.function(mo_transform)) {
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
  }
  x$`.mo`[x$`.mo` %in% c(NA, "UNKNOWN")] <- "(??)"

  # get antimicrobials
  ab_trycatch <- tryCatch(colnames(suppressWarnings(x[, antimicrobials, drop = FALSE])), error = function(e) NULL)
  if (is.null(ab_trycatch)) {
    # try with tidyverse
    ab_trycatch <- tryCatch(colnames(dplyr::select(x, {{ antimicrobials }})), error = function(e) NULL)
  }
  if (is.null(ab_trycatch)) {
    stop_ifnot(is.character(suppressMessages(antimicrobials)), "{.arg antimicrobials} must be an antimicrobial selector, or a character vector.")
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
        S_values <- c("S", "WT")
        if (isTRUE(combine_SI)) {
          S_values <- c(S_values, "SDD", "I")
        }
        other_values <- setdiff(c("S", "SDD", "I", "R", "WT", "NWT", "NS"), S_values)
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
        "Detected column name(s) containing the '+' character, which conflicts with the expected syntax in {.help [{.fun antibiogram}](AMR::antibiogram)}: the '+' is used to combine separate antimicrobial drug columns (e.g., \"AMP+GEN\").\n\n",
        "To avoid incorrectly guessing which antimicrobials this represents, {.arg ab_transform} was automatically set to {.code NULL}.\n\n",
        "If this is unintended, please rename the column(s) to avoid using '+' in the name, or set {.code ab_transform = NULL} explicitly to suppress this message."
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

  out$n_susceptible <- out$S + out$WT
  if (isTRUE(combine_SI)) {
    out$n_susceptible <- out$n_susceptible + out$I + out$SDD
  }
  if (all(out$n_tested < minimum, na.rm = TRUE) && wisca == FALSE) {
    warning_("All combinations had less than {.arg minimum} = ", minimum, " results, returning an empty antibiogram")
    return(as_original_data_class(data.frame(), class(x), extra_class = "antibiogram"))
  } else if (any(out$n_tested < minimum, na.rm = TRUE)) {
    mins <- sum(out$n_tested < minimum, na.rm = TRUE)
    if (wisca == FALSE) {
      out <- out %pm>%
        subset(n_tested >= minimum)
      if (isTRUE(info) && mins > 0) {
        message_("NOTE: ", mins, " combinations had less than {.arg minimum} = ", minimum, " results and were ignored")
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

  # WISCA START ----
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
    wisca_draws <- list()
    wisca_components <- list()

    # quantile probabilities are constant across all groups
    probs <- if (interval_side == "two-tailed") {
      c((1 - conf_interval) / 2, 1 - (1 - conf_interval) / 2)
    } else if (interval_side == "left") {
      c(0, conf_interval)
    } else {
      c(1 - conf_interval, 1)
    }

    unique_groups <- as.character(unique(wisca_parameters$group))

    use_parallel_wisca <- isTRUE(parallel) && n_workers > 1L && length(unique_groups) > 0L

    if (use_parallel_wisca) {
      ## WISCA parallel ----
      if (isTRUE(info)) {
        message_("Running WISCA in parallel mode using ", n_workers, " workers...", as_note = FALSE, appendLF = FALSE)
      }
      # chunks_per_group gives ~n_workers total jobs so all workers stay busy
      # even when the number of regimens is smaller than n_workers
      chunks_per_group <- max(1L, ceiling(n_workers / length(unique_groups)))
      chunk_sizes <- diff(c(0L, round(seq_len(chunks_per_group) * simulations / chunks_per_group)))

      params_g_lookup <- list()

      # precompute priors per group and build (group, chunk) job list
      jobs <- unlist(lapply(unique_groups, function(g) {
        params_g <- wisca_parameters[wisca_parameters$group == g, , drop = FALSE]
        if (sum(params_g$n_tested, na.rm = TRUE) == 0L) {
          return(NULL)
        }
        # store for later reassembly
        params_g_lookup[[g]] <<- params_g
        priors_g <- create_wisca_priors(params_g, sep = sep)
        lapply(seq_along(chunk_sizes), function(ch) {
          list(group = g, priors = priors_g, n_sims = chunk_sizes[ch])
        })
      }), recursive = FALSE)
      jobs <- Filter(Negate(is.null), jobs)

      flat <- future.apply::future_lapply(jobs, function(job) {
        n_p <- length(job$priors$gamma_posterior)
        n_s <- job$n_sims
        inc_mat <- matrix(NA_real_, nrow = n_s, ncol = n_p)
        susc_mat <- matrix(NA_real_, nrow = n_s, ncol = n_p)
        cov_vec <- numeric(n_s)
        for (i in seq_len(n_s)) {
          inc_raw <- stats::rgamma(n_p, shape = job$priors$gamma_posterior, scale = 1)
          inc_norm <- inc_raw / sum(inc_raw)
          susc <- stats::rbeta(n_p,
            shape1 = job$priors$beta_posterior_1,
            shape2 = job$priors$beta_posterior_2
          )
          inc_mat[i, ] <- inc_norm
          susc_mat[i, ] <- susc
          cov_vec[i] <- sum(inc_norm * susc)
        }
        list(coverage = cov_vec, incidence = inc_mat, susceptibility = susc_mat)
      }, future.seed = TRUE)

      # reassemble per group: concatenate chunks, then summarise
      for (g in unique_groups) {
        g_idx <- vapply(jobs, function(j) identical(j$group, g), logical(1))
        if (!any(g_idx)) next
        chunks <- flat[g_idx]
        sims <- unlist(lapply(chunks, `[[`, "coverage"), use.names = FALSE)
        inc_combined <- do.call(rbind, lapply(chunks, `[[`, "incidence"))
        susc_combined <- do.call(rbind, lapply(chunks, `[[`, "susceptibility"))
        colnames(inc_combined) <- as.character(params_g_lookup[[g]]$mo)
        colnames(susc_combined) <- as.character(params_g_lookup[[g]]$mo)
        wisca_draws[[g]] <- sims
        wisca_components[[g]] <- list(incidence = inc_combined, susceptibility = susc_combined)
        out_wisca$coverage[out_wisca$group == g] <- mean(sims)
        ci_vals <- unname(stats::quantile(sims, probs = probs))
        out_wisca$lower_ci[out_wisca$group == g] <- ci_vals[1]
        out_wisca$upper_ci[out_wisca$group == g] <- ci_vals[2]
      }

      if (isTRUE(info)) message_(font_green_bg("\u00a0DONE\u00a0"), as_note = FALSE)
    } else {
      ## WISCA sequential ----
      progress <- progress_ticker(
        n = length(unique_groups) * simulations,
        n_min = 25,
        print = info,
        title = paste("Calculating WISCA for", length(unique_groups), "regimens")
      )
      on.exit(close(progress), add = TRUE)

      for (group in unique_groups) {
        params_current <- wisca_parameters[wisca_parameters$group == group, , drop = FALSE]
        if (sum(params_current$n_tested, na.rm = TRUE) == 0) next
        priors_current <- create_wisca_priors(params_current, sep = sep)
        # replace the vapply block in the sequential branch with:
        n_pathogens_g <- length(priors_current$gamma_posterior)
        sim_coverage <- numeric(simulations)
        sim_incidence <- matrix(NA_real_, nrow = simulations, ncol = n_pathogens_g)
        sim_susceptibility <- matrix(NA_real_, nrow = simulations, ncol = n_pathogens_g)
        colnames(sim_incidence) <- as.character(params_current$mo)
        colnames(sim_susceptibility) <- as.character(params_current$mo)

        for (i in seq_len(simulations)) {
          progress$tick()
          inc_raw <- stats::rgamma(n_pathogens_g, shape = priors_current$gamma_posterior, scale = 1)
          inc_norm <- inc_raw / sum(inc_raw)
          susc <- stats::rbeta(n_pathogens_g,
            shape1 = priors_current$beta_posterior_1,
            shape2 = priors_current$beta_posterior_2
          )
          sim_incidence[i, ] <- inc_norm
          sim_susceptibility[i, ] <- susc
          sim_coverage[i] <- sum(inc_norm * susc)
        }

        wisca_draws[[group]] <- sim_coverage
        wisca_components[[group]] <- list(
          incidence = sim_incidence,
          susceptibility = sim_susceptibility
        )
        out_wisca$coverage[out_wisca$group == group] <- mean(sim_coverage)
        ci_vals <- unname(stats::quantile(sim_coverage, probs = probs))
        out_wisca$lower_ci[out_wisca$group == group] <- ci_vals[1]
        out_wisca$upper_ci[out_wisca$group == group] <- ci_vals[2]
      }
      close(progress)
      if (isTRUE(info) && simulations >= 500 && length(unique_groups) >= 3) {
        suggest <- ifelse(.Platform$OS.type == "windows" || in_rstudio(),
          "plan(multisession)",
          "plan(multicore)"
        )
        if (requireNamespace("future.apply", quietly = TRUE)) {
          message_("Running in sequential mode. To speed up WISCA, set a parallel {.help [{.fun future::plan}](future::plan)} such as {.code ", suggest, "} and use {.code parallel = TRUE}.")
        } else {
          message_("Running in sequential mode. To speed up WISCA, install the {.pkg future.apply} package and then set {.code parallel = TRUE}.")
        }
      }
    }

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
    message_("Using WISCA with a {.arg formatting_type} that includes the denominator is not useful")
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

  if (wisca) {
    names(wisca_draws) <- out$ab[!is.na(out$out_value)]
    names(wisca_components) <- out$ab[!is.na(out$out_value)]
  }
  if (!is.null(attr(x, "antibiogram_groups", exact = TRUE))) {
    new_df <- as.data.frame(new_df)
    grps <- attr(x, "antibiogram_groups", exact = TRUE)
    parts <- strsplit(new_df[[1]], "||", fixed = TRUE)
    new_cols <- do.call(rbind, parts)
    colnames(new_cols) <- grps
    new_df <- cbind(as.data.frame(new_cols), new_df[-1])
  }

  out <- structure(as_original_data_class(new_df, class(x), extra_class = "antibiogram"),
    has_syndromic_group = has_syndromic_group,
    antibiogram_groups = attr(x, "antibiogram_groups", exact = TRUE),
    combine_SI = combine_SI,
    wisca = wisca,
    conf_interval = conf_interval,
    simulations = if (isFALSE(wisca)) NULL else simulations,
    formatting_type = formatting_type,
    sep = sep,
    wisca_parameters = if (isFALSE(wisca)) NULL else as_original_data_class(wisca_parameters, class(x)),
    long_numeric = as_original_data_class(long_numeric, class(x)),
    wisca_draws = if (isFALSE(wisca)) NULL else wisca_draws,
    wisca_components = if (isFALSE(wisca)) NULL else wisca_components
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
                                   parallel = FALSE,
                                   ...) {
  stop_ifnot(is.null(mo_transform), "{.arg mo_transform} must not be set if creating an antibiogram using a grouped tibble. The groups will become the variables over which the antimicrobials are calculated, which could include the pathogen information (though not necessary). Nonetheless, this makes {.arg mo_transform} redundant.", call = FALSE)
  stop_ifnot(is.null(syndromic_group), "{.arg syndromic_group} must not be set if creating an antibiogram using a grouped tibble. The groups will become the variables over which the antimicrobials are calculated, making {.arg syndromic_group} redundant.", call = FALSE)

  groups <- attributes(x)$groups
  group_cols <- intersect(colnames(groups), colnames(x))
  # paste group together, will be split later in antibiogram.default()
  x$.group <- do.call(paste, c(x[group_cols], list(sep = "||")))
  attr(x, "antibiogram_groups") <- group_cols

  antibiogram.default(x,
    antimicrobials = antimicrobials,
    mo_transform = mo_transform,
    ab_transform = ab_transform,
    syndromic_group = ".group",
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
    info = info,
    parallel = parallel,
    ...
  )
}

create_wisca_priors <- function(data, sep) {
  pathogens <- unique(data$mo)
  n_pathogens <- length(pathogens)

  # Dirichlet prior (gamma parameters)
  gamma_prior <- rep(1, times = n_pathogens)
  multinomial_obs <- data$n_total
  gamma_posterior <- gamma_prior + multinomial_obs

  # Beta priors: Jeffreys prior Beta(0.5, 0.5) by default (Bielicki et al., 2016)
  beta_prior_alpha <- rep(0.5, n_pathogens)
  beta_prior_beta <- rep(0.5, n_pathogens)

  # strongly informative Beta(1, 9999) for intrinsically resistant bug-drug pairs (Bielicki et al., 2016)
  is_intrinsic <- vapply(
    FUN.VALUE = logical(1),
    seq_len(nrow(data)),
    function(i) {
      # split by " + ", or wherever `sep` is set to
      ab_components <- as.ab(trimws(strsplit(as.character(data$ab[i]), trimws(sep), fixed = TRUE)[[1]]))
      ab_components <- ab_components[!is.na(ab_components)]
      length(ab_components) > 0 &&
        all(vapply(
          FUN.VALUE = logical(1),
          ab_components,
          function(ab) any(AMR::intrinsic_resistant$mo == data$mo[i] & AMR::intrinsic_resistant$ab == ab)
        ))
    }
  )
  beta_prior_alpha[is_intrinsic] <- 1
  beta_prior_beta[is_intrinsic] <- 9999

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
  stop_ifnot(isTRUE(attributes(wisca_model)$wisca), "This function only applies to WISCA models. Use {.help [{.fun wisca}](AMR::wisca)} or {.help [{.fun antibiogram}](AMR::antibiogram)} (with {.code wisca = TRUE}) to create a WISCA model.")
  attributes(wisca_model)$wisca_parameters
}

# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(pillar::tbl_sum, antibiogram)
tbl_sum.antibiogram <- function(x, ...) {
  dims <- NextMethod()
  if (isTRUE(attributes(x)$wisca)) {
    dims <- c(dims,
      Type = "Weighted-Incidence Syndromic Combination Antibiogram (WISCA)",
      "Cred. interval" = paste0(attributes(x)$conf_interval * 100, "%"),
      Simulations = paste0(attributes(x)$simulations, " per stratum")
    )
  } else {
    type <- ifelse(any(attributes(x)$long_numeric$ab %like% "/", na.rm = TRUE),
      "Combination Antibiogram",
      ifelse(colnames(attributes(x)$long_numeric)[1] == "syndromic_group",
        "Syndromic Antibiogram",
        "Traditional Antibiogram"
      )
    )
    dims <- c(dims, Type = type)
    if (isTRUE(attributes(x)$formatting_type >= 13)) {
      dims <- c(dims, "Conf. interval" = paste0(attributes(x)$conf_interval * 100, "%"))
    }
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
  if (isTRUE(attributes(x)$wisca)) {
    c(
      footer,
      font_subtle(format_inline_("# Use {.fn ggplot2::autoplot} or base R {.fn plot} to create a plot of this antibiogram,\n# and use {.help [{.fn wisca_plot}](AMR::antibiogram)} to assess the simulation outcomes.\n# Or, use it directly in R Markdown or Quarto, see {.help [{.fn antibiogram}](AMR::antibiogram)}."))
    )
  } else {
    c(
      footer,
      font_subtle(format_inline_("# Use {.fn ggplot2::autoplot} or base R {.fn plot} to create a plot of this antibiogram,\n# or use it directly in R Markdown or Quarto, see {.help [{.fn antibiogram}](AMR::antibiogram)}."))
    )
  }
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
#' @param geom The plotting style for the point estimate. One of `"pointrange"` (default), `"point"`, `"col"`/`"bar"`, or `"errorbar"`. `"pointrange"` is recommended for coverage data: bars imply a meaningful baseline at zero, which coverage estimates rarely have.
#' @param ci Logical, whether to draw the credible/confidence interval. Defaults to `TRUE`. Ignored (forced `TRUE`) when `geom = "pointrange"` or `"errorbar"`, since the interval is intrinsic to those geoms.
#' @param sort Logical, whether to order regimens by coverage. Defaults to `TRUE`. When faceted (per pathogen) or grouped (syndromic), ordering is applied within each panel/group.
#' @param flip Logical, whether to draw regimens on the y-axis (horizontal). Defaults to `NULL`, which flips automatically when any regimen label exceeds 20 characters (long combination names read poorly on the x-axis). Set `TRUE`/`FALSE` to override.
#' @param caption Text to show as caption, will explain non-inferiority for WISCA models.
autoplot.antibiogram <- function(object,
                                 geom = c("pointrange", "point", "col", "bar", "errorbar"),
                                 ci = TRUE,
                                 sort = TRUE,
                                 flip = NULL,
                                 caption = NULL,
                                 ...) {
  geom <- match.arg(geom)
  if (geom == "bar") geom <- "col"
  meet_criteria(ci, allow_class = "logical", has_length = 1)
  meet_criteria(sort, allow_class = "logical", has_length = 1)
  meet_criteria(flip, allow_class = "logical", has_length = 1, allow_NULL = TRUE)
  meet_criteria(caption, allow_class = "logical", has_length = 1, allow_NULL = TRUE, allow_NA = TRUE)

  df <- attributes(object)$long_numeric
  combine_SI <- isTRUE(attributes(object)$combine_SI)
  is_wisca <- isTRUE(attributes(object)$wisca)

  if (!"mo" %in% colnames(df)) {
    df$mo <- ""
  }
  groups <- colnames(df)[seq_len(which(colnames(df) %in% c("mo", "ab"))[1] - 1)]
  group_name <- paste(groups, collapse = "/")
  if (length(groups) > 1) {
    df$syndromic_group <- apply(df[groups], 1, function(x) {
      paste(stats::na.omit(x), collapse = " / ")
    })
  } else if ("syndromic_group" %in% colnames(df)) {
    if (is.null(attributes(object)$antibiogram_groups)) {
      # translated value of "Syndromic group"
      group_name <- colnames(object)[1]
    } else {
      # multiple groups, created with a grouped tibble or when syndromic_group was length >1
      group_name <- paste0(attributes(object)$antibiogram_groups, collapse = " / ")
      df$syndromic_group <- gsub("||", " / ", df$syndromic_group, fixed = TRUE)
    }
  }
  has_syndromic <- "syndromic_group" %in% colnames(df)
  has_facet <- !all(as.character(df$mo) == "", na.rm = TRUE)

  # coverage on the percentage scale
  df$.coverage <- df$coverage * 100
  df$.lower <- df$lower_ci * 100
  df$.upper <- df$upper_ci * 100

  # decide orientation: auto-flip when labels are long
  if (is.null(flip)) {
    flip <- max(nchar(as.character(df$ab)), na.rm = TRUE) > 20
  }

  # ordering by coverage, applied within facet/group so each panel ranks correctly
  if (isTRUE(sort)) {
    split_keys <- interaction(
      if (has_facet) as.character(df$mo) else rep("", nrow(df)),
      if (has_syndromic) df$syndromic_group else rep("", nrow(df)),
      drop = TRUE
    )
    # build a within-group rank, then a global ordered factor whose level order
    # respects that rank; reorder_within-style without the tidytext dependency
    ord <- order(split_keys, df$.coverage)
    df <- df[ord, , drop = FALSE]
    df$ab <- factor(df$ab, levels = unique(df$ab[order(split_keys[ord], df$.coverage[ord])]))
    # note: with multiple facets the level order is a compromise (one global
    # axis), acceptable because each facet shows its own subset in coverage order
  }

  fill_var <- if (has_syndromic) "syndromic_group" else NULL

  out <- ggplot2::ggplot(
    df,
    mapping = ggplot2::aes(
      x = ab,
      y = .coverage,
      fill = if (has_syndromic) syndromic_group else NULL,
      colour = if (has_syndromic) syndromic_group else NULL
    )
  )

  dodge <- ggplot2::position_dodge2(preserve = "single", width = 0.6)

  if (geom == "col") {
    out <- out + ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single"))
    if (isTRUE(ci)) {
      out <- out + ggplot2::geom_errorbar(
        mapping = ggplot2::aes(ymin = .lower, ymax = .upper),
        position = ggplot2::position_dodge2(preserve = "single", width = 1),
        width = 0.7
      )
    }
  } else if (geom == "point") {
    out <- out + ggplot2::geom_point(position = dodge, size = 2)
    if (isTRUE(ci)) {
      out <- out + ggplot2::geom_errorbar(
        mapping = ggplot2::aes(ymin = .lower, ymax = .upper),
        position = dodge, width = 0.4
      )
    }
  } else if (geom == "errorbar") {
    out <- out + ggplot2::geom_errorbar(
      mapping = ggplot2::aes(ymin = .lower, ymax = .upper),
      position = dodge, width = 0.4
    )
  } else {
    # pointrange (default)
    out <- out + ggplot2::geom_pointrange(
      mapping = ggplot2::aes(ymin = .lower, ymax = .upper),
      position = dodge, size = 0.5
    )
  }

  if (is.null(caption)) {
    if (is_wisca) {
      out <- out + ggplot2::labs(caption = "Overlapping credible intervals:\nclinically non-inferior (Bielicki 2020)")
    }
  } else if (!caption %in% c(FALSE, NA)) {
    out <- out + ggplot2::labs(caption = caption)
  }

  out <- out +
    ggplot2::labs(
      y = ifelse(combine_SI, "%SI", "%S"),
      x = NULL,
      fill = if (has_syndromic) group_name else NULL,
      colour = if (has_syndromic) group_name else NULL
    )

  if (isTRUE(flip)) {
    out <- out + ggplot2::coord_flip()
  }

  if (has_facet) {
    out <- out + ggplot2::facet_wrap("mo")
  }

  out
}

#' @param wisca_plot_type Either `"susceptibility_incidence"` (default) or `"posterior_coverage"`.
#' @param ... Currently unused.
#' @rdname antibiogram
#' @export
wisca_plot <- function(wisca_model,
                       wisca_plot_type = c("susceptibility_incidence", "posterior_coverage"),
                       ...) {
  stop_ifnot_installed("ggplot2")
  stop_ifnot(
    isTRUE(attributes(wisca_model)$wisca),
    "This function only applies to WISCA models."
  )
  meet_criteria(wisca_plot_type, allow_class = "character", has_length = 1, is_in = c("susceptibility_incidence", "posterior_coverage"))
  wisca_plot_type <- match.arg(wisca_plot_type)

  sep <- attributes(wisca_model)$sep %||% " + "

  if (wisca_plot_type == "posterior_coverage") {
    plot_wisca_posterior_coverage(wisca_model, sep = sep)
  } else {
    plot_wisca_susceptibility_incidence(wisca_model, sep = sep)
  }
}

# ---- posterior_coverage ----
plot_wisca_posterior_coverage <- function(wisca_model, sep) {
  draws <- attributes(wisca_model)$wisca_draws
  stop_if(
    is.null(draws),
    "No simulation draws found. Re-run {.fun wisca} with the latest AMR version to retain draws."
  )

  if (!is.null(sep)) {
    names(draws) <- gsub(sep, paste0(trimws(sep, which = "right"), "\n"), names(draws), fixed = TRUE)
  }

  df <- do.call(rbind, lapply(names(draws), function(nm) {
    data.frame(regimen = nm, coverage = draws[[nm]] * 100, stringsAsFactors = FALSE)
  }))

  medians <- tapply(df$coverage, df$regimen, stats::median)
  df$regimen <- factor(df$regimen, levels = names(sort(medians, decreasing = TRUE)))

  ggplot2::ggplot(df, ggplot2::aes(x = coverage, fill = regimen, colour = regimen)) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.7) +
    ggplot2::scale_y_continuous(n.breaks = 5, expand = ggplot2::expansion(c(0, 0.05))) +
    ggplot2::scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      n.breaks = 5,
      limits = c(NA, 100)
    ) +
    ggplot2::labs(
      title = "WISCA",
      subtitle = "Posteriors coverage",
      x = "Coverage",
      y = "Relative likelihood",
      fill = translate_AMR("Regimen", language = get_AMR_locale()),
      colour = translate_AMR("Regimen", language = get_AMR_locale())
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.key.spacing.y = ggplot2::unit(0.25, "lines"),
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
}

# ---- susceptibility_incidence scatter, faceted by regimen ----
plot_wisca_susceptibility_incidence <- function(wisca_model, sep) {
  components <- attributes(wisca_model)$wisca_components
  stop_if(
    is.null(components),
    "No simulation components found. Re-run {.fun wisca} with the latest AMR version to retain draws."
  )

  df_list <- lapply(names(components), function(g) {
    comp <- components[[g]]
    n_sims <- nrow(comp$incidence)
    n_path <- ncol(comp$incidence)
    mo_names <- colnames(comp$incidence)

    reg_label <- g
    if (!is.null(sep)) {
      reg_label <- gsub(sep, paste0(trimws(sep, which = "right"), "\n"), g, fixed = TRUE)
    }

    data.frame(
      regimen = rep(reg_label, n_sims * n_path),
      pathogen = rep(mo_names, each = n_sims),
      incidence = as.vector(comp$incidence) * 100,
      susceptibility = as.vector(comp$susceptibility) * 100,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, df_list)
  df$pathogen <- mo_shortname(df$pathogen, keep_synonyms = TRUE, info = FALSE)

  # order pathogens by median incidence across all regimens
  med_inc <- tapply(df$incidence, df$pathogen, stats::median)
  df$pathogen <- factor(df$pathogen, levels = names(sort(med_inc, decreasing = TRUE)))

  # order regimens by median coverage (from wisca_draws)
  draws <- attributes(wisca_model)$wisca_draws
  if (!is.null(draws)) {
    med_cov <- vapply(names(draws), function(g) stats::median(draws[[g]]), double(1))
    reg_labels <- unique(df$regimen)
    # match order: draws names -> display labels
    draw_order <- names(sort(med_cov, decreasing = TRUE))
    label_order <- vapply(draw_order, function(g) {
      if (!is.null(sep)) gsub(sep, paste0(trimws(sep, which = "right"), "\n"), g, fixed = TRUE) else g
    }, character(1))
    label_order <- unique(label_order[label_order %in% reg_labels])
    df$regimen <- factor(df$regimen, levels = label_order)
  }

  # retrieve coverage + CI
  coverage <- attributes(wisca_model)$long_numeric
  if (!is.null(sep)) {
    coverage$ab <- gsub(sep, paste0(trimws(sep, which = "right"), "\n"), coverage$ab, fixed = TRUE)
  }
  df$coverage <- coverage$coverage[match(df$regimen, coverage$ab)] * 100
  df$lower_ci <- coverage$lower_ci[match(df$regimen, coverage$ab)] * 100
  df$upper_ci <- coverage$upper_ci[match(df$regimen, coverage$ab)] * 100
  ci_df <- df[!duplicated(df$regimen), c("regimen", "coverage", "lower_ci", "upper_ci"), drop = FALSE]

  ggplot2::ggplot(df, ggplot2::aes(x = susceptibility, y = incidence, colour = pathogen)) +
    ggplot2::geom_rect(
      data = ci_df,
      ggplot2::aes(xmin = lower_ci, xmax = upper_ci, ymin = -Inf, ymax = Inf),
      fill = "grey50", alpha = 0.15, colour = NA,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(
      data = ci_df,
      ggplot2::aes(xintercept = coverage),
      linewidth = 0.5,
      linetype = 2,
      colour = "grey50",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(size = 0.5, alpha = 0.2, shape = 16) +
    ggplot2::facet_wrap(~regimen) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      n.breaks = 5
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(0, 100)
    ) +
    ggplot2::labs(
      title = "WISCA",
      subtitle = "Susceptibility vs. incidence weight",
      x = translate_AMR("Susceptibility", language = get_AMR_locale()),
      y = translate_AMR("Incidence weight (normalised)", language = get_AMR_locale()),
      colour = translate_AMR("Pathogen", language = get_AMR_locale()),
      caption = paste(attributes(wisca_model)$simulations, "Monte Carlo simulations")
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(alpha = 1, size = 3)
      )
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.key.spacing.y = ggplot2::unit(0.25, "lines"),
      legend.text = ggplot2::element_text(face = "italic"),
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
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
