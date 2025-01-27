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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Generate Traditional, Combination, Syndromic, or WISCA Antibiograms
#'
#' @description
#' Create detailed antibiograms with options for traditional, combination, syndromic, and Bayesian WISCA methods.
#' 
#' Adhering to previously described approaches (see *Source*) and especially the Bayesian WISCA model (Weighted-Incidence Syndromic Combination Antibiogram) by Bielicki *et al.*, these functions provides flexible output formats including plots and tables, ideal for integration with R Markdown and Quarto reports.
#' @param x a [data.frame] containing at least a column with microorganisms and columns with antimicrobial results (class 'sir', see [as.sir()])
#' @param antibiotics vector of any antimicrobial name or code (will be evaluated with [as.ab()], column name of `x`, or (any combinations of) [antimicrobial selectors][antimicrobial_class_selectors] such as [aminoglycosides()] or [carbapenems()]. For combination antibiograms, this can also be set to values separated with `"+"`, such as "TZP+TOB" or "cipro + genta", given that columns resembling such antimicrobials exist in `x`. See *Examples*.
#' @param mo_transform a character to transform microorganism input - must be `"name"`, `"shortname"` (default), `"gramstain"`, or one of the column names of the [microorganisms] data set: `r vector_or(colnames(microorganisms), sort = FALSE, quotes = TRUE)`. Can also be `NULL` to not transform the input.
#' @param ab_transform a character to transform antimicrobial input - must be one of the column names of the [antibiotics] data set (defaults to `"name"`): `r vector_or(colnames(antibiotics), sort = FALSE, quotes = TRUE)`. Can also be `NULL` to not transform the input.
#' @param syndromic_group a column name of `x`, or values calculated to split rows of `x`, e.g. by using [ifelse()] or [`case_when()`][dplyr::case_when()]. See *Examples*.
#' @param add_total_n a [logical] to indicate whether total available numbers per pathogen should be added to the table (default is `TRUE`). This will add the lowest and highest number of available isolates per antimicrobial (e.g, if for *E. coli* 200 isolates are available for ciprofloxacin and 150 for amoxicillin, the returned number will be "150-200").
#' @param only_all_tested (for combination antibiograms): a [logical] to indicate that isolates must be tested for all antimicrobials, see *Details*
#' @param digits number of digits to use for rounding the susceptibility percentage
#' @param formatting_type numeric value (1–22 for WISCA, 1-12 for non-WISCA) indicating how the 'cells' of the antibiogram table should be formatted. See *Details* > *Formatting Type* for a list of options.
#' @param col_mo column name of the names or codes of the microorganisms (see [as.mo()]) - the default is the first column of class [`mo`]. Values will be coerced using [as.mo()].
#' @param language language to translate text, which defaults to the system language (see [get_AMR_locale()])
#' @param minimum the minimum allowed number of available (tested) isolates. Any isolate count lower than `minimum` will return `NA` with a warning. The default number of `30` isolates is advised by the Clinical and Laboratory Standards Institute (CLSI) as best practice, see *Source*.
#' @param combine_SI a [logical] to indicate whether all susceptibility should be determined by results of either S, SDD, or I, instead of only S (default is `TRUE`)
#' @param sep a separating character for antimicrobial columns in combination antibiograms
#' @param wisca a [logical] to indicate whether a Weighted-Incidence Syndromic Combination Antibiogram (WISCA) must be generated (default is `FALSE`). This will use a Bayesian hierarchical model to estimate regimen coverage probabilities using Montecarlo simulations. Set `simulations` to adjust.
#' @param simulations (for WISCA) a numerical value to set the number of Montecarlo simulations
#' @param conf_interval (for WISCA) a numerical value to set confidence interval (default is `0.95`)
#' @param interval_side (for WISCA) the side of the confidence interval, either `"two-tailed"` (default), `"left"` or `"right"`
#' @param info 	a [logical] to indicate info should be printed - the default is `TRUE` only in interactive mode
#' @param object an [antibiogram()] object
#' @param ... when used in [R Markdown or Quarto][knitr::kable()]: arguments passed on to [knitr::kable()] (otherwise, has no use)
#' @details This function returns a table with values between 0 and 100 for *susceptibility*, not resistance.
#' 
#' **Remember that you should filter your data to let it contain only first isolates!** This is needed to exclude duplicates and to reduce selection bias. Use [first_isolate()] to determine them in your data set with one of the four available algorithms.
#' 
#' For estimating antimicrobial coverage, especially when creating a WISCA, the outcome might become more reliable by only including the top *n* species encountered in the data. You can filter on this top *n* using [top_n_microorganisms()]. For example, use `top_n_microorganisms(your_data, n = 10)` as a pre-processing step to only include the top 10 species in the data.
#' 
#' The numeric values of an antibiogram are stored in a long format as the [attribute] `long_numeric`. You can retrieve them using `attributes(x)$long_numeric`, where `x` is the outcome of [antibiogram()] or [wisca()]. This is ideal for e.g. advanced plotting.
#' 
#' ### Formatting Type
#' 
#' The formatting of the 'cells' of the table can be set with the argument `formatting_type`. In these examples, `5` is the susceptibility percentage (for WISCA: `4-6` indicates the confidence level), `15` the numerator, and `300` the denominator:
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
#' 10. 5% (15/300) - **default for non-WISCA**
#' 11. 5 (N=15/300)
#' 12. 5% (N=15/300)
#'     
#'     Additional options for WISCA (using `antibiogram(..., wisca = TRUE)` or `wisca()`):
#' 13. 5 (4-6)
#' 14. 5% (4-6%)
#' 15. 5 (4-6,300)
#' 16. 5% (4-6%,300)
#' 17. 5 (4-6,N=300)
#' 18. 5% (4-6%,N=300) - **default for WISCA**
#' 19. 5 (4-6,15/300)
#' 20. 5% (4-6%,15/300)
#' 21. 5 (4-6,N=15/300)
#' 22. 5% (4-6%,N=15/300)
#' 
#' The default is `18` for WISCA and `10` for non-WISCA, which can be set globally with the package option [`AMR_antibiogram_formatting_type`][AMR-options], e.g. `options(AMR_antibiogram_formatting_type = 5)`.
#' 
#' Set `digits` (defaults to `0`) to alter the rounding of the susceptibility percentages.
#'
#' ### Antibiogram Types
#'
#' There are various antibiogram types, as summarised by Klinker *et al.* (2021, \doi{10.1177/20499361211011373}), and they are all supported by [antibiogram()].
#' 
#' **Use WISCA whenever possible**, since it provides more precise coverage estimates by accounting for pathogen incidence and antimicrobial susceptibility, as has been shown by Bielicki *et al.* (2020, \doi{10.1001.jamanetworkopen.2019.21124}). See the section *Why Use WISCA?* on this page.
#'
#' 1. **Traditional Antibiogram**
#'
#'    Case example: Susceptibility of *Pseudomonas aeruginosa* to piperacillin/tazobactam (TZP)
#'
#'    Code example:
#'
#'    ```r
#'    antibiogram(your_data,
#'                antibiotics = "TZP")
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
#'                antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"))
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
#'                antibiotics = penicillins(),
#'                syndromic_group = "ward")
#'    ```
#'
#' 4. **Weighted-Incidence Syndromic Combination Antibiogram (WISCA)**
#' 
#'    WISCA can be applied to any antibiogram, see the section *Why Use WISCA?* on this page for more information.
#'
#'    Code example:
#'
#'    ```r
#'    antibiogram(your_data,
#'                antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"),
#'                wisca = TRUE)
#'                
#'    # this is equal to:
#'    wisca(your_data,
#'          antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"))
#'    ```
#'    
#'    WISCA uses a sophisticated Bayesian decision model to combine both local and pooled antimicrobial resistance data. This approach not only evaluates local patterns but can also draw on multi-centre datasets to improve regimen accuracy, even in low-incidence infections like paediatric bloodstream infections (BSIs).
#'    
#' Grouped [tibbles][tibble::tibble] can also be used to calculate susceptibilities over various groups.
#' 
#' Code example:
#'
#' ```r
#' your_data %>%
#'   group_by(has_sepsis, is_neonate, sex) %>%
#'   wisca(antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"))
#' ```
#'
#' ### Inclusion in Combination Antibiogram and Syndromic Antibiogram
#'
#' Note that for types 2 and 3 (Combination Antibiogram and Syndromic Antibiogram), it is important to realise that susceptibility can be calculated in two ways, which can be set with the `only_all_tested` argument (default is `FALSE`). See this example for two antimicrobials, Drug A and Drug B, about how [antibiogram()] works to calculate the %SI:
#'
#' ```
#' --------------------------------------------------------------------
#'                     only_all_tested = FALSE  only_all_tested = TRUE
#'                     -----------------------  -----------------------
#'  Drug A    Drug B   include as  include as   include as  include as
#'                     numerator   denominator  numerator   denominator
#' --------  --------  ----------  -----------  ----------  -----------
#'  S or I    S or I       X            X            X            X
#'    R       S or I       X            X            X            X
#'   <NA>     S or I       X            X            -            -
#'  S or I      R          X            X            X            X
#'    R         R          -            X            -            X
#'   <NA>       R          -            -            -            -
#'  S or I     <NA>        X            X            -            -
#'    R        <NA>        -            -            -            -
#'   <NA>      <NA>        -            -            -            -
#' --------------------------------------------------------------------
#' ```
#' 
#' ### Plotting
#' 
#' All types of antibiograms as listed above can be plotted (using [ggplot2::autoplot()] or base \R's [plot()] and [barplot()]).
#' 
#' THe outcome of [antibiogram()] can also be used directly in R Markdown / Quarto (i.e., `knitr`) for reports. In this case, [knitr::kable()] will be applied automatically and microorganism names will even be printed in italics at default (see argument `italicise`).
#' 
#' You can also use functions from specific 'table reporting' packages to transform the output of [antibiogram()] to your needs, e.g. with `flextable::as_flextable()` or `gt::gt()`.
#'
#' @section Why Use WISCA?:
#' WISCA, as outlined by Barbieri *et al.* (\doi{10.1186/s13756-021-00939-2}), stands for 
#' Weighted-Incidence Syndromic Combination Antibiogram, which estimates the probability 
#' of adequate empirical antimicrobial regimen coverage for specific infection syndromes. 
#' This method leverages a Bayesian hierarchical logistic regression framework with random 
#' effects for pathogens and regimens, enabling robust estimates in the presence of sparse 
#' data. 
#'
#' The Bayesian model assumes conjugate priors for parameter estimation. For example, the 
#' coverage probability \ifelse{latex}{\deqn{$theta$}}{$theta$} for a given antimicrobial regimen 
#' is modeled using a Beta distribution as a prior: 
#'
#' \ifelse{latex}{\deqn{$theta$ \sim \text{Beta}($alpha$_0, $beta$_0)}}{
#' \ifelse{html}{\figure{beta_prior.png}{options: width="300" alt="Beta prior"}}{$theta$ ~ Beta($alpha$_0, $beta$_0)}}
#'
#' where \eqn{$alpha$_0} and \eqn{$beta$_0} represent prior successes and failures, respectively, 
#' informed by expert knowledge or weakly informative priors (e.g., \eqn{$alpha$_0 = 1, $beta$_0 = 1}).
#' 
#' The likelihood function is constructed based on observed data, where the number of covered 
#' cases for a regimen follows a binomial distribution:
#'
#' \ifelse{latex}{\deqn{y \sim \text{Binomial}(n, $theta$)}}{
#' \ifelse{html}{\figure{binomial_likelihood.png}{options: width="300" alt="Binomial likelihood"}}{y ~ Binomial(n, $theta$)}}
#'
#' Posterior parameter estimates are obtained by combining the prior and likelihood using 
#' Bayes' theorem. The posterior distribution of \eqn{$theta$} is also a Beta distribution:
#'
#' \ifelse{latex}{\deqn{$theta$ | y \sim \text{Beta}($alpha$_0 + y, $beta$_0 + n - y)}}{
#' \ifelse{html}{\figure{posterior_beta.png}{options: width="300" alt="Beta posterior"}}{$theta$ | y ~ Beta($alpha$_0 + y, $beta$_0 + n - y)}}
#'
#' For hierarchical modeling, pathogen-level effects (e.g., differences in resistance 
#' patterns) and regimen-level effects are modelled using Gaussian priors on log-odds. 
#' This hierarchical structure ensures partial pooling of estimates across groups, 
#' improving stability in strata with small sample sizes. The model is implemented using 
#' Hamiltonian Monte Carlo (HMC) sampling.
#'
#' Stratified results are provided based on covariates such as age, sex, and clinical 
#' complexity (e.g., prior antimicrobial treatments or renal/urological comorbidities). 
#' For example, posterior odds ratios (ORs) are derived to quantify the effect of these 
#' covariates on coverage probabilities:
#'
#' \ifelse{latex}{\deqn{\text{OR}_{\text{covariate}} = \frac{\exp($beta$_{\text{covariate}})}{\exp($beta$_0)}}}{
#' \ifelse{html}{\figure{odds_ratio.png}{options: width="300" alt="Odds ratio formula"}}{OR_covariate = exp(beta_covariate) / exp(beta_0)}}
#'
#' By combining empirical data with prior knowledge, WISCA overcomes the limitations 
#' of traditional combination antibiograms, offering disease-specific, patient-stratified 
#' estimates with robust uncertainty quantification. This tool is invaluable for antimicrobial 
#' stewardship programs and empirical treatment guideline refinement.
#' 
#' @source
#' * Bielicki JA *et al.* (2016). **Selecting appropriate empirical antibiotic regimens for paediatric bloodstream infections: application of a Bayesian decision model to local and pooled antimicrobial resistance surveillance data** *Journal of Antimicrobial Chemotherapy* 71(3); \doi{10.1093/jac/dkv397}
#' * Bielicki JA *et al.* (2020). **Evaluation of the coverage of 3 antibiotic regimens for neonatal sepsis in the hospital setting across Asian countries** *JAMA Netw Open.* 3(2):e1921124; \doi{10.1001.jamanetworkopen.2019.21124}
#' * Klinker KP *et al.* (2021). **Antimicrobial stewardship and antibiograms: importance of moving beyond traditional antibiograms**. *Therapeutic Advances in Infectious Disease*, May 5;8:20499361211011373; \doi{10.1177/20499361211011373}
#' * Barbieri E *et al.* (2021). **Development of a Weighted-Incidence Syndromic Combination Antibiogram (WISCA) to guide the choice of the empiric antibiotic treatment for urinary tract infection in paediatric patients: a Bayesian approach** *Antimicrobial Resistance & Infection Control* May 1;10(1):74; \doi{10.1186/s13756-021-00939-2}
#' * **M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 5th Edition**, 2022, *Clinical and Laboratory Standards Institute (CLSI)*. <https://clsi.org/standards/products/microbiology/documents/m39/>.
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
#'   antibiotics = c(aminoglycosides(), carbapenems())
#' )
#'
#' antibiogram(example_isolates,
#'   antibiotics = aminoglycosides(),
#'   ab_transform = "atc",
#'   mo_transform = "gramstain"
#' )
#'
#' antibiogram(example_isolates,
#'   antibiotics = carbapenems(),
#'   ab_transform = "name",
#'   mo_transform = "name"
#' )
#'
#'
#' # Combined antibiogram -------------------------------------------------
#'
#' # combined antibiotics yield higher empiric coverage
#' antibiogram(example_isolates,
#'   antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"),
#'   mo_transform = "gramstain"
#' )
#'
#' # names of antibiotics do not need to resemble columns exactly:
#' antibiogram(example_isolates,
#'   antibiotics = c("Cipro", "cipro + genta"),
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
#'   antibiotics = c(aminoglycosides(), carbapenems()),
#'   syndromic_group = "ward"
#' )
#'
#' # now define a data set with only E. coli
#' ex1 <- example_isolates[which(mo_genus() == "Escherichia"), ]
#'
#' # with a custom language, though this will be determined automatically
#' # (i.e., this table will be in Spanish on Spanish systems)
#' antibiogram(ex1,
#'   antibiotics = aminoglycosides(),
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
#' # can be used for any of the above types - just add `wisca = TRUE`
#' antibiogram(example_isolates,
#'   antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"),
#'   mo_transform = "gramstain",
#'   wisca = TRUE
#' )
#'
#'
#' # Print the output for R Markdown / Quarto -----------------------------
#'
#' ureido <- antibiogram(example_isolates,
#'   antibiotics = ureidopenicillins(),
#'   ab_transform = "name"
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
#'   antibiotics = c("AMC", "CIP", "TZP", "TZP+TOB"),
#'   mo_transform = "gramstain",
#'   wisca = TRUE
#' )
#' ab2 <- antibiogram(example_isolates,
#'   antibiotics = c("AMC", "CIP", "TZP", "TZP+TOB"),
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
                        antibiotics = where(is.sir),
                        mo_transform = "shortname",
                        ab_transform = "name",
                        syndromic_group = NULL,
                        add_total_n = FALSE,
                        only_all_tested = FALSE,
                        digits = 0,
                        formatting_type = getOption("AMR_antibiogram_formatting_type", ifelse(wisca, 18, 10)),
                        col_mo = NULL,
                        language = get_AMR_locale(),
                        minimum = 30,
                        combine_SI = TRUE,
                        sep = " + ",
                        wisca = FALSE,
                        simulations = 1000,
                        conf_interval = 0.95,
                        interval_side = "two-tailed",
                        info = interactive()) {
  UseMethod("antibiogram")
}

#' @method antibiogram default
#' @export
antibiogram.default <- function(x,
                                antibiotics = where(is.sir),
                                mo_transform = "shortname",
                                ab_transform = "name",
                                syndromic_group = NULL,
                                add_total_n = FALSE,
                                only_all_tested = FALSE,
                                digits = 0,
                                formatting_type = getOption("AMR_antibiogram_formatting_type", ifelse(wisca, 18, 10)),
                                col_mo = NULL,
                                language = get_AMR_locale(),
                                minimum = 30,
                                combine_SI = TRUE,
                                sep = " + ",
                                wisca = FALSE,
                                simulations = 1000,
                                conf_interval = 0.95,
                                interval_side = "two-tailed",
                                info = interactive()) {
  meet_criteria(x, allow_class = "data.frame")
  x <- ascertain_sir_classes(x, "x")
  if (!is.function(mo_transform)) {
    meet_criteria(mo_transform, allow_class = "character", has_length = 1, is_in = c("name", "shortname", "gramstain", colnames(AMR::microorganisms)), allow_NULL = TRUE)
  }
  if (!is.function(ab_transform)) {
    meet_criteria(ab_transform, allow_class = "character", has_length = 1, is_in = colnames(AMR::antibiotics), allow_NULL = TRUE)
  }
  meet_criteria(syndromic_group, allow_class = "character", allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(add_total_n, allow_class = "logical", has_length = 1)
  meet_criteria(only_all_tested, allow_class = "logical", has_length = 1)
  meet_criteria(digits, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE)
  meet_criteria(wisca, allow_class = "logical", has_length = 1)
  meet_criteria(formatting_type, allow_class = c("numeric", "integer"), has_length = 1, is_in = if (wisca == TRUE) c(1:22) else c(1:12))
  meet_criteria(col_mo, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  language <- validate_language(language)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(sep, allow_class = "character", has_length = 1)
  meet_criteria(simulations, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE, is_positive = TRUE)
  meet_criteria(conf_interval, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE, is_positive = TRUE)
  meet_criteria(interval_side, allow_class = "character", has_length = 1, is_in = c("two-tailed", "left", "right"))
  meet_criteria(info, allow_class = "logical", has_length = 1)
  
  # try to find columns based on type
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = interactive())
    stop_if(is.null(col_mo), "`col_mo` must be set")
  }
  # transform MOs
  x$`.mo` <- x[, col_mo, drop = TRUE]
  if (is.null(mo_transform)) {
    # leave as is
  } else if (is.function(mo_transform)) {
    x$`.mo` <- mo_transform(x$`.mo`)
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
  
  # get antibiotics
  ab_trycatch <- tryCatch(colnames(suppressWarnings(x[, antibiotics, drop = FALSE])), error = function(e) NULL)
  if (is.null(ab_trycatch)) {
    stop_ifnot(is.character(suppressMessages(antibiotics)), "`antibiotics` must be an antimicrobial selector, or a character vector.")
    antibiotics.bak <- antibiotics
    # split antibiotics on separator and make it a list
    antibiotics <- strsplit(gsub(" ", "", antibiotics), "+", fixed = TRUE)
    # get available antibiotics in data set
    df_ab <- get_column_abx(x, verbose = FALSE, info = FALSE)
    # get antibiotics from user
    user_ab <- suppressMessages(suppressWarnings(lapply(antibiotics, as.ab, flag_multiple_results = FALSE, info = FALSE)))
    non_existing <- character(0)
    user_ab <- lapply(user_ab, function(x) {
      out <- unname(df_ab[match(x, names(df_ab))])
      non_existing <<- c(non_existing, x[is.na(out) & !is.na(x)])
      # remove non-existing columns
      out[!is.na(out)]
    })
    user_ab <- user_ab[unlist(lapply(user_ab, length)) > 0]
    
    if (length(non_existing) > 0) {
      warning_("The following antibiotics were not available and ignored: ", vector_and(ab_name(non_existing, language = NULL, tolower = TRUE), quotes = FALSE))
    }
    
    # make list unique
    antibiotics <- unique(user_ab)
    # go through list to set AMR in combinations
    for (i in seq_len(length(antibiotics))) {
      abx <- antibiotics[[i]]
      for (ab in abx) {
        # make sure they are SIR columns
        x[, ab] <- as.sir(x[, ab, drop = TRUE])
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
        other_values <- setdiff(c("S", "SDD", "I", "R", "NI"), S_values)
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
      antibiotics[[i]] <- new_colname
    }
    antibiotics <- unlist(antibiotics)
  } else {
    antibiotics <- ab_trycatch
  }
  
  if (isTRUE(has_syndromic_group)) {
    out <- x %pm>%
      pm_select(.syndromic_group, .mo, antibiotics) %pm>%
      pm_group_by(.syndromic_group)
  } else {
    out <- x %pm>%
      pm_select(.mo, antibiotics)
  }
  
  
  # get numbers of S, I, R (per group)
  out <- out %pm>%
    bug_drug_combinations(
      col_mo = ".mo",
      FUN = function(x) x,
      include_n_rows = TRUE
    )
  counts <- out
  
  
  if (wisca == TRUE) {
    # WISCA ----
    
    # set up progress bar
    progress <- progress_ticker(n = NROW(out[which(out$total > 0), , drop = FALSE]),
                                n_min = 10,
                                print = info,
                                title = "Calculating beta/gamma parameters for WISCA")
    on.exit(close(progress))
    
    out$percentage = NA_real_
    out$lower = NA_real_
    out$upper = NA_real_
    
    for (i in seq_len(NROW(out))) {
      if (out$total[i] == 0) {
        next
      }
      progress$tick()
      out_current <- out[i, , drop = FALSE]
      priors <- calculate_priors(out_current, combine_SI = combine_SI)
      
      # Monte Carlo simulation
      coverage_simulations <- replicate(simulations, {
        # simulate pathogen incidence
        # = Dirichlet (Gamma) parameters
        simulated_incidence <- stats::rgamma(
          n = length(priors$gamma_posterior),
          shape = priors$gamma_posterior,
          rate = 1 # Scale = 1 for gamma
        )
        # normalise
        simulated_incidence <- simulated_incidence / sum(simulated_incidence)
        
        # simulate susceptibility
        # = Beta parameters
        simulated_susceptibility <- stats::rbeta(
          n = length(priors$beta_posterior_1),
          shape1 = priors$beta_posterior_1,
          shape2 = priors$beta_posterior_2
        )
        sum(simulated_incidence * simulated_susceptibility)
      })
      
      # calculate coverage statistics
      coverage_mean <- mean(coverage_simulations)
      if (interval_side == "two-tailed") {
        probs <- c((1 - conf_interval) / 2, 1 - (1 - conf_interval) / 2)
      } else if (interval_side == "left") {
        probs <- c(0, conf_interval)
      } else if (interval_side == "right") {
        probs <- c(1 - conf_interval, 1)
      }
      coverage_ci <- unname(stats::quantile(coverage_simulations, probs = probs))
      
      out$percentage[i] <- coverage_mean
      out$lower[i] <- coverage_ci[1]
      out$upper[i] <- coverage_ci[2]
    }
    # remove progress bar from console
    close(progress)
  }
  
  if (isTRUE(combine_SI)) {
    out$numerator <- out$S + out$I + out$SDD
  } else {
    out$numerator <- out$S
  }
  if (all(out$total < minimum, na.rm = TRUE) && wisca == FALSE) {
    warning_("All combinations had less than `minimum = ", minimum, "` results, returning an empty antibiogram")
    return(as_original_data_class(data.frame(), class(out), extra_class = "antibiogram"))
  } else if (any(out$total < minimum, na.rm = TRUE)) {
    out <- out %pm>%
      # also for WISCA, refrain from anything below 15 isolates:
      subset(total > 15)
    mins <- sum(out$total < minimum, na.rm = TRUE)
    if (wisca == FALSE) {
      out <- out %pm>%
        subset(total >= minimum)
      if (isTRUE(info) && mins > 0) {
        message_("NOTE: ", mins, " combinations had less than `minimum = ", minimum, "` results and were ignored", add_fn = font_red)
      }
    } else if (isTRUE(info)) {
      warning_("Number of tested isolates per regimen should exceed ", minimum, ". Coverage estimates will be inaccurate for ", mins, " regimen", ifelse(mins == 1, "", "s"), ".", call = FALSE)
    }
  }
  
  if (NROW(out) == 0) {
    return(as_original_data_class(data.frame(), class(out), extra_class = "antibiogram"))
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
  
  if (wisca == TRUE) {
    long_numeric <- out %pm>%
      pm_summarise(percentage = percentage,
                   lower = lower,
                   upper = upper,
                   numerator = numerator,
                   total = total)
  } else {
    long_numeric <- out %pm>%
      pm_summarise(percentage = numerator / total,
                   numerator = numerator,
                   total = total)  
  }
  
  out$digits <- digits # since pm_sumarise() cannot work with an object outside the current frame
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
  if (formatting_type == 1) out <- out %pm>% pm_summarise(out_value = round((numerator / total) * 100, digits = digits))
  if (formatting_type == 2) out <- out %pm>% pm_summarise(out_value = numerator)
  if (formatting_type == 3) out <- out %pm>% pm_summarise(out_value = total)
  if (formatting_type == 4) out <- out %pm>% pm_summarise(out_value = paste0(numerator, "/", total))
  if (formatting_type == 5) out <- out %pm>% pm_summarise(out_value = paste0(round((numerator / total) * 100, digits = digits), " (", total, ")"))
  if (formatting_type == 6) out <- out %pm>% pm_summarise(out_value = paste0(round((numerator / total) * 100, digits = digits), "% (", total, ")"))
  if (formatting_type == 7) out <- out %pm>% pm_summarise(out_value = paste0(round((numerator / total) * 100, digits = digits), " (N=", total, ")"))
  if (formatting_type == 8) out <- out %pm>% pm_summarise(out_value = paste0(round((numerator / total) * 100, digits = digits), "% (N=", total, ")"))
  if (formatting_type == 9) out <- out %pm>% pm_summarise(out_value = paste0(round((numerator / total) * 100, digits = digits), " (", numerator, "/", total, ")"))
  if (formatting_type == 10) out <- out %pm>% pm_summarise(out_value = paste0(round((numerator / total) * 100, digits = digits), "% (", numerator, "/", total, ")"))
  if (formatting_type == 11) out <- out %pm>% pm_summarise(out_value = paste0(round((numerator / total) * 100, digits = digits), " (N=", numerator, "/", total, ")"))
  if (formatting_type == 12) out <- out %pm>% pm_summarise(out_value = paste0(round((numerator / total) * 100, digits = digits), "% (N=", numerator, "/", total, ")"))
  if (formatting_type == 13) out <- out %pm>% pm_summarise(out_value = paste0(round(percentage * 100, digits = digits), " (", round(lower * 100, digits = digits), "-", round(upper * 100, digits = digits), ")"))
  if (formatting_type == 14) out <- out %pm>% pm_summarise(out_value = paste0(round(percentage * 100, digits = digits), "% (", round(lower * 100, digits = digits), "-", round(upper * 100, digits = digits), "%)"))
  if (formatting_type == 15) out <- out %pm>% pm_summarise(out_value = paste0(round(percentage * 100, digits = digits), " (", round(lower * 100, digits = digits), "-", round(upper * 100, digits = digits), ",", total, ")"))
  if (formatting_type == 16) out <- out %pm>% pm_summarise(out_value = paste0(round(percentage * 100, digits = digits), "% (", round(lower * 100, digits = digits), "-", round(upper * 100, digits = digits), "%,", total, ")"))
  if (formatting_type == 17) out <- out %pm>% pm_summarise(out_value = paste0(round(percentage * 100, digits = digits), " (", round(lower * 100, digits = digits), "-", round(upper * 100, digits = digits), ",N=", total, ")"))
  if (formatting_type == 18) out <- out %pm>% pm_summarise(out_value = paste0(round(percentage * 100, digits = digits), "% (", round(lower * 100, digits = digits), "-", round(upper * 100, digits = digits), "%,N=", total, ")"))
  
  # transform names of antibiotics
  ab_naming_function <- function(x, t, l, s) {
    x <- strsplit(x, s, fixed = TRUE)
    out <- character(length = length(x))
    for (i in seq_len(length(x))) {
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
    object <- object %pm>%
      # an unclassed data.frame is required for stats::reshape()
      as.data.frame(stringsAsFactors = FALSE) %pm>%
      stats::reshape(direction = "wide", idvar = "mo", timevar = "ab", v.names = "out_value")
    colnames(object) <- gsub("^out_value?[.]", "", colnames(object))
    return(object)
  }
  
  # ungroup for long -> wide transformation
  attr(out, "pm_groups") <- NULL
  attr(out, "groups") <- NULL
  class(out) <- class(out)[!class(out) %in% c("grouped_df", "grouped_data")]
  
  if (isTRUE(has_syndromic_group)) {
    grps <- unique(out$syndromic_group)
    for (i in seq_len(length(grps))) {
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
    # sort rows
    new_df <- new_df %pm>% pm_arrange(mo, syndromic_group)
    # sort columns
    new_df <- new_df[, c("syndromic_group", "mo", sort(colnames(new_df)[!colnames(new_df) %in% c("syndromic_group", "mo")])), drop = FALSE]
    colnames(new_df)[1:2] <- translate_AMR(c("Syndromic Group", "Pathogen"), language = language)
  } else {
    new_df <- long_to_wide(out)
    # sort rows
    new_df <- new_df %pm>% pm_arrange(mo)
    # sort columns
    new_df <- new_df[, c("mo", sort(colnames(new_df)[colnames(new_df) != "mo"])), drop = FALSE]
    colnames(new_df)[1] <- translate_AMR("Pathogen", language = language)
  }
  
  # add total N if indicated
  if (isTRUE(add_total_n)) {
    if (isTRUE(has_syndromic_group)) {
      n_per_mo <- counts %pm>%
        pm_group_by(mo, .syndromic_group) %pm>%
        pm_summarise(paste0(min(total, na.rm = TRUE), "-", max(total, na.rm = TRUE)))
      colnames(n_per_mo) <- c("mo", "syn", "count")
      count_group <- n_per_mo$count[match(paste(new_df[[2]], new_df[[1]]), paste(n_per_mo$mo, n_per_mo$syn))]
      edit_col <- 2
    } else {
      n_per_mo <- counts %pm>%
        pm_group_by(mo) %pm>%
        pm_summarise(paste0(min(total, na.rm = TRUE), "-", max(total, na.rm = TRUE)))
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
  
  out <- as_original_data_class(new_df, class(x), extra_class = "antibiogram")
  rownames(out) <- NULL
  rownames(long_numeric) <- NULL
  
  structure(out,
            has_syndromic_group = has_syndromic_group,
            combine_SI = combine_SI,
            wisca = wisca,
            conf_interval = conf_interval,
            long_numeric = as_original_data_class(long_numeric, class(out))
  )
}

#' @method antibiogram grouped_df
#' @export
antibiogram.grouped_df <- function(x,
                                   antibiotics = where(is.sir),
                                   mo_transform = function (...) "no_mo",
                                   ab_transform = "name",
                                   syndromic_group = NULL,
                                   add_total_n = FALSE,
                                   only_all_tested = FALSE,
                                   digits = 0,
                                   formatting_type = getOption("AMR_antibiogram_formatting_type", ifelse(wisca, 18, 10)),
                                   col_mo = NULL,
                                   language = get_AMR_locale(),
                                   minimum = 30,
                                   combine_SI = TRUE,
                                   sep = " + ",
                                   wisca = FALSE,
                                   simulations = 1000,
                                   conf_interval = 0.95,
                                   interval_side = "two-tailed",
                                   info = interactive()) {
  stop_ifnot(is.null(syndromic_group), "`syndromic_group` must not be set if creating an antibiogram using a grouped tibble. The groups will become the variables over which the antimicrobials are calculated, making `syndromic_groups` redundant.", call = FALSE)
  groups <- attributes(x)$groups
  n_groups <- NROW(groups)
  progress <- progress_ticker(n = n_groups,
                              n_min = 5,
                              print = info,
                              title = paste("Calculating AMR for", n_groups, "groups"))
  on.exit(close(progress))
  
  for (i in seq_len(n_groups)) {
    if (i > 1) progress$tick()
    rows <- unlist(groups[i, ]$.rows)
    if (length(rows) == 0) {
      next
    }
    
    new_out <- antibiogram(as.data.frame(x)[rows, , drop = FALSE],
                           antibiotics = antibiotics,
                           mo_transform = function(x) "no_mo",
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
                           wisca = wisca,
                           simulations = simulations,
                           conf_interval = conf_interval,
                           interval_side = interval_side,
                           info = i == 1 && info == TRUE)
    new_long_numeric <- attributes(new_out)$long_numeric
    
    if (i == 1) progress$tick()
    
    if (NROW(new_out) == 0) {
      next
    }
    
    # remove first column 'Pathogen' (in whatever language)
    new_out <- new_out[, -1, drop = FALSE]
    new_long_numeric <- new_long_numeric[, -1, drop = FALSE]
    
    # add group names to data set
    for (col in rev(seq_len(NCOL(groups) - 1))) {
      col_name <- colnames(groups)[col]
      col_value <- groups[i, col, drop = TRUE]
      new_out[, col_name] <- col_value
      new_out <- new_out[, c(col_name, setdiff(names(new_out), col_name))] # set place to 1st col
      new_long_numeric[, col_name] <- col_value
      new_long_numeric <- new_long_numeric[, c(col_name, setdiff(names(new_long_numeric), col_name))] # set place to 1st col
    }
    
    if (i == 1) {
      # the first go
      out <- new_out
      long_numeric <- new_long_numeric
    } else {
      out <- rbind_AMR(out, new_out)
      long_numeric <- rbind_AMR(long_numeric, new_long_numeric)
    }
  }

  close(progress)
  
  out <- structure(as_original_data_class(out, class(x), extra_class = "antibiogram"),
                   has_syndromic_group = FALSE,
                   combine_SI = isTRUE(combine_SI),
                   wisca = isTRUE(wisca),
                   conf_interval = conf_interval,
                   long_numeric = as_original_data_class(long_numeric, class(x)))
}

#' @export
#' @rdname antibiogram
wisca <- function(x,
                  antibiotics = where(is.sir),
                  mo_transform = "shortname",
                  ab_transform = "name",
                  syndromic_group = NULL,
                  add_total_n = FALSE,
                  only_all_tested = FALSE,
                  digits = 0,
                  formatting_type = getOption("AMR_antibiogram_formatting_type", 18),
                  col_mo = NULL,
                  language = get_AMR_locale(),
                  minimum = 30,
                  combine_SI = TRUE,
                  sep = " + ",
                  simulations = 1000,
                  info = interactive()) {
  antibiogram(x = x,
              antibiotics = antibiotics,
              mo_transform = mo_transform,
              ab_transform = ab_transform,
              syndromic_group = syndromic_group,
              add_total_n = add_total_n,
              only_all_tested = only_all_tested,
              digits = digits,
              formatting_type = formatting_type,
              col_mo = col_mo,
              language = language,
              minimum = minimum,
              combine_SI = combine_SI,
              sep = sep,
              wisca = TRUE,
              simulations = simulations,
              info = info)
}

calculate_priors <- function(data, combine_SI = TRUE) {
  # Ensure data has required columns
  stopifnot(all(c("mo", "total_rows", "total", "S") %in% colnames(data)))
  if (combine_SI == TRUE && "I" %in% colnames(data)) {
    data$S <- data$S + data$I
  }
  if (combine_SI == TRUE && "SDD" %in% colnames(data)) {
    data$S <- data$S + data$SDD
  }
  
  # Pathogen incidence (Dirichlet distribution)
  gamma_prior <- rep(1, length(unique(data$mo)))   # Dirichlet prior
  gamma_posterior <- gamma_prior + data$total_rows # Posterior parameters
  
  # Regimen susceptibility (Beta distribution)
  beta_prior <- rep(1, length(unique(data$mo))) # Beta prior
  r <- data$S                                   # Number of pathogens tested susceptible
  n <- data$total                               # Total tested
  beta_posterior_1 <- beta_prior + r            # Posterior alpha
  beta_posterior_2 <- beta_prior + (n - r)      # Posterior beta
  
  # Return parameters as a list
  list(
    gamma_posterior = gamma_posterior,
    beta_posterior_1 = beta_posterior_1,
    beta_posterior_2 = beta_posterior_2
  )
}

# will be exported in R/zzz.R
tbl_sum.antibiogram <- function(x, ...) {
  dims <- paste(format(NROW(x), big.mark = ","), AMR_env$cross_icon, format(NCOL(x), big.mark = ","))
  if (isTRUE(attributes(x)$wisca)) {
    names(dims) <- paste0("An Antibiogram (WISCA / ", attributes(x)$conf_interval * 100, "% CI)")
  } else {
    names(dims) <- "An Antibiogram (non-WISCA)"
  }
  dims
}

# will be exported in R/zzz.R
tbl_format_footer.antibiogram <- function(x, ...) {
  footer <- NextMethod()
  if (NROW(x) == 0) {
    return(footer)
  }
  c(footer, font_subtle(paste0("# Use `plot()` or `ggplot2::autoplot()` to create a plot of this antibiogram,\n",
                               "# or use it directly in R Markdown or ",
                               font_url("https://quarto.org", "Quarto"), ", see ", word_wrap("?antibiogram"))))
}

#' @export
#' @rdname antibiogram
plot.antibiogram <- function(x, ...) {
  df <- attributes(x)$long_numeric
  if (!"mo" %in% colnames(df)) {
    stop_("Plotting antibiograms using `plot()` is only possible if they were not created using dplyr groups. See `?antibiogram` for how to retrieve numeric values in a long format for advanced plotting.",
          call = FALSE)
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
    df_sub <- df[df$mo == mo, , drop = FALSE]
    
    bp <- barplot(
      height = df_sub$percentage * 100,
      xlab = NULL,
      ylab = ifelse(isTRUE(attributes(x)$combine_SI), "%SI", "%S"),
      names.arg = df_sub$ab,
      col = "#aaaaaa",
      beside = TRUE,
      main = mo,
      legend = NULL
    )
    
    if (isTRUE(attributes(x)$wisca)) {
      lower <- df_sub$lower * 100
      upper <- df_sub$upper * 100
      arrows(
        x0 = bp, y0 = lower,  # Start of error bar (lower bound)
        x1 = bp, y1 = upper,  # End of error bar (upper bound)
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
# will be exported using s3_register() in R/zzz.R
autoplot.antibiogram <- function(object, ...) {
  df <- attributes(object)$long_numeric
  if (!"mo" %in% colnames(df)) {
    stop_("Plotting antibiograms using `autoplot()` is only possible if they were not created using dplyr groups. See `?antibiogram` for how to retrieve numeric values in a long format for advanced plotting.",
          call = FALSE)
  }
  out <- ggplot2::ggplot(df,
                         mapping = ggplot2::aes(
                           x = ab,
                           y = percentage * 100,
                           fill = if ("syndromic_group" %in% colnames(df)) {
                             syndromic_group
                           } else {
                             NULL
                           }
                         )) +
    ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::facet_wrap("mo") +
    ggplot2::labs(
      y = ifelse(isTRUE(attributes(object)$combine_SI), "%SI", "%S"),
      x = NULL,
      fill = if ("syndromic_group" %in% colnames(df)) {
        colnames(object)[1]
      } else {
        NULL
      }
    )
  if (isTRUE(attributes(object)$wisca)) {
    out <- out + 
      ggplot2::geom_errorbar(mapping = ggplot2::aes(ymin = lower * 100, ymax = upper * 100),
                             position = ggplot2::position_dodge2(preserve = "single"),
                             width = 0.5)
  }
  out
}

# will be exported in zzz.R
#' @method knit_print antibiogram
#' @param italicise a [logical] to indicate whether the microorganism names in the [knitr][knitr::kable()] table should be made italic, using [italicise_taxonomy()].
#' @param na character to use for showing `NA` values
#' @rdname antibiogram
knit_print.antibiogram <- function(x, italicise = TRUE, na = getOption("knitr.kable.NA", default = ""), ...) {
  stop_ifnot_installed("knitr")
  meet_criteria(italicise, allow_class = "logical", has_length = 1)
  meet_criteria(na, allow_class = "character", has_length = 1, allow_NA = TRUE)
  
  if (isTRUE(italicise) && "mo" %in% colnames(attributes(x)$long_numeric)) {
    # make all microorganism names italic, according to nomenclature
    names_col <- ifelse(isTRUE(attributes(x)$has_syndromic_group), 2, 1)
    x[[names_col]] <- italicise_taxonomy(x[[names_col]], type = "markdown")
  }
  
  old_option <- getOption("knitr.kable.NA")
  options(knitr.kable.NA = na)
  on.exit(options(knitr.kable.NA = old_option))
  
  out <- paste(c("", "", knitr::kable(x, ..., output = FALSE)), collapse = "\n")
  knitr::asis_output(out)
}
