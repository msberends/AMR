# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# doi:10.18637/jss.v104.i03                                            #
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

#' Generate Antibiogram: Traditional, Combined, Syndromic, or Weighted-Incidence Syndromic Combination (WISCA)
#'
#' Generate an antibiogram, and communicate the results in plots or tables. These functions follow the logic of Klinker *et al.* and Barbieri *et al.* (see *Source*), and allow reporting in e.g. R Markdown and Quarto as well.
#' @param x a [data.frame] containing at least a column with microorganisms and columns with antibiotic results (class 'sir', see [as.sir()])
#' @param antibiotics vector of column names, or (any combinations of) [antibiotic selectors][antibiotic_class_selectors] such as [aminoglycosides()] or [carbapenems()]. For combination antibiograms, this can also be column names separated with `"+"`, such as "TZP+TOB" given that the data set contains columns "TZP" and "TOB". See *Examples*.
#' @param mo_transform a character to transform microorganism input - must be "name", "shortname", "gramstain", or one of the column names of the [microorganisms] data set: `r vector_or(colnames(microorganisms), sort = FALSE, quotes = TRUE)`. Can also be `NULL` to not transform the input.
#' @param ab_transform a character to transform antibiotic input - must be one of the column names of the [antibiotics] data set: `r vector_or(colnames(antibiotics), sort = FALSE, quotes = TRUE)`. Can also be `NULL` to not transform the input.
#' @param syndromic_group a column name of `x`, or values calculated to split rows of `x`, e.g. by using [ifelse()] or [`case_when()`][dplyr::case_when()]. See *Examples*.
#' @param add_total_n a [logical] to indicate whether total available numbers per pathogen should be added to the table (defaults to `TRUE`). This will add the lowest and highest number of available isolate per antibiotic (e.g, if for *E. coli* 200 isolates are available for ciprofloxacin and 150 for amoxicillin, the returned number will be "150-200").
#' @param only_all_tested (for combination antibiograms): a [logical] to indicate that isolates must be tested for all antibiotics, see *Details*
#' @param digits number of digits to use for rounding
#' @param col_mo column name of the names or codes of the microorganisms (see [as.mo()]), defaults to the first column of class [`mo`]. Values will be coerced using [as.mo()].
#' @param language language to translate text, which defaults to the system language (see [get_AMR_locale()])
#' @param minimum the minimum allowed number of available (tested) isolates. Any isolate count lower than `minimum` will return `NA` with a warning. The default number of `30` isolates is advised by the Clinical and Laboratory Standards Institute (CLSI) as best practice, see *Source*.
#' @param combine_SI a [logical] to indicate whether all susceptibility should be determined by results of either S or I, instead of only S (defaults to `TRUE`)
#' @param sep a separating character for antibiotic columns in combination antibiograms
#' @param object an [antibiogram()] object
#' @param ... method extensions
#' @details This function returns a table with values between 0 and 100 for *susceptibility*, not resistance.
#'
#' **Remember that you should filter your data to let it contain only first isolates!** This is needed to exclude duplicates and to reduce selection bias. Use [first_isolate()] to determine them in your data set with one of the four available algorithms.
#'
#' There are four antibiogram types, as proposed by Klinker *et al.* (2021, \doi{10.1177/20499361211011373}), and they are all supported by [antibiogram()]:
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
#'    Case example: Susceptibility of *Pseudomonas aeruginosa* to TZP among respiratory specimens (obtained among ICU patients only) for male patients age >=65 years with heart failure
#'
#'    Code example:
#'
#'    ```r
#'    library(dplyr)
#'    your_data %>% 
#'      filter(ward == "ICU" & specimen_type == "Respiratory") %>% 
#'      antibiogram(antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"),
#'                  syndromic_group = ifelse(.$age >= 65 &
#'                                             .$gender == "Male" &
#'                                             .$condition == "Heart Disease",
#'                                           "Study Group", "Control Group"))
#'    ```
#'
#' All types of antibiograms can be generated with the functions as described on this page, and can be plotted (using [ggplot2::autoplot()] or base \R [plot()]/[barplot()]) or printed into R Markdown / Quarto formats for reports using `print()`. Use functions from specific 'table reporting' packages to transform the output of [antibiogram()] to your needs, e.g. `flextable::as_flextable()` or `gt::gt()`.
#'
#' Note that for combination antibiograms, it is important to realise that susceptibility can be calculated in two ways, which can be set with the `only_all_tested` argument (defaults to `FALSE`). See this example for two antibiotics, Drug A and Drug B, about how [antibiogram()] works to calculate the %SI:
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
#' @source
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
#' antibiogram(example_isolates,
#'   antibiotics = c("TZP", "TZP+TOB"),
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
#' # Weighted-incidence syndromic combination antibiogram (WISCA) ---------
#'
#' # the data set could contain a filter for e.g. respiratory specimens/ICU
#' antibiogram(example_isolates,
#'   antibiotics = c("AMC", "AMC+CIP", "TZP", "TZP+TOB"),
#'   mo_transform = "gramstain",
#'   minimum = 10, # this should be >=30, but now just as example
#'   syndromic_group = ifelse(example_isolates$age >= 65 &
#'     example_isolates$gender == "M",
#'   "WISCA Group 1", "WISCA Group 2"
#'   )
#' )
#'
#'
#' # Generate plots with ggplot2 or base R --------------------------------
#'
#' ab1 <- antibiogram(example_isolates,
#'   antibiotics = c("AMC", "CIP", "TZP", "TZP+TOB"),
#'   mo_transform = "gramstain"
#' )
#' ab2 <- antibiogram(example_isolates,
#'   antibiotics = c("AMC", "CIP", "TZP", "TZP+TOB"),
#'   mo_transform = "gramstain",
#'   syndromic_group = "ward"
#' )
#'
#' plot(ab1)
#'
#' if (requireNamespace("ggplot2")) {
#'   ggplot2::autoplot(ab1)
#' }
#'
#' plot(ab2)
#'
#' if (requireNamespace("ggplot2")) {
#'   ggplot2::autoplot(ab2)
#' }
#' }
antibiogram <- function(x,
                        antibiotics = where(is.sir),
                        mo_transform = "shortname",
                        ab_transform = NULL,
                        syndromic_group = NULL,
                        add_total_n = TRUE,
                        only_all_tested = FALSE,
                        digits = 0,
                        col_mo = NULL,
                        language = get_AMR_locale(),
                        minimum = 30,
                        combine_SI = TRUE,
                        sep = " + ") {
  meet_criteria(x, allow_class = "data.frame", contains_column_class = "sir")
  meet_criteria(mo_transform, allow_class = "character", has_length = 1, is_in = c("name", "shortname", "gramstain", colnames(AMR::microorganisms)), allow_NULL = TRUE)
  meet_criteria(ab_transform, allow_class = "character", has_length = 1, is_in = colnames(AMR::antibiotics), allow_NULL = TRUE)
  meet_criteria(syndromic_group, allow_class = "character", allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(add_total_n, allow_class = "logical", has_length = 1)
  meet_criteria(only_all_tested, allow_class = "logical", has_length = 1)
  meet_criteria(digits, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE)
  meet_criteria(col_mo, allow_class = "character", has_length = 1, allow_NULL = TRUE, is_in = colnames(x))
  language <- validate_language(language)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(sep, allow_class = "character", has_length = 1)

  # try to find columns based on type
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = interactive())
    stop_if(is.null(col_mo), "`col_mo` must be set")
  }
  # transform MOs
  x$`.mo` <- x[, col_mo, drop = TRUE]
  if (is.null(mo_transform)) {
    # leave as is
  } else if (mo_transform == "gramstain") {
    x$`.mo` <- mo_gramstain(x$`.mo`, language = language)
  } else if (mo_transform == "shortname") {
    x$`.mo` <- mo_shortname(x$`.mo`, language = language)
  } else if (mo_transform == "name") {
    x$`.mo` <- mo_name(x$`.mo`, language = language)
  } else {
    x$`.mo` <- mo_property(x$`.mo`, language = language)
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
  if (tryCatch(is.character(antibiotics), error = function(e) FALSE)) {
    antibiotics <- strsplit(gsub(" ", "", antibiotics), "+", fixed = TRUE)
    non_existing <- unlist(antibiotics)[!unlist(antibiotics) %in% colnames(x)]
    if (length(non_existing) > 0) {
      warning_("The following antibiotics were not available and ignored: ", vector_and(non_existing, sort = FALSE))
      antibiotics <- lapply(antibiotics, function(ab) ab[!ab %in% non_existing])
    }
    # make list unique
    antibiotics <- unique(antibiotics)
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
          S_values <- c("S", "I")
        } else {
          S_values <- "S"
        }
        other_values <- setdiff(c("S", "I", "R"), S_values)
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
    antibiotics <- colnames(suppressWarnings(x[, antibiotics, drop = FALSE]))
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
      FUN = function(x) x
    )
  counts <- out
  
  if (isTRUE(combine_SI)) {
    out$numerator <- out$S + out$I
  } else {
    out$numerator <- out$S
  }
  if (any(out$total < minimum, na.rm = TRUE)) {
    message_("NOTE: ", sum(out$total < minimum, na.rm = TRUE), " combinations had less than `minimum = ", minimum, "` results and were ignored", add_fn = font_red)
    out <- out %pm>%
      subset(total >= minimum)
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
  out <- out %pm>%
    pm_summarise(SI = numerator / total)

  # transform names of antibiotics
  ab_naming_function <- function(x, t, l, s) {
    x <- strsplit(x, s, fixed = TRUE)
    out <- character(length = length(x))
    for (i in seq_len(length(x))) {
      a <- x[[i]]
      if (is.null(t)) {
        # leave as is
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

  # transform long to wide
  long_to_wide <- function(object, digs) {
    object$SI <- round(object$SI * 100, digits = digs)
    object <- object %pm>%
      # an unclassed data.frame is required for stats::reshape()
      as.data.frame(stringsAsFactors = FALSE) %pm>%
      stats::reshape(direction = "wide", idvar = "mo", timevar = "ab", v.names = "SI")
    colnames(object) <- gsub("^SI?[.]", "", colnames(object))
    return(object)
  }

  # ungroup for long -> wide transformation
  attr(out, "pm_groups") <- NULL
  attr(out, "groups") <- NULL
  class(out) <- class(out)[!class(out) %in% c("grouped_df", "grouped_data")]
  long <- out

  if (isTRUE(has_syndromic_group)) {
    grps <- unique(out$syndromic_group)
    for (i in seq_len(length(grps))) {
      grp <- grps[i]
      if (i == 1) {
        new_df <- long_to_wide(out[which(out$syndromic_group == grp), , drop = FALSE], digs = digits)
      } else {
        new_df <- rbind2(
          new_df,
          long_to_wide(out[which(out$syndromic_group == grp), , drop = FALSE], digs = digits)
        )
      }
    }
    # sort rows
    new_df <- new_df %pm>% pm_arrange(mo, syndromic_group)
    # sort columns
    new_df <- new_df[, c("syndromic_group", "mo", sort(colnames(new_df)[!colnames(new_df) %in% c("syndromic_group", "mo")])), drop = FALSE]
    colnames(new_df)[1:2] <- translate_AMR(c("Syndromic Group", "Pathogen"), language = language)
  } else {
    new_df <- long_to_wide(out, digs = digits)
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
    new_df[[edit_col]] <- paste0(new_df[[edit_col]], " (", count_group, ")")
    colnames(new_df)[edit_col] <- paste(colnames(new_df)[edit_col], "(N min-max)")
  }

  structure(as_original_data_class(new_df, class(x), extra_class = "antibiogram"),
    long = long,
    combine_SI = combine_SI
  )
}

#' @export
#' @rdname antibiogram
plot.antibiogram <- function(x, ...) {
  df <- attributes(x)$long
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

    barplot(
      height = df_sub$SI * 100,
      xlab = NULL,
      ylab = ifelse(isTRUE(attributes(x)$combine_SI), "%SI", "%S"),
      names.arg = df_sub$ab,
      col = "#aaaaaa",
      beside = TRUE,
      main = mo,
      legend = NULL
    )
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
  df <- attributes(object)$long
  ggplot2::ggplot(df) +
    ggplot2::geom_col(
      ggplot2::aes(
        x = ab,
        y = SI * 100,
        fill = if ("syndromic_group" %in% colnames(df)) {
          syndromic_group
        } else {
          NULL
        }
      ),
      position = ggplot2::position_dodge2(preserve = "single")
    ) +
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
}

#' @export
#' @param as_kable a [logical] to indicate whether the printing should be done using [knitr::kable()] (which is the default in non-interactive sessions)
#' @details Printing the antibiogram in non-interactive sessions will be done by [knitr::kable()], with support for [all their implemented formats][knitr::kable()], such as "markdown". The knitr format will be automatically determined if printed inside a knitr document (LaTeX, HTML, etc.).
#' @rdname antibiogram
print.antibiogram <- function(x, as_kable = !interactive(), ...) {
  meet_criteria(as_kable, allow_class = "logical", has_length = 1)
  
  kable <- import_fn("kable", "knitr", error_on_fail = FALSE)
  if (!is.null(kable) &&
      isTRUE(as_kable) &&
      # be sure not to run kable in pkgdown for our website generation
      !identical(Sys.getenv("IN_PKGDOWN"), "true")) {
    kable(x, ...)
  
  } else {
    # remove 'antibiogram' class and print with default method
    class(x) <- class(x)[class(x) != "antibiogram"]
    print(x, ...)
  }
}
