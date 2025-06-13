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

#' Random MIC Values/Disk Zones/SIR Generation
#'
#' These functions can be used for generating random MIC values and disk diffusion diameters, for AMR data analysis practice. By providing a microorganism and antimicrobial drug, the generated results will reflect reality as much as possible.
#' @param size Desired size of the returned vector. If used in a [data.frame] call or `dplyr` verb, will get the current (group) size if left blank.
#' @param mo Any [character] that can be coerced to a valid microorganism code with [as.mo()]. Can be the same length as `size`.
#' @param ab Any [character] that can be coerced to a valid antimicrobial drug code with [as.ab()].
#' @param prob_SIR A vector of length 3: the probabilities for "S" (1st value), "I" (2nd value) and "R" (3rd value).
#' @param skew Direction of skew for MIC or disk values, either `"right"` or `"left"`. A left-skewed distribution has the majority of the data on the right.
#' @param severity Skew severity; higher values will increase the skewedness. Default is `2`; use `0` to prevent skewedness.
#' @param ... Ignored, only in place to allow future extensions.
#' @details
#' Internally, MIC and disk zone values are sampled based on clinical breakpoints defined in the [clinical_breakpoints] data set. To create specific generated values per bug or drug, set the `mo` and/or `ab` argument. The MICs are sampled on a log2 scale and disks linearly, using weighted probabilities. The weights are based on the `skew` and `severity` arguments:
#' * `skew = "right"` places more emphasis on lower MIC or higher disk values.
#' * `skew = "left"` places more emphasis on higher MIC or lower disk values.
#' * `severity` controls the exponential bias applied.
#' @return class `mic` for [random_mic()] (see [as.mic()]) and class `disk` for [random_disk()] (see [as.disk()])
#' @name random
#' @rdname random
#' @export
#' @examples
#' random_mic(25)
#' random_disk(25)
#' random_sir(25)
#'
#' # add more skewedness, make more realistic by setting a bug and/or drug:
#' disks <- random_disk(100, severity = 2, mo = "Escherichia coli", ab = "CIP")
#' plot(disks)
#' # `plot()` and `ggplot2::autoplot()` allow for coloured bars if `mo` and `ab` are set
#' plot(disks, mo = "Escherichia coli", ab = "CIP", guideline = "CLSI 2025")
#'
#' \donttest{
#' random_mic(25, "Klebsiella pneumoniae") # range 0.0625-64
#' random_mic(25, "Klebsiella pneumoniae", "meropenem") # range 0.0625-16
#' random_mic(25, "Streptococcus pneumoniae", "meropenem") # range 0.0625-4
#'
#' random_disk(25, "Klebsiella pneumoniae") # range 8-50
#' random_disk(25, "Klebsiella pneumoniae", "ampicillin") # range 11-17
#' random_disk(25, "Streptococcus pneumoniae", "ampicillin") # range 12-27
#' }
random_mic <- function(size = NULL, mo = NULL, ab = NULL, skew = "right", severity = 1, ...) {
  meet_criteria(size, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = TRUE, allow_NULL = TRUE)
  meet_criteria(mo, allow_class = "character", has_length = c(1, size), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(skew, allow_class = "character", is_in = c("right", "left"), has_length = 1)
  meet_criteria(severity, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)

  if (is.null(size)) {
    size <- NROW(get_current_data(arg_name = "size", call = -3))
  }
  if (length(mo) > 1) {
    out <- rep(NA_mic_, length(size))
    p <- progress_ticker(n = length(unique(mo)), n_min = 10, title = "Generating random MIC values")
    for (mo_ in unique(mo)) {
      p$tick()
      out[which(mo == mo_)] <- random_exec("MIC", size = sum(mo == mo_), mo = mo_, ab = ab, skew = skew, severity = severity)
    }
    out <- as.mic(out, keep_operators = "none")
    if (stats::runif(1) > 0.5 && length(unique(out)) > 1) {
      out[out == min(out)] <- paste0("<=", out[out == min(out)])
    }
    if (stats::runif(1) > 0.5 && length(unique(out)) > 1) {
      out[out == max(out)] <- paste0(">=", out[out == max(out)])
    }
    return(out)
  } else {
    random_exec("MIC", size = size, mo = mo, ab = ab, skew = skew, severity = severity)
  }
}

#' @rdname random
#' @export
random_disk <- function(size = NULL, mo = NULL, ab = NULL, skew = "left", severity = 1, ...) {
  meet_criteria(size, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = TRUE, allow_NULL = TRUE)
  meet_criteria(mo, allow_class = "character", has_length = c(1, size), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(skew, allow_class = "character", is_in = c("right", "left"), has_length = 1)
  meet_criteria(severity, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)

  if (is.null(size)) {
    size <- NROW(get_current_data(arg_name = "size", call = -3))
  }
  if (length(mo) > 1) {
    out <- rep(NA_mic_, length(size))
    p <- progress_ticker(n = length(unique(mo)), n_min = 10, title = "Generating random MIC values")
    for (mo_ in unique(mo)) {
      p$tick()
      out[which(mo == mo_)] <- random_exec("DISK", size = sum(mo == mo_), mo = mo_, ab = ab, skew = skew, severity = severity)
    }
    out <- as.disk(out)
    return(out)
  } else {
    random_exec("DISK", size = size, mo = mo, ab = ab, skew = skew, severity = severity)
  }
}

#' @rdname random
#' @export
random_sir <- function(size = NULL, prob_SIR = c(0.33, 0.33, 0.33), ...) {
  meet_criteria(size, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = TRUE, allow_NULL = TRUE)
  meet_criteria(prob_SIR, allow_class = c("numeric", "integer"), has_length = 3)
  if (is.null(size)) {
    size <- NROW(get_current_data(arg_name = "size", call = -3))
  }
  sample(as.sir(c("S", "I", "R")), size = size, replace = TRUE, prob = prob_SIR)
}


random_exec <- function(method_type, size, mo = NULL, ab = NULL, skew = "right", severity = 1) {
  df <- AMR::clinical_breakpoints %pm>% subset(method == method_type & type == "human")

  if (!is.null(mo)) {
    mo_coerced <- as.mo(mo, info = FALSE)
    mo_include <- c(mo_coerced, as.mo(mo_genus(mo_coerced)), as.mo(mo_family(mo_coerced)), as.mo(mo_order(mo_coerced)))
    df_new <- df %pm>% subset(mo %in% mo_include)
    if (nrow(df_new) > 0) df <- df_new
  }

  if (!is.null(ab)) {
    ab_coerced <- as.ab(ab)
    df_new <- df %pm>% subset(ab %in% ab_coerced)
    if (nrow(df_new) > 0) df <- df_new
  }

  if (method_type == "MIC") {
    lowest_mic <- min(df$breakpoint_S, na.rm = TRUE)
    lowest_mic <- log2(lowest_mic) + sample(c(-3:2), 1)
    lowest_mic <- 2^lowest_mic
    highest_mic <- max(df$breakpoint_R, na.rm = TRUE)
    highest_mic <- log2(highest_mic) + sample(c(-3:1), 1)
    highest_mic <- max(lowest_mic * 2, 2^highest_mic)

    out <- skewed_values(COMMON_MIC_VALUES, size = size, min = lowest_mic, max = highest_mic, skew = skew, severity = severity)
    if (stats::runif(1) > 0.5 && length(unique(out)) > 1) {
      out[out == min(out)] <- paste0("<=", out[out == min(out)])
    }
    if (stats::runif(1) > 0.5 && length(unique(out)) > 1) {
      out[out == max(out)] <- paste0(">=", out[out == max(out)])
    }
    return(as.mic(out))
  } else if (method_type == "DISK") {
    disk_range <- seq(
      from = floor(min(df$breakpoint_R[!is.na(df$breakpoint_R)], na.rm = TRUE) / 1.25),
      to = ceiling(max(df$breakpoint_S[df$breakpoint_S != 50], na.rm = TRUE) * 1.25),
      by = 1
    )
    disk_range <- disk_range[disk_range >= 6 & disk_range <= 50]
    out <- skewed_values(disk_range, size = size, min = min(disk_range), max = max(disk_range), skew = skew, severity = severity)
    return(as.disk(out))
  }
}

skewed_values <- function(values, size, min, max, skew = c("right", "left"), severity = 1) {
  skew <- match.arg(skew)
  range_vals <- values[values >= min & values <= max]
  if (length(range_vals) < 2) range_vals <- values
  ranks <- seq_along(range_vals)
  weights <- switch(skew,
    right = rev(ranks)^severity,
    left = ranks^severity
  )
  weights <- weights / sum(weights)
  sample(range_vals, size = size, replace = TRUE, prob = weights)
}
