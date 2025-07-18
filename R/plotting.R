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

#' Plotting Helpers for AMR Data Analysis
#'
#' @description
#' Functions to plot classes `sir`, `mic` and `disk`, with support for base \R and `ggplot2`.
#'
#' Especially the `scale_*_mic()` functions are relevant wrappers to plot MIC values for `ggplot2`. They allows custom MIC ranges and to plot intermediate log2 levels for missing MIC values.
#' @param x,object Values created with [as.mic()], [as.disk()] or [as.sir()] (or their `random_*` variants, such as [random_mic()]).
#' @param mo Any (vector of) text that can be coerced to a valid microorganism code with [as.mo()].
#' @param ab Any (vector of) text that can be coerced to a valid antimicrobial drug code with [as.ab()].
#' @param guideline Interpretation guideline to use - the default is the latest included EUCAST guideline, see *Details*.
#' @param main,title Title of the plot.
#' @param xlab,ylab Axis title.
#' @param colours_SIR Colours to use for filling in the bars, must be a vector of three values (in the order S, I and R). The default colours are colour-blind friendly.
#' @param language Language to be used to translate 'Susceptible', 'Increased exposure'/'Intermediate' and 'Resistant' - the default is system language (see [get_AMR_locale()]) and can be overwritten by setting the package option [`AMR_locale`][AMR-options], e.g. `options(AMR_locale = "de")`, see [translate]. Use `language = NULL` to prevent translation.
#' @param expand A [logical] to indicate whether the range on the x axis should be expanded between the lowest and highest value. For MIC values, intermediate values will be factors of 2 starting from the highest MIC value. For disk diameters, the whole diameter range will be filled.
#' @param aesthetics Aesthetics to apply the colours to - the default is "fill" but can also be (a combination of) "alpha", "colour", "fill", "linetype", "shape" or "size".
#' @param eucast_I A [logical] to indicate whether the 'I' must be interpreted as "Susceptible, under increased exposure". Will be `TRUE` if the default [AMR interpretation guideline][as.sir()] is set to EUCAST (which is the default). With `FALSE`, it will be interpreted as "Intermediate".
#' @inheritParams as.sir
#' @param mic_range A manual range to rescale the MIC values (using [rescale_mic()]), e.g., `mic_range = c(0.001, 32)`. Use `NA` to prevent rescaling on one side, e.g., `mic_range = c(NA, 32)`. **Note:** This rescales values but does not filter them - use the ggplot2 `limits` argument separately to exclude values from the plot.
#' @inheritParams as.mic
#' @inheritParams ggplot_sir
#' @inheritParams proportion
#' @details
#' ### The `scale_*_mic()` Functions
#'
#' The functions [scale_x_mic()], [scale_y_mic()], [scale_colour_mic()], and [scale_fill_mic()] functions allow to plot the [mic][as.mic()] class (MIC values) on a continuous, logarithmic scale. They also allow to rescale the MIC range with an 'inside' or 'outside' range if required, and retain the operators in MIC values (such as `>=`) if desired. Missing intermediate log2 levels will be plotted too.
#'
#' ### The `scale_*_sir()` Functions
#'
#' The functions [scale_x_sir()], [scale_colour_sir()], and [scale_fill_sir()] functions allow to plot the [sir][as.sir()] class in the right order (`r paste(levels(NA_sir_), collapse = " < ")`). At default, they translate the S/I/R values to an interpretative text ("Susceptible", "Resistant", etc.) in any of the `r length(AMR:::LANGUAGES_SUPPORTED)` supported languages (use `language = NULL` to keep S/I/R). Also, except for [scale_x_sir()], they set colour-blind friendly colours to the `colour` and `fill` aesthetics.
#'
#' ### Additional `ggplot2` Functions
#'
#' This package contains more functions that extend the `ggplot2` package, to help in visualising AMR data results. All these functions are internally used by [ggplot_sir()] too.
#'
#' * [facet_sir()] creates 2d plots (at default based on S/I/R) using [ggplot2::facet_wrap()].
#' * [scale_y_percent()] transforms the y axis to a 0 to 100% range using [ggplot2::scale_y_continuous()].
#' * [scale_sir_colours()] allows to set colours to any aesthetic, even for `shape` or `linetype`.
#' * [theme_sir()] is a [ggplot2 theme][ggplot2::theme()] with minimal distraction.
#' * [labels_sir_count()] print datalabels on the bars with percentage and number of isolates, using [ggplot2::geom_text()].
#'
#' The interpretation of "I" will be named "Increased exposure" for all EUCAST guidelines since 2019, and will be named "Intermediate" in all other cases.
#'
#' For interpreting MIC values as well as disk diffusion diameters, the default guideline is `r AMR::clinical_breakpoints$guideline[1]`, unless the package option [`AMR_guideline`][AMR-options] is set. See [as.sir()] for more information.
#' @name plot
#' @rdname plot
#' @return The `autoplot()` functions return a [`ggplot`][ggplot2::ggplot()] model that is extendible with any `ggplot2` function.
#' @param ... Arguments passed on to methods.
#' @examples
#' some_mic_values <- random_mic(size = 100)
#' some_disk_values <- random_disk(size = 100, mo = "Escherichia coli", ab = "cipro")
#' some_sir_values <- random_sir(50, prob_SIR = c(0.55, 0.05, 0.30))
#'
#' \donttest{
#' # Plotting using ggplot2's autoplot() for MIC, disk, and SIR -----------
#' if (require("ggplot2")) {
#'   autoplot(some_mic_values)
#' }
#' if (require("ggplot2")) {
#'   # when providing the microorganism and antibiotic, colours will show interpretations:
#'   autoplot(some_mic_values, mo = "Escherichia coli", ab = "cipro")
#' }
#' if (require("ggplot2")) {
#'   autoplot(some_mic_values, mo = "Staph aureus", ab = "Ceftaroline", guideline = "CLSI")
#' }
#'
#' if (require("ggplot2")) {
#'   # support for 27 languages, various guidelines, and many options
#'   autoplot(some_disk_values,
#'     mo = "Escherichia coli", ab = "cipro",
#'     guideline = "CLSI 2024", language = "no",
#'     title = "Disk diffusion from the North"
#'   )
#' }
#'
#'
#' # Plotting using scale_x_mic() -----------------------------------------
#' if (require("ggplot2")) {
#'   mic_plot <- ggplot(
#'     data.frame(
#'       mics = as.mic(c(0.25, "<=4", 4, 8, 32, ">=32")),
#'       counts = c(1, 1, 2, 2, 3, 3)
#'     ),
#'     aes(mics, counts)
#'   ) +
#'     geom_col()
#'   mic_plot +
#'     labs(title = "without scale_x_mic()")
#' }
#' if (require("ggplot2")) {
#'   mic_plot +
#'     scale_x_mic() +
#'     labs(title = "with scale_x_mic()")
#' }
#' if (require("ggplot2")) {
#'   mic_plot +
#'     scale_x_mic(keep_operators = "all") +
#'     labs(title = "with scale_x_mic() keeping all operators")
#' }
#' if (require("ggplot2")) {
#'   mic_plot +
#'     scale_x_mic(mic_range = c(1, 16)) +
#'     labs(title = "with scale_x_mic() using a manual 'within' range")
#' }
#' if (require("ggplot2")) {
#'   mic_plot +
#'     scale_x_mic(mic_range = c(0.032, 256)) +
#'     labs(title = "with scale_x_mic() using a manual 'outside' range")
#' }
#'
#'
#' # Plotting using scale_y_mic() -----------------------------------------
#' some_groups <- sample(LETTERS[1:5], 20, replace = TRUE)
#'
#' if (require("ggplot2")) {
#'   ggplot(
#'     data.frame(
#'       mic = some_mic_values,
#'       group = some_groups
#'     ),
#'     aes(group, mic)
#'   ) +
#'     geom_boxplot() +
#'     geom_violin(linetype = 2, colour = "grey30", fill = NA) +
#'     scale_y_mic()
#' }
#' if (require("ggplot2")) {
#'   ggplot(
#'     data.frame(
#'       mic = some_mic_values,
#'       group = some_groups
#'     ),
#'     aes(group, mic)
#'   ) +
#'     geom_boxplot() +
#'     geom_violin(linetype = 2, colour = "grey30", fill = NA) +
#'     scale_y_mic(mic_range = c(NA, 0.25))
#' }
#'
#'
#' # Plotting using scale_x_sir() -----------------------------------------
#' if (require("ggplot2")) {
#'   ggplot(
#'     data.frame(
#'       x = c("I", "R", "S"),
#'       y = c(45, 323, 573)
#'     ),
#'     aes(x, y)
#'   ) +
#'     geom_col() +
#'     scale_x_sir()
#' }
#'
#'
#' # Plotting using scale_y_mic() and scale_colour_sir() ------------------
#' if (require("ggplot2")) {
#'   plain <- ggplot(
#'     data.frame(
#'       mic = some_mic_values,
#'       group = some_groups,
#'       sir = as.sir(some_mic_values,
#'         mo = "E. coli",
#'         ab = "cipro"
#'       )
#'     ),
#'     aes(x = group, y = mic, colour = sir)
#'   ) +
#'     theme_minimal() +
#'     geom_boxplot(fill = NA, colour = "grey30") +
#'     geom_jitter(width = 0.25)
#'
#'   plain
#' }
#' if (require("ggplot2")) {
#'   # and now with our MIC and SIR scale functions:
#'   plain +
#'     scale_y_mic() +
#'     scale_colour_sir()
#' }
#' if (require("ggplot2")) {
#'   plain +
#'     scale_y_mic(mic_range = c(0.005, 32), name = "Our MICs!") +
#'     scale_colour_sir(
#'       language = "pt",
#'       name = "Support in 27 languages"
#'     )
#' }
#' }
#'
#' # Plotting using base R's plot() ---------------------------------------
#'
#' plot(some_mic_values)
#' # when providing the microorganism and antibiotic, colours will show interpretations:
#' plot(some_mic_values, mo = "S. aureus", ab = "ampicillin")
#'
#' plot(some_disk_values)
#' plot(some_disk_values, mo = "Escherichia coli", ab = "cipro")
#' plot(some_disk_values, mo = "Escherichia coli", ab = "cipro", language = "nl")
#'
#' plot(some_sir_values)
NULL

create_scale_mic <- function(aest, keep_operators, mic_range = NULL, ...) {
  ggplot_fn <- getExportedValue(paste0("scale_", aest, "_continuous"),
    ns = asNamespace("ggplot2")
  )
  args <- list(...)
  breaks_set <- args$breaks
  limits_set <- args$limits

  # do not take these arguments into account, as they will be overwritten and seem to allow weird behaviour if set anyway
  args[c("aesthetics", "trans", "transform", "transform_df", "breaks", "labels", "limits")] <- NULL
  scale <- do.call(ggplot_fn, args)
  scale$mic_breaks_set <- breaks_set
  scale$mic_limits_set <- limits_set

  scale$transform <- function(x) {
    as.double(rescale_mic(x = as.double(as.mic(x)), keep_operators = keep_operators, mic_range = mic_range, as.mic = TRUE))
  }
  scale$transform_df <- function(self, df) {
    if (!aest %in% colnames(df)) {
      # support for geom_hline(), geom_vline(), etc
      other_x <- c("xintercept", "xmin", "xmax", "xend", "width")
      other_y <- c("yintercept", "ymin", "ymax", "yend", "height")
      if (any(other_y %in% colnames(df))) {
        aest_val <- intersect(other_y, colnames(df))[1]
      } else if (any(other_x %in% colnames(df))) {
        aest_val <- intersect(other_x, colnames(df))[1]
      } else {
        stop_("No support for plotting df with `scale_", aest, "_mic()` with columns ", vector_and(colnames(df), sort = FALSE))
      }
      out <- rescale_mic(x = as.double(as.mic(df[[aest_val]])), keep_operators = "none", mic_range = NULL, as.mic = TRUE)
      if (!is.null(self$mic_values_rescaled) && any(out < min(self$mic_values_rescaled, na.rm = TRUE) | out > max(self$mic_values_rescaled, na.rm = TRUE), na.rm = TRUE)) {
        warning_("The value for `", aest_val, "` is outside the plotted MIC range, consider using/updating the `mic_range` argument in `scale_", aest, "_mic()`.")
      }
      df[[aest_val]] <- log2(as.double(out))
    } else {
      self$mic_values_rescaled <- rescale_mic(x = as.double(as.mic(df[[aest]])), keep_operators = keep_operators, mic_range = mic_range, as.mic = TRUE)
      # create new breaks and labels here
      lims <- range(self$mic_values_rescaled, na.rm = TRUE)
      # support inner and outer 'mic_range' settings (e.g., the data ranges 0.5-8 and 'mic_range' is set to 0.025-32)
      if (!is.null(mic_range) && !is.na(mic_range[1]) && !is.na(lims[1]) && mic_range[1] < lims[1]) {
        lims[1] <- mic_range[1]
      }
      if (!is.null(mic_range) && !is.na(mic_range[2]) && !is.na(lims[2]) && mic_range[2] > lims[2]) {
        lims[2] <- mic_range[2]
      }
      ind_min <- which(COMMON_MIC_VALUES <= lims[1])[which.min(abs(COMMON_MIC_VALUES[COMMON_MIC_VALUES <= lims[1]] - lims[1]))] # Closest index where COMMON_MIC_VALUES <= lims[1]
      ind_max <- which(COMMON_MIC_VALUES >= lims[2])[which.min(abs(COMMON_MIC_VALUES[COMMON_MIC_VALUES >= lims[2]] - lims[2]))] # Closest index where COMMON_MIC_VALUES >= lims[2]

      self$mic_values_levels <- as.mic(COMMON_MIC_VALUES[ind_min:ind_max])

      if (keep_operators %in% c("edges", "all") && length(unique(self$mic_values_levels)) > 1) {
        self$mic_values_levels[1] <- paste0("<=", self$mic_values_levels[1])
        self$mic_values_levels[length(self$mic_values_levels)] <- paste0(">=", self$mic_values_levels[length(self$mic_values_levels)])
      }

      self$mic_values_log <- log2(as.double(self$mic_values_rescaled))
      if (aest == "y" && "group" %in% colnames(df) && "x" %in% colnames(df)) {
        df$group <- as.integer(factor(df$x))
      }
      df[[aest]] <- self$mic_values_log
    }
    df
  }

  scale$breaks <- function(..., self) {
    if (!is.null(self$mic_breaks_set)) {
      if (is.function(self$mic_breaks_set)) {
        self$mic_breaks_set(...)
      } else {
        log2(as.mic(self$mic_breaks_set))
      }
    } else {
      log2(as.mic(self$mic_values_levels))
    }
  }
  scale$labels <- function(..., self) {
    if (is.null(self$mic_breaks_set)) {
      self$mic_values_levels
    } else {
      breaks <- tryCatch(scale$breaks(), error = function(e) NULL)
      if (!is.null(breaks)) {
        # for when breaks are set by the user
        2^breaks
      } else {
        self$mic_values_levels
      }
    }
  }

  scale$limits <- function(x, ..., self) {
    if (!is.null(self$mic_limits_set)) {
      if (is.function(self$mic_limits_set)) {
        self$mic_limits_set(...)
      } else {
        log2(as.mic(self$mic_limits_set))
      }
    } else {
      rng <- range(log2(as.mic(self$mic_values_levels)))
      # add 0.5 extra space
      rng <- c(rng[1] - 0.5, rng[2] + 0.5)
      if (!is.na(x[1]) && x[1] == 0) {
        # scale that start at 0 must remain so, e.g. in case of geom_col()
        rng[1] <- 0
      }
      rng
    }
  }

  scale
}

#' @export
#' @rdname plot
scale_x_mic <- function(keep_operators = "edges", mic_range = NULL, ...) {
  meet_criteria(keep_operators, allow_class = c("character", "logical"), is_in = c("all", "none", "edges", FALSE, TRUE), has_length = 1)
  meet_criteria(mic_range, allow_class = c("numeric", "integer", "logical", "mic"), has_length = 2, allow_NA = TRUE, allow_NULL = TRUE)
  create_scale_mic("x", keep_operators = keep_operators, mic_range = mic_range, ...)
}

#' @export
#' @rdname plot
scale_y_mic <- function(keep_operators = "edges", mic_range = NULL, ...) {
  meet_criteria(keep_operators, allow_class = c("character", "logical"), is_in = c("all", "none", "edges", FALSE, TRUE), has_length = 1)
  meet_criteria(mic_range, allow_class = c("numeric", "integer", "logical", "mic"), has_length = 2, allow_NA = TRUE, allow_NULL = TRUE)
  create_scale_mic("y", keep_operators = keep_operators, mic_range = mic_range, ...)
}

#' @export
#' @rdname plot
scale_colour_mic <- function(keep_operators = "edges", mic_range = NULL, ...) {
  meet_criteria(keep_operators, allow_class = c("character", "logical"), is_in = c("all", "none", "edges", FALSE, TRUE), has_length = 1)
  meet_criteria(mic_range, allow_class = c("numeric", "integer", "logical", "mic"), has_length = 2, allow_NA = TRUE, allow_NULL = TRUE)
  create_scale_mic("colour", keep_operators = keep_operators, mic_range = mic_range, ...)
}

#' @export
#' @rdname plot
#' @usage NULL
scale_color_mic <- scale_colour_mic

#' @export
#' @rdname plot
scale_fill_mic <- function(keep_operators = "edges", mic_range = NULL, ...) {
  meet_criteria(keep_operators, allow_class = c("character", "logical"), is_in = c("all", "none", "edges", FALSE, TRUE), has_length = 1)
  meet_criteria(mic_range, allow_class = c("numeric", "integer", "logical", "mic"), has_length = 2, allow_NA = TRUE, allow_NULL = TRUE)
  create_scale_mic("fill", keep_operators = keep_operators, mic_range = mic_range, ...)
}

create_scale_sir <- function(aesthetics, colours_SIR, language, eucast_I, ...) {
  args <- list(...)
  args[c("value", "labels", "limits")] <- NULL

  colours_SIR <- expand_SIR_colours(colours_SIR)

  if (identical(aesthetics, "x")) {
    ggplot_fn <- ggplot2::scale_x_discrete
  } else {
    ggplot_fn <- ggplot2::scale_discrete_manual
    args <- c(
      args,
      list(
        aesthetics = aesthetics,
        values = c(
          S = colours_SIR[1],
          SDD = colours_SIR[2],
          I = colours_SIR[3],
          R = colours_SIR[4],
          NI = "grey30"
        )
      )
    )
  }
  scale <- do.call(ggplot_fn, args)

  scale$labels <- function(x) {
    stop_ifnot(all(x %in% c(levels(NA_sir_), NA)),
      "Apply `scale_", aesthetics[1], "_sir()` to a variable of class 'sir', see `?as.sir`.",
      call = FALSE
    )
    x <- as.character(as.sir(x))
    if (!is.null(language)) {
      x[x == "S"] <- "(S) Susceptible"
      x[x == "SDD"] <- "(SDD) Susceptible dose-dependent"
      if (eucast_I == TRUE) {
        x[x == "I"] <- "(I) Susceptible, incr. exp."
      } else {
        x[x == "I"] <- "(I) Intermediate"
      }
      x[x == "R"] <- "(R) Resistant"
      x[x == "NI"] <- "(NI) Non-interpretable"
      x <- translate_AMR(x, language = language)
    }
    x
  }
  scale$limits <- function(x, ...) {
    # force SIR in the right order
    as.character(sort(factor(x, levels = levels(NA_sir_))))
  }

  scale
}

#' @rdname plot
#' @export
scale_x_sir <- function(colours_SIR = c(
                          S = "#3CAEA3",
                          SDD = "#8FD6C4",
                          I = "#F6D55C",
                          R = "#ED553B"
                        ),
                        language = get_AMR_locale(),
                        eucast_I = getOption("AMR_guideline", "EUCAST") == "EUCAST",
                        ...) {
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(eucast_I, allow_class = "logical", has_length = 1)
  create_scale_sir(aesthetics = "x", colours_SIR = colours_SIR, language = language, eucast_I = eucast_I)
}

#' @rdname plot
#' @export
scale_colour_sir <- function(colours_SIR = c(
                               S = "#3CAEA3",
                               SDD = "#8FD6C4",
                               I = "#F6D55C",
                               R = "#ED553B"
                             ),
                             language = get_AMR_locale(),
                             eucast_I = getOption("AMR_guideline", "EUCAST") == "EUCAST",
                             ...) {
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(eucast_I, allow_class = "logical", has_length = 1)
  args <- list(...)
  args$colours_SIR <- colours_SIR
  args$language <- language
  args$eucast_I <- eucast_I
  if (is.null(args$aesthetics)) {
    args$aesthetics <- "colour"
  }
  do.call(create_scale_sir, args)
}

#' @export
#' @rdname plot
#' @usage NULL
scale_color_sir <- scale_colour_sir

#' @rdname plot
#' @export
scale_fill_sir <- function(colours_SIR = c(
                             S = "#3CAEA3",
                             SDD = "#8FD6C4",
                             I = "#F6D55C",
                             R = "#ED553B"
                           ),
                           language = get_AMR_locale(),
                           eucast_I = getOption("AMR_guideline", "EUCAST") == "EUCAST",
                           ...) {
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(eucast_I, allow_class = "logical", has_length = 1)
  args <- list(...)
  args$colours_SIR <- colours_SIR
  args$language <- language
  args$eucast_I <- eucast_I
  if (is.null(args$aesthetics)) {
    args$aesthetics <- "fill"
  }
  do.call(create_scale_sir, args)
}

#' @method plot mic
#' @importFrom graphics barplot axis mtext legend
#' @export
#' @rdname plot
plot.mic <- function(x,
                     mo = NULL,
                     ab = NULL,
                     guideline = getOption("AMR_guideline", "EUCAST"),
                     main = deparse(substitute(x)),
                     ylab = translate_AMR("Frequency", language = language),
                     xlab = translate_AMR("Minimum Inhibitory Concentration (mg/L)", language = language),
                     colours_SIR = c(
                       S = "#3CAEA3",
                       SDD = "#8FD6C4",
                       I = "#F6D55C",
                       R = "#ED553B"
                     ),
                     language = get_AMR_locale(),
                     expand = TRUE,
                     include_PKPD = getOption("AMR_include_PKPD", TRUE),
                     breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                     ...) {
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  x <- as.mic(x) # make sure that currently implemented MIC levels are used
  main <- gsub(" +", " ", paste0(main, collapse = " "))
  colours_SIR <- expand_SIR_colours(colours_SIR)

  x <- plotrange_as_table(x, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    fn = as.mic,
    language = language,
    method = "MIC",
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    ...
  )
  barplot(x,
    col = cols_sub$cols,
    main = main,
    ylim = c(0, max(x) * ifelse(any(colours_SIR %in% cols_sub$cols), 1.1, 1)),
    ylab = ylab,
    xlab = xlab,
    axes = FALSE
  )
  axis(2, seq(0, max(x)))
  if (!is.null(cols_sub$sub)) {
    mtext(side = 3, line = 0.5, adj = 0.5, cex = 0.75, cols_sub$sub)
  }

  if (any(colours_SIR %in% cols_sub$cols)) {
    legend_txt <- character(0)
    legend_col <- character(0)
    if (any(cols_sub$cols == colours_SIR[1] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "(S) Susceptible")
      legend_col <- colours_SIR[1]
    }
    if (any(cols_sub$cols == colours_SIR[2] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "(SDD) Susceptible dose-dependent")
      legend_col <- c(legend_col, colours_SIR[2])
    }
    if (any(cols_sub$cols == colours_SIR[3] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, paste("(I)", plot_name_of_I(cols_sub$guideline)))
      legend_col <- c(legend_col, colours_SIR[3])
    }
    if (any(cols_sub$cols == colours_SIR[4] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "(R) Resistant")
      legend_col <- c(legend_col, colours_SIR[4])
    }

    legend("top",
      x.intersp = 0.5,
      legend = translate_into_language(legend_txt, language = language),
      fill = legend_col,
      horiz = TRUE,
      cex = 0.75,
      box.lwd = 0,
      box.col = "#FFFFFF55",
      bg = "#FFFFFF55"
    )
  }
}

#' @method barplot mic
#' @export
#' @noRd
barplot.mic <- function(height,
                        mo = NULL,
                        ab = NULL,
                        guideline = getOption("AMR_guideline", "EUCAST"),
                        main = deparse(substitute(height)),
                        ylab = translate_AMR("Frequency", language = language),
                        xlab = translate_AMR("Minimum Inhibitory Concentration (mg/L)", language = language),
                        colours_SIR = c(
                          S = "#3CAEA3",
                          SDD = "#8FD6C4",
                          I = "#F6D55C",
                          R = "#ED553B"
                        ),
                        language = get_AMR_locale(),
                        expand = TRUE,
                        ...) {
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  main <- gsub(" +", " ", paste0(main, collapse = " "))

  height <- as.mic(height) # make sure that currently implemented MIC levels are used

  plot(
    x = height,
    main = main,
    ylab = ylab,
    xlab = xlab,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    ...
  )
}

#' @method autoplot mic
#' @rdname plot
# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(ggplot2::autoplot, mic)
autoplot.mic <- function(object,
                         mo = NULL,
                         ab = NULL,
                         guideline = getOption("AMR_guideline", "EUCAST"),
                         title = deparse(substitute(object)),
                         ylab = translate_AMR("Frequency", language = language),
                         xlab = translate_AMR("Minimum Inhibitory Concentration (mg/L)", language = language),
                         colours_SIR = c(
                           S = "#3CAEA3",
                           SDD = "#8FD6C4",
                           I = "#F6D55C",
                           R = "#ED553B"
                         ),
                         language = get_AMR_locale(),
                         expand = TRUE,
                         include_PKPD = getOption("AMR_include_PKPD", TRUE),
                         breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                         ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(title, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  if ("main" %in% names(list(...))) {
    title <- list(...)$main
  }
  if (!is.null(title)) {
    title <- gsub(" +", " ", paste0(title, collapse = " "))
  }

  colours_SIR <- expand_SIR_colours(colours_SIR)

  object <- as.mic(object) # make sure that currently implemented MIC levels are used
  x <- plotrange_as_table(object, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    fn = as.mic,
    language = language,
    method = "MIC",
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    ...
  )
  df <- as.data.frame(x, stringsAsFactors = TRUE)
  colnames(df) <- c("mic", "count")
  df$cols <- cols_sub$cols
  df$cols[df$cols == colours_SIR[1]] <- "(S) Susceptible"
  df$cols[df$cols == colours_SIR[2]] <- "(SDD) Susceptible dose-dependent"
  df$cols[df$cols == colours_SIR[3]] <- paste("(I)", plot_name_of_I(cols_sub$guideline))
  df$cols[df$cols == colours_SIR[4]] <- "(R) Resistant"
  df$cols <- factor(translate_into_language(df$cols, language = language),
    levels = translate_into_language(
      c(
        "(S) Susceptible",
        "(SDD) Susceptible dose-dependent",
        paste("(I)", plot_name_of_I(cols_sub$guideline)),
        "(R) Resistant"
      ),
      language = language
    ),
    ordered = TRUE
  )
  p <- ggplot2::ggplot(df)

  if (any(colours_SIR %in% cols_sub$cols)) {
    vals <- c(
      "(S) Susceptible" = colours_SIR[1],
      "(SDD) Susceptible dose-dependent" = colours_SIR[2],
      "(I) Susceptible, incr. exp." = colours_SIR[3],
      "(I) Intermediate" = colours_SIR[3],
      "(R) Resistant" = colours_SIR[4],
      "(NI) Non-interpretable" = "grey30"
    )
    names(vals) <- translate_into_language(names(vals), language = language)
    p <- p +
      ggplot2::geom_col(ggplot2::aes(x = mic, y = count, fill = cols)) +
      # limits = force is needed because of a ggplot2 >= 3.3.4 bug (#4511)
      ggplot2::scale_fill_manual(
        values = vals,
        name = NULL,
        limits = force
      )
  } else {
    p <- p +
      ggplot2::geom_col(ggplot2::aes(x = mic, y = count))
  }

  p +
    ggplot2::labs(title = title, x = xlab, y = ylab, subtitle = cols_sub$sub)
}

#' @method fortify mic
#' @noRd
# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(ggplot2::fortify, mic)
fortify.mic <- function(object, ...) {
  object <- as.mic(object) # make sure that currently implemented MIC levels are used
  stats::setNames(
    as.data.frame(plotrange_as_table(object, expand = FALSE)),
    c("x", "y")
  )
}


#' @method plot disk
#' @export
#' @importFrom graphics barplot axis mtext legend
#' @rdname plot
plot.disk <- function(x,
                      main = deparse(substitute(x)),
                      ylab = translate_AMR("Frequency", language = language),
                      xlab = translate_AMR("Disk diffusion diameter (mm)", language = language),
                      mo = NULL,
                      ab = NULL,
                      guideline = getOption("AMR_guideline", "EUCAST"),
                      colours_SIR = c(
                        S = "#3CAEA3",
                        SDD = "#8FD6C4",
                        I = "#F6D55C",
                        R = "#ED553B"
                      ),
                      language = get_AMR_locale(),
                      expand = TRUE,
                      include_PKPD = getOption("AMR_include_PKPD", TRUE),
                      breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                      ...) {
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  main <- gsub(" +", " ", paste0(main, collapse = " "))
  colours_SIR <- expand_SIR_colours(colours_SIR)

  x <- plotrange_as_table(x, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    fn = as.disk,
    language = language,
    method = "disk",
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    ...
  )

  barplot(x,
    col = cols_sub$cols,
    main = main,
    ylim = c(0, max(x) * ifelse(any(colours_SIR %in% cols_sub$cols), 1.1, 1)),
    ylab = ylab,
    xlab = xlab,
    axes = FALSE
  )
  axis(2, seq(0, max(x)))
  if (!is.null(cols_sub$sub)) {
    mtext(side = 3, line = 0.5, adj = 0.5, cex = 0.75, cols_sub$sub)
  }

  if (any(colours_SIR %in% cols_sub$cols)) {
    legend_txt <- character(0)
    legend_col <- character(0)
    if (any(cols_sub$cols == colours_SIR[4] & cols_sub$count > 0)) {
      legend_txt <- "(R) Resistant"
      legend_col <- colours_SIR[4]
    }
    if (any(cols_sub$cols == colours_SIR[3] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, paste("(I)", plot_name_of_I(cols_sub$guideline)))
      legend_col <- c(legend_col, colours_SIR[3])
    }
    if (any(cols_sub$cols == colours_SIR[2] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "(SDD) Susceptible dose-dependent")
      legend_col <- c(legend_col, colours_SIR[2])
    }
    if (any(cols_sub$cols == colours_SIR[1] & cols_sub$count > 0)) {
      legend_txt <- c(legend_txt, "(S) Susceptible")
      legend_col <- c(legend_col, colours_SIR[1])
    }
    legend("top",
      x.intersp = 0.5,
      legend = translate_into_language(legend_txt, language = language),
      fill = legend_col,
      horiz = TRUE,
      cex = 0.75,
      box.lwd = 0,
      box.col = "#FFFFFF55",
      bg = "#FFFFFF55"
    )
  }
}

#' @method barplot disk
#' @export
#' @noRd
barplot.disk <- function(height,
                         main = deparse(substitute(height)),
                         ylab = translate_AMR("Frequency", language = language),
                         xlab = translate_AMR("Disk diffusion diameter (mm)", language = language),
                         mo = NULL,
                         ab = NULL,
                         guideline = getOption("AMR_guideline", "EUCAST"),
                         colours_SIR = c(
                           S = "#3CAEA3",
                           SDD = "#8FD6C4",
                           I = "#F6D55C",
                           R = "#ED553B"
                         ),
                         language = get_AMR_locale(),
                         expand = TRUE,
                         ...) {
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  main <- gsub(" +", " ", paste0(main, collapse = " "))

  plot(
    x = height,
    main = main,
    ylab = ylab,
    xlab = xlab,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    ...
  )
}

#' @method autoplot disk
#' @rdname plot
# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(ggplot2::autoplot, disk)
autoplot.disk <- function(object,
                          mo = NULL,
                          ab = NULL,
                          title = deparse(substitute(object)),
                          ylab = translate_AMR("Frequency", language = language),
                          xlab = translate_AMR("Disk diffusion diameter (mm)", language = language),
                          guideline = getOption("AMR_guideline", "EUCAST"),
                          colours_SIR = c(
                            S = "#3CAEA3",
                            SDD = "#8FD6C4",
                            I = "#F6D55C",
                            R = "#ED553B"
                          ),
                          language = get_AMR_locale(),
                          expand = TRUE,
                          include_PKPD = getOption("AMR_include_PKPD", TRUE),
                          breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                          ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(title, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(mo, allow_class = c("mo", "character"), allow_NULL = TRUE)
  meet_criteria(ab, allow_class = c("ab", "character"), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  if ("main" %in% names(list(...))) {
    title <- list(...)$main
  }
  if (!is.null(title)) {
    title <- gsub(" +", " ", paste0(title, collapse = " "))
  }

  colours_SIR <- expand_SIR_colours(colours_SIR)

  x <- plotrange_as_table(object, expand = expand)
  cols_sub <- plot_colours_subtitle_guideline(
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    colours_SIR = colours_SIR,
    fn = as.disk,
    language = language,
    method = "disk",
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    ...
  )
  df <- as.data.frame(x, stringsAsFactors = TRUE)
  colnames(df) <- c("disk", "count")
  df$cols <- cols_sub$cols
  df$cols[df$cols == colours_SIR[1]] <- "(S) Susceptible"
  df$cols[df$cols == colours_SIR[2]] <- "(SDD) Susceptible dose-dependent"
  df$cols[df$cols == colours_SIR[3]] <- paste("(I)", plot_name_of_I(cols_sub$guideline))
  df$cols[df$cols == colours_SIR[4]] <- "(R) Resistant"
  df$cols <- factor(translate_into_language(df$cols, language = language),
    levels = translate_into_language(
      c(
        "(S) Susceptible",
        paste("(I)", plot_name_of_I(cols_sub$guideline)),
        "(R) Resistant"
      ),
      language = language
    ),
    ordered = TRUE
  )
  p <- ggplot2::ggplot(df)

  if (any(colours_SIR %in% cols_sub$cols)) {
    vals <- c(
      "(S) Susceptible" = colours_SIR[1],
      "(SDD) Susceptible dose-dependent" = colours_SIR[2],
      "(I) Susceptible, incr. exp." = colours_SIR[3],
      "(I) Intermediate" = colours_SIR[3],
      "(R) Resistant" = colours_SIR[4],
      "(NI) Non-interpretable" = "grey30"
    )
    names(vals) <- translate_into_language(names(vals), language = language)
    p <- p +
      ggplot2::geom_col(ggplot2::aes(x = disk, y = count, fill = cols)) +
      # limits = force is needed because of a ggplot2 >= 3.3.4 bug (#4511)
      ggplot2::scale_fill_manual(
        values = vals,
        name = NULL,
        limits = force
      )
  } else {
    p <- p +
      ggplot2::geom_col(ggplot2::aes(x = disk, y = count))
  }

  p +
    ggplot2::labs(title = title, x = xlab, y = ylab, subtitle = cols_sub$sub)
}

#' @method fortify disk
#' @noRd
# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(ggplot2::fortify, disk)
fortify.disk <- function(object, ...) {
  stats::setNames(
    as.data.frame(plotrange_as_table(object, expand = FALSE)),
    c("x", "y")
  )
}

#' @method plot sir
#' @export
#' @importFrom graphics plot text axis
#' @rdname plot
plot.sir <- function(x,
                     ylab = translate_AMR("Percentage", language = language),
                     xlab = translate_AMR("Antimicrobial Interpretation", language = language),
                     main = deparse(substitute(x)),
                     language = get_AMR_locale(),
                     ...) {
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)

  data <- as.data.frame(table(x), stringsAsFactors = FALSE)
  colnames(data) <- c("x", "n")
  data$s <- round((data$n / sum(data$n)) * 100, 1)

  if (!"S" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "S", n = 0, s = 0, stringsAsFactors = FALSE))
  }
  if (!"SDD" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "SDD", n = 0, s = 0, stringsAsFactors = FALSE))
  }
  if (!"I" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "I", n = 0, s = 0, stringsAsFactors = FALSE))
  }
  if (!"R" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "R", n = 0, s = 0, stringsAsFactors = FALSE))
  }
  if (!"NI" %in% data$x) {
    data <- rbind_AMR(data, data.frame(x = "NI", n = 0, s = 0, stringsAsFactors = FALSE))
  }

  data <- data[!(data$n == 0 & data$x %in% c("SDD", "I", "NI")), , drop = FALSE]
  data$x <- factor(data$x, levels = intersect(unique(data$x), c("S", "SDD", "I", "R", "NI")), ordered = TRUE)

  ymax <- pm_if_else(max(data$s) > 95, 105, 100)

  plot(
    x = data$x,
    y = data$s,
    lwd = 2,
    ylim = c(0, ymax),
    ylab = ylab,
    xlab = xlab,
    main = main,
    axes = FALSE
  )
  # x axis
  axis(side = 1, at = 1:pm_n_distinct(data$x), labels = levels(data$x), lwd = 0)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, 100, 5))

  text(
    x = data$x,
    y = data$s + 4,
    labels = paste0(data$s, "% (n = ", data$n, ")")
  )
}


#' @method barplot sir
#' @importFrom graphics barplot axis
#' @export
#' @noRd
barplot.sir <- function(height,
                        main = deparse(substitute(height)),
                        xlab = translate_AMR("Antimicrobial Interpretation", language = language),
                        ylab = translate_AMR("Frequency", language = language),
                        colours_SIR = c(
                          S = "#3CAEA3",
                          SDD = "#8FD6C4",
                          I = "#F6D55C",
                          R = "#ED553B"
                        ),
                        language = get_AMR_locale(),
                        expand = TRUE,
                        ...) {
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(main, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))
  language <- validate_language(language)
  meet_criteria(expand, allow_class = "logical", has_length = 1)

  colours_SIR <- expand_SIR_colours(colours_SIR)

  # add SDD and N to colours
  colours_SIR <- c(colours_SIR, "grey30")
  main <- gsub(" +", " ", paste0(main, collapse = " "))

  x <- table(height)
  # remove missing I, SDD, and N
  colours_SIR <- colours_SIR[!(names(x) %in% c("SDD", "I", "NI") & x == 0)]
  x <- x[!(names(x) %in% c("SDD", "I", "NI") & x == 0)]
  # plot it
  barplot(x,
    col = colours_SIR,
    xlab = xlab,
    main = main,
    ylab = ylab,
    axes = FALSE
  )
  axis(2, seq(0, max(x)))
}

#' @method autoplot sir
#' @rdname plot
# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(ggplot2::autoplot, sir)
autoplot.sir <- function(object,
                         title = deparse(substitute(object)),
                         xlab = translate_AMR("Antimicrobial Interpretation", language = language),
                         ylab = translate_AMR("Frequency", language = language),
                         colours_SIR = c(
                           S = "#3CAEA3",
                           SDD = "#8FD6C4",
                           I = "#F6D55C",
                           R = "#ED553B"
                         ),
                         language = get_AMR_locale(),
                         ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(title, allow_class = "character", allow_NULL = TRUE)
  meet_criteria(ylab, allow_class = "character", has_length = 1)
  meet_criteria(xlab, allow_class = "character", has_length = 1)
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))

  if ("main" %in% names(list(...))) {
    title <- list(...)$main
  }
  if (!is.null(title)) {
    title <- gsub(" +", " ", paste0(title, collapse = " "))
  }

  colours_SIR <- expand_SIR_colours(colours_SIR)

  df <- as.data.frame(table(object), stringsAsFactors = TRUE)
  colnames(df) <- c("x", "n")
  df <- df[!(df$n == 0 & df$x %in% c("SDD", "I", "NI")), , drop = FALSE]
  ggplot2::ggplot(df) +
    ggplot2::geom_col(ggplot2::aes(x = x, y = n, fill = x)) +
    # limits = force is needed because of a ggplot2 >= 3.3.4 bug (#4511)
    ggplot2::scale_fill_manual(
      values = c(
        "S" = colours_SIR[1],
        "SDD" = colours_SIR[2],
        "I" = colours_SIR[3],
        "R" = colours_SIR[4],
        "NI" = "grey30"
      ),
      limits = force
    ) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme(legend.position = "none")
}

#' @method fortify sir
#' @noRd
# this prevents the requirement for putting the dependency in Imports:
#' @rawNamespace if(getRversion() >= "3.0.0") S3method(ggplot2::fortify, sir)
fortify.sir <- function(object, ...) {
  stats::setNames(
    as.data.frame(table(object)),
    c("x", "y")
  )
}

plotrange_as_table <- function(x, expand, keep_operators = "all", mic_range = NULL) {
  x <- x[!is.na(x)]
  if (is.mic(x)) {
    x <- as.mic(x, keep_operators = keep_operators)
    if (expand == TRUE) {
      # expand range for MIC by adding common intermediate factors levels
      if (!is.null(mic_range) && !all(is.na(mic_range))) {
        # base on mic_range
        `%na_or%` <- function(x, y) if (is.na(x)) y else x
        extra_range <- COMMON_MIC_VALUES[COMMON_MIC_VALUES >= (mic_range[1] %na_or% min(x, na.rm = TRUE)) & COMMON_MIC_VALUES <= (mic_range[2] %na_or% max(x, na.rm = TRUE))]
      } else {
        # base on x
        extra_range <- COMMON_MIC_VALUES[COMMON_MIC_VALUES > min(x, na.rm = TRUE) & COMMON_MIC_VALUES < max(x, na.rm = TRUE)]
      }
      # remove the ones that are in 25% range of user values
      extra_range <- extra_range[!vapply(FUN.VALUE = logical(1), extra_range, function(r) any(abs(r - x) / x < 0.25, na.rm = TRUE))]
      nms <- extra_range
      extra_range <- rep(0, length(extra_range))
      names(extra_range) <- nms
      x <- table(droplevels(x, as.mic = FALSE))
      extra_range <- extra_range[!names(extra_range) %in% names(x) & names(extra_range) %in% VALID_MIC_LEVELS]
      x <- as.table(c(x, extra_range))
    } else {
      x <- table(droplevels(x, as.mic = FALSE))
    }
    x <- x[order(as.double(as.mic(names(x))))]
  } else if (is.disk(x)) {
    if (expand == TRUE) {
      # expand range for disks from lowest to highest so all mm's in between also print
      extra_range <- rep(0, max(x) - min(x) - 1)
      names(extra_range) <- seq(min(x) + 1, max(x) - 1)
      x <- table(x)
      extra_range <- extra_range[!names(extra_range) %in% names(x)]
      x <- as.table(c(x, extra_range))
    } else {
      x <- table(x)
    }
    x <- x[order(as.double(names(x)))]
  }
  as.table(x)
}

ggplot2_get_from_dots <- function(arg, default, ...) {
  dots <- list(...)
  if (!arg %in% names(dots)) {
    default
  } else {
    dots[[arg]]
  }
}

plot_name_of_I <- function(guideline) {
  if (guideline %unlike% "CLSI" && as.double(gsub("[^0-9]+", "", guideline)) >= 2019) {
    # interpretation since 2019
    "Susceptible, incr. exp."
  } else {
    # interpretation until 2019
    "Intermediate"
  }
}

plot_colours_subtitle_guideline <- function(x, mo, ab, guideline, colours_SIR, fn, language, method, breakpoint_type, include_PKPD, ...) {
  stop_if(length(x) == 0, "no observations to plot", call = FALSE)

  guideline <- get_guideline(guideline, AMR::clinical_breakpoints)

  # store previous interpretations to backup
  sir_history <- AMR_env$sir_interpretation_history
  # and clear previous interpretations
  AMR_env$sir_interpretation_history <- AMR_env$sir_interpretation_history[0, , drop = FALSE]

  if (!is.null(mo) && !is.null(ab)) {
    # interpret and give colour based on MIC values
    mo <- as.mo(mo)
    moname <- mo_name(mo, language = language)
    ab <- as.ab(ab)
    abname <- ab_name(ab, language = language)

    sir <- suppressWarnings(suppressMessages(as.sir(fn(names(x)), mo = mo, ab = ab, guideline = guideline, include_screening = FALSE, include_PKPD = include_PKPD, breakpoint_type = breakpoint_type, ...)))
    guideline_txt <- guideline
    if (all(is.na(sir))) {
      sir_screening <- suppressWarnings(suppressMessages(as.sir(fn(names(x)), mo = mo, ab = ab, guideline = guideline, include_screening = TRUE, include_PKPD = include_PKPD, breakpoint_type = breakpoint_type, ...)))
      if (!all(is.na(sir_screening))) {
        message_(
          "Only ", guideline, " ", method, " interpretations found for ",
          ab_name(ab, language = NULL, tolower = TRUE), " in ", italicise(moname), " for screening"
        )
        sir <- sir_screening
        guideline_txt <- paste0("(Screen, ", guideline_txt, ")")
      } else {
        message_(
          "No ", guideline, " ", method, " interpretations found for ",
          ab_name(ab, language = NULL, tolower = TRUE), " in ", italicise(moname)
        )
        guideline_txt <- paste0("(", guideline_txt, ")")
      }
    } else {
      if (isTRUE(list(...)$uti)) {
        guideline_txt <- paste("UTIs,", guideline_txt)
      }
      ref_tbl <- paste0('"', unique(AMR_env$sir_interpretation_history$ref_table), '"', collapse = "/")
      guideline_txt <- paste0("(", guideline_txt, ": ", ref_tbl, ")")
    }
    cols <- character(length = length(sir))
    cols[is.na(sir)] <- "#BEBEBE"
    cols[sir == "S"] <- colours_SIR[1]
    cols[sir == "SDD"] <- colours_SIR[2]
    cols[sir == "I"] <- colours_SIR[3]
    cols[sir == "R"] <- colours_SIR[4]
    cols[sir == "NI"] <- "grey30"
    sub <- bquote(.(abname) ~ "-" ~ italic(.(moname)) ~ .(guideline_txt))
  } else {
    cols <- "#BEBEBE"
    sub <- NULL
  }

  # restore previous interpretations to backup
  AMR_env$sir_interpretation_history <- sir_history

  list(cols = cols, count = as.double(x), sub = sub, guideline = guideline)
}

#' @rdname plot
#' @export
facet_sir <- function(facet = c("interpretation", "antibiotic"), nrow = NULL) {
  facet <- facet[1]
  stop_ifnot_installed("ggplot2")
  meet_criteria(facet, allow_class = "character", has_length = 1)
  meet_criteria(nrow, allow_class = c("numeric", "integer"), has_length = 1, allow_NULL = TRUE, is_positive = TRUE, is_finite = TRUE)

  facet_deparse <- deparse(substitute(facet))
  if (facet_deparse != "facet") {
    facet <- facet_deparse
  }
  if (facet %like% '".*"') {
    facet <- substr(facet, 2, nchar(facet) - 1)
  }

  if (tolower(facet) %in% tolower(c("SIR", "sir", "interpretations", "result"))) {
    facet <- "interpretation"
  } else if (tolower(facet) %in% tolower(c("ab", "abx", "antimicrobials"))) {
    facet <- "antibiotic"
  }

  ggplot2::facet_wrap(facets = facet, scales = "free_x", nrow = nrow)
}

#' @rdname plot
#' @export
scale_y_percent <- function(breaks = function(x) seq(0, max(x, na.rm = TRUE), 0.1), limits = c(0, NA)) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(breaks, allow_class = c("numeric", "integer", "function"))
  meet_criteria(limits, allow_class = c("numeric", "integer"), has_length = 2, allow_NULL = TRUE, allow_NA = TRUE)

  if (!is.function(breaks) && all(breaks[breaks != 0] > 1)) {
    breaks <- breaks / 100
  }
  ggplot2::scale_y_continuous(
    breaks = breaks,
    labels = if (is.function(breaks)) function(x) percentage(breaks(x)) else percentage(breaks),
    limits = limits
  )
}

#' @rdname plot
#' @export
scale_sir_colours <- function(...,
                              aesthetics,
                              colours_SIR = c(
                                S = "#3CAEA3",
                                SDD = "#8FD6C4",
                                I = "#F6D55C",
                                R = "#ED553B"
                              )) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(aesthetics, allow_class = "character", is_in = c("alpha", "colour", "color", "fill", "linetype", "shape", "size"))
  meet_criteria(colours_SIR, allow_class = "character", has_length = c(1, 3, 4))

  if ("fill" %in% aesthetics && message_not_thrown_before("scale_sir_colours", "fill", entire_session = TRUE)) {
    warning_("Using `scale_sir_colours()` for the `fill` aesthetic has been superseded by `scale_fill_sir()`, please use that instead. This warning will be shown once per session.")
  }
  if (any(c("colour", "color") %in% aesthetics) && message_not_thrown_before("scale_sir_colours", "colour", entire_session = TRUE)) {
    warning_("Using `scale_sir_colours()` for the `colour` aesthetic has been superseded by `scale_colour_sir()`, please use that instead. This warning will be shown once per session.")
  }

  if ("colours" %in% names(list(...))) {
    colours_SIR <- list(...)$colours
  }

  colours_SIR <- expand_SIR_colours(colours_SIR, unname = FALSE)

  # behaviour when coming from ggplot_sir()
  if ("colours" %in% names(list(...))) {
    # limits = force is needed in ggplot2 3.3.4 and 3.3.5, see here;
    # https://github.com/tidyverse/ggplot2/issues/4511#issuecomment-866185530
    return(ggplot2::scale_fill_manual(values = colours_SIR, limits = force, aesthetics = aesthetics))
  }
  if (identical(unlist(list(...)), FALSE)) {
    return(invisible())
  }

  colours_SIR <- unname(colours_SIR)

  names_susceptible <- c("S", "SI", "IS", "S+I", "I+S", "susceptible", "Susceptible")
  names_susceptible_dose_dep <- c("SDD", "susceptible dose-dependent", "Susceptible dose-dependent")
  names_incr_exposure <- c(
    "I", "intermediate", "increased exposure", "incr. exposure",
    "Increased exposure", "Incr. exposure", "Susceptible, incr. exp."
  )
  names_resistant <- c("R", "IR", "RI", "R+I", "I+R", "resistant", "Resistant")

  susceptible <- rep(colours_SIR[1], length(names_susceptible))
  names(susceptible) <- names_susceptible
  susceptible_dose_dep <- rep(colours_SIR[2], length(names_susceptible_dose_dep))
  names(susceptible_dose_dep) <- names_susceptible_dose_dep
  incr_exposure <- rep(colours_SIR[3], length(names_incr_exposure))
  names(incr_exposure) <- names_incr_exposure
  resistant <- rep(colours_SIR[4], length(names_resistant))
  names(resistant) <- names_resistant

  original_cols <- c(susceptible, susceptible_dose_dep, incr_exposure, resistant)
  dots <- c(...)
  # replace S, SDD, I, R as colours: scale_sir_colours(mydatavalue = "S")
  dots[dots == "S"] <- colours_SIR[1]
  dots[dots == "SDD"] <- colours_SIR[2]
  dots[dots == "I"] <- colours_SIR[3]
  dots[dots == "R"] <- colours_SIR[4]
  cols <- replace(original_cols, names(dots), dots)
  # limits = force is needed in ggplot2 3.3.4 and 3.3.5, see here;
  # https://github.com/tidyverse/ggplot2/issues/4511#issuecomment-866185530
  ggplot2::scale_discrete_manual(aesthetics = aesthetics, values = cols, limits = force)
}

#' @export
#' @rdname plot
#' @usage NULL
scale_sir_colors <- scale_sir_colours

#' @rdname plot
#' @export
theme_sir <- function() {
  stop_ifnot_installed("ggplot2")
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey75"),
      # center title and subtitle
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
}

#' @rdname plot
#' @export
labels_sir_count <- function(position = NULL,
                             x = "antibiotic",
                             translate_ab = "name",
                             minimum = 30,
                             language = get_AMR_locale(),
                             combine_SI = TRUE,
                             datalabels.size = 3,
                             datalabels.colour = "grey15") {
  stop_ifnot_installed("ggplot2")
  meet_criteria(position, allow_class = "character", has_length = 1, is_in = c("fill", "stack", "dodge"), allow_NULL = TRUE)
  meet_criteria(x, allow_class = "character", has_length = 1)
  meet_criteria(translate_ab, allow_class = c("character", "logical"), has_length = 1, allow_NA = TRUE)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  language <- validate_language(language)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(datalabels.size, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = TRUE)
  meet_criteria(datalabels.colour, allow_class = "character", has_length = 1)

  if (is.null(position)) {
    position <- "fill"
  }
  if (identical(position, "fill")) {
    position <- ggplot2::position_fill(vjust = 0.5, reverse = TRUE)
  }

  x_name <- x
  ggplot2::geom_text(
    mapping = utils::modifyList(ggplot2::aes(), list(label = str2lang("lbl"), x = str2lang(x), y = str2lang("value"))),
    position = position,
    inherit.aes = FALSE,
    size = datalabels.size,
    colour = datalabels.colour,
    lineheight = 0.75,
    data = function(x) {
      transformed <- sir_df(
        data = x,
        translate_ab = translate_ab,
        combine_SI = combine_SI,
        minimum = minimum,
        language = language
      )
      transformed$gr <- transformed[, x_name, drop = TRUE]
      transformed %pm>%
        pm_group_by(gr) %pm>%
        pm_mutate(lbl = paste0("n=", isolates)) %pm>%
        pm_ungroup() %pm>%
        pm_select(-gr)
    }
  )
}

expand_SIR_colours <- function(colours_SIR, unname = TRUE) {
  sir_order <- c("S", "SDD", "I", "R")

  if (is.null(names(colours_SIR))) {
    if (length(colours_SIR) == 1) {
      colours_SIR <- rep(colours_SIR, 4)
    } else if (length(colours_SIR) == 3) {
      # old method for AMR < 3.0.1 which allowed for 3 colours
      # fill in green for SDD as extra colour
      colours_SIR <- c(colours_SIR[1], colours_SIR[1], colours_SIR[2], colours_SIR[3])
    }
    names(colours_SIR) <- sir_order
  } else {
    # named input: match and reorder
    stop_ifnot(
      all(names(colours_SIR) %in% sir_order),
      "Unknown names in `colours_SIR`. Expected any of: ", vector_or(sir_order, quotes = FALSE, sort = FALSE), "."
    )
    colours_SIR <- colours_SIR[sir_order]
  }

  if (unname) {
    colours_SIR <- unname(colours_SIR)
  }

  return(colours_SIR)
}
