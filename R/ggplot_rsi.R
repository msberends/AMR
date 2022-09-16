# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       #
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

#' AMR Plots with `ggplot2`
#'
#' Use these functions to create bar plots for AMR data analysis. All functions rely on [ggplot2][ggplot2::ggplot()] functions.
#' @param data a [data.frame] with column(s) of class [`rsi`] (see [as.rsi()])
#' @param position position adjustment of bars, either `"fill"`, `"stack"` or `"dodge"`
#' @param x variable to show on x axis, either `"antibiotic"` (default) or `"interpretation"` or a grouping variable
#' @param fill variable to categorise using the plots legend, either `"antibiotic"` (default) or `"interpretation"` or a grouping variable
#' @param breaks a [numeric] vector of positions
#' @param limits a [numeric] vector of length two providing limits of the scale, use `NA` to refer to the existing minimum or maximum
#' @param facet variable to split plots by, either `"interpretation"` (default) or `"antibiotic"` or a grouping variable
#' @inheritParams proportion
#' @param nrow (when using `facet`) number of rows
#' @param colours a named vactor with colour to be used for filling. The default colours are colour-blind friendly.
#' @param aesthetics aesthetics to apply the colours to, defaults to "fill" but can also be (a combination of) "alpha", "colour", "fill", "linetype", "shape" or "size"
#' @param datalabels show datalabels using [labels_rsi_count()]
#' @param datalabels.size size of the datalabels
#' @param datalabels.colour colour of the datalabels
#' @param title text to show as title of the plot
#' @param subtitle text to show as subtitle of the plot
#' @param caption text to show as caption of the plot
#' @param x.title text to show as x axis description
#' @param y.title text to show as y axis description
#' @param ... other arguments passed on to [geom_rsi()] or, in case of [scale_rsi_colours()], named values to set colours. The default colours are colour-blind friendly, while maintaining the convention that e.g. 'susceptible' should be green and 'resistant' should be red. See *Examples*.
#' @details At default, the names of antibiotics will be shown on the plots using [ab_name()]. This can be set with the `translate_ab` argument. See [count_df()].
#'
#' ## The Functions
#' [geom_rsi()] will take any variable from the data that has an [`rsi`] class (created with [as.rsi()]) using [rsi_df()] and will plot bars with the percentage R, I and S. The default behaviour is to have the bars stacked and to have the different antibiotics on the x axis.
#'
#' [facet_rsi()] creates 2d plots (at default based on S/I/R) using [ggplot2::facet_wrap()].
#'
#' [scale_y_percent()] transforms the y axis to a 0 to 100% range using [ggplot2::scale_y_continuous()].
#'
#' [scale_rsi_colours()] sets colours to the bars (green for S, yellow for I, and red for R). with multilingual support. The default colours are colour-blind friendly, while maintaining the convention that e.g. 'susceptible' should be green and 'resistant' should be red.
#'
#' [theme_rsi()] is a [ggplot2 theme][[ggplot2::theme()] with minimal distraction.
#'
#' [labels_rsi_count()] print datalabels on the bars with percentage and amount of isolates using [ggplot2::geom_text()].
#'
#' [ggplot_rsi()] is a wrapper around all above functions that uses data as first input. This makes it possible to use this function after a pipe (`%>%`). See *Examples*.
#' @rdname ggplot_rsi
#' @export
#' @examples
#' \donttest{
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # get antimicrobial results for drugs against a UTI:
#'   ggplot(example_isolates %>% select(AMX, NIT, FOS, TMP, CIP)) +
#'     geom_rsi()
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # prettify the plot using some additional functions:
#'   df <- example_isolates %>% select(AMX, NIT, FOS, TMP, CIP)
#'   ggplot(df) +
#'     geom_rsi() +
#'     scale_y_percent() +
#'     scale_rsi_colours() +
#'     labels_rsi_count() +
#'     theme_rsi()
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # or better yet, simplify this using the wrapper function - a single command:
#'   example_isolates %>%
#'     select(AMX, NIT, FOS, TMP, CIP) %>%
#'     ggplot_rsi()
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # get only proportions and no counts:
#'   example_isolates %>%
#'     select(AMX, NIT, FOS, TMP, CIP) %>%
#'     ggplot_rsi(datalabels = FALSE)
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # add other ggplot2 arguments as you like:
#'   example_isolates %>%
#'     select(AMX, NIT, FOS, TMP, CIP) %>%
#'     ggplot_rsi(
#'       width = 0.5,
#'       colour = "black",
#'       size = 1,
#'       linetype = 2,
#'       alpha = 0.25
#'     )
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # you can alter the colours with colour names:
#'   example_isolates %>%
#'     select(AMX) %>%
#'     ggplot_rsi(colours = c(SI = "yellow"))
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # but you can also use the built-in colour-blind friendly colours for
#'   # your plots, where "S" is green, "I" is yellow and "R" is red:
#'   data.frame(
#'     x = c("Value1", "Value2", "Value3"),
#'     y = c(1, 2, 3),
#'     z = c("Value4", "Value5", "Value6")
#'   ) %>%
#'     ggplot() +
#'     geom_col(aes(x = x, y = y, fill = z)) +
#'     scale_rsi_colours(Value4 = "S", Value5 = "I", Value6 = "R")
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # resistance of ciprofloxacine per age group
#'   example_isolates %>%
#'     mutate(first_isolate = first_isolate()) %>%
#'     filter(
#'       first_isolate == TRUE,
#'       mo == as.mo("Escherichia coli")
#'     ) %>%
#'     # age_groups() is also a function in this AMR package:
#'     group_by(age_group = age_groups(age)) %>%
#'     select(age_group, CIP) %>%
#'     ggplot_rsi(x = "age_group")
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # a shorter version which also adjusts data label colours:
#'   example_isolates %>%
#'     select(AMX, NIT, FOS, TMP, CIP) %>%
#'     ggplot_rsi(colours = FALSE)
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'
#'   # it also supports groups (don't forget to use the group var on `x` or `facet`):
#'   example_isolates %>%
#'     filter(mo_is_gram_negative(), ward != "Outpatient") %>%
#'     # select only UTI-specific drugs
#'     select(ward, AMX, NIT, FOS, TMP, CIP) %>%
#'     group_by(ward) %>%
#'     ggplot_rsi(
#'       x = "ward",
#'       facet = "antibiotic",
#'       nrow = 1,
#'       title = "AMR of Anti-UTI Drugs Per Ward",
#'       x.title = "Ward",
#'       datalabels = FALSE
#'     )
#' }
#' }
ggplot_rsi <- function(data,
                       position = NULL,
                       x = "antibiotic",
                       fill = "interpretation",
                       # params = list(),
                       facet = NULL,
                       breaks = seq(0, 1, 0.1),
                       limits = NULL,
                       translate_ab = "name",
                       combine_SI = TRUE,
                       combine_IR = FALSE,
                       minimum = 30,
                       language = get_AMR_locale(),
                       nrow = NULL,
                       colours = c(
                         S = "#3CAEA3",
                         SI = "#3CAEA3",
                         I = "#F6D55C",
                         IR = "#ED553B",
                         R = "#ED553B"
                       ),
                       datalabels = TRUE,
                       datalabels.size = 2.5,
                       datalabels.colour = "grey15",
                       title = NULL,
                       subtitle = NULL,
                       caption = NULL,
                       x.title = "Antimicrobial",
                       y.title = "Proportion",
                       ...) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(data, allow_class = "data.frame", contains_column_class = "rsi")
  meet_criteria(position, allow_class = "character", has_length = 1, is_in = c("fill", "stack", "dodge"), allow_NULL = TRUE)
  meet_criteria(x, allow_class = "character", has_length = 1)
  meet_criteria(fill, allow_class = "character", has_length = 1)
  meet_criteria(facet, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(breaks, allow_class = c("numeric", "integer"))
  meet_criteria(limits, allow_class = c("numeric", "integer"), has_length = 2, allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(translate_ab, allow_class = c("character", "logical"), has_length = 1, allow_NA = TRUE)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(combine_IR, allow_class = "logical", has_length = 1)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE)
  language <- validate_language(language)
  meet_criteria(nrow, allow_class = c("numeric", "integer"), has_length = 1, allow_NULL = TRUE, is_positive = TRUE, is_finite = TRUE)
  meet_criteria(colours, allow_class = c("character", "logical"))
  meet_criteria(datalabels, allow_class = "logical", has_length = 1)
  meet_criteria(datalabels.size, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = TRUE)
  meet_criteria(datalabels.colour, allow_class = "character", has_length = 1)
  meet_criteria(title, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(subtitle, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(caption, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(x.title, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(y.title, allow_class = "character", has_length = 1, allow_NULL = TRUE)

  # we work with aes_string later on
  x_deparse <- deparse(substitute(x))
  if (x_deparse != "x") {
    x <- x_deparse
  }
  if (x %like% '".*"') {
    x <- substr(x, 2, nchar(x) - 1)
  }
  facet_deparse <- deparse(substitute(facet))
  if (facet_deparse != "facet") {
    facet <- facet_deparse
  }
  if (facet %like% '".*"') {
    facet <- substr(facet, 2, nchar(facet) - 1)
  }
  if (facet %in% c("NULL", "")) {
    facet <- NULL
  }

  if (is.null(position)) {
    position <- "fill"
  }

  p <- ggplot2::ggplot(data = data) +
    geom_rsi(
      position = position, x = x, fill = fill, translate_ab = translate_ab,
      minimum = minimum, language = language,
      combine_SI = combine_SI, combine_IR = combine_IR, ...
    ) +
    theme_rsi()

  if (fill == "interpretation") {
    p <- p + scale_rsi_colours(colours = colours)
  }

  if (identical(position, "fill")) {
    # proportions, so use y scale with percentage
    p <- p + scale_y_percent(breaks = breaks, limits = limits)
  }

  if (datalabels == TRUE) {
    p <- p + labels_rsi_count(
      position = position,
      x = x,
      translate_ab = translate_ab,
      minimum = minimum,
      language = language,
      combine_SI = combine_SI,
      combine_IR = combine_IR,
      datalabels.size = datalabels.size,
      datalabels.colour = datalabels.colour
    )
  }

  if (!is.null(facet)) {
    p <- p + facet_rsi(facet = facet, nrow = nrow)
  }

  p <- p + ggplot2::labs(
    title = title,
    subtitle = subtitle,
    caption = caption,
    x = x.title,
    y = y.title
  )

  p
}

#' @rdname ggplot_rsi
#' @export
geom_rsi <- function(position = NULL,
                     x = c("antibiotic", "interpretation"),
                     fill = "interpretation",
                     translate_ab = "name",
                     minimum = 30,
                     language = get_AMR_locale(),
                     combine_SI = TRUE,
                     combine_IR = FALSE,
                     ...) {
  x <- x[1]
  stop_ifnot_installed("ggplot2")
  stop_if(is.data.frame(position), "`position` is invalid. Did you accidentally use '%>%' instead of '+'?")
  meet_criteria(position, allow_class = "character", has_length = 1, is_in = c("fill", "stack", "dodge"), allow_NULL = TRUE)
  meet_criteria(x, allow_class = "character", has_length = 1)
  meet_criteria(fill, allow_class = "character", has_length = 1)
  meet_criteria(translate_ab, allow_class = c("character", "logical"), has_length = 1, allow_NA = TRUE)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE)
  language <- validate_language(language)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(combine_IR, allow_class = "logical", has_length = 1)

  y <- "value"
  if (missing(position) || is.null(position)) {
    position <- "fill"
  }

  if (identical(position, "fill")) {
    position <- ggplot2::position_fill(vjust = 0.5, reverse = TRUE)
  }

  # we work with aes_string later on
  x_deparse <- deparse(substitute(x))
  if (x_deparse != "x") {
    x <- x_deparse
  }
  if (x %like% '".*"') {
    x <- substr(x, 2, nchar(x) - 1)
  }

  if (tolower(x) %in% tolower(c("ab", "abx", "antibiotics"))) {
    x <- "antibiotic"
  } else if (tolower(x) %in% tolower(c("SIR", "RSI", "interpretations", "result"))) {
    x <- "interpretation"
  }

  ggplot2::geom_col(
    data = function(x) {
      rsi_df(
        data = x,
        translate_ab = translate_ab,
        language = language,
        minimum = minimum,
        combine_SI = combine_SI,
        combine_IR = combine_IR
      )
    },
    mapping = ggplot2::aes_string(x = x, y = y, fill = fill),
    position = position,
    ...
  )
}

#' @rdname ggplot_rsi
#' @export
facet_rsi <- function(facet = c("interpretation", "antibiotic"), nrow = NULL) {
  facet <- facet[1]
  stop_ifnot_installed("ggplot2")
  meet_criteria(facet, allow_class = "character", has_length = 1)
  meet_criteria(nrow, allow_class = c("numeric", "integer"), has_length = 1, allow_NULL = TRUE, is_positive = TRUE, is_finite = TRUE)

  # we work with aes_string later on
  facet_deparse <- deparse(substitute(facet))
  if (facet_deparse != "facet") {
    facet <- facet_deparse
  }
  if (facet %like% '".*"') {
    facet <- substr(facet, 2, nchar(facet) - 1)
  }

  if (tolower(facet) %in% tolower(c("SIR", "RSI", "interpretations", "result"))) {
    facet <- "interpretation"
  } else if (tolower(facet) %in% tolower(c("ab", "abx", "antibiotics"))) {
    facet <- "antibiotic"
  }

  ggplot2::facet_wrap(facets = facet, scales = "free_x", nrow = nrow)
}

#' @rdname ggplot_rsi
#' @export
scale_y_percent <- function(breaks = seq(0, 1, 0.1), limits = NULL) {
  stop_ifnot_installed("ggplot2")
  meet_criteria(breaks, allow_class = c("numeric", "integer"))
  meet_criteria(limits, allow_class = c("numeric", "integer"), has_length = 2, allow_NULL = TRUE, allow_NA = TRUE)

  if (all(breaks[breaks != 0] > 1)) {
    breaks <- breaks / 100
  }
  ggplot2::scale_y_continuous(
    breaks = breaks,
    labels = percentage(breaks),
    limits = limits
  )
}

#' @rdname ggplot_rsi
#' @export
scale_rsi_colours <- function(...,
                              aesthetics = "fill") {
  stop_ifnot_installed("ggplot2")
  meet_criteria(aesthetics, allow_class = "character", is_in = c("alpha", "colour", "color", "fill", "linetype", "shape", "size"))
  # behaviour until AMR pkg v1.5.0 and also when coming from ggplot_rsi()
  if ("colours" %in% names(list(...))) {
    original_cols <- c(
      S = "#3CAEA3",
      SI = "#3CAEA3",
      I = "#F6D55C",
      IR = "#ED553B",
      R = "#ED553B"
    )
    colours <- replace(original_cols, names(list(...)$colours), list(...)$colours)
    # limits = force is needed in ggplot2 3.3.4 and 3.3.5, see here;
    # https://github.com/tidyverse/ggplot2/issues/4511#issuecomment-866185530
    return(ggplot2::scale_fill_manual(values = colours, limits = force))
  }
  if (identical(unlist(list(...)), FALSE)) {
    return(invisible())
  }

  names_susceptible <- c(
    "S", "SI", "IS", "S+I", "I+S", "susceptible", "Susceptible",
    unique(TRANSLATIONS[which(TRANSLATIONS$pattern == "Susceptible"),
      "replacement",
      drop = TRUE
    ])
  )
  names_incr_exposure <- c(
    "I", "intermediate", "increased exposure", "incr. exposure",
    "Increased exposure", "Incr. exposure", "Susceptible, incr. exp.",
    unique(TRANSLATIONS[which(TRANSLATIONS$pattern == "Intermediate"),
      "replacement",
      drop = TRUE
    ]),
    unique(TRANSLATIONS[which(TRANSLATIONS$pattern == "Susceptible, incr. exp."),
      "replacement",
      drop = TRUE
    ])
  )
  names_resistant <- c(
    "R", "IR", "RI", "R+I", "I+R", "resistant", "Resistant",
    unique(TRANSLATIONS[which(TRANSLATIONS$pattern == "Resistant"),
      "replacement",
      drop = TRUE
    ])
  )

  susceptible <- rep("#3CAEA3", length(names_susceptible))
  names(susceptible) <- names_susceptible
  incr_exposure <- rep("#F6D55C", length(names_incr_exposure))
  names(incr_exposure) <- names_incr_exposure
  resistant <- rep("#ED553B", length(names_resistant))
  names(resistant) <- names_resistant

  original_cols <- c(susceptible, incr_exposure, resistant)
  dots <- c(...)
  # replace S, I, R as colours: scale_rsi_colours(mydatavalue = "S")
  dots[dots == "S"] <- "#3CAEA3"
  dots[dots == "I"] <- "#F6D55C"
  dots[dots == "R"] <- "#ED553B"
  cols <- replace(original_cols, names(dots), dots)
  # limits = force is needed in ggplot2 3.3.4 and 3.3.5, see here;
  # https://github.com/tidyverse/ggplot2/issues/4511#issuecomment-866185530
  ggplot2::scale_discrete_manual(aesthetics = aesthetics, values = cols, limits = force)
}

#' @rdname ggplot_rsi
#' @export
theme_rsi <- function() {
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

#' @rdname ggplot_rsi
#' @export
labels_rsi_count <- function(position = NULL,
                             x = "antibiotic",
                             translate_ab = "name",
                             minimum = 30,
                             language = get_AMR_locale(),
                             combine_SI = TRUE,
                             combine_IR = FALSE,
                             datalabels.size = 3,
                             datalabels.colour = "grey15") {
  stop_ifnot_installed("ggplot2")
  meet_criteria(position, allow_class = "character", has_length = 1, is_in = c("fill", "stack", "dodge"), allow_NULL = TRUE)
  meet_criteria(x, allow_class = "character", has_length = 1)
  meet_criteria(translate_ab, allow_class = c("character", "logical"), has_length = 1, allow_NA = TRUE)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_finite = TRUE)
  language <- validate_language(language)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(combine_IR, allow_class = "logical", has_length = 1)
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
    mapping = ggplot2::aes_string(
      label = "lbl",
      x = x,
      y = "value"
    ),
    position = position,
    inherit.aes = FALSE,
    size = datalabels.size,
    colour = datalabels.colour,
    lineheight = 0.75,
    data = function(x) {
      transformed <- rsi_df(
        data = x,
        translate_ab = translate_ab,
        combine_SI = combine_SI,
        combine_IR = combine_IR,
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
