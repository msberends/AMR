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

#' AMR Plots with `ggplot2`
#'
#' Use these functions to create bar plots for AMR data analysis. All functions rely on [ggplot2][ggplot2::ggplot()] functions.
#' @param data a [data.frame] with column(s) of class [`sir`] (see [as.sir()])
#' @param position position adjustment of bars, either `"fill"`, `"stack"` or `"dodge"`
#' @param x variable to show on x axis, either `"antibiotic"` (default) or `"interpretation"` or a grouping variable
#' @param fill variable to categorise using the plots legend, either `"antibiotic"` (default) or `"interpretation"` or a grouping variable
#' @param breaks a [numeric] vector of positions
#' @param limits a [numeric] vector of length two providing limits of the scale, use `NA` to refer to the existing minimum or maximum
#' @param facet variable to split plots by, either `"interpretation"` (default) or `"antibiotic"` or a grouping variable
#' @inheritParams proportion
#' @param nrow (when using `facet`) number of rows
#' @param colours a named vactor with colour to be used for filling. The default colours are colour-blind friendly.
#' @param datalabels show datalabels using [labels_sir_count()]
#' @param datalabels.size size of the datalabels
#' @param datalabels.colour colour of the datalabels
#' @param title text to show as title of the plot
#' @param subtitle text to show as subtitle of the plot
#' @param caption text to show as caption of the plot
#' @param x.title text to show as x axis description
#' @param y.title text to show as y axis description
#' @param ... other arguments passed on to [geom_sir()] or, in case of [scale_sir_colours()], named values to set colours. The default colours are colour-blind friendly, while maintaining the convention that e.g. 'susceptible' should be green and 'resistant' should be red. See *Examples*.
#' @details At default, the names of antibiotics will be shown on the plots using [ab_name()]. This can be set with the `translate_ab` argument. See [count_df()].
#'
#' [geom_sir()] will take any variable from the data that has an [`sir`] class (created with [as.sir()]) using [sir_df()] and will plot bars with the percentage S, I, and R. The default behaviour is to have the bars stacked and to have the different antibiotics on the x axis.
#' 
#' Additional functions include:
#'
#' * [facet_sir()] creates 2d plots (at default based on S/I/R) using [ggplot2::facet_wrap()].
#' * [scale_y_percent()] transforms the y axis to a 0 to 100% range using [ggplot2::scale_y_continuous()].
#' * [scale_sir_colours()] sets colours to the bars (green for S, yellow for I, and red for R). with multilingual support. The default colours are colour-blind friendly, while maintaining the convention that e.g. 'susceptible' should be green and 'resistant' should be red.
#' * [theme_sir()] is a [ggplot2 theme][[ggplot2::theme()] with minimal distraction.
#' * [labels_sir_count()] print datalabels on the bars with percentage and amount of isolates using [ggplot2::geom_text()].
#'
#' [ggplot_sir()] is a wrapper around all above functions that uses data as first input. This makes it possible to use this function after a pipe (`%>%`). See *Examples*.
#' @rdname ggplot_sir
#' @export
#' @examples
#' \donttest{
#' if (require("ggplot2") && require("dplyr")) {
#'   # get antimicrobial results for drugs against a UTI:
#'   ggplot(example_isolates %>% select(AMX, NIT, FOS, TMP, CIP)) +
#'     geom_sir()
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'   # prettify the plot using some additional functions:
#'   df <- example_isolates %>% select(AMX, NIT, FOS, TMP, CIP)
#'   ggplot(df) +
#'     geom_sir() +
#'     scale_y_percent() +
#'     scale_sir_colours() +
#'     labels_sir_count() +
#'     theme_sir()
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'   # or better yet, simplify this using the wrapper function - a single command:
#'   example_isolates %>%
#'     select(AMX, NIT, FOS, TMP, CIP) %>%
#'     ggplot_sir()
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'   # get only proportions and no counts:
#'   example_isolates %>%
#'     select(AMX, NIT, FOS, TMP, CIP) %>%
#'     ggplot_sir(datalabels = FALSE)
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'   # add other ggplot2 arguments as you like:
#'   example_isolates %>%
#'     select(AMX, NIT, FOS, TMP, CIP) %>%
#'     ggplot_sir(
#'       width = 0.5,
#'       colour = "black",
#'       size = 1,
#'       linetype = 2,
#'       alpha = 0.25
#'     )
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'   # you can alter the colours with colour names:
#'   example_isolates %>%
#'     select(AMX) %>%
#'     ggplot_sir(colours = c(SI = "yellow"))
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'   # but you can also use the built-in colour-blind friendly colours for
#'   # your plots, where "S" is green, "I" is yellow and "R" is red:
#'   data.frame(
#'     x = c("Value1", "Value2", "Value3"),
#'     y = c(1, 2, 3),
#'     z = c("Value4", "Value5", "Value6")
#'   ) %>%
#'     ggplot() +
#'     geom_col(aes(x = x, y = y, fill = z)) +
#'     scale_sir_colours(Value4 = "S", Value5 = "I", Value6 = "R")
#' }
#' if (require("ggplot2") && require("dplyr")) {
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
#'     ggplot_sir(x = "age_group")
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'   # a shorter version which also adjusts data label colours:
#'   example_isolates %>%
#'     select(AMX, NIT, FOS, TMP, CIP) %>%
#'     ggplot_sir(colours = FALSE)
#' }
#' if (require("ggplot2") && require("dplyr")) {
#'   # it also supports groups (don't forget to use the group var on `x` or `facet`):
#'   example_isolates %>%
#'     filter(mo_is_gram_negative(), ward != "Outpatient") %>%
#'     # select only UTI-specific drugs
#'     select(ward, AMX, NIT, FOS, TMP, CIP) %>%
#'     group_by(ward) %>%
#'     ggplot_sir(
#'       x = "ward",
#'       facet = "antibiotic",
#'       nrow = 1,
#'       title = "AMR of Anti-UTI Drugs Per Ward",
#'       x.title = "Ward",
#'       datalabels = FALSE
#'     )
#' }
#' }
ggplot_sir <- function(data,
                       position = NULL,
                       x = "antibiotic",
                       fill = "interpretation",
                       # params = list(),
                       facet = NULL,
                       breaks = seq(0, 1, 0.1),
                       limits = NULL,
                       translate_ab = "name",
                       combine_SI = TRUE,
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
  meet_criteria(data, allow_class = "data.frame")
  data <- ascertain_sir_classes(data, "data")
  meet_criteria(position, allow_class = "character", has_length = 1, is_in = c("fill", "stack", "dodge"), allow_NULL = TRUE)
  meet_criteria(x, allow_class = "character", has_length = 1)
  meet_criteria(fill, allow_class = "character", has_length = 1)
  meet_criteria(facet, allow_class = "character", has_length = 1, allow_NULL = TRUE)
  meet_criteria(breaks, allow_class = c("numeric", "integer"))
  meet_criteria(limits, allow_class = c("numeric", "integer"), has_length = 2, allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(translate_ab, allow_class = c("character", "logical"), has_length = 1, allow_NA = TRUE)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
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
    geom_sir(
      position = position, x = x, fill = fill, translate_ab = translate_ab,
      minimum = minimum, language = language,
      combine_SI = combine_SI, ...
    ) +
    theme_sir()

  if (fill == "interpretation") {
    p <- p + scale_sir_colours(colours = colours)
  }

  if (identical(position, "fill")) {
    # proportions, so use y scale with percentage
    p <- p + scale_y_percent(breaks = breaks, limits = limits)
  }

  if (datalabels == TRUE) {
    p <- p + labels_sir_count(
      position = position,
      x = x,
      translate_ab = translate_ab,
      minimum = minimum,
      language = language,
      combine_SI = combine_SI,
      datalabels.size = datalabels.size,
      datalabels.colour = datalabels.colour
    )
  }

  if (!is.null(facet)) {
    p <- p + facet_sir(facet = facet, nrow = nrow)
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

#' @rdname ggplot_sir
#' @export
geom_sir <- function(position = NULL,
                     x = c("antibiotic", "interpretation"),
                     fill = "interpretation",
                     translate_ab = "name",
                     minimum = 30,
                     language = get_AMR_locale(),
                     combine_SI = TRUE,
                     ...) {
  x <- x[1]
  stop_ifnot_installed("ggplot2")
  stop_if(is.data.frame(position), "`position` is invalid. Did you accidentally use '%>%' instead of '+'?")
  meet_criteria(position, allow_class = "character", has_length = 1, is_in = c("fill", "stack", "dodge"), allow_NULL = TRUE)
  meet_criteria(x, allow_class = "character", has_length = 1)
  meet_criteria(fill, allow_class = "character", has_length = 1)
  meet_criteria(translate_ab, allow_class = c("character", "logical"), has_length = 1, allow_NA = TRUE)
  meet_criteria(minimum, allow_class = c("numeric", "integer"), has_length = 1, is_positive_or_zero = TRUE, is_finite = TRUE)
  language <- validate_language(language)
  meet_criteria(combine_SI, allow_class = "logical", has_length = 1)

  y <- "value"
  if (missing(position) || is.null(position)) {
    position <- "fill"
  }

  if (identical(position, "fill")) {
    position <- ggplot2::position_fill(vjust = 0.5, reverse = TRUE)
  }
  
  x_deparse <- deparse(substitute(x))
  if (x_deparse != "x") {
    x <- x_deparse
  }
  if (x %like% '".*"') {
    x <- substr(x, 2, nchar(x) - 1)
  }

  if (tolower(x) %in% tolower(c("ab", "abx", "antibiotics"))) {
    x <- "antibiotic"
  } else if (tolower(x) %in% tolower(c("SIR", "sir", "interpretations", "result"))) {
    x <- "interpretation"
  }
  
  ggplot2::geom_col(
    data = function(x) {
      sir_df(
        data = x,
        translate_ab = translate_ab,
        language = language,
        minimum = minimum,
        combine_SI = combine_SI
      )
    },
    mapping = utils::modifyList(ggplot2::aes(), list(x = str2lang(x), y = str2lang(y), fill = str2lang(fill))),
    position = position,
    ...
  )
}
