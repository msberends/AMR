# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

#' AMR plots with \code{ggplot2}
#'
#' Use these functions to create bar plots for antimicrobial resistance analysis. All functions rely on internal \code{\link[ggplot2]{ggplot}2} functions.
#' @param data a \code{data.frame} with column(s) of class \code{"rsi"} (see \code{\link{as.rsi}})
#' @param position position adjustment of bars, either \code{"fill"} (default when \code{fun} is \code{\link{count_df}}), \code{"stack"} (default when \code{fun} is \code{\link{portion_df}}) or \code{"dodge"}
#' @param x variable to show on x axis, either \code{"Antibiotic"} (default) or \code{"Interpretation"} or a grouping variable
#' @param fill variable to categorise using the plots legend, either \code{"Antibiotic"} (default) or \code{"Interpretation"} or a grouping variable
#' @param breaks numeric vector of positions
#' @param limits numeric vector of length two providing limits of the scale, use \code{NA} to refer to the existing minimum or maximum
#' @param facet variable to split plots by, either \code{"Interpretation"} (default) or \code{"Antibiotic"} or a grouping variable
#' @param fun function to transform \code{data}, either \code{\link{count_df}} (default) or \code{\link{portion_df}}
#' @inheritParams portion
#' @param nrow (when using \code{facet}) number of rows
#' @param colours a named vector with colours for the bars. The names must be one or more of: S, SI, I, IR, R or be \code{FALSE} to use default \code{ggplot2} colours.
#' @param datalabels show datalabels using \code{labels_rsi_count}, will only be shown when \code{fun = count_df}
#' @param datalabels.size size of the datalabels
#' @param datalabels.colour colour of the datalabels
#' @param title text to show as title of the plot
#' @param subtitle text to show as subtitle of the plot
#' @param caption text to show as caption of the plot
#' @param x.title text to show as x axis description
#' @param y.title text to show as y axis description
#' @param ... other parameters passed on to \code{geom_rsi}
#' @details At default, the names of antibiotics will be shown on the plots using \code{\link{ab_name}}. This can be set with the \code{translate_ab} parameter. See \code{\link{count_df}}.
#'
#' \strong{The functions}\cr
#' \code{geom_rsi} will take any variable from the data that has an \code{rsi} class (created with \code{\link{as.rsi}}) using \code{fun} (\code{\link{count_df}} at default, can also be \code{\link{portion_df}}) and will plot bars with the percentage R, I and S. The default behaviour is to have the bars stacked and to have the different antibiotics on the x axis.
#'
#' \code{facet_rsi} creates 2d plots (at default based on S/I/R) using \code{\link[ggplot2]{facet_wrap}}.
#'
#' \code{scale_y_percent} transforms the y axis to a 0 to 100\% range using \code{\link[ggplot2]{scale_continuous}}.
#'
#' \code{scale_rsi_colours} sets colours to the bars: pastel blue for S, pastel turquoise for I and pastel red for R, using \code{\link[ggplot2]{scale_brewer}}.
#'
#' \code{theme_rsi} is a \code{ggplot \link[ggplot2]{theme}} with minimal distraction.
#'
#' \code{labels_rsi_count} print datalabels on the bars with percentage and amount of isolates using \code{\link[ggplot2]{geom_text}}
#'
#' \code{ggplot_rsi} is a wrapper around all above functions that uses data as first input. This makes it possible to use this function after a pipe (\code{\%>\%}). See Examples.
#' @rdname ggplot_rsi
#' @importFrom utils installed.packages
#' @export
#' @inheritSection AMR Read more on our website!
#' @examples
#' library(dplyr)
#' library(ggplot2)
#'
#' # get antimicrobial results for drugs against a UTI:
#' ggplot(septic_patients %>% select(AMX, NIT, FOS, TMP, CIP)) +
#'   geom_rsi()
#'
#' # prettify the plot using some additional functions:
#' df <- septic_patients %>% select(AMX, NIT, FOS, TMP, CIP)
#' ggplot(df) +
#'   geom_rsi() +
#'   scale_y_percent() +
#'   scale_rsi_colours() +
#'   labels_rsi_count() +
#'   theme_rsi()
#'
#' # or better yet, simplify this using the wrapper function - a single command:
#' septic_patients %>%
#'   select(AMX, NIT, FOS, TMP, CIP) %>%
#'   ggplot_rsi()
#'
#' # get only portions and no counts:
#' septic_patients %>%
#'   select(AMX, NIT, FOS, TMP, CIP) %>%
#'   ggplot_rsi(fun = portion_df)
#'
#' # add other ggplot2 parameters as you like:
#' septic_patients %>%
#'   select(AMX, NIT, FOS, TMP, CIP) %>%
#'   ggplot_rsi(width = 0.5,
#'              colour = "black",
#'              size = 1,
#'              linetype = 2,
#'              alpha = 0.25)
#'
#' septic_patients %>%
#'   select(AMX) %>%
#'   ggplot_rsi(colours = c(SI = "yellow"))
#'
#' # resistance of ciprofloxacine per age group
#' septic_patients %>%
#'   mutate(first_isolate = first_isolate(.)) %>%
#'   filter(first_isolate == TRUE,
#'          mo == as.mo("E. coli")) %>%
#'   # `age_group` is also a function of this package:
#'   group_by(age_group = age_groups(age)) %>%
#'   select(age_group,
#'          CIP) %>%
#'   ggplot_rsi(x = "age_group")
#' \donttest{
#'
#' # for colourblind mode, use divergent colours from the viridis package:
#' septic_patients %>%
#'   select(AMX, NIT, FOS, TMP, CIP) %>%
#'   ggplot_rsi() + scale_fill_viridis_d()
#' # a shorter version which also adjusts data label colours:
#' septic_patients %>%
#'   select(AMX, NIT, FOS, TMP, CIP) %>%
#'   ggplot_rsi(colours = FALSE)
#'
#'
#' # it also supports groups (don't forget to use the group var on `x` or `facet`):
#' septic_patients %>%
#'   select(hospital_id, AMX, NIT, FOS, TMP, CIP) %>%
#'   group_by(hospital_id) %>%
#'   ggplot_rsi(x = "hospital_id",
#'              facet = "Antibiotic",
#'              nrow = 1,
#'              title = "AMR of Anti-UTI Drugs Per Hospital",
#'              x.title = "Hospital",
#'              datalabels = FALSE)
#'
#' # genuine analysis: check 3 most prevalent microorganisms
#' septic_patients %>%
#'   # create new bacterial ID's, with all CoNS under the same group (Becker et al.)
#'   mutate(mo = as.mo(mo, Becker = TRUE)) %>%
#'   # filter on top three bacterial ID's
#'   filter(mo %in% top_freq(freq(.$mo), 3)) %>%
#'   # filter on first isolates
#'   filter_first_isolate() %>%
#'   # get short MO names (like "E. coli")
#'   mutate(bug = mo_shortname(mo, Becker = TRUE)) %>%
#'   # select this short name and some antiseptic drugs
#'   select(bug, CXM, GEN, CIP) %>%
#'   # group by MO
#'   group_by(bug) %>%
#'   # plot the thing, putting MOs on the facet
#'   ggplot_rsi(x = "Antibiotic",
#'              facet = "bug",
#'              translate_ab = FALSE,
#'              nrow = 1,
#'              title = "AMR of Top Three Microorganisms In Blood Culture Isolates",
#'              subtitle = expression(paste("Only First Isolates, CoNS grouped according to Becker ",
#'                                          italic("et al."), " (2014)")),
#'              x.title = "Antibiotic (EARS-Net code)")
#' }
ggplot_rsi <- function(data,
                       position = NULL,
                       x = "Antibiotic",
                       fill = "Interpretation",
                       # params = list(),
                       facet = NULL,
                       breaks = seq(0, 1, 0.1),
                       limits = NULL,
                       translate_ab = "name",
                       combine_SI = TRUE,
                       combine_IR = FALSE,
                       language = get_locale(),
                       fun = count_df,
                       nrow = NULL,
                       colours = c(S = "#61a8ff",
                                   SI = "#61a8ff",
                                   I = "#61f7ff",
                                   IR = "#ff6961",
                                   R = "#ff6961"),
                       datalabels = TRUE,
                       datalabels.size = 2.5,
                       datalabels.colour = "gray15",
                       title = NULL,
                       subtitle = NULL,
                       caption = NULL,
                       x.title = NULL,
                       y.title = NULL,
                       ...) {

  stopifnot_installed_package("ggplot2")

  fun_name <- deparse(substitute(fun))
  if (!fun_name %in% c("portion_df", "count_df")) {
    stop("`fun` must be portion_df or count_df")
  }

  x <- x[1]
  facet <- facet[1]

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
    geom_rsi(position = position, x = x, fill = fill, translate_ab = translate_ab,
             fun = fun, combine_SI = combine_SI, combine_IR = combine_IR, ...) +
    theme_rsi()

  if (fill == "Interpretation") {
    # set RSI colours
    if (isFALSE(colours) & missing(datalabels.colour)) {
      # set datalabel colour to middle gray
      datalabels.colour <- "gray50"
    }
    p <- p + scale_rsi_colours(colours = colours)
  }

  if (fun_name == "portion_df"
      | (fun_name == "count_df" & identical(position, "fill"))) {
    # portions, so use y scale with percentage
    p <- p + scale_y_percent(breaks = breaks, limits = limits)
  }

  if (fun_name == "count_df" & datalabels == TRUE) {
    p <- p + labels_rsi_count(position = position,
                              x = x,
                              translate_ab = translate_ab,
                              combine_SI = combine_SI,
                              combine_IR = combine_IR,
                              datalabels.size = datalabels.size,
                              datalabels.colour = datalabels.colour)
  }

  if (!is.null(facet)) {
    p <- p + facet_rsi(facet = facet, nrow = nrow)
  }

  p <- p + ggplot2::labs(title = title,
                         subtitle = subtitle,
                         caption = caption,
                         x = x.title,
                         y = y.title)

  p
}

#' @rdname ggplot_rsi
#' @export
geom_rsi <- function(position = NULL,
                     x = c("Antibiotic", "Interpretation"),
                     fill = "Interpretation",
                     translate_ab = "name",
                     language = get_locale(),
                     combine_SI = TRUE,
                     combine_IR = FALSE,
                     fun = count_df,
                     ...)  {

  stopifnot_installed_package("ggplot2")

  if (is.data.frame(position)) {
    stop("`position` is invalid. Did you accidentally use '%>%' instead of '+'?", call. = FALSE)
  }

  fun_name <- deparse(substitute(fun))
  if (!fun_name %in% c("portion_df", "count_df", "fun")) {
    stop("`fun` must be portion_df or count_df")
  }
  y <- "Value"
  if (identical(fun, count_df)) {
    if (missing(position) | is.null(position)) {
      position <- "fill"
    }
  } else {
    if (missing(position) | is.null(position)) {
      position <- "stack"
    }
  }

  if (identical(position, "fill")) {
    position <- ggplot2::position_fill(vjust = 0.5, reverse = TRUE)
  }

  x <- x[1]

  # we work with aes_string later on
  x_deparse <- deparse(substitute(x))
  if (x_deparse != "x") {
    x <- x_deparse
  }
  if (x %like% '".*"') {
    x <- substr(x, 2, nchar(x) - 1)
  }

  if (tolower(x) %in% tolower(c('ab', 'antibiotic', 'abx', 'antibiotics'))) {
    x <- "Antibiotic"
  } else if (tolower(x) %in% tolower(c('SIR', 'RSI', 'interpretation', 'interpretations', 'result'))) {
    x <- "Interpretation"
  }

  ggplot2::layer(geom = "bar", stat = "identity", position = position,
                 mapping = ggplot2::aes_string(x = x, y = y, fill = fill),
                 params = list(...), data = function(x) {
                   fun(data = x,
                       translate_ab = translate_ab,
                       language = language,
                       combine_SI = combine_SI,
                       combine_IR = combine_IR)
                 })

}

#' @rdname ggplot_rsi
#' @export
facet_rsi <- function(facet = c("Interpretation", "Antibiotic"), nrow = NULL) {

  stopifnot_installed_package("ggplot2")

  facet <- facet[1]

  # we work with aes_string later on
  facet_deparse <- deparse(substitute(facet))
  if (facet_deparse != "facet") {
    facet <- facet_deparse
  }
  if (facet %like% '".*"') {
    facet <- substr(facet, 2, nchar(facet) - 1)
  }

  if (tolower(facet) %in% tolower(c('SIR', 'RSI', 'interpretation', 'interpretations', 'result'))) {
    facet <- "Interpretation"
  } else if (tolower(facet) %in% tolower(c('ab', 'antibiotic', 'abx', 'antibiotics'))) {
    facet <- "Antibiotic"
  }

  ggplot2::facet_wrap(facets = facet, scales = "free_x", nrow = nrow)
}

#' @rdname ggplot_rsi
#' @export
scale_y_percent <- function(breaks = seq(0, 1, 0.1), limits = NULL) {
  stopifnot_installed_package("ggplot2")

  if (all(breaks[breaks != 0] > 1)) {
    breaks <- breaks / 100
  }
  ggplot2::scale_y_continuous(breaks = breaks,
                              labels = percent(breaks),
                              limits = limits)
}

#' @rdname ggplot_rsi
#' @export
scale_rsi_colours <- function(colours = c(S = "#61a8ff",
                                          SI = "#61a8ff",
                                          I = "#61f7ff",
                                          IR = "#ff6961",
                                          R = "#ff6961")) {
  stopifnot_installed_package("ggplot2")
  #ggplot2::scale_fill_brewer(palette = "RdYlGn")
  #ggplot2::scale_fill_manual(values = c("#b22222", "#ae9c20", "#7cfc00"))

  if (!identical(colours, FALSE)) {
    original_cols <- c(S = "#61a8ff",
                       SI = "#61a8ff",
                       I = "#61f7ff",
                       IR = "#ff6961",
                       R = "#ff6961")
    colours <- replace(original_cols, names(colours), colours)
    ggplot2::scale_fill_manual(values = colours)
  }
}

#' @rdname ggplot_rsi
#' @export
theme_rsi <- function() {
  stopifnot_installed_package("ggplot2")
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey75"),
                   # center title and subtitle
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5))
}

#' @rdname ggplot_rsi
#' @importFrom dplyr mutate %>% group_by_at
#' @export
labels_rsi_count <- function(position = NULL,
                             x = "Antibiotic",
                             translate_ab = "name",
                             combine_SI = TRUE,
                             combine_IR = FALSE,
                             datalabels.size = 3,
                             datalabels.colour = "gray15") {
  stopifnot_installed_package("ggplot2")
  if (is.null(position)) {
    position <- "fill"
  }
  if (identical(position, "fill")) {
    position <- ggplot2::position_fill(vjust = 0.5, reverse = TRUE)
  }
  x_name <- x
  ggplot2::geom_text(mapping = ggplot2::aes_string(label = "lbl",
                                                   x = x,
                                                   y = "Value"),
                     position = position,
                     inherit.aes = FALSE,
                     size = datalabels.size,
                     colour = datalabels.colour,
                     lineheight = 0.75,
                     data = function(x) {
                       # labels are only shown when function is count_df,
                       # so no need parameterise it here
                       count_df(data = x,
                                translate_ab = translate_ab,
                                combine_SI = combine_SI,
                                combine_IR = combine_IR) %>%
                         group_by_at(x_name) %>%
                         mutate(lbl = paste0(percent(Value / sum(Value, na.rm = TRUE), force_zero = TRUE),
                                             "\n(n=", Value, ")"))
                     })
}
