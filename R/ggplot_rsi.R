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
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

#' AMR bar plots with \code{ggplot}
#'
#' Use these functions to create bar plots for antimicrobial resistance analysis. All functions rely on internal \code{\link[ggplot2]{ggplot}} functions.
#' @param data a \code{data.frame} with column(s) of class \code{"rsi"} (see \code{\link{as.rsi}})
#' @param position position adjustment of bars, either \code{"fill"} (default when \code{fun} is \code{\link{count_df}}), \code{"stack"} (default when \code{fun} is \code{\link{portion_df}}) or \code{"dodge"}
#' @param x variable to show on x axis, either \code{"Antibiotic"} (default) or \code{"Interpretation"} or a grouping variable
#' @param fill variable to categorise using the plots legend, either \code{"Antibiotic"} (default) or \code{"Interpretation"} or a grouping variable
#' @param breaks numeric vector of positions
#' @param limits numeric vector of length two providing limits of the scale, use \code{NA} to refer to the existing minimum or maximum
#' @param facet variable to split plots by, either \code{"Interpretation"} (default) or \code{"Antibiotic"} or a grouping variable
#' @param translate_ab a column name of the \code{\link{antibiotics}} data set to translate the antibiotic abbreviations into, using \code{\link{abname}}. Default behaviour is to translate to official names according to the WHO. Use \code{translate_ab = FALSE} to disable translation.
#' @param fun function to transform \code{data}, either \code{\link{count_df}} (default) or \code{\link{portion_df}}
#' @param nrow (when using \code{facet}) number of rows
#' @param datalabels show datalabels using \code{labels_rsi_count}, will at default only be shown when \code{fun = count_df}
#' @param datalabels.size size of the datalabels
#' @param datalabels.colour colour of the datalabels
#' @param ... other parameters passed on to \code{geom_rsi}
#' @details At default, the names of antibiotics will be shown on the plots using \code{\link{abname}}. This can be set with the option \code{get_antibiotic_names} (a logical value), so change it e.g. to \code{FALSE} with \code{options(get_antibiotic_names = FALSE)}.
#'
#' \strong{The functions}\cr
#' \code{geom_rsi} will take any variable from the data that has an \code{rsi} class (created with \code{\link{as.rsi}}) using \code{fun} (\code{\link{count_df}} at default, can also be \code{\link{portion_df}}) and will plot bars with the percentage R, I and S. The default behaviour is to have the bars stacked and to have the different antibiotics on the x axis.
#'
#' \code{facet_rsi} creates 2d plots (at default based on S/I/R) using \code{\link[ggplot2]{facet_wrap}}.
#'
#' \code{scale_y_percent} transforms the y axis to a 0 to 100\% range using \code{\link[ggplot2]{scale_continuous}}.
#'
#' \code{scale_rsi_colours} sets colours to the bars: green for S, yellow for I and red for R, using \code{\link[ggplot2]{scale_brewer}}.
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
#' ggplot(septic_patients %>% select(amox, nitr, fosf, trim, cipr)) +
#'   geom_rsi()
#'
#' # prettify the plot using some additional functions:
#' df <- septic_patients[, c("amox", "nitr", "fosf", "trim", "cipr")]
#' ggplot(df) +
#'   geom_rsi() +
#'   scale_y_percent() +
#'   scale_rsi_colours() +
#'   labels_rsi_count() +
#'   theme_rsi()
#'
#' # or better yet, simplify this using the wrapper function - a single command:
#' septic_patients %>%
#'   select(amox, nitr, fosf, trim, cipr) %>%
#'   ggplot_rsi()
#'
#' # get only portions and no counts:
#' septic_patients %>%
#'   select(amox, nitr, fosf, trim, cipr) %>%
#'   ggplot_rsi(fun = portion_df)
#'
#' # add other ggplot2 parameters as you like:
#' septic_patients %>%
#'   select(amox, nitr, fosf, trim, cipr) %>%
#'   ggplot_rsi(width = 0.5,
#'              colour = "black",
#'              size = 1,
#'              linetype = 2,
#'              alpha = 0.25)
#'
#' # resistance of ciprofloxacine per age group
#' septic_patients %>%
#'   mutate(first_isolate = first_isolate(.)) %>%
#'   filter(first_isolate == TRUE,
#'          mo == as.mo("E. coli")) %>%
#'   # `age_group` is also a function of this package:
#'   group_by(age_group = age_groups(age)) %>%
#'   select(age_group,
#'          cipr) %>%
#'   ggplot_rsi(x = "age_group")
#' \donttest{
#'
#' # for colourblind mode, use divergent colours from the viridis package:
#' septic_patients %>%
#'   select(amox, nitr, fosf, trim, cipr) %>%
#'   ggplot_rsi() + scale_fill_viridis_d()
#'
#'
#' # it also supports groups (don't forget to use the group var on `x` or `facet`):
#' septic_patients %>%
#'   select(hospital_id, amox, nitr, fosf, trim, cipr) %>%
#'   group_by(hospital_id) %>%
#'   ggplot_rsi(x = hospital_id,
#'              facet = Antibiotic,
#'              nrow = 1) +
#'   labs(title = "AMR of Anti-UTI Drugs Per Hospital",
#'        x = "Hospital")
#'
#' # genuine analysis: check 2 most prevalent microorganisms
#' septic_patients %>%
#'   # create new bacterial ID's, with all CoNS under the same group (Becker et al.)
#'   mutate(mo = as.mo(mo, Becker = TRUE)) %>%
#'   # filter on top three bacterial ID's
#'   filter(mo %in% top_freq(freq(.$mo), 3)) %>%
#'   # determine first isolates
#'   mutate(first_isolate = first_isolate(.,
#'                                        col_date = "date",
#'                                        col_patient_id = "patient_id",
#'                                        col_mo = "mo")) %>%
#'   # filter on first isolates
#'   filter(first_isolate == TRUE) %>%
#'   # get short MO names (like "E. coli")
#'   mutate(mo = mo_shortname(mo, Becker = TRUE)) %>%
#'   # select this short name and some antiseptic drugs
#'   select(mo, cfur, gent, cipr) %>%
#'   # group by MO
#'   group_by(mo) %>%
#'   # plot the thing, putting MOs on the facet
#'   ggplot_rsi(x = Antibiotic,
#'              facet = mo,
#'              translate_ab = FALSE,
#'              nrow = 1) +
#'   labs(title = "AMR of Top Three Microorganisms In Blood Culture Isolates",
#'        subtitle = "Only First Isolates, CoNS grouped according to Becker et al. (2014)",
#'        x = "Microorganisms")
#' }
ggplot_rsi <- function(data,
                       position = NULL,
                       x = "Antibiotic",
                       fill = "Interpretation",
                       # params = list(),
                       facet = NULL,
                       breaks = seq(0, 1, 0.1),
                       limits = NULL,
                       translate_ab = "official",
                       fun = count_df,
                       nrow = NULL,
                       datalabels = TRUE,
                       datalabels.size = 3,
                       datalabels.colour = "grey15",
                       ...) {

  if (!"ggplot2" %in% rownames(installed.packages())) {
    stop('this function requires the ggplot2 package.', call. = FALSE)
  }

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

  p <- ggplot2::ggplot(data = data) +
    geom_rsi(position = position, x = x, fill = fill, translate_ab = translate_ab, fun = fun, ...) +
    theme_rsi()

  if (fill == "Interpretation") {
    # set RSI colours
    p <- p + scale_rsi_colours()
  }
  if (is.null(position)) {
    position <- "fill"
  }
  if (fun_name == "portion_df"
      | (fun_name == "count_df" & position == "fill")) {
    # portions, so use y scale with percentage
    p <- p + scale_y_percent(breaks = breaks, limits = limits)
  }

  if (fun_name == "count_df" & datalabels == TRUE) {
    p <- p + labels_rsi_count(position = position,
                              x = x,
                              datalabels.size = datalabels.size,
                              datalabels.colour = datalabels.colour)
  }

  if (!is.null(facet)) {
    p <- p + facet_rsi(facet = facet, nrow = nrow)
  }

  p
}

#' @rdname ggplot_rsi
#' @export
geom_rsi <- function(position = NULL,
                     x = c("Antibiotic", "Interpretation"),
                     fill = "Interpretation",
                     translate_ab = "official",
                     fun = count_df,
                     ...)  {

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

  options(get_antibiotic_names = translate_ab)

  ggplot2::layer(geom = "bar", stat = "identity", position = position,
                 mapping = ggplot2::aes_string(x = x, y = y, fill = fill),
                 data = fun, params = list(...))

}

#' @rdname ggplot_rsi
#' @export
facet_rsi <- function(facet = c("Interpretation", "Antibiotic"), nrow = NULL) {

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
  if (all(breaks[breaks != 0] > 1)) {
    breaks <- breaks / 100
  }
  ggplot2::scale_y_continuous(breaks = breaks,
                              labels = percent(breaks),
                              limits = limits)
}

#' @rdname ggplot_rsi
#' @export
scale_rsi_colours <- function() {
  ggplot2::scale_fill_brewer(palette = "RdYlGn")
  #ggplot2::scale_fill_gradient2(low = "#d5613e", mid = "#ae5ac0", high = "#7daf44")
}

#' @rdname ggplot_rsi
#' @export
theme_rsi <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey75"))
}

#' @rdname ggplot_rsi
#' @export
labels_rsi_count <- function(position = NULL,
                             x = "Antibiotic",
                             datalabels.size = 3,
                             datalabels.colour = "grey15") {
  if (is.null(position)) {
    position <- "fill"
  }
  if (position == "fill") {
    position <- ggplot2::position_fill(vjust = 0.5)
  }
  ggplot2::geom_text(mapping = ggplot2::aes_string(label = "lbl",
                                                   x = x,
                                                   y = "Value"),
                     position = position,
                     data = getlbls,
                     inherit.aes = FALSE,
                     size = datalabels.size,
                     colour = datalabels.colour)
}

#' @importFrom dplyr %>% group_by mutate
getlbls <- function(data) {
  data %>%
    count_df() %>%
    group_by(Antibiotic) %>%
    mutate(lbl = paste0(percent(Value / sum(Value, na.rm = TRUE), force_zero = TRUE),
                        " (n=", Value, ")")) %>%
    mutate(lbl = ifelse(lbl == "0.0% (n=0)", "", lbl))
}
