# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

#' AMR bar plots with \code{ggplot}
#'
#' Use these functions to create bar plots for antimicrobial resistance analysis. All functions rely on internal \code{\link[ggplot2]{ggplot}} functions.
#' @param data a \code{data.frame} with column(s) of class \code{"rsi"} (see \code{\link{as.rsi}})
#' @param position position adjustment of bars, either \code{"stack"} (default) or \code{"dodge"}
#' @param x variable to show on x axis, either \code{"Antibiotic"} (default) or \code{"Interpretation"}
#' @param fill variable to categorise using the plots legend
#' @param facet variable to split plots by, either \code{"Interpretation"} (default) or \code{"Antibiotic"}
#' @details At default, the names of antibiotics will be shown on the plots using \code{\link{abname}}. This can be set with the option \code{get_antibiotic_names} (a logical value), so change it e.g. to \code{FALSE} with \code{options(get_antibiotic_names = FALSE)}.
#'
#' \strong{The functions}\cr
#' \code{geom_rsi} will take any variable from the data that has an \code{rsi} class (created with \code{\link{as.rsi}}) using \code{\link{portion_df}} and will plot bars with the percentage R, I and S. The default behaviour is to have the bars stacked and to have the different antibiotics on the x axis.
#'
#' \code{facet_rsi} creates 2d plots (at default based on S/I/R) using \code{\link[ggplot2]{facet_wrap}}.
#'
#' \code{scale_y_percent} transforms the y axis to a 0 to 100% range.
#'
#' \code{scale_rsi_colours} sets colours to the bars: green for S, yellow for I and red for R.
#'
#' \code{theme_rsi} is a \code{\link[ggplot2]{theme}} with minimal distraction.
#'
#' \code{ggplot_rsi} is a wrapper around all above functions that uses data as first input. This makes it possible to use this function after a pipe (\code{\%>\%}). See Examples.
#' @rdname ggplot_rsi
#' @export
#' @examples
#' library(dplyr)
#' library(ggplot2)
#'
#' # get antimicrobial results for drugs against a UTI:
#' ggplot(septic_patients %>% select(amox, nitr, fosf, trim, cipr)) +
#'   geom_rsi()
#'
#' # prettify it using some additional functions
#' df <- septic_patients[, c("amox", "nitr", "fosf", "trim", "cipr")]
#' ggplot(df) +
#'   geom_rsi(x = "Interpretation") +
#'   facet_rsi(facet = "Antibiotic") +
#'   scale_y_percent() +
#'   scale_rsi_colours() +
#'   theme_rsi()
#'
#' # or better yet, simplify this using the wrapper function - a single command:
#' septic_patients %>%
#'   select(amox, nitr, fosf, trim, cipr) %>%
#'   ggplot_rsi()
#'
#' septic_patients %>%
#'   select(amox, nitr, fosf, trim, cipr) %>%
#'   ggplot_rsi(x = "Interpretation", facet = "Antibiotic")
#'
#' # it also supports groups (don't forget to use facet on the group):
#' septic_patients %>%
#'   select(hospital_id, amox, cipr) %>%
#'   group_by(hospital_id) %>%
#'   ggplot_rsi() +
#'   facet_grid("hospital_id") +
#'   labs(title = "AMR of Amoxicillin And Ciprofloxacine Per Hospital")
ggplot_rsi <- function(data,
                       position = "stack",
                       x = "Antibiotic",
                       fill = "Interpretation",
                       facet = NULL) {

  if (!"ggplot2" %in% rownames(installed.packages())) {
    stop('this function requires the ggplot2 package.', call. = FALSE)
  }

  p <- ggplot2::ggplot(data = data) +
    geom_rsi(position = position, x = x, fill = fill) +
    scale_y_percent() +
    theme_rsi()

  if (fill == "Interpretation") {
    # set RSI colours
    p <- p + scale_rsi_colours()
  }

  if (!is.null(facet)) {
    p <- p + facet_rsi(facet = facet)
  }

  p
}

#' @rdname ggplot_rsi
#' @export
geom_rsi <- function(position = "stack", x = c("Antibiotic", "Interpretation"), fill = "Interpretation")  {

  x <- x[1]
  if (!x %in% c("Antibiotic", "Interpretation")) {
    stop("`x` must be 'Antibiotic' or 'Interpretation'")
  }

  ggplot2::layer(geom = "bar", stat = "identity", position = position,
                 mapping = ggplot2::aes_string(x = x, y = "Percentage", fill = fill),
                 data = AMR::portion_df, params = list())

}

#' @rdname ggplot_rsi
#' @export
facet_rsi <- function(facet = c("Interpretation", "Antibiotic")) {
  facet <- facet[1]
  if (!facet %in% c("Antibiotic", "Interpretation")) {
    stop("`facet` must be 'Antibiotic' or 'Interpretation'")
  }
  ggplot2::facet_wrap(facets = facet, scales = "free")
}

#' @rdname ggplot_rsi
#' @export
scale_y_percent <- function() {
  ggplot2::scale_y_continuous(name = "Percentage",
                              breaks = seq(0, 1, 0.1),
                              limits = c(0, 1),
                              labels = percent(seq(0, 1, 0.1)))
}

#' @rdname ggplot_rsi
#' @export
scale_rsi_colours <- function() {
  ggplot2::scale_fill_brewer(palette = "RdYlGn")
}

#' @rdname ggplot_rsi
#' @export
theme_rsi <- function() {
  theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "grey75"))
}
