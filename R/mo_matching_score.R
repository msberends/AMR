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

#' Calculate the Matching Score for Microorganisms
#'
#' This algorithm is used by [as.mo()] and all the [`mo_*`][mo_property()] functions to determine the most probable match of taxonomic records based on user input.
#' @author Dr. Matthijs Berends, 2018
#' @param x Any user input value(s)
#' @param n A full taxonomic name, that exists in [`microorganisms$fullname`][microorganisms]
#' @note This algorithm was originally described in: Berends MS *et al.* (2022). **AMR: An R Package for Working with Antimicrobial Resistance Data**. *Journal of Statistical Software*, 104(3), 1-31; \doi{10.18637/jss.v104.i03}.
#'
#' Later, the work of Bartlett A *et al.* about bacterial pathogens infecting humans (2022, \doi{10.1099/mic.0.001269}) was incorporated.
#' @section Matching Score for Microorganisms:
#' With ambiguous user input in [as.mo()] and all the [`mo_*`][mo_property()] functions, the returned results are chosen based on their matching score using [mo_matching_score()]. This matching score \eqn{m}, is calculated as:
#'
#' \ifelse{latex}{\deqn{m_{(x, n)} = \frac{l_{n} - 0.5 \cdot \min \begin{cases}l_{n} \\ \textrm{lev}(x, n)\end{cases}}{l_{n} \cdot p_{n} \cdot k_{n}}}}{
#'
#' \ifelse{html}{\figure{mo_matching_score.png}{options: width="300" alt="mo matching score"}}{m(x, n) = ( l_n * min(l_n, lev(x, n) ) ) / ( l_n * p_n * k_n )}}
#'
#' where:
#'
#' * \eqn{x} is the user input;
#' * \eqn{n} is a taxonomic name (genus, species, and subspecies);
#' * \eqn{l_n} is the length of \eqn{n};
#' * \eqn{lev} is the [Levenshtein distance function](https://en.wikipedia.org/wiki/Levenshtein_distance) (counting any insertion as 1, and any deletion or substitution as 2) that is needed to change \eqn{x} into \eqn{n};
#' * \eqn{p_n} is the human pathogenic prevalence group of \eqn{n}, as described below;
#' * \eqn{k_n} is the taxonomic kingdom of \eqn{n}, set as Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5.
#'
#' The grouping into human pathogenic prevalence \eqn{p} is based on recent work from Bartlett *et al.* (2022, \doi{10.1099/mic.0.001269}) who extensively studied medical-scientific literature to categorise all bacterial species into these groups:
#'
#' - **Established**, if a taxonomic species has infected at least three persons in three or more references. These records have `prevalence = 1.0` in the [microorganisms] data set;
#' - **Putative**, if a taxonomic species has fewer than three known cases. These records have `prevalence = 1.25` in the [microorganisms] data set.
#'
#' Furthermore,
#'
#' - Any genus present in the **established** list also has `prevalence = 1.0` in the [microorganisms] data set;
#' - Any other genus present in the **putative** list has `prevalence = 1.25` in the [microorganisms] data set;
#' - Any other species or subspecies of which the genus is present in the two aforementioned groups, has `prevalence = 1.5` in the [microorganisms] data set;
#' - Any *non-bacterial* genus, species or subspecies of which the genus is present in the following list, has `prevalence = 1.5` in the [microorganisms] data set: `r vector_or(MO_PREVALENT_GENERA, quotes = "*")`;
#' - All other records have `prevalence = 2.0` in the [microorganisms] data set.
#'
#' When calculating the matching score, all characters in \eqn{x} and \eqn{n} are ignored that are other than A-Z, a-z, 0-9, spaces and parentheses.
#'
#' All matches are sorted descending on their matching score and for all user input values, the top match will be returned. This will lead to the effect that e.g., `"E. coli"` will return the microbial ID of *Escherichia coli* (\eqn{m = `r round(mo_matching_score("E. coli", "Escherichia coli"), 3)`}, a highly prevalent microorganism found in humans) and not *Entamoeba coli* (\eqn{m = `r round(mo_matching_score("E. coli", "Entamoeba coli"), 3)`}, a less prevalent microorganism in humans), although the latter would alphabetically come first.
#' @export
#' @inheritSection AMR Reference Data Publicly Available
#' @examples
#' mo_reset_session()
#'
#' as.mo("E. coli")
#' mo_uncertainties()
#'
#' mo_matching_score(
#'   x = "E. coli",
#'   n = c("Escherichia coli", "Entamoeba coli")
#' )
mo_matching_score <- function(x, n) {
  meet_criteria(x, allow_class = c("character", "data.frame", "list"))
  meet_criteria(n, allow_class = "character")

  add_MO_lookup_to_AMR_env()

  x <- parse_and_convert(x)
  # no dots and other non-whitespace characters
  x <- gsub("[^a-zA-Z0-9 \\(\\)]+", "", x)

  # only keep one space
  x <- gsub(" +", " ", x)

  # force a capital letter, so this conversion will not count as a substitution
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))

  # n is always a taxonomically valid full name
  if (length(n) == 1) {
    n <- rep(n, length(x))
  }
  if (length(x) == 1) {
    x <- rep(x, length(n))
  }

  # length of fullname
  l_n <- nchar(n)
  lev <- double(length = length(x))
  l_n.lev <- double(length = length(x))
  # get Levenshtein distance
  lev <- unlist(Map(f = function(a, b) {
    as.double(utils::adist(a, b,
      ignore.case = FALSE,
      fixed = TRUE,
      costs = c(insertions = 1, deletions = 2, substitutions = 2),
      counts = FALSE
    ))
  }, x, n, USE.NAMES = FALSE))

  l_n.lev[l_n < lev] <- l_n[l_n < lev]
  l_n.lev[lev < l_n] <- lev[lev < l_n]
  l_n.lev[lev == l_n] <- lev[lev == l_n]

  # human pathogenic prevalence (1 to 3), see ?as.mo
  p_n <- AMR_env$MO_lookup[match(n, AMR_env$MO_lookup$fullname), "prevalence", drop = TRUE]
  # kingdom index (Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5)
  k_n <- AMR_env$MO_lookup[match(n, AMR_env$MO_lookup$fullname), "kingdom_index", drop = TRUE]

  # matching score:
  (l_n - 0.5 * l_n.lev) / (l_n * p_n * k_n)
}
