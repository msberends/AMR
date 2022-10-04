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

#' Calculate the Matching Score for Microorganisms
#'
#' This algorithm is used by [as.mo()] and all the [`mo_*`][mo_property()] functions to determine the most probable match of taxonomic records based on user input.
#' @author Dr. Matthijs Berends
#' @param x Any user input value(s)
#' @param n A full taxonomic name, that exists in [`microorganisms$fullname`][microorganisms]
#' @note This algorithm was described in: Berends MS *et al.* (2022). **AMR: An R Package for Working with Antimicrobial Resistance Data**. *Journal of Statistical Software*, 104(3), 1-31; \doi{10.18637/jss.v104.i03}.
#' @section Matching Score for Microorganisms:
#' With ambiguous user input in [as.mo()] and all the [`mo_*`][mo_property()] functions, the returned results are chosen based on their matching score using [mo_matching_score()]. This matching score \eqn{m}, is calculated as:
#'
#' \ifelse{latex}{\deqn{m_{(x, n)} = \frac{l_{n} - 0.5 \cdot \min \begin{cases}l_{n} \\ \textrm{lev}(x, n)\end{cases}}{l_{n} \cdot p_{n} \cdot k_{n}}}}{\ifelse{html}{\figure{mo_matching_score.png}{options: width="300" alt="mo matching score"}}{m(x, n) = ( l_n * min(l_n, lev(x, n) ) ) / ( l_n * p_n * k_n )}}
#'
#' where:
#'
#' * \ifelse{html}{\out{<i>x</i> is the user input;}}{\eqn{x} is the user input;}
#' * \ifelse{html}{\out{<i>n</i> is a taxonomic name (genus, species, and subspecies);}}{\eqn{n} is a taxonomic name (genus, species, and subspecies);}
#' * \ifelse{html}{\out{<i>l<sub>n</sub></i> is the length of <i>n</i>;}}{l_n is the length of \eqn{n};}
#' * \ifelse{html}{\out{<i>lev</i> is the <a href="https://en.wikipedia.org/wiki/Levenshtein_distance">Levenshtein distance function</a> (counting any insertion as 1, and any deletion or substitution as 2) that is needed to change <i>x</i> into <i>n</i>;}}{lev is the Levenshtein distance function (counting any insertion as 1, and any deletion or substitution as 2) that is needed to change \eqn{x} into \eqn{n};}
#' * \ifelse{html}{\out{<i>p<sub>n</sub></i> is the human pathogenic prevalence group of <i>n</i>, as described below;}}{p_n is the human pathogenic prevalence group of \eqn{n}, as described below;}
#' * \ifelse{html}{\out{<i>k<sub>n</sub></i> is the taxonomic kingdom of <i>n</i>, set as Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5.}}{l_n is the taxonomic kingdom of \eqn{n}, set as Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5.}
#'
#' The grouping into human pathogenic prevalence (\eqn{p}) is based on experience from several microbiological laboratories in the Netherlands in conjunction with international reports on pathogen prevalence:
#'
#' **Group 1** (most prevalent microorganisms) consists of all microorganisms where the taxonomic class is Gammaproteobacteria or where the taxonomic genus is *Enterococcus*, *Staphylococcus* or *Streptococcus*. This group consequently contains all common Gram-negative bacteria, such as *Pseudomonas* and *Legionella* and all species within the order Enterobacterales.
#'
#' **Group 2** consists of all microorganisms where the taxonomic phylum is Proteobacteria, Firmicutes, Actinobacteria or Sarcomastigophora, or where the taxonomic genus is `r vector_or(MO_PREVALENT_GENERA, quotes = "*")`.
#'
#' **Group 3** consists of all other microorganisms.
#'
#' All characters in \eqn{x} and \eqn{n} are ignored that are other than A-Z, a-z, 0-9, spaces and parentheses.
#'
#' All matches are sorted descending on their matching score and for all user input values, the top match will be returned. This will lead to the effect that e.g., `"E. coli"` will return the microbial ID of *Escherichia coli* (\eqn{m = `r round(mo_matching_score("E. coli", "Escherichia coli"), 3)`}, a highly prevalent microorganism found in humans) and not *Entamoeba coli* (\eqn{m = `r round(mo_matching_score("E. coli", "Entamoeba coli"), 3)`}, a less prevalent microorganism in humans), although the latter would alphabetically come first.
#' @export
#' @inheritSection AMR Reference Data Publicly Available
#' @examples
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

  x <- parse_and_convert(x)
  # no dots and other non-whitespace characters
  x <- gsub("[^a-zA-Z0-9 \\(\\)]+", "", x)

  # only keep one space
  x <- gsub(" +", " ", x)
  
  # start with a capital letter
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
  lev <- unlist(Map(f = function(a, b) {
    as.double(utils::adist(a, b, 
                           ignore.case = FALSE,
                           fixed = TRUE,
                           costs = c(insertions = 1, deletions = 2, substitutions = 2),
                           counts = FALSE))
  }, x, n, USE.NAMES = FALSE))
  
  l_n.lev[l_n < lev] <- l_n[l_n < lev]
  l_n.lev[lev < l_n] <- lev[lev < l_n]
  l_n.lev[lev == l_n] <- lev[lev == l_n]

  # human pathogenic prevalence (1 to 3), see ?as.mo
  p_n <- MO_lookup[match(n, MO_lookup$fullname), "prevalence", drop = TRUE]
  # kingdom index (Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5)
  k_n <- MO_lookup[match(n, MO_lookup$fullname), "kingdom_index", drop = TRUE]

  # matching score:
  (l_n - 0.5 * l_n.lev) / (l_n * p_n * k_n)
}
