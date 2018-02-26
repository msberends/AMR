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

#' Dataset with 420 antibiotics
#'
#' A dataset containing all antibiotics with a J0 code, with their DDD's.
#' @format A data.frame with 420 observations and 12 variables:
#' \describe{
#'   \item{\code{atc}}{ATC code, like \code{J01CR02}}
#'   \item{\code{molis}}{MOLIS code, like \code{amcl}}
#'   \item{\code{umcg}}{UMCG code, like \code{AMCL}}
#'   \item{\code{official}}{Official name by the WHO, like \code{"amoxicillin and enzyme inhibitor"}}
#'   \item{\code{official_nl}}{Official name in the Netherlands, like \code{"Amoxicilline met enzymremmer"}}
#'   \item{\code{trivial_nl}}{Trivial name in Dutch, like \code{"Amoxicilline/clavulaanzuur"}}
#'   \item{\code{oral_ddd}}{Daily Defined Dose (DDD) according to the WHO, oral treatment}
#'   \item{\code{oral_units}}{Units of \code{ddd_units}}
#'   \item{\code{iv_ddd}}{Daily Defined Dose (DDD) according to the WHO, parenteral treatment}
#'   \item{\code{iv_units}}{Units of \code{iv_ddd}}
#'   \item{\code{atc_group1}}{ATC group in Dutch, like \code{"Macroliden, lincosamiden en streptograminen"}}
#'   \item{\code{atc_group2}}{Subgroup of \code{atc_group1} in Dutch, like \code{"Macroliden"}}
#' }
#' @source MOLIS (LIS of Certe) - \url{https://www.certe.nl} \cr \cr GLIMS (LIS of UMCG) - \url{https://www.umcg.nl} \cr \cr World Health Organization - \url{https://www.whocc.no/atc_ddd_index/}
#' @seealso \code{\link{bactlist}}
# todo:
# ablist <- ablist %>% mutate(useful_gramnegative = if_else(atc_group2 == 'Tetracyclines', FALSE, TRUE))
# ablist <- ablist %>% mutate(useful_gramnegative = if_else(atc_group2 %like% 'Glycopept', FALSE, useful_gramnegative))
# Tbl1 Enterobacteriaceae are also intrinsically resistant to benzylpenicillin, glycopeptides, fusidic acid, macrolides (with some exceptions1), lincosamides, streptogramins, rifampicin, daptomycin and linezolid.
# Tbl2 Non-fermentative Gram-negative bacteria are also generally intrinsically resistant to benzylpenicillin, first and second generation cephalosporins, glycopeptides, fusidic acid, macrolides, lincosamides, streptogramins, rifampicin, daptomycin and linezolid
# Tbl3 Gram-negative bacteria other than Enterobacteriaceae and non-fermentative Gram-negative bacteria listed are also intrinsically resistant to glycopeptides, lincosamides, daptomycin and linezolid.
"ablist"

#' Dataset with ~2500 microorganisms
#'
#' A dataset containing all microorganisms of MOLIS. MO codes of the UMCG can be looked up using \code{\link{bactlist.umcg}}.
#' @format A data.frame with 2507 observations and 10 variables:
#' \describe{
#'   \item{\code{bactid}}{ID of microorganism}
#'   \item{\code{bactsys}}{Bactsyscode of microorganism}
#'   \item{\code{family}}{Family name of microorganism}
#'   \item{\code{genus}}{Genus name of microorganism, like \code{"Echerichia"}}
#'   \item{\code{species}}{Species name of microorganism, like \code{"coli"}}
#'   \item{\code{subspecies}}{Subspecies name of bio-/serovar of microorganism, like \code{"EHEC"}}
#'   \item{\code{fullname}}{Full name, like \code{"Echerichia coli (EHEC)"}}
#'   \item{\code{type}}{Type of microorganism in Dutch, like \code{"Bacterie"} and \code{"Schimmel/gist"}}
#'   \item{\code{gramstain}}{Gram of microorganism in Dutch, like \code{"Negatieve staven"}}
#'   \item{\code{aerobic}}{Type aerobe/anaerobe of bacteria}
#' }
#' @source MOLIS (LIS of Certe) - \url{https://www.certe.nl}
#' @seealso \code{\link{ablist}} \code{\link{bactlist.umcg}}
"bactlist"

#' Translation table for UMCG with ~1100 microorganisms
#'
#' A dataset containing all bacteria codes of UMCG MMB. These codes can be joined to data with an ID from \code{\link{bactlist}$bactid}, using \code{\link{left_join_bactlist}}.
#' @format A data.frame with 1090 observations and 2 variables:
#' \describe{
#'   \item{\code{mocode}}{Code of microorganism according to UMCG MMB}
#'   \item{\code{bactid}}{Code of microorganism in \code{\link{bactlist}}}
#' }
#' @source MOLIS (LIS of Certe) - \url{https://www.certe.nl} \cr \cr GLIMS (LIS of UMCG) - \url{https://www.umcg.nl}
#' @seealso \code{\link{bactlist}}
"bactlist.umcg"
