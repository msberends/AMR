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
#' A dataset containing all antibiotics with a J0 code, with their DDD's. Properties were downloaded from the WHO, see Source.
#' @format A data.frame with 420 observations and 16 variables:
#' \describe{
#'   \item{\code{atc}}{ATC code, like \code{J01CR02}}
#'   \item{\code{molis}}{MOLIS code, like \code{amcl}}
#'   \item{\code{umcg}}{UMCG code, like \code{AMCL}}
#'   \item{\code{official}}{Official name by the WHO, like \code{"Amoxicillin and enzyme inhibitor"}}
#'   \item{\code{official_nl}}{Official name in the Netherlands, like \code{"Amoxicilline met enzymremmer"}}
#'   \item{\code{trivial_nl}}{Trivial name in Dutch, like \code{"Amoxicilline/clavulaanzuur"}}
#'   \item{\code{oral_ddd}}{Defined Daily Dose (DDD), oral treatment}
#'   \item{\code{oral_units}}{Units of \code{ddd_units}}
#'   \item{\code{iv_ddd}}{Defined Daily Dose (DDD), parenteral treatment}
#'   \item{\code{iv_units}}{Units of \code{iv_ddd}}
#'   \item{\code{atc_group1}}{ATC group, like \code{"Macrolides, lincosamides and streptogramins"}}
#'   \item{\code{atc_group2}}{Subgroup of \code{atc_group1}, like \code{"Macrolides"}}
#'   \item{\code{atc_group1_nl}}{ATC group in Dutch, like \code{"Macroliden, lincosamiden en streptograminen"}}
#'   \item{\code{atc_group2_nl}}{Subgroup of \code{atc_group1} in Dutch, like \code{"Macroliden"}}
#'   \item{\code{useful_gramnegative}}{\code{FALSE} if not useful according to EUCAST, \code{NA} otherwise (see Source)}
#'   \item{\code{useful_grampositive}}{\code{FALSE} if not useful according to EUCAST, \code{NA} otherwise (see Source)}
#' }
#' @source - World Health Organization: \url{https://www.whocc.no/atc_ddd_index/} \cr - EUCAST - Expert rules intrinsic exceptional V3.1 \cr - MOLIS (LIS of Certe): \url{https://www.certe.nl} \cr - GLIMS (LIS of UMCG): \url{https://www.umcg.nl}
#' @seealso \code{\link{microorganisms}}
# last two columns created with:
# antibiotics %>%
#   mutate(useful_gramnegative = 
#            if_else(
#              atc_group1 %like% '(fusidic|glycopeptide|macrolide|lincosamide|daptomycin|linezolid)' |
#                atc_group2 %like% '(fusidic|glycopeptide|macrolide|lincosamide|daptomycin|linezolid)' |
#                official %like% '(fusidic|glycopeptide|macrolide|lincosamide|daptomycin|linezolid)',
#              FALSE,
#              NA
#            ),
#          useful_grampositive =
#            if_else(
#              atc_group1 %like% '(aztreonam|temocillin|polymyxin|colistin|nalidixic)' |
#                atc_group2 %like% '(aztreonam|temocillin|polymyxin|colistin|nalidixic)' |
#                official %like% '(aztreonam|temocillin|polymyxin|colistin|nalidixic)',
#              FALSE,
#              NA
#            )
#   )
"antibiotics"

#' Dataset with ~2500 microorganisms
#'
#' A dataset containing 2500 microorganisms. MO codes of the UMCG can be looked up using \code{\link{microorganisms.umcg}}.
#' @format A data.frame with 2507 observations and 12 variables:
#' \describe{
#'   \item{\code{bactid}}{ID of microorganism}
#'   \item{\code{bactsys}}{Bactsyscode of microorganism}
#'   \item{\code{family}}{Family name of microorganism}
#'   \item{\code{genus}}{Genus name of microorganism, like \code{"Echerichia"}}
#'   \item{\code{species}}{Species name of microorganism, like \code{"coli"}}
#'   \item{\code{subspecies}}{Subspecies name of bio-/serovar of microorganism, like \code{"EHEC"}}
#'   \item{\code{fullname}}{Full name, like \code{"Echerichia coli (EHEC)"}}
#'   \item{\code{type}}{Type of microorganism, like \code{"Bacteria"} and \code{"Fungus/yeast"}}
#'   \item{\code{gramstain}}{Gram of microorganism, like \code{"Negative rods"}}
#'   \item{\code{aerobic}}{Logical whether bacteria is aerobic}
#'   \item{\code{type_nl}}{Type of microorganism in Dutch, like \code{"Bacterie"} and \code{"Schimmel/gist"}}
#'   \item{\code{gramstain_nl}}{Gram of microorganism in Dutch, like \code{"Negatieve staven"}}
#' }
#' @source MOLIS (LIS of Certe) - \url{https://www.certe.nl}
#' @seealso \code{\link{guess_bactid}} \code{\link{antibiotics}} \code{\link{microorganisms.umcg}}
"microorganisms"

#' Translation table for UMCG with ~1100 microorganisms
#'
#' A dataset containing all bacteria codes of UMCG MMB. These codes can be joined to data with an ID from \code{\link{microorganisms}$bactid} (using \code{\link{left_join_microorganisms}}). GLIMS codes can also be translated to valid \code{bactid}'s with \code{\link{guess_bactid}}.
#' @format A data.frame with 1090 observations and 2 variables:
#' \describe{
#'   \item{\code{mocode}}{Code of microorganism according to UMCG MMB}
#'   \item{\code{bactid}}{Code of microorganism in \code{\link{microorganisms}}}
#' }
#' @source MOLIS (LIS of Certe) - \url{https://www.certe.nl} \cr \cr GLIMS (LIS of UMCG) - \url{https://www.umcg.nl}
#' @seealso \code{\link{guess_bactid}} \code{\link{microorganisms}}
"microorganisms.umcg"

#' Dataset with 2000 blood culture isolates of septic patients
#'
#' An anonymised dataset containing 2000 microbial blood culture isolates with their antibiogram of septic patients found in 5 different hospitals in the Netherlands, between 2001 and 2017. This data.frame can be used to practice AMR analysis. For examples, press F1.
#' @format A data.frame with 2000 observations and 47 variables:
#' \describe{
#'   \item{\code{date}}{date of receipt at the laboratory}
#'   \item{\code{hospital_id}}{ID of the hospital}
#'   \item{\code{ward_icu}}{logical to determine if ward is an intensive care unit}
#'   \item{\code{ward_clinical}}{logical to determine if ward is a regular clinical ward}
#'   \item{\code{ward_outpatient}}{logical to determine if ward is an outpatient clinic}
#'   \item{\code{age}}{age of the patient}
#'   \item{\code{sex}}{sex of the patient}
#'   \item{\code{patient_id}}{ID of the patient, first 10 characters of an SHA hash containing irretrievable information}
#'   \item{\code{bactid}}{ID of microorganism, see \code{\link{microorganisms}}}
#'   \item{\code{peni:mupi}}{38 different antibiotics with class \code{rsi} (see \code{\link{as.rsi}}); these column names occur in \code{\link{antibiotics}} and can be translated with \code{\link{abname}}}
#' }
#' @source MOLIS (LIS of Certe) - \url{https://www.certe.nl}
#' @examples
#' # ----------- #
#' # PREPARATION #
#' # ----------- #
#' 
#' # Save this example dataset to an object, so we can edit it:
#' my_data <- septic_patients
#' 
#' # load the dplyr package to make data science A LOT easier
#' library(dplyr)
#' 
#' # Add first isolates to our dataset:
#' my_data <- my_data %>% 
#'   mutate(first_isolates = first_isolate(my_data, date, patient_id, bactid))
#' 
#' # -------- #
#' # ANALYSIS #
#' # -------- #
#' 
#' # 1. Get the amoxicillin resistance percentages 
#' #    of E. coli, divided by hospital:
#' 
#' my_data %>%
#'   filter(bactid == "ESCCOL",
#'          first_isolates == TRUE) %>% 
#'   group_by(hospital_id) %>% 
#'   summarise(n = n(),
#'             amoxicillin_resistance = rsi(amox))
#'   
#'   
#' # 2. Get the amoxicillin/clavulanic acid resistance 
#' #    percentages of E. coli, trend over the years:
#' 
#' my_data %>% 
#'   filter(bactid == guess_bactid("E. coli"),
#'          first_isolates == TRUE) %>% 
#'   group_by(year = format(date, "%Y")) %>% 
#'   summarise(n = n(),
#'             amoxclav_resistance = rsi(amcl, minimum = 20))
"septic_patients"
