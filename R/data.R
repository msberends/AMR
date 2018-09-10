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

#' Data set with 423 antibiotics
#'
#' A data set containing all antibiotics with a J0 code and some other antimicrobial agents, with their DDDs. Except for trade names and abbreviations, all properties were downloaded from the WHO, see Source.
#' @format A \code{\link{tibble}} with 423 observations and 18 variables:
#' \describe{
#'   \item{\code{atc}}{ATC code, like \code{J01CR02}}
#'   \item{\code{certe}}{Certe code, like \code{amcl}}
#'   \item{\code{umcg}}{UMCG code, like \code{AMCL}}
#'   \item{\code{abbr}}{Abbreviation as used by many countries, used internally by \code{\link{as.atc}}}
#'   \item{\code{official}}{Official name by the WHO, like \code{"Amoxicillin and beta-lactamase inhibitor"}}
#'   \item{\code{official_nl}}{Official name in the Netherlands, like \code{"Amoxicilline met enzymremmer"}}
#'   \item{\code{trivial_nl}}{Trivial name in Dutch, like \code{"Amoxicilline/clavulaanzuur"}}
#'   \item{\code{trade_name}}{Trade name as used by many countries (a total of 294), used internally by \code{\link{as.atc}}}
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
# use this later to further fill AMR::antibiotics
# drug <- "Ciprofloxacin"
# url <- xml2::read_html(paste0("https://www.ncbi.nlm.nih.gov/pccompound?term=", drug)) %>%
#   html_nodes(".rslt") %>%
#   .[[1]] %>%
#   html_nodes(".title a") %>%
#   html_attr("href") %>%
#   gsub("/compound/", "/rest/pug_view/data/compound/", ., fixed = TRUE) %>%
#   paste0("/XML/?response_type=display")
# synonyms <- url %>%
#   read_xml() %>%
#   xml_contents() %>% .[[6]] %>%
#   xml_contents() %>% .[[8]] %>%
#   xml_contents() %>% .[[3]] %>%
#   xml_contents() %>% .[[3]] %>%
#   xml_contents() %>%
#   paste() %>%
#   .[. %like% "StringValueList"] %>%
#   gsub("[</]+StringValueList[>]", "", .)

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
#
# ADD NEW TRADE NAMES FROM OTHER DATAFRAME
# antibiotics_add_to_property <- function(ab_df, atc, property, value) {
#   if (length(atc) > 1L) {
#     stop("only one atc at a time")
#   }
#   if (!property %in% c("abbr", "trade_name")) {
#     stop("only possible for abbr and trade_name")
#   }
#
#   value <- gsub(ab_df[which(ab_df$atc == atc),] %>% pull("official"), "", value, fixed = TRUE)
#   value <- gsub("||", "|", value, fixed = TRUE)
#   value <- gsub("[äáàâ]", "a", value)
#   value <- gsub("[ëéèê]", "e", value)
#   value <- gsub("[ïíìî]", "i", value)
#   value <- gsub("[öóòô]", "o", value)
#   value <- gsub("[üúùû]", "u", value)
#   if (!atc %in% ab_df$atc) {
#     message("SKIPPING - UNKNOWN ATC: ", atc)
#   }
#   if (is.na(value)) {
#     message("SKIPPING - VALUE MISSES: ", atc)
#   }
#   if (atc %in% ab_df$atc & !is.na(value)) {
#     current <- ab_df[which(ab_df$atc == atc),] %>% pull(property)
#     if (!is.na(current)) {
#       value <- paste(current, value, sep = "|")
#     }
#     value <- strsplit(value, "|", fixed = TRUE) %>% unlist() %>% unique() %>% paste(collapse = "|")
#     value <- gsub("||", "|", value, fixed = TRUE)
#     # print(value)
#     ab_df[which(ab_df$atc == atc), property] <- value
#     message("Added ", value, " to ", ab_official(atc), " (", atc, ", ", ab_certe(atc), ")")
#   }
#   ab_df
# }
#
"antibiotics"

#' Data set with human pathogenic microorganisms
#'
#' A data set containing 2,630 (potential) human pathogenic microorganisms. MO codes can be looked up using \code{\link{guess_mo}}.
#' @format A \code{\link{tibble}} with 2,630 observations and 10 variables:
#' \describe{
#'   \item{\code{mo}}{ID of microorganism}
#'   \item{\code{bactsys}}{Bactsyscode of microorganism}
#'   \item{\code{family}}{Family name of microorganism}
#'   \item{\code{genus}}{Genus name of microorganism, like \code{"Echerichia"}}
#'   \item{\code{species}}{Species name of microorganism, like \code{"coli"}}
#'   \item{\code{subspecies}}{Subspecies name of bio-/serovar of microorganism, like \code{"EHEC"}}
#'   \item{\code{fullname}}{Full name, like \code{"Echerichia coli (EHEC)"}}
#'   \item{\code{aerobic}}{Logical whether bacteria is aerobic}
#'   \item{\code{type}}{Type of microorganism, like \code{"Bacteria"} and \code{"Fungus/yeast"}}
#'   \item{\code{gramstain}}{Gram of microorganism, like \code{"Negative rods"}}
#' }
#  source MOLIS (LIS of Certe) - \url{https://www.certe.nl}
# new <- microorganisms %>% filter(genus == "Bacteroides") %>% .[1,]
# new[1, 'mo'] <- "DIAPNU"
# new[1, 'bactsys'] <- "DIAPNU"
# new[1, 'family'] <- "Veillonellaceae"
# new[1, 'genus'] <- "Dialister"
# new[1, 'species'] <- "pneumosintes"
# new[1, 'subspecies'] <- NA
# new[1, 'fullname'] <- paste(new[1, 'genus'], new[1, 'species'])
# microorganisms <- microorganisms %>% bind_rows(new) %>% arrange(mo)
#' @seealso \code{\link{guess_mo}} \code{\link{antibiotics}} \code{\link{microorganisms.umcg}}
"microorganisms"

#' Translation table for UMCG with ~1,100 microorganisms
#'
#' A data set containing all bacteria codes of UMCG MMB. These codes can be joined to data with an ID from \code{\link{microorganisms}$mo} (using \code{\link{left_join_microorganisms}}). GLIMS codes can also be translated to valid \code{MO}s with \code{\link{guess_mo}}.
#' @format A \code{\link{tibble}} with 1,095 observations and 2 variables:
#' \describe{
#'   \item{\code{umcg}}{Code of microorganism according to UMCG MMB}
#'   \item{\code{mo}}{Code of microorganism in \code{\link{microorganisms}}}
#' }
# source MOLIS (LIS of Certe) - \url{https://www.certe.nl} \cr \cr GLIMS (LIS of UMCG) - \url{https://www.umcg.nl}
#' @seealso \code{\link{guess_mo}} \code{\link{microorganisms}}
"microorganisms.umcg"

#' Data set with 2000 blood culture isolates of septic patients
#'
#' An anonymised data set containing 2,000 microbial blood culture isolates with their full antibiograms found in septic patients in 4 different hospitals in the Netherlands, between 2001 and 2017. It is true, genuine data. This \code{data.frame} can be used to practice AMR analysis. For examples, press F1.
#' @format A \code{\link{tibble}} with 2,000 observations and 49 variables:
#' \describe{
#'   \item{\code{date}}{date of receipt at the laboratory}
#'   \item{\code{hospital_id}}{ID of the hospital, from A to D}
#'   \item{\code{ward_icu}}{logical to determine if ward is an intensive care unit}
#'   \item{\code{ward_clinical}}{logical to determine if ward is a regular clinical ward}
#'   \item{\code{ward_outpatient}}{logical to determine if ward is an outpatient clinic}
#'   \item{\code{age}}{age of the patient}
#'   \item{\code{sex}}{sex of the patient}
#'   \item{\code{patient_id}}{ID of the patient, first 10 characters of an SHA hash containing irretrievable information}
#'   \item{\code{mo}}{ID of microorganism, see \code{\link{microorganisms}}}
#'   \item{\code{peni:rifa}}{40 different antibiotics with class \code{rsi} (see \code{\link{as.rsi}}); these column names occur in \code{\link{antibiotics}} data set and can be translated with \code{\link{abname}}}
#' }
# source MOLIS (LIS of Certe) - \url{https://www.certe.nl}
#' @examples
#' # ----------- #
#' # PREPARATION #
#' # ----------- #
#'
#' # Save this example data set to an object, so we can edit it:
#' my_data <- septic_patients
#'
#' # load the dplyr package to make data science A LOT easier
#' library(dplyr)
#'
#' # Add first isolates to our data set:
#' my_data <- my_data %>%
#'   mutate(first_isolates = first_isolate(my_data, "date", "patient_id", "mo"))
#'
#' # -------- #
#' # ANALYSIS #
#' # -------- #
#'
#' # 1. Get the amoxicillin resistance percentages (p)
#' #     and numbers (n) of E. coli, divided by hospital:
#'
#' my_data %>%
#'   filter(mo == guess_mo("E. coli"),
#'          first_isolates == TRUE) %>%
#'   group_by(hospital_id) %>%
#'   summarise(n = n_rsi(amox),
#'             p = portion_IR(amox))
#'
#'
#' # 2. Get the amoxicillin/clavulanic acid resistance
#' #    percentages of E. coli, trend over the years:
#'
#' my_data %>%
#'   filter(mo == guess_mo("E. coli"),
#'          first_isolates == TRUE) %>%
#'   group_by(year = format(date, "%Y")) %>%
#'   summarise(n = n_rsi(amcl),
#'             p = portion_IR(amcl, minimum = 20))
"septic_patients"
