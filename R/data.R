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

#' Data set with ~500 antibiotics
#'
#' A data set containing all antibiotics with a J0 code and some other antimicrobial agents, with their DDDs. Except for trade names and abbreviations, all properties were downloaded from the WHO, see Source.
#' @format A \code{\link{data.frame}} with 488 observations and 17 variables:
#' \describe{
#'   \item{\code{atc}}{ATC code (Anatomical Therapeutic Chemical), like \code{J01CR02}}
#'   \item{\code{ears_net}}{EARS-Net code (European Antimicrobial Resistance Surveillance Network), like \code{AMC}}
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
#'   \item{\code{useful_gramnegative}}{\code{FALSE} if not useful according to EUCAST, \code{NA} otherwise (see Source)}
#'   \item{\code{useful_grampositive}}{\code{FALSE} if not useful according to EUCAST, \code{NA} otherwise (see Source)}
#' }
#' @source World Health Organization (WHO) Collaborating Centre for Drug Statistics Methodology: \url{https://www.whocc.no/atc_ddd_index/}
#'
#' Table antibiotic coding EARSS (from WHONET 5.3): \url{http://www.madsonline.dk/Tutorials/landskoder_antibiotika_WM.pdf}
#'
#' EUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes Tables. Version 3.1, 2016: \url{http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf}
#'
#' European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER: \url{http://ec.europa.eu/health/documents/community-register/html/atc.htm}
#' @inheritSection WHOCC WHOCC
#' @inheritSection AMR Read more on our website!
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

#' Data set with ~60,000 microorganisms
#'
#' A data set containing the microbial taxonomy of six kingdoms from the Catalogue of Life. MO codes can be looked up using \code{\link{as.mo}}.
#' @inheritSection catalogue_of_life Catalogue of Life
#' @format A \code{\link{data.frame}} with 56,672 observations and 14 variables:
#' \describe{
#'   \item{\code{mo}}{ID of microorganism as used by this package}
#'   \item{\code{col_id}}{Catalogue of Life ID}
#'   \item{\code{fullname}}{Full name, like \code{"Echerichia coli"}}
#'   \item{\code{kingdom}}{Taxonomic kingdom of the microorganism}
#'   \item{\code{phylum}}{Taxonomic phylum of the microorganism}
#'   \item{\code{class}}{Taxonomic class of the microorganism}
#'   \item{\code{order}}{Taxonomic order of the microorganism}
#'   \item{\code{family}}{Taxonomic family of the microorganism}
#'   \item{\code{genus}}{Taxonomic genus of the microorganism}
#'   \item{\code{species}}{Taxonomic species of the microorganism}
#'   \item{\code{subspecies}}{Taxonomic subspecies of the microorganism}
#'   \item{\code{rank}}{Taxonomic rank of the microorganism, like \code{"species"} or \code{"genus"}}
#'   \item{\code{ref}}{Author(s) and year of concerning scientific publication}
#'   \item{\code{species_id}}{ID of the species as used by the Catalogue of Life}
#' }
#' @source Catalogue of Life: Annual Checklist (public online database), \url{www.catalogueoflife.org}.
#' @details Manually added were:
#' \itemize{
#'   \item{9 species of \emph{Streptococcus} (beta haemolytic groups A, B, C, D, F, G, H, K and unspecified)}
#'   \item{2 species of \emph{Staphylococcus} (coagulase-negative [CoNS] and coagulase-positive [CoPS])}
#'   \item{2 other undefined (unknown Gram negatives and unknown Gram positives)}
#' }
#' @inheritSection AMR Read more on our website!
#' @seealso \code{\link{as.mo}}, \code{\link{mo_property}}, \code{\link{microorganisms.codes}}
"microorganisms"

catalogue_of_life <- list(
  version = "Catalogue of Life: 2018 Annual Checklist",
  url = "http://www.catalogueoflife.org/annual-checklist/2018"
)

#' Version info of included Catalogue of Life
#' @seealso \code{\link{microorganisms}}
#' @inheritSection catalogue_of_life Catalogue of Life
#' @export
catalogue_of_life_version <- function() {
  list(version = catalogue_of_life$version,
       url = catalogue_of_life$url,
       no_of_species = nrow(AMR::microorganisms),
       no_of_synonyms = nrow(AMR::microorganisms.old))
}

#' Data set with previously accepted taxonomic names
#'
#' A data set containing old (previously valid or accepted) taxonomic names according to the Catalogue of Life. This data set is used internally by \code{\link{as.mo}}.
#' @inheritSection catalogue_of_life Catalogue of Life
#' @format A \code{\link{data.frame}} with 14,506 observations and 4 variables:
#' \describe{
#'   \item{\code{col_id}}{Catalogue of Life ID}
#'   \item{\code{tsn_new}}{New Catalogue of Life ID}
#'   \item{\code{fullname}}{Old taxonomic name of the microorganism}
#'   \item{\code{ref}}{Author(s) and year of concerning scientific publication}
#' }
#' @source [3] Catalogue of Life: Annual Checklist (public online database), \url{www.catalogueoflife.org}.
#' @inheritSection AMR Read more on our website!
#' @seealso \code{\link{as.mo}} \code{\link{mo_property}} \code{\link{microorganisms}}
"microorganisms.old"

#' Translation table for microorganism codes
#'
#' A data set containing commonly used codes for microorganisms, from laboratory systems and WHONET. Define your own with \code{\link{set_mo_source}}.
#' @format A \code{\link{data.frame}} with 4,731 observations and 2 variables:
#' \describe{
#'   \item{\code{certe}}{Commonly used code of a microorganism}
#'   \item{\code{mo}}{ID of the microorganism in the \code{\link{microorganisms}} data set}
#' }
#' @inheritSection catalogue_of_life Catalogue of Life
#' @inheritSection AMR Read more on our website!
#' @seealso \code{\link{as.mo}} \code{\link{microorganisms}}
"microorganisms.codes"

#' Data set with 2,000 blood culture isolates from septic patients
#'
#' An anonymised data set containing 2,000 microbial blood culture isolates with their full antibiograms found in septic patients in 4 different hospitals in the Netherlands, between 2001 and 2017. It is true, genuine data. This \code{data.frame} can be used to practice AMR analysis. For examples, please read \href{https://msberends.gitlab.io/AMR/articles/AMR.html}{the tutorial on our website}.
#' @format A \code{\link{data.frame}} with 2,000 observations and 49 variables:
#' \describe{
#'   \item{\code{date}}{date of receipt at the laboratory}
#'   \item{\code{hospital_id}}{ID of the hospital, from A to D}
#'   \item{\code{ward_icu}}{logical to determine if ward is an intensive care unit}
#'   \item{\code{ward_clinical}}{logical to determine if ward is a regular clinical ward}
#'   \item{\code{ward_outpatient}}{logical to determine if ward is an outpatient clinic}
#'   \item{\code{age}}{age of the patient}
#'   \item{\code{gender}}{gender of the patient}
#'   \item{\code{patient_id}}{ID of the patient, first 10 characters of an SHA hash containing irretrievable information}
#'   \item{\code{mo}}{ID of microorganism created with \code{\link{as.mo}}, see also \code{\link{microorganisms}}}
#'   \item{\code{peni:rifa}}{40 different antibiotics with class \code{rsi} (see \code{\link{as.rsi}}); these column names occur in \code{\link{antibiotics}} data set and can be translated with \code{\link{abname}}}
#' }
#' @inheritSection AMR Read more on our website!
"septic_patients"

#' Data set with 500 isolates - WHONET example
#'
#' This example data set has the exact same structure as an export file from WHONET. Such files can be used with this package, as this example data set shows. The data itself was based on our \code{\link{septic_patients}} data set.
#' @format A \code{\link{data.frame}} with 500 observations and 53 variables:
#' \describe{
#'   \item{\code{Identification number}}{ID of the sample}
#'   \item{\code{Specimen number}}{ID of the specimen}
#'   \item{\code{Organism}}{Name of the microorganism. Before analysis, you should transform this to a valid microbial class, using \code{\link{as.mo}}.}
#'   \item{\code{Country}}{Country of origin}
#'   \item{\code{Laboratory}}{Name of laboratory}
#'   \item{\code{Last name}}{Last name of patient}
#'   \item{\code{First name}}{Initial of patient}
#'   \item{\code{Sex}}{Gender of patient}
#'   \item{\code{Age}}{Age of patient}
#'   \item{\code{Age category}}{Age group, can also be looked up using \code{\link{age_groups}}}
#'   \item{\code{Date of admission}}{Date of hospital admission}
#'   \item{\code{Specimen date}}{Date when specimen was received at laboratory}
#'   \item{\code{Specimen type}}{Specimen type or group}
#'   \item{\code{Specimen type (Numeric)}}{Translation of \code{"Specimen type"}}
#'   \item{\code{Reason}}{Reason of request with Differential Diagnosis}
#'   \item{\code{Isolate number}}{ID of isolate}
#'   \item{\code{Organism type}}{Type of microorganism, can also be looked up using \code{\link{mo_type}}}
#'   \item{\code{Serotype}}{Serotype of microorganism}
#'   \item{\code{Beta-lactamase}}{Microorganism produces beta-lactamase?}
#'   \item{\code{ESBL}}{Microorganism produces extended spectrum beta-lactamase?}
#'   \item{\code{Carbapenemase}}{Microorganism produces carbapenemase?}
#'   \item{\code{MRSA screening test}}{Microorganism is possible MRSA?}
#'   \item{\code{Inducible clindamycin resistance}}{Clindamycin can be induced?}
#'   \item{\code{Comment}}{Other comments}
#'   \item{\code{Date of data entry}}{Date this data was entered in WHONET}
#'   \item{\code{AMP_ND10:CIP_EE}}{27 different antibiotics. You can lookup the abbreviatons in the \code{\link{antibiotics}} data set, or use e.g. \code{\link{atc_name}("AMP")} to get the official name immediately. Before analysis, you should transform this to a valid antibiotic class, using \code{\link{as.rsi}}.}
#' }
#' @inheritSection AMR Read more on our website!
"WHONET"
