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

#' Data set with ~450 antibiotics
#'
#' A data set containing all antibiotics. Use \code{\link{as.ab}} or one of the \code{\link{ab_property}} functions to retrieve values from this data set. Three identifiers are included in this data set: an antibiotic ID (\code{ab}, primarily used in this package) as defined by WHONET/EARS-Net, an ATC code (\code{atc}) as defined by the WHO, and a Compound ID (\code{cid}) as found in PubChem. Other properties in this data set are derived from one or more of these codes.
#' @format A \code{\link{data.frame}} with 453 observations and 13 variables:
#' \describe{
#'   \item{\code{ab}}{Antibiotic ID as used in this package (like \code{AMC}), using the official EARS-Net (European Antimicrobial Resistance Surveillance Network) codes where available}
#'   \item{\code{atc}}{ATC code (Anatomical Therapeutic Chemical) as defined by the WHOCC, like \code{J01CR02}}
#'   \item{\code{cid}}{Compound ID as found in PubChem}
#'   \item{\code{name}}{Official name as used by WHONET/EARS-Net or the WHO}
#'   \item{\code{group}}{A short and concise group name, based on WHONET and WHOCC definitions}
#'   \item{\code{atc_group1}}{Official pharmacological subgroup (3rd level ATC code) as defined by the WHOCC, like \code{"Macrolides, lincosamides and streptogramins"}}
#'   \item{\code{atc_group2}}{Official chemical subgroup (4th level ATC code) as defined by the WHOCC, like \code{"Macrolides"}}
#'   \item{\code{abbr}}{List of abbreviations as used in many countries, also for antibiotic susceptibility testing (AST)}
#'   \item{\code{synonyms}}{Synonyms (often trade names) of a drug, as found in PubChem based on their compound ID}
#'   \item{\code{oral_ddd}}{Defined Daily Dose (DDD), oral treatment}
#'   \item{\code{oral_units}}{Units of \code{ddd_units}}
#'   \item{\code{iv_ddd}}{Defined Daily Dose (DDD), parenteral treatment}
#'   \item{\code{iv_units}}{Units of \code{iv_ddd}}
#' }
#' @details Properties that are based on an ATC code are only available when an ATC is available. These properties are: \code{atc_group1}, \code{atc_group2}, \code{oral_ddd}, \code{oral_units}, \code{iv_ddd} and \code{iv_units}
#'
#' Synonyms (i.e. trade names) are derived from the Compound ID (\code{cid}) and consequently only available where a CID is available.
#' @source World Health Organization (WHO) Collaborating Centre for Drug Statistics Methodology (WHOCC): \url{https://www.whocc.no/atc_ddd_index/}
#'
#' WHONET 2019 software: \url{http://www.whonet.org/software.html}
#'
#' European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER: \url{http://ec.europa.eu/health/documents/community-register/html/atc.htm}
#' @inheritSection WHOCC WHOCC
#' @inheritSection AMR Read more on our website!
#' @seealso \code{\link{microorganisms}}
"antibiotics"

#' Data set with ~70,000 microorganisms
#'
#' A data set containing the microbial taxonomy of six kingdoms from the Catalogue of Life. MO codes can be looked up using \code{\link{as.mo}}.
#' @inheritSection catalogue_of_life Catalogue of Life
#' @format A \code{\link{data.frame}} with 69,454 observations and 16 variables:
#' \describe{
#'   \item{\code{mo}}{ID of microorganism as used by this package}
#'   \item{\code{col_id}}{Catalogue of Life ID}
#'   \item{\code{fullname}}{Full name, like \code{"Escherichia coli"}}
#'   \item{\code{kingdom}, \code{phylum}, \code{class}, \code{order}, \code{family}, \code{genus}, \code{species}, \code{subspecies}}{Taxonomic rank of the microorganism}
#'   \item{\code{rank}}{Text of the taxonomic rank of the microorganism, like \code{"species"} or \code{"genus"}}
#'   \item{\code{ref}}{Author(s) and year of concerning scientific publication}
#'   \item{\code{species_id}}{ID of the species as used by the Catalogue of Life}
#'   \item{\code{source}}{Either "CoL", "DSMZ" (see Source) or "manually added"}
#'   \item{\code{prevalence}}{Prevalence of the microorganism, see \code{?as.mo}}
#' }
#' @details Manually added were:
#' \itemize{
#'   \item{11 entries of \emph{Streptococcus} (beta-haemolytic: groups A, B, C, D, F, G, H, K and unspecified; other: viridans, milleri)}
#'   \item{2 entries of \emph{Staphylococcus} (coagulase-negative [CoNS] and coagulase-positive [CoPS])}
#'   \item{3 entries of \emph{Trichomonas} (\emph{Trichomonas vaginalis}, and its family and genus)}
#'   \item{5 other 'undefined' entries (unknown, unknown Gram negatives, unknown Gram positives, unknown yeast and unknown fungus)}
#'   \item{9,460 species from the DSMZ (Deutsche Sammlung von Mikroorganismen und Zellkulturen) since the DSMZ contain the latest taxonomic information based on recent publications}
#' }
#' @section About the records from DSMZ (see source):
#' Names of prokaryotes are defined as being validly published by the International Code of Nomenclature of Bacteria. Validly published are all names which are included in the Approved Lists of Bacterial Names and the names subsequently published in the International Journal of Systematic Bacteriology (IJSB) and, from January 2000, in the International Journal of Systematic and Evolutionary Microbiology (IJSEM) as original articles or in the validation lists.
#'
#' From: \url{https://www.dsmz.de/support/bacterial-nomenclature-up-to-date-downloads/readme.html}
#' @source Catalogue of Life: Annual Checklist (public online taxonomic database), \url{http://www.catalogueoflife.org} (check included annual version with \code{\link{catalogue_of_life_version}()}).
#'
#' Leibniz Institute DSMZ-German Collection of Microorganisms and Cell Cultures, Germany, Prokaryotic Nomenclature Up-to-Date, \url{http://www.dsmz.de/bacterial-diversity/prokaryotic-nomenclature-up-to-date} (check included version with \code{\link{catalogue_of_life_version}()}).
#' @inheritSection AMR Read more on our website!
#' @seealso \code{\link{as.mo}}, \code{\link{mo_property}}, \code{\link{microorganisms.codes}}
"microorganisms"

catalogue_of_life <- list(
  year = 2018,
  version = "Catalogue of Life: {year} Annual Checklist",
  url_CoL = "http://www.catalogueoflife.org/annual-checklist/{year}/",
  url_DSMZ = "https://www.dsmz.de/services/online-tools/prokaryotic-nomenclature-up-to-date/prokaryotic-nomenclature-up-to-date/genus-search",
  yearmonth_DSMZ = "August 2019"
)

#' Data set with previously accepted taxonomic names
#'
#' A data set containing old (previously valid or accepted) taxonomic names according to the Catalogue of Life. This data set is used internally by \code{\link{as.mo}}.
#' @inheritSection catalogue_of_life Catalogue of Life
#' @format A \code{\link{data.frame}} with 24,246 observations and 5 variables:
#' \describe{
#'   \item{\code{col_id}}{Catalogue of Life ID that was originally given}
#'   \item{\code{col_id_new}}{New Catalogue of Life ID that responds to an entry in the \code{\link{microorganisms}} data set}
#'   \item{\code{fullname}}{Old full taxonomic name of the microorganism}
#'   \item{\code{ref}}{Author(s) and year of concerning scientific publication}
#'   \item{\code{prevalence}}{Prevalence of the microorganism, see \code{?as.mo}}
#' }
#' @source Catalogue of Life: Annual Checklist (public online taxonomic database), \url{http://www.catalogueoflife.org} (check included annual version with \code{\link{catalogue_of_life_version}()}).
#' @inheritSection AMR Read more on our website!
#' @seealso \code{\link{as.mo}} \code{\link{mo_property}} \code{\link{microorganisms}}
"microorganisms.old"

#' Translation table for common microorganism codes
#'
#' A data set containing commonly used codes for microorganisms, from laboratory systems and WHONET. Define your own with \code{\link{set_mo_source}}.
#' @format A \code{\link{data.frame}} with 4,927 observations and 2 variables:
#' \describe{
#'   \item{\code{code}}{Commonly used code of a microorganism}
#'   \item{\code{mo}}{ID of the microorganism in the \code{\link{microorganisms}} data set}
#' }
#' @inheritSection catalogue_of_life Catalogue of Life
#' @inheritSection AMR Read more on our website!
#' @seealso \code{\link{as.mo}} \code{\link{microorganisms}}
"microorganisms.codes"

#' Data set with 2,000 blood culture isolates
#'
#' An anonymised data set containing 2,000 microbial blood culture isolates with their full antibiograms found 4 different hospitals in the Netherlands, between 2001 and 2017. This \code{data.frame} can be used to practice AMR analysis. For examples, please read \href{https://msberends.gitlab.io/AMR/articles/AMR.html}{the tutorial on our website}.
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
#'   \item{\code{PEN:RIF}}{40 different antibiotics with class \code{rsi} (see \code{\link{as.rsi}}); these column names occur in \code{\link{antibiotics}} data set and can be translated with \code{\link{ab_name}}}
#' }
#' @inheritSection AMR Read more on our website!
"example_isolates"

#' Data set with 500 isolates - WHONET example
#'
#' This example data set has the exact same structure as an export file from WHONET. Such files can be used with this package, as this example data set shows. The data itself was based on our \code{\link{example_isolates}} data set.
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
#'   \item{\code{AMP_ND10:CIP_EE}}{27 different antibiotics. You can lookup the abbreviatons in the \code{\link{antibiotics}} data set, or use e.g. \code{\link{ab_name}("AMP")} to get the official name immediately. Before analysis, you should transform this to a valid antibiotic class, using \code{\link{as.rsi}}.}
#' }
#' @inheritSection AMR Read more on our website!
"WHONET"

#' Data set for RSI interpretation
#'
#' Data set to interpret MIC and disk diffusion to RSI values. Included guidelines are CLSI (2011-2019) and EUCAST (2011-2019). Use \code{\link{as.rsi}} to transform MICs or disks measurements to RSI values.
#' @format A \code{\link{data.frame}} with 11,559 observations and 9 variables:
#' \describe{
#'   \item{\code{guideline}}{Name of the guideline}
#'   \item{\code{mo}}{Microbial ID, see \code{\link{as.mo}}}
#'   \item{\code{ab}}{Antibiotic ID, see \code{\link{as.ab}}}
#'   \item{\code{ref_tbl}}{Info about where the guideline rule can be found}
#'   \item{\code{S_mic}}{Lowest MIC value that leads to "S"}
#'   \item{\code{R_mic}}{Highest MIC value that leads to "R"}
#'   \item{\code{dose_disk}}{Dose of the used disk diffusion method}
#'   \item{\code{S_disk}}{Lowest number of millimeters that leads to "S"}
#'   \item{\code{R_disk}}{Highest number of millimeters that leads to "R"}
#' }
#' @inheritSection AMR Read more on our website!
"rsi_translation"

# transforms data set to data.frame with only ASCII values, to comply with CRAN policies
dataset_UTF8_to_ASCII <- function(df) {
  trans <- function(vect) {
    iconv(vect, from = "UTF-8", to = "ASCII//TRANSLIT")
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (i in 1:NCOL(df)) {
    col <- df[, i]
    if (is.list(col)) {
      for (j in 1:length(col)) {
        col[[j]] <- trans(col[[j]])
      }
      df[, i] <- list(col)
    } else {
      if (is.factor(col)) {
        levels(col) <- trans(levels(col))
      } else if (is.character(col)) {
        col <- trans(col)
      } else {
        col
      }
      df[, i] <- col
    }
  }
  df
}
