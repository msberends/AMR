# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
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

format_included_data_number <- function(data) {
  if (is.data.frame(data)) {
    n <- nrow(data)
  } else {
    n <- length(unique(data))
  }
  if (n > 10000) {
    rounder <- -3 # round on thousands
  } else if (n > 1000) {
    rounder <- -2 # round on hundreds
  } else {
    rounder <- -1 # round on tens
  }
  paste0("~", format(round(n, rounder), decimal.mark = ".", big.mark = ","))
}

#' The Catalogue of Life
#'
#' This package contains the complete taxonomic tree (last updated: `r CATALOGUE_OF_LIFE$yearmonth_LPSN`) of almost all microorganisms from the authoritative and comprehensive Catalogue of Life (CoL), supplemented with data from the List of Prokaryotic names with Standing in Nomenclature (LPSN).
#' @section Catalogue of Life:
#' \if{html}{\figure{logo_col.png}{options: height="40" style=margin-bottom:"5"} \cr}
#' This package contains the complete taxonomic tree of almost all microorganisms (`r format_included_data_number(microorganisms)` species) from the authoritative and comprehensive Catalogue of Life (CoL, <http://www.catalogueoflife.org>). The CoL is the most comprehensive and authoritative global index of species currently available. Nonetheless, we supplemented the CoL data with data from the List of Prokaryotic names with Standing in Nomenclature (LPSN, [lpsn.dsmz.de](https://lpsn.dsmz.de)). This supplementation is needed until the [CoL+ project](https://github.com/CatalogueOfLife/general) is finished, which we await.
#'
#' [Click here][catalogue_of_life] for more information about the included taxa. Check which versions of the CoL and LPSN were included in this package with [catalogue_of_life_version()].
#' @section Included Taxa:
#' Included are:
#' - All `r format_included_data_number(microorganisms[which(microorganisms$kingdom %in% c("Archeae", "Bacteria", "Chromista", "Protozoa")), ])` (sub)species from the kingdoms of Archaea, Bacteria, Chromista and Protozoa
#' - All `r format_included_data_number(microorganisms[which(microorganisms$kingdom == "Fungi" & microorganisms$order %in% c("Eurotiales", "Microascales", "Mucorales", "Onygenales", "Pneumocystales", "Saccharomycetales", "Schizosaccharomycetales", "Tremellales")), ])` (sub)species from these orders of the kingdom of Fungi: Eurotiales, Microascales, Mucorales, Onygenales, Pneumocystales, Saccharomycetales, Schizosaccharomycetales and Tremellales, as well as `r format_included_data_number(microorganisms[which(microorganisms$kingdom == "Fungi" & !microorganisms$order %in% c("Eurotiales", "Microascales", "Mucorales", "Onygenales", "Pneumocystales", "Saccharomycetales", "Schizosaccharomycetales", "Tremellales")), ])` other fungal (sub)species. The kingdom of Fungi is a very large taxon with almost 300,000 different (sub)species, of which most are not microbial (but rather macroscopic, like mushrooms). Because of this, not all fungi fit the scope of this package and including everything would tremendously slow down our algorithms too. By only including the aforementioned taxonomic orders, the most relevant fungi are covered (such as all species of *Aspergillus*, *Candida*, *Cryptococcus*, *Histplasma*, *Pneumocystis*, *Saccharomyces* and *Trichophyton*).
#' - All `r format_included_data_number(microorganisms[which(microorganisms$kingdom == "Animalia"), ])` (sub)species from `r format_included_data_number(microorganisms[which(microorganisms$kingdom == "Animalia"), "genus"])` other relevant genera from the kingdom of Animalia (such as *Strongyloides* and *Taenia*)
#' - All `r format_included_data_number(microorganisms.old)` previously accepted names of all included (sub)species (these were taxonomically renamed)
#' - The complete taxonomic tree of all included (sub)species: from kingdom to subspecies
#' - The responsible author(s) and year of scientific publication
#'
#' The Catalogue of Life (<http://www.catalogueoflife.org>) is the most comprehensive and authoritative global index of species currently available. It holds essential information on the names, relationships and distributions of over 1.9 million species. The Catalogue of Life is used to support the major biodiversity and conservation information services such as the Global Biodiversity Information Facility (GBIF), Encyclopedia of Life (EoL) and the International Union for Conservation of Nature Red List. It is recognised by the Convention on Biological Diversity as a significant component of the Global Taxonomy Initiative and a contribution to Target 1 of the Global Strategy for Plant Conservation.
#'
#' The syntax used to transform the original data to a cleansed \R format, can be found here: <https://github.com/msberends/AMR/blob/main/data-raw/reproduction_of_microorganisms.R>.
#' @inheritSection AMR Read more on Our Website!
#' @name catalogue_of_life
#' @rdname catalogue_of_life
#' @seealso Data set [microorganisms] for the actual data. \cr
#' Function [as.mo()] to use the data for intelligent determination of microorganisms.
#' @examples
#' # Get version info of included data set
#' catalogue_of_life_version()
#'
#'
#' # Get a note when a species was renamed
#' mo_shortname("Chlamydophila psittaci")
#' # Note: 'Chlamydophila psittaci' (Everett et al., 1999) was renamed back to
#' #       'Chlamydia psittaci' (Page, 1968)
#' #> [1] "C. psittaci"
#'
#' # Get any property from the entire taxonomic tree for all included species
#' mo_class("E. coli")
#' #> [1] "Gammaproteobacteria"
#'
#' mo_family("E. coli")
#' #> [1] "Enterobacteriaceae"
#'
#' mo_gramstain("E. coli") # based on kingdom and phylum, see ?mo_gramstain
#' #> [1] "Gram-negative"
#'
#' mo_ref("E. coli")
#' #> [1] "Castellani et al., 1919"
#'
#' # Do not get mistaken - this package is about microorganisms
#' mo_kingdom("C. elegans")
#' #> [1] "Fungi"                 # Fungi?!
#' mo_name("C. elegans")
#' #> [1] "Cladosporium elegans"  # Because a microorganism was found
NULL

#' Version info of included Catalogue of Life
#'
#' This function returns information about the included data from the Catalogue of Life.
#' @seealso [microorganisms]
#' @details For LPSN, see [microorganisms].
#' @return a [list], which prints in pretty format
#' @inheritSection catalogue_of_life Catalogue of Life
#' @inheritSection AMR Read more on Our Website!
#' @export
catalogue_of_life_version <- function() {
  
  check_dataset_integrity()
  
  # see the `CATALOGUE_OF_LIFE` list in R/globals.R
  lst <- list(CoL =
                list(version = gsub("{year}", CATALOGUE_OF_LIFE$year, CATALOGUE_OF_LIFE$version, fixed = TRUE),
                     url = gsub("{year}", CATALOGUE_OF_LIFE$year, CATALOGUE_OF_LIFE$url_CoL, fixed = TRUE),
                     n = nrow(pm_filter(microorganisms, source == "CoL"))),
              LPSN =
                list(version = "List of Prokaryotic names with Standing in Nomenclature",
                     url = CATALOGUE_OF_LIFE$url_LPSN,
                     yearmonth = CATALOGUE_OF_LIFE$yearmonth_LPSN,
                     n = nrow(pm_filter(microorganisms, source == "LPSN"))),
              total_included =
                list(
                  n_total_species = nrow(microorganisms),
                  n_total_synonyms = nrow(microorganisms.old)))
  
  set_clean_class(lst,
                  new_class = c("catalogue_of_life_version", "list"))
}

#' @method print catalogue_of_life_version
#' @export
#' @noRd
print.catalogue_of_life_version <- function(x, ...) {
  cat(paste0(font_bold("Included in this AMR package (v", utils::packageDescription("AMR")$Version, ") are:\n\n", collapse = ""),
             font_underline(x$CoL$version), "\n",
             "  Available at: ", font_blue(x$CoL$url), "\n",
             "  Number of included microbial species: ", format(x$CoL$n, big.mark = ","), "\n",
             font_underline(paste0(x$LPSN$version, " (",
                                   x$LPSN$yearmonth, ")")), "\n",
             "  Available at: ", font_blue(x$LPSN$url), "\n",
             "  Number of included bacterial species: ", format(x$LPSN$n, big.mark = ","), "\n\n",
             "=> Total number of species included:  ", format(x$total_included$n_total_species, big.mark = ","), "\n",
             "=> Total number of synonyms included: ", format(x$total_included$n_total_synonyms, big.mark = ","), "\n\n",
             "See for more info ", font_grey_bg("`?microorganisms`"), " and ", font_grey_bg("`?catalogue_of_life`"), ".\n"))
}
