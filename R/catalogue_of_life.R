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

#' The Catalogue of Life
#'
#' This package contains the complete taxonomic tree of almost all microorganisms from the authoritative and comprehensive Catalogue of Life.
#' @section Catalogue of Life:
#' \if{html}{\figure{logo_col.png}{options: height=40px style=margin-bottom:5px} \cr}
#' This package contains the complete taxonomic tree of almost all microorganisms (~70,000 species) from the authoritative and comprehensive Catalogue of Life (\url{http://www.catalogueoflife.org}). The Catalogue of Life is the most comprehensive and authoritative global index of species currently available.
#'
#' \link[=catalogue_of_life]{Click here} for more information about the included taxa. Check which version of the Catalogue of Life was included in this package with \code{\link{catalogue_of_life_version}()}.
#' @section Included taxa:
#' Included are:
#' \itemize{
#'   \item{All ~61,000 (sub)species from the kingdoms of Archaea, Bacteria, Chromista and Protozoa}
#'   \item{All ~8,500 (sub)species from these orders of the kingdom of Fungi: Eurotiales, Microascales, Mucorales, Onygenales, Pneumocystales, Saccharomycetales, Schizosaccharomycetales and Tremellales. The kingdom of Fungi is a very large taxon with almost 300,000 different (sub)species, of which most are not microbial (but rather macroscopic, like mushrooms). Because of this, not all fungi fit the scope of this package and including everything would tremendously slow down our algorithms too. By only including the aforementioned taxonomic orders, the most relevant fungi are covered (like all species of \emph{Aspergillus}, \emph{Candida}, \emph{Cryptococcus}, \emph{Histplasma}, \emph{Pneumocystis}, \emph{Saccharomyces} and \emph{Trichophyton}).}
#'   \item{All ~150 (sub)species from ~100 other relevant genera from the kingdom of Animalia (like \emph{Strongyloides} and \emph{Taenia})}
#'   \item{All ~23,000 previously accepted names of all included (sub)species (these were taxonomically renamed)}
#'   \item{The complete taxonomic tree of all included (sub)species: from kingdom to subspecies}
#'   \item{The responsible author(s) and year of scientific publication}
#' }
#'
#' The Catalogue of Life (\url{http://www.catalogueoflife.org}) is the most comprehensive and authoritative global index of species currently available. It holds essential information on the names, relationships and distributions of over 1.9 million species. The Catalogue of Life is used to support the major biodiversity and conservation information services such as the Global Biodiversity Information Facility (GBIF), Encyclopedia of Life (EoL) and the International Union for Conservation of Nature Red List. It is recognised by the Convention on Biological Diversity as a significant component of the Global Taxonomy Initiative and a contribution to Target 1 of the Global Strategy for Plant Conservation.
#'
#' The syntax used to transform the original data to a cleansed R format, can be found here: \url{https://gitlab.com/msberends/AMR/blob/master/data-raw/reproduction_of_microorganisms.R}.
#' @inheritSection AMR Read more on our website!
#' @name catalogue_of_life
#' @rdname catalogue_of_life
#' @seealso Data set \code{\link{microorganisms}} for the actual data. \cr
#' Function \code{\link{as.mo}()} to use the data for intelligent determination of microorganisms.
#' @examples
#' # Get version info of included data set
#' catalogue_of_life_version()
#'
#'
#' # Get a note when a species was renamed
#' mo_shortname("Chlamydia psittaci")
#' # Note: 'Chlamydia psittaci' (Page, 1968) was renamed
#' #       'Chlamydophila psittaci' (Everett et al., 1999)
#' # [1] "C. psittaci"
#'
#' # Get any property from the entire taxonomic tree for all included species
#' mo_class("E. coli")
#' # [1] "Gammaproteobacteria"
#'
#' mo_family("E. coli")
#' # [1] "Enterobacteriaceae"
#'
#' mo_gramstain("E. coli") # based on kingdom and phylum, see ?mo_gramstain
#' # [1] "Gram negative"
#'
#' mo_ref("E. coli")
#' # [1] "Castellani et al., 1919"
#'
#' # Do not get mistaken - the package only includes microorganisms
#' mo_kingdom("C. elegans")
#' # [1] "Bacteria"                        # Bacteria?!
#' mo_name("C. elegans")
#' # [1] "Chroococcus limneticus elegans"  # Because a microorganism was found
NULL

#' Version info of included Catalogue of Life
#'
#' This function returns information about the included data from the Catalogue of Life.
#' @seealso \code{\link{microorganisms}}
#' @details For DSMZ, see \code{?microorganisms}.
#' @return a \code{list}, which prints in pretty format
#' @inheritSection catalogue_of_life Catalogue of Life
#' @inheritSection AMR Read more on our website!
#' @importFrom crayon bold underline
#' @importFrom dplyr filter
#' @export
#' @examples
#' library(dplyr)
#' microorganisms %>% freq(kingdom)
#' microorganisms %>% group_by(kingdom) %>% freq(phylum, nmax = NULL)
catalogue_of_life_version <- function() {
  # see the `catalogue_of_life` list in R/data.R
  lst <- list(catalogue_of_life =
                list(version = gsub("{year}", catalogue_of_life$year, catalogue_of_life$version, fixed = TRUE),
                     url = gsub("{year}", catalogue_of_life$year, catalogue_of_life$url_CoL, fixed = TRUE),
                     n = nrow(filter(AMR::microorganisms, source == "CoL"))),
              deutsche_sammlung_von_mikroorganismen_und_zellkulturen =
                list(version = "Prokaryotic Nomenclature Up-to-Date from DSMZ",
                     url = catalogue_of_life$url_DSMZ,
                     yearmonth = catalogue_of_life$yearmonth_DSMZ,
                     n = nrow(filter(AMR::microorganisms, source == "DSMZ"))),
              total_included =
                list(
                  n_total_species = nrow(AMR::microorganisms),
                  n_total_synonyms = nrow(AMR::microorganisms.old)))

  structure(.Data = lst,
            class = c("catalogue_of_life_version", "list"))
}

#' @exportMethod print.catalogue_of_life_version
#' @export
#' @noRd
print.catalogue_of_life_version <- function(x, ...) {
  lst <- x
  cat(paste0(bold("Included in this AMR package are:\n\n"),
             underline(lst$catalogue_of_life$version), "\n",
             "  Available at: ", lst$catalogue_of_life$url, "\n",
             "  Number of included species: ", format(lst$catalogue_of_life$n, big.mark = ","), "\n",
             underline(paste0(lst$deutsche_sammlung_von_mikroorganismen_und_zellkulturen$version, " (",
                              lst$deutsche_sammlung_von_mikroorganismen_und_zellkulturen$yearmonth, ")")), "\n",
             "  Available at: ", lst$deutsche_sammlung_von_mikroorganismen_und_zellkulturen$url, "\n",
             "  Number of included species: ", format(lst$deutsche_sammlung_von_mikroorganismen_und_zellkulturen$n, big.mark = ","), "\n\n",
             "=> Total number of species included:  ", format(lst$total_included$n_total_species, big.mark = ","), "\n",
             "=> Total number of synonyms included: ", format(lst$total_included$n_total_synonyms, big.mark = ","), "\n\n",
             "See for more info ?microorganisms and ?catalogue_of_life.\n"))
}
