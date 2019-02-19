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

#' The Catalogue of Life
#'
#' This package contains the complete taxonomic tree of almost all microorganisms from the authoritative and comprehensive Catalogue of Life.
#' @section Catalogue of Life:
#' \if{html}{\figure{logo_col.png}{options: height=60px style=margin-bottom:5px} \cr}
#' This package contains the complete taxonomic tree of almost all microorganisms from the authoritative and comprehensive Catalogue of Life (\url{http://www.catalogueoflife.org}). This data is updated annually - check the included version with \code{\link{catalogue_of_life_version}}.
#'
#' Included are:
#' \itemize{
#'   \item{All ~55,000 (sub)species from the kingdoms of Archaea, Bacteria, Protozoa and Viruses}
#'   \item{All ~3,000 (sub)species from these orders of the kingdom of Fungi: Eurotiales, Onygenales, Pneumocystales, Saccharomycetales and Schizosaccharomycetales. The kingdom of Fungi is a very large taxon with almost 300,000 different species, of which most are not microbial. Including everything tremendously slows down our algortihms, and not all fungi fit the scope of this package. By only including the aforementioned taxonomic orders, the most relevant species are covered (like genera \emph{Aspergillus}, \emph{Candida}, \emph{Pneumocystis}, \emph{Saccharomyces} and \emph{Trichophyton}).}
#'   \item{All ~15,000 previously accepted names of (sub)species that have been taxonomically renamed}
#'   \item{The complete taxonomic tree of all included (sub)species: from kingdom to subspecies}
#'   \item{The responsible author(s) and year of scientific publication}
#' }
#'
#' The Catalogue of Life (\url{http://www.catalogueoflife.org}) is the most comprehensive and authoritative global index of species currently available. It holds essential information on the names, relationships and distributions of over 1.6 million species. The Catalogue of Life is used to support the major biodiversity and conservation information services such as the Global Biodiversity Information Facility (GBIF), Encyclopedia of Life (EoL) and the International Union for Conservation of Nature Red List. It is recognised by the Convention on Biological Diversity as a significant component of the Global Taxonomy Initiative and a contribution to Target 1 of the Global Strategy for Plant Conservation.
#'
#' The syntax used to transform the original data to a cleansed R format, can be found here: \url{https://gitlab.com/msberends/AMR/blob/master/reproduction_of_microorganisms.R}.
#' @inheritSection AMR Read more on our website!
#' @name catalogue_of_life
#' @rdname catalogue_of_life
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
#' mo_phylum("C. elegans")
#' # [1] "Cyanobacteria"                   # Bacteria?!
#' mo_fullname("C. elegans")
#' # [1] "Chroococcus limneticus elegans"  # Because a microorganism was found
NULL
