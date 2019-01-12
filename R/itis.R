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

#' ITIS: Integrated Taxonomic Information System
#'
#' All taxonomic names of all microorganisms are included in this package, using the authoritative Integrated Taxonomic Information System (ITIS).
#' @section ITIS:
#' \if{html}{\figure{itis_logo.jpg}{options: height=60px style=margin-bottom:5px} \cr}
#' This package contains the \strong{complete microbial taxonomic data} (with all nine taxonomic ranks - from kingdom to subspecies) from the publicly available Integrated Taxonomic Information System (ITIS, \url{https://www.itis.gov}).
#'
#' All ~20,000 (sub)species from \strong{the taxonomic kingdoms Bacteria, Fungi and Protozoa are included in this package}, as well as all ~2,500 previously accepted names known to ITIS. Furthermore, the responsible authors and year of publication are available. This allows users to use authoritative taxonomic information for their data analysis on any microorganism, not only human pathogens. It also helps to quickly determine the Gram stain of bacteria, since all bacteria are classified into subkingdom Negibacteria or Posibacteria.
#'
#' ITIS is a partnership of U.S., Canadian, and Mexican agencies and taxonomic specialists [3].
#' @inheritSection AMR Read more on our website!
#' @name ITIS
#' @rdname ITIS
#' @examples
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
#' mo_subkingdom("E. coli")
#' # [1] "Negibacteria"
#'
#' mo_gramstain("E. coli") # based on subkingdom
#' # [1] "Gram negative"
#'
#' mo_ref("E. coli")
#' # [1] "Castellani and Chalmers, 1919"
#'
#' # Do not get mistaken - the package only includes microorganisms
#' mo_phylum("C. elegans")
#' # [1] "Cyanobacteria"                   # Bacteria?!
#' mo_fullname("C. elegans")
#' # [1] "Chroococcus limneticus elegans"  # Because a microorganism was found
NULL
