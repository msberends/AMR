# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This package is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This R package is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License version 2.0 for more details.             #
# ==================================================================== #

install_if_needed <- function(package_to_install) {
  package_path <- find.package(package_to_install, quiet = TRUE)

  if(length(package_path) == 0){
    # Only install if not present
    install.packages(package_to_install)
  }
}

ci_setup <- function() {
  install_if_needed("packrat")
  packrat::restore()
}

ci_check <- function() {
  install_if_needed("devtools")
  devtools::check()
}

ci_coverage <- function() {
  install_if_needed("covr")
  cc <- covr::package_coverage(type = c("tests", "examples"))
  covr::codecov(coverage = cc, token = "50ffa0aa-fee0-4f8b-a11d-8c7edc6d32ca")
  cat("Code coverage:", covr::percent_coverage(cc))
}

ci_pages <- function() {
  install_if_needed("pkgdown")
  pkgdown::build_site(examples = FALSE, override = list(destination = "public"))
}
