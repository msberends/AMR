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

install_if_needed <- function(pkg, repos, quiet) {
  package_path <- find.package(pkg, quiet = quiet)
  if (length(package_path) == 0) {
    message("NOTE: pkg ", pkg, " missing, installing...")
    install.packages(pkg, repos = repos, quiet = quiet)
  }
}

gl_update_pkg_all <- function(repos = "https://cran.rstudio.com",
                              quiet = TRUE,
                              install_pkgdown = FALSE) {
  # update existing
  update.packages(ask = FALSE, repos = repos, quiet = quiet)

  install_if_needed(pkg = "devtools", repos = repos, quiet = quiet)
  if (install_pkgdown == TRUE) {
    install_if_needed(pkg = "pkgdown", repos = repos, quiet = quiet)
  }

  devtools::install_dev_deps(repos = repos, quiet = quiet, upgrade = TRUE)

  cat("\nINSTALLED:\n\n")
  instld <- as.data.frame(installed.packages())
  rownames(instld) <- NULL
  print(instld[, c("Package", "Version")])

  return(invisible(TRUE))

  # which ones are needed now?
  # pkg_needed <-

  # if (length(pkg_needed) > 0) {
  #  # install them
  #  for (i %in% 1:length(pkg_needed)) {
  #    install_if_needed(pkg = pkg_needed[i], repos = repos, quiet = quiet)
  #  }
  # }
}
