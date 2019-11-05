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

install_if_needed <- function(pkg, repos, quiet) {
  package_path <- find.package(pkg, quiet = quiet)
  if (length(package_path) == 0) {
    message("NOTE: pkg ", pkg, " missing, installing...")
    install.packages(pkg, repos = repos, quiet = quiet)
  }
}

gl_update_pkg_all <- function(repos = "https://cran.rstudio.com",
                              quiet = TRUE,
                              install_pkgdown = FALSE,
                              install_lintr = FALSE) {
  # update existing
  update.packages(ask = FALSE, repos = repos, quiet = quiet)

  install_if_needed(pkg = "devtools", repos = repos, quiet = quiet)
  if (install_pkgdown == TRUE) {
    install_if_needed(pkg = "pkgdown", repos = repos, quiet = quiet)
  }
  if (install_lintr == TRUE) {
    install_if_needed(pkg = "lintr", repos = repos, quiet = quiet)
  }
  install_if_needed(pkg = "cleaner", repos = repos, quiet = quiet)
  
  devtools::install_dev_deps(repos = repos, quiet = quiet, upgrade = TRUE)

  cat("INSTALLED:\n")
  instld <- as.data.frame(installed.packages())
  rownames(instld) <- NULL
  print(instld[, c("Package", "Version")])

  return(invisible(TRUE))
}
