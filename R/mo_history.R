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

# print successful as.mo coercions to file, not uncertain ones
#' @importFrom dplyr %>% filter
set_mo_history <- function(x, mo, force = FALSE) {
  file_location <- base::path.expand('~/.Rhistory_mo')
  if ((base::interactive() & mo != "UNKNOWN") | force == TRUE) {
    mo_hist <- read_mo_history(force = force)
    if (NROW(mo_hist[base::which(mo_hist$x == x & mo_hist$package_version == utils::packageVersion("AMR")),]) == 0) {
      base::write(x = c(x, mo, base::as.character(utils::packageVersion("AMR"))),
                  file = file_location,
                  ncolumns = 3,
                  append = TRUE,
                  sep = "\t")
    }
  }
  return(base::invisible())
}

get_mo_history <- function(x, force = FALSE) {
  file_read <- read_mo_history(force = force)
  if (base::is.null(file_read)) {
    NA
  } else {
    data.frame(x, stringsAsFactors = FALSE) %>%
      left_join(file_read, by = "x") %>%
      pull(mo)
  }
}

read_mo_history <- function(force = FALSE) {
  file_location <- base::path.expand('~/.Rhistory_mo')
  if (!base::file.exists(file_location) | (!base::interactive() & force == FALSE)) {
    return(NULL)
  }
  file_read <- utils::read.table(file = file_location,
                                 header = FALSE,
                                 sep = "\t",
                                 col.names = c("x", "mo", "package_version"),
                                 stringsAsFactors = FALSE)
  # Below: filter on current package version.
  # Future fullnames may even be replaced by new taxonomic names, so new versions of
  # the Catalogue of Life must not lead to data corruption.
  file_read[base::which(file_read$package_version == utils::packageVersion("AMR")), c("x", "mo")]
}

#' @rdname as.mo
#' @export
clean_mo_history <- function() {
  file_location <- base::path.expand('~/.Rhistory_mo')
  if (base::file.exists(file_location)) {
    base::unlink(file_location)
  }
}

