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
#' @importFrom dplyr %>% distinct filter
set_mo_history <- function(x, mo, force = FALSE) {
  file_location <- base::path.expand('~/.Rhistory_mo')
  if (base::interactive() | force == TRUE) {
    mo_hist <- read_mo_history(force = force)
    df <- data.frame(x, mo, stringsAsFactors = FALSE) %>%
      distinct(x, .keep_all = TRUE) %>%
      filter(!is.na(x) & !is.na(mo))
    if (nrow(df) == 0) {
      return(base::invisible())
    }
    x <- toupper(df$x)
    mo <- df$mo
    for (i in 1:length(x)) {
      # save package version too, as both the as.mo() algorithm and the reference data set may change
      if (NROW(mo_hist[base::which(mo_hist$x == x[i] & mo_hist$package_version == utils::packageVersion("AMR")),]) == 0) {
        base::write(x = c(x[i], mo[i], base::as.character(utils::packageVersion("AMR"))),
                    file = file_location,
                    ncolumns = 3,
                    append = TRUE,
                    sep = "\t")
      }
    }
  }
  return(base::invisible())
}

get_mo_history <- function(x, force = FALSE) {
  file_read <- read_mo_history(force = force)
  if (base::is.null(file_read)) {
    NA
  } else {
    data.frame(x = toupper(x), stringsAsFactors = FALSE) %>%
      left_join(file_read, by = "x") %>%
      pull(mo)
  }
}

#' @importFrom dplyr %>% filter distinct
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
  # Even current fullnames may be replaced by new taxonomic names, so new versions of
  # the Catalogue of Life must not lead to data corruption.
  file_read %>%
    filter(package_version == utils::packageVersion("AMR")) %>%
    distinct(x, mo, .keep_all = TRUE)
}

#' @rdname as.mo
#' @export
clean_mo_history <- function() {
  file_location <- base::path.expand('~/.Rhistory_mo')
  if (base::file.exists(file_location)) {
    base::unlink(file_location)
  }
}

