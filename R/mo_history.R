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

# print successful as.mo coercions to AMR environment
#' @importFrom dplyr %>% distinct filter
set_mo_history <- function(x, mo, uncertainty_level, force = FALSE) {
  if (base::interactive() | force == TRUE) {
    mo_hist <- read_mo_history(uncertainty_level = uncertainty_level, force = force)
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
      if (NROW(mo_hist[base::which(mo_hist$x == x[i] &
                                   mo_hist$uncertainty_level >= uncertainty_level &
                                   mo_hist$package_v == utils::packageVersion("AMR")),]) == 0) {
        tryCatch(
          assign(x = "mo_history",
                 value = rbind(mo_hist,
                               data.frame(
                                 x = x[i],
                                 mo = mo[i],
                                 uncertainty_level = uncertainty_level,
                                 package_v = base::as.character(utils::packageVersion("AMR")),
                                 stringsAsFactors = FALSE)),
                 envir = asNamespace("AMR")),
          error = function(e) invisible())
      }
    }
  }
  return(base::invisible())
}

get_mo_history <- function(x, uncertainty_level, force = FALSE) {
  history <- read_mo_history(uncertainty_level = uncertainty_level, force = force)
  if (base::is.null(history)) {
    NA
  } else {
    data.frame(x = toupper(x), stringsAsFactors = FALSE) %>%
      left_join(history, by = "x") %>%
      pull(mo)
  }
}

#' @importFrom dplyr %>% filter distinct
read_mo_history <- function(uncertainty_level = 2, force = FALSE, unfiltered = FALSE) {
  if ((!base::interactive() & force == FALSE)) {
    return(NULL)
  }
  uncertainty_level_param <- uncertainty_level

  history <- tryCatch(get("mo_history", envir = asNamespace("AMR")),
                        error = function(e) NULL)
  if (is.null(history)) {
    return(NULL)
  }
  # Below: filter on current package version.
  # Even current fullnames may be replaced by new taxonomic names, so new versions of
  # the Catalogue of Life must not lead to data corruption.

  if (unfiltered == FALSE) {
    history <- history %>%
      filter(package_v == as.character(utils::packageVersion("AMR")),
             # only take unknowns if uncertainty_level_param is higher
             ((mo == "UNKNOWN" & uncertainty_level_param == uncertainty_level) |
                (mo != "UNKNOWN" & uncertainty_level_param >= uncertainty_level))) %>%
      arrange(desc(uncertainty_level)) %>%
      distinct(x, mo, .keep_all = TRUE)
  }

  if (nrow(history) == 0) {
    NULL
  } else {
    history
  }
}

#' @rdname as.mo
#' @importFrom crayon red
#' @importFrom utils menu
#' @export
clean_mo_history <- function(...) {
  if (!is.null(read_mo_history())) {
    if (interactive() & !isTRUE(list(...)$force)) {
      q <- menu(title = paste("This will remove all",
                              format(nrow(read_mo_history(999, unfiltered = TRUE)), big.mark = ","),
                              "microbial IDs determined previously in this session. Are you sure?"),
                choices = c("Yes", "No"),
                graphics = FALSE)
      if (q != 1) {
        return(invisible())
      }
    }
    tryCatch(
      assign(x = "mo_history",
             value = NULL,
             envir = asNamespace("AMR")),
      error = function(e) invisible())
    cat(red("History removed."))
  }
}

