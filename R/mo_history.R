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

mo_history_file <- file.path(file.path(system.file(package = "AMR"), "mo_history"), "mo_history.csv")

# print successful as.mo coercions to a options entry
#' @importFrom dplyr %>% distinct filter
set_mo_history <- function(x, mo, uncertainty_level, force = FALSE, disable = FALSE) {
  if (isTRUE(disable)) {
    return(base::invisible())
  }
  
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
        # # Not using the file system:        
        # tryCatch(options(mo_remembered_results = rbind(mo_hist,
        #                                                data.frame(
        #                                                  x = x[i],
        #                                                  mo = mo[i],
        #                                                  uncertainty_level = uncertainty_level,
        #                                                  package_v = base::as.character(utils::packageVersion("AMR")),
        #                                                  stringsAsFactors = FALSE))),
        #          error = function(e) base::invisible())
        # # don't remember more than 1,000 different input values
        # if (tryCatch(nrow(getOption("mo_remembered_results")), error = function(e) 1001) > 1000) {
        #   return(base::invisible())
        # }
        if (is.null(mo_hist)) {
          message(blue(paste0("NOTE: results are saved to ", mo_history_file, ".")))
        }
        tryCatch(write.csv(rbind(mo_hist,
                                 data.frame(
                                   x = x[i],
                                   mo = mo[i],
                                   uncertainty_level = uncertainty_level,
                                   package_v = base::as.character(utils::packageVersion("AMR")),
                                   stringsAsFactors = FALSE)),
                           file = mo_history_file, row.names = FALSE),
                 error = function(e) base::invisible())
      }
    }
  }
  return(base::invisible())
}

get_mo_history <- function(x, uncertainty_level, force = FALSE, disable = FALSE) {
  if (isTRUE(disable)) {
    return(to_class_mo(NA))
  }

  history <- read_mo_history(uncertainty_level = uncertainty_level, force = force)
  if (base::is.null(history)) {
    result <- NA
  } else {
    result <- data.frame(x = toupper(x), stringsAsFactors = FALSE) %>%
      left_join(history, by = "x") %>%
      pull(mo)
  }
  to_class_mo(result)
}

#' @importFrom dplyr %>% filter distinct
read_mo_history <- function(uncertainty_level = 2, force = FALSE, unfiltered = FALSE, disable = FALSE) {
  if (isTRUE(disable)) {
    return(NULL)
  }

  if ((!base::interactive() & force == FALSE)) {
    return(NULL)
  }
  uncertainty_level_param <- uncertainty_level
  
  # # Not using the file system:
  # history <- tryCatch(getOption("mo_remembered_results"),
  #                     error = function(e) NULL)
  history <- tryCatch(read.csv(mo_history_file, stringsAsFactors = FALSE),
                      warning = function(w) invisible(),
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
clear_mo_history <- function(...) {
  if (!is.null(read_mo_history())) {
    if (interactive() & !isTRUE(list(...)$force)) {
      q <- menu(title = paste("This will clear all",
                              format(nrow(read_mo_history(999, unfiltered = TRUE)), big.mark = ","),
                              "previously determined microbial IDs. Are you sure?"),
                choices = c("Yes", "No"),
                graphics = FALSE)
      if (q != 1) {
        return(invisible())
      }
    }
    # # Not using the file system:
    # success <- tryCatch(options(mo_remembered_results = NULL),
    #                     error = function(e) FALSE)
    success <- create_blank_mo_history()
    if (!isFALSE(success)) {
      cat(red(paste("File", mo_history_file, "cleared.")))
    }
  }
}

create_blank_mo_history <- function() {
  tryCatch(
    write.csv(x = data.frame(x = character(0),
                             mo = character(0),
                             uncertainty_level = integer(0),
                             package_v = character(0),
                             stringsAsFactors = FALSE),
              file = mo_history_file),
    warning = function(w) invisible(),
    error = function(e) TRUE)
}
