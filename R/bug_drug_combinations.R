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

#' Determine bug-drug combinations
#' 
#' Determine antimicrobial resistance (AMR) of all bug-drug combinations in your data set where at least 30 (default) isolates are available per species. Use \code{format} on the result to prettify it to a printable format, see Examples.
#' @inheritParams eucast_rules
#' @inheritParams rsi_df
#' @importFrom dplyr rename
#' @importFrom tidyr spread
#' @importFrom clean freq
#' @export
#' @source \strong{M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition}, 2014, \emph{Clinical and Laboratory Standards Institute (CLSI)}. \url{https://clsi.org/standards/products/microbiology/documents/m39/}.
#' @inheritSection AMR Read more on our website!
#' @examples 
#' \donttest{
#' x <- bug_drug_combinations(septic_patients)
#' x
#' format(x)
#' }
bug_drug_combinations <- function(x, col_mo = NULL, minimum = 30) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }
  
  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo")
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }
  
  x <- x %>%
    mutate(col_mo = x %>% pull(col_mo)) %>% 
    filter(mo %in% (clean::freq(mo) %>%
                      filter(count >= minimum) %>% 
                      pull(item))) %>%
    group_by(mo) %>% 
    AMR::rsi_df(translate_ab = FALSE, combine_SI = FALSE) %>% 
    select(-value) %>%
    spread(interpretation, isolates) %>%
    mutate(total = S + I + R) %>%
    filter(total >= minimum) %>% 
    rename(ab = antibiotic)
  
  structure(.Data = x, class = c("bugdrug", class(x)))
}

#' @importFrom dplyr everything rename
#' @importFrom tidyr spread
#' @exportMethod format.bugdrug
#' @export
format.bugdrug <- function(x, combine_SI = TRUE, add_ab_group = TRUE, ...) {
  if (combine_SI == TRUE) {
    x$isolates <- x$R
  } else {
    x$isolates <- x$R + x$I
  }
  y <- x %>%
    mutate(mo = mo_name(mo),
           txt = paste0(percent(isolates / total, force_zero = TRUE), 
                        " (", trimws(format(isolates, big.mark = ",")), "/", 
                        trimws(format(total, big.mark = ",")), ")")) %>% 
    select(ab, mo, txt) %>% 
    spread(mo, txt) %>%
    mutate_all(~ifelse(is.na(.), "", .)) %>% 
    mutate(ab = paste0(ab_name(ab), " (", as.ab(ab), ", ", ab_atc(ab), ")"),
           ab_group = ab_group(ab)) %>% 
    select(ab_group, ab, everything()) %>% 
    arrange(ab_group, ab) %>% 
    mutate(ab_group = ifelse(ab_group != lag(ab_group) | is.na(lag(ab_group)), ab_group, ""))
  
  if (add_ab_group == FALSE) {
    y <- y %>% select(-ab_group)
  }
  y <- y %>% rename("Group" = ab_group,
                    "Antibiotic" = ab)
  y
}

#' @exportMethod print.bugdrug
#' @export
#' @importFrom crayon blue
print.bugdrug <- function(x, ...) {
  print(as.data.frame(x, stringsAsFactors = FALSE))
  message(blue("NOTE: Use 'format()' on this result to get a format that is ready for export or printing."))
}
