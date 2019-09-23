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
#' @param combine_IR logical to indicate whether values R and I should be summed
#' @param add_ab_group logical to indicate where the group of the antimicrobials must be included as a first column
#' @param FUN the function to call on the \code{mo} column to transform the microorganism IDs, defaults to \code{\link{mo_shortname}} 
#' @param ... argumments passed on to \code{FUN}
#' @inheritParams rsi_df
#' @inheritParams base::formatC
#' @importFrom dplyr rename
#' @importFrom tidyr spread
#' @importFrom clean freq
#' @details The function \code{format} calculates the resistance per bug-drug combination. Use \code{combine_IR = FALSE} (default) to test R vs. S+I and \code{combine_IR = TRUE} to test R+I vs. S. 
#' 
#' The language of the output can be overwritten with \code{options(AMR_locale)}, please see \link{translate}.
#' @export
#' @rdname bug_drug_combinations
#' @source \strong{M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition}, 2014, \emph{Clinical and Laboratory Standards Institute (CLSI)}. \url{https://clsi.org/standards/products/microbiology/documents/m39/}.
#' @inheritSection AMR Read more on our website!
#' @examples 
#' \donttest{
#' x <- bug_drug_combinations(example_isolates)
#' x
#' format(x)
#' 
#' # Use FUN to change to transformation of microorganism codes
#' x <- bug_drug_combinations(example_isolates, 
#'                            FUN = mo_gramstain)
#'                            
#' x <- bug_drug_combinations(example_isolates,
#'                            FUN = function(x) ifelse(x == "B_ESCHR_COLI",
#'                                                     "E. coli",
#'                                                     "Others"))
#' }
bug_drug_combinations <- function(x, 
                                  col_mo = NULL, 
                                  minimum = 30,
                                  FUN = mo_shortname,
                                  ...) {
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
    mutate(mo = x %>% pull(col_mo) %>% FUN(...)) %>% 
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
  
  structure(.Data = x, class = c("bug_drug_combinations", class(x)))
}

#' @importFrom dplyr everything rename
#' @importFrom tidyr spread
#' @exportMethod format.bug_drug_combinations
#' @export
#' @rdname bug_drug_combinations
format.bug_drug_combinations <- function(x, 
                                         combine_IR = FALSE, 
                                         add_ab_group = TRUE,
                                         decimal.mark = getOption("OutDec"),
                                         big.mark = ifelse(decimal.mark == ",", ".", ",")) {
  if (combine_IR == FALSE) {
    x$isolates <- x$R
  } else {
    x$isolates <- x$R + x$I
  }
  y <- x %>%
    mutate(txt = paste0(percent(isolates / total, force_zero = TRUE, decimal.mark = decimal.mark, big.mark = big.mark), 
                        " (", trimws(format(isolates, big.mark = big.mark)), "/", 
                        trimws(format(total, big.mark = big.mark)), ")")) %>% 
    select(ab, mo, txt) %>% 
    spread(mo, txt) %>%
    mutate_all(~ifelse(is.na(.), "", .)) %>% 
    mutate(ab_group = ab_group(ab),
           ab = paste0(ab_name(ab), " (", as.ab(ab), ", ", ab_atc(ab), ")")) %>% 
    select(ab_group, ab, everything()) %>% 
    arrange(ab_group, ab) %>% 
    mutate(ab_group = ifelse(ab_group != lag(ab_group) | is.na(lag(ab_group)), ab_group, ""))
  
  if (add_ab_group == FALSE) {
    y <- y %>% select(-ab_group) %>% rename("Antibiotic" = ab)
    colnames(y)[1] <- translate_AMR(colnames(y)[1], language = get_locale(), only_unknown = FALSE)
  } else {
    y <- y %>% rename("Group" = ab_group,
                      "Antibiotic" = ab)
    colnames(y)[1:2] <- translate_AMR(colnames(y)[1:2], language = get_locale(), only_unknown = FALSE)
  }
  y
}

#' @exportMethod print.bug_drug_combinations
#' @export
#' @importFrom crayon blue
print.bug_drug_combinations <- function(x, ...) {
  print(as.data.frame(x, stringsAsFactors = FALSE))
  message(blue("NOTE: Use 'format()' on this result to get a format that is ready for export or printing."))
}
