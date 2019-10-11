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
#' Determine antimicrobial resistance (AMR) of all bug-drug combinations in your data set where at least 30 (default) isolates are available per species. Use \code{format} on the result to prettify it to a publicable/printable format, see Examples.
#' @inheritParams eucast_rules
#' @param combine_IR logical to indicate whether values R and I should be summed
#' @param add_ab_group logical to indicate where the group of the antimicrobials must be included as a first column
#' @param remove_intrinsic_resistant logical to indicate that rows with 100\% resistance for all tested antimicrobials must be removed from the table
#' @param FUN the function to call on the \code{mo} column to transform the microorganism IDs, defaults to \code{\link{mo_shortname}} 
#' @param translate_ab a character of length 1 containing column names of the \code{\link{antibiotics}} data set
#' @param ... arguments passed on to \code{FUN}
#' @inheritParams rsi_df
#' @inheritParams base::formatC
#' @importFrom dplyr %>% rename group_by select mutate filter summarise ungroup
#' @importFrom tidyr spread
# @importFrom clean freq percentage
#' @details The function \code{format} calculates the resistance per bug-drug combination. Use \code{combine_IR = FALSE} (default) to test R vs. S+I and \code{combine_IR = TRUE} to test R+I vs. S. 
#' 
#' The language of the output can be overwritten with \code{options(AMR_locale)}, please see \link{translate}.
#' @export
#' @rdname bug_drug_combinations
#' @return The function \code{bug_drug_combinations} returns a \code{data.frame} with columns "mo", "ab", "S", "I", "R" and "total".
#' @source \strong{M39 Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data, 4th Edition}, 2014, \emph{Clinical and Laboratory Standards Institute (CLSI)}. \url{https://clsi.org/standards/products/microbiology/documents/m39/}.
#' @inheritSection AMR Read more on our website!
#' @examples 
#' \donttest{
#' x <- bug_drug_combinations(example_isolates)
#' x
#' format(x, translate_ab = "name (atc)")
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
    as.data.frame(stringsAsFactors = FALSE) %>% 
    mutate(mo = x %>% 
             pull(col_mo) %>% 
             FUN(...)) %>% 
    group_by(mo) %>% 
    select_if(is.rsi) %>% 
    gather("ab", "value", -mo) %>% 
    group_by(mo, ab) %>% 
    summarise(S = sum(value == "S", na.rm = TRUE),
              I = sum(value == "I", na.rm = TRUE),
              R = sum(value == "R", na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(total = S + I + R) %>%
    as.data.frame(stringsAsFactors = FALSE)

  structure(.Data = x, class = c("bug_drug_combinations", class(x)))
}

#' @importFrom dplyr everything rename %>% ungroup group_by summarise mutate_all arrange everything lag
#' @importFrom tidyr spread
#' @exportMethod format.bug_drug_combinations
#' @export
#' @rdname bug_drug_combinations
format.bug_drug_combinations <- function(x,
                                         translate_ab = "name (ab, atc)",
                                         language = get_locale(),
                                         minimum = 30,
                                         combine_SI = TRUE,
                                         combine_IR = FALSE,
                                         add_ab_group = TRUE,
                                         remove_intrinsic_resistant = FALSE,
                                         decimal.mark = getOption("OutDec"),
                                         big.mark = ifelse(decimal.mark == ",", ".", ","),
                                         ...) {
  x <- x %>% filter(total >= minimum)
  
  if (remove_intrinsic_resistant == TRUE) {
    x <- x %>% filter(R != total)
  }
  if (combine_SI == TRUE | combine_IR == FALSE) {
    x$isolates <- x$R
  } else {
    x$isolates <- x$R + x$I
  }
  
  give_ab_name <- function(ab, format, language) {
    format <- tolower(format)
    ab_txt <- rep(format, length(ab))
    for (i in seq_len(length(ab_txt))) {
      ab_txt[i] <- gsub("ab", ab[i], ab_txt[i])
      ab_txt[i] <- gsub("cid", ab_cid(ab[i]), ab_txt[i])
      ab_txt[i] <- gsub("group", ab_group(ab[i], language = language), ab_txt[i])
      ab_txt[i] <- gsub("atc_group1", ab_atc_group1(ab[i], language = language), ab_txt[i])
      ab_txt[i] <- gsub("atc_group2", ab_atc_group2(ab[i], language = language), ab_txt[i])
      ab_txt[i] <- gsub("atc", ab_atc(ab[i]), ab_txt[i])
      ab_txt[i] <- gsub("name", ab_name(ab[i], language = language), ab_txt[i])
      ab_txt[i]
    }
    ab_txt
  }

  y <- x %>%
    mutate(ab = as.ab(ab),
           ab_txt = give_ab_name(ab = ab, format = translate_ab, language = language)) %>% 
    group_by(ab, ab_txt, mo) %>% 
    summarise(isolates = sum(isolates, na.rm = TRUE),
              total = sum(total, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(txt = paste0(percentage(isolates / total, decimal.mark = decimal.mark, big.mark = big.mark), 
                        " (", trimws(format(isolates, big.mark = big.mark)), "/", 
                        trimws(format(total, big.mark = big.mark)), ")")) %>% 
    select(ab, ab_txt, mo, txt) %>% 
    spread(mo, txt) %>%
    mutate_all(~ifelse(is.na(.), "", .)) %>% 
    mutate(ab_group = ab_group(ab, language = language),
           ab_txt) %>% 
    select(ab_group, ab_txt, everything(), -ab) %>% 
    arrange(ab_group, ab_txt) %>% 
    mutate(ab_group = ifelse(ab_group != lag(ab_group) | is.na(lag(ab_group)), ab_group, ""))
  
  if (add_ab_group == FALSE) {
    y <- y %>% select(-ab_group) %>% rename("Drug" = ab_txt)
    colnames(y)[1] <- translate_AMR(colnames(y)[1], language = get_locale(), only_unknown = FALSE)
  } else {
    y <- y %>% rename("Group" = ab_group,
                      "Drug" = ab_txt)
    colnames(y)[1:2] <- translate_AMR(colnames(y)[1:2], language = get_locale(), only_unknown = FALSE)
  }
  y
}

#' @exportMethod print.bug_drug_combinations
#' @export
#' @importFrom crayon blue
print.bug_drug_combinations <- function(x, ...) {
  print(as.data.frame(x, stringsAsFactors = FALSE))
  message(blue("NOTE: Use 'format()' on this result to get a publicable/printable format."))
}
