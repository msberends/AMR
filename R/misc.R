# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

#' Pattern Matching
#'
#' Convenience function to compare a vector with a pattern, like \code{\link[base]{grep}}. It always returns a \code{logical} vector and is always case-insensitive.
#' @inheritParams base::grep
#' @return A \code{logical} vector
#' @name like
#' @rdname like
#' @export
#' @source Inherited from the \href{https://github.com/Rdatatable/data.table/blob/master/R/like.R}{\code{like} function from the \code{data.table} package}, but made it case insensitive at default.
#' @examples
#' library(dplyr)
#' # get unique occurences of bacteria whose name start with 'Ent'
#' septic_patients %>%
#'   left_join_microorganisms() %>%
#'   filter(fullname %like% '^Ent') %>%
#'   pull(fullname) %>%
#'   unique()
"%like%" <- function(x, pattern) {
  if (length(pattern) > 1) {
    pattern <- pattern[1]
    warning('only the first element of argument `pattern` used for `%like%`', call. = FALSE)
  }
  if (is.factor(x)) {
    as.integer(x) %in% base::grep(pattern, levels(x), ignore.case = TRUE)
  } else {
    base::grepl(pattern, x, ignore.case = TRUE)
  }
}

# No export, no Rd
percent <- function(x, round = 1, force_zero = FALSE, ...) {
  val <- base::round(x * 100, digits = round)
  if (force_zero == TRUE & any(val == as.integer(val) & !is.na(val))) {
    val[val == as.integer(val)] <- paste0(val[val == as.integer(val)], ".", strrep(0, round))
  }
  pct <- base::paste0(val, "%")
  pct[pct == "NA%"] <- NA_character_
  pct
}

check_available_columns <- function(tbl, col.list, info = TRUE) {
  # check columns
  col.list <- col.list[!is.na(col.list)]
  names(col.list) <- col.list
  col.list.bak <- col.list
  # are they available as upper case or lower case then?
  for (i in 1:length(col.list)) {
    if (toupper(col.list[i]) %in% colnames(tbl)) {
      col.list[i] <- toupper(col.list[i])
    } else if (tolower(col.list[i]) %in% colnames(tbl)) {
      col.list[i] <- tolower(col.list[i])
    } else if (!col.list[i] %in% colnames(tbl)) {
      col.list[i] <- NA
    }
  }
  if (!all(col.list %in% colnames(tbl))) {
    if (info == TRUE) {
      warning('These columns do not exist and will be ignored: ',
              col.list.bak[!(col.list %in% colnames(tbl))] %>% toString(),
              '.\nTHIS MAY STRONGLY INFLUENCE THE OUTCOME.',
              immediate. = TRUE,
              call. = FALSE)
    }
  }
  col.list
}

# Coefficient of variation (CV)
cv <- function(x, na.rm = TRUE) {
  cv.x <- sd(x, na.rm = na.rm) / abs(mean(x, na.rm = na.rm))
  cv.x
}

# Coefficient of dispersion, or coefficient of quartile variation (CQV).
# (Bonett et al., 2006: Confidence interval for a coefficient of quartile variation).
cqv <- function(x, na.rm = TRUE) {
  cqv.x <-
    (quantile(x, 0.75, na.rm = na.rm, type = 6) - quantile(x, 0.25, na.rm = na.rm, type = 6)) /
    (quantile(x, 0.75, na.rm = na.rm, type = 6) + quantile(x, 0.25, na.rm = na.rm, type = 6))
  unname(cqv.x)
}

# show bytes as kB/MB/GB
# size_humanreadable(123456) # 121 kB
# size_humanreadable(12345678) # 11.8 MB
size_humanreadable <- function(bytes, decimals = 1) {
  bytes <- bytes %>% as.double()
  # Adapted from:
  # http://jeffreysambells.com/2012/10/25/human-readable-filesize-php
  size <- c('B','kB','MB','GB','TB','PB','EB','ZB','YB')
  factor <- floor((nchar(bytes) - 1) / 3)
  # added slight improvement; no decimals for B and kB:
  decimals <- rep(decimals, length(bytes))
  decimals[size[factor + 1] %in% c('B', 'kB')] <- 0

  out <- paste(sprintf(paste0("%.", decimals, "f"), bytes / (1024 ^ factor)), size[factor + 1])
  out
}
