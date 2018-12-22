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

# No export, no Rd
addin_insert_in <- function() {
  rstudioapi::insertText(" %in% ")
}

# No export, no Rd
addin_insert_like <- function() {
  rstudioapi::insertText(" %like% ")
}

# No export, no Rd
# works exactly like round(), but rounds `round(0.55, 1)` as 0.6
round2 <- function(x, digits = 0) {
  # https://stackoverflow.com/a/12688836/4575331
  (trunc((abs(x) * 10 ^ digits) + 0.5) / 10 ^ digits) * sign(x)
}

# No export, no Rd
percent <- function(x, round = 1, force_zero = FALSE, decimal.mark = getOption("OutDec"), ...) {

  decimal.mark.options <- getOption("OutDec")
  options(OutDec = ".")

  val <- round2(x, round + 2) # round up 0.5
  val <- round(x = val * 100, digits = round) # remove floating point error

  if (force_zero == TRUE) {
    if (any(val == as.integer(val) & !is.na(val))) {
      # add zeroes to all integers
      val[val == as.integer(as.character(val))] <- paste0(val[val == as.integer(val)], ".", strrep(0, round))
    }
    # add extra zeroes if needed
    val_decimals <- nchar(gsub(".*[.](.*)", "\\1", as.character(val)))
    val[val_decimals < round] <- paste0(val[val_decimals < round], strrep(0, max(0, round - val_decimals)))
  }
  pct <- base::paste0(val, "%")
  pct[pct %in% c("NA%", "NaN%")] <- NA_character_
  if (decimal.mark != ".") {
    pct <- gsub(".", decimal.mark, pct, fixed = TRUE)
  }
  options(OutDec = decimal.mark.options)
  pct
}

check_available_columns <- function(tbl, col.list, info = TRUE) {
  # check columns
  col.list <- col.list[!is.na(col.list) & !is.null(col.list)]
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
  stats::sd(x, na.rm = na.rm) / base::abs(base::mean(x, na.rm = na.rm))
}

# Coefficient of dispersion, or coefficient of quartile variation (CQV).
# (Bonett et al., 2006: Confidence interval for a coefficient of quartile variation).
cqv <- function(x, na.rm = TRUE) {
  fives <- stats::fivenum(x, na.rm = na.rm)
  (fives[4] - fives[2]) / (fives[4] + fives[2])
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
