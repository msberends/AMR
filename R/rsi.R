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

#' Class 'rsi'
#'
#' Interpret MIC values according to EUCAST or CLSI, or clean up existing RSI values. This transforms the input to a new class \code{rsi}, which is an ordered factor with levels \code{S < I < R}. Invalid antimicrobial interpretations will be translated as \code{NA} with a warning.
#' @rdname as.rsi
#' @param x vector of values (for class \code{mic}: an MIC value in mg/L, for class \code{disk}: a disk diffusion radius in millimeters)
#' @param mo a microorganism code, generated with \code{\link{as.mo}}
#' @param ab an antibiotic code, generated with \code{\link{as.ab}}
#' @inheritParams first_isolate
#' @param guideline defaults to the latest included EUCAST guideline, run \code{unique(AMR::rsi_translation$guideline)} for all options
#' @param threshold maximum fraction of \code{x} that is allowed to fail transformation, see Examples
#' @param ... parameters passed on to methods
#' @details Run \code{unique(AMR::rsi_translation$guideline)} for a list of all supported guidelines.
#'
#' After using \code{as.rsi}, you can use \code{\link{eucast_rules}} to (1) apply inferred susceptibility and resistance based on results of other antibiotics and (2) apply intrinsic resistance based on taxonomic properties of a microorganism.
#'
#' The function \code{is.rsi.eligible} returns \code{TRUE} when a columns contains at most 5\% invalid antimicrobial interpretations (not S and/or I and/or R), and \code{FALSE} otherwise. The threshold of 5\% can be set with the \code{threshold} parameter.
#' @return Ordered factor with new class \code{rsi}
#' @keywords rsi
#' @export
#' @importFrom dplyr %>%
#' @seealso \code{\link{as.mic}}
#' @inheritSection AMR Read more on our website!
#' @examples
#' rsi_data <- as.rsi(c(rep("S", 474), rep("I", 36), rep("R", 370)))
#' rsi_data <- as.rsi(c(rep("S", 474), rep("I", 36), rep("R", 370), "A", "B", "C"))
#' is.rsi(rsi_data)
#'
#' # this can also coerce combined MIC/RSI values:
#' as.rsi("<= 0.002; S") # will return S
#'
#' # interpret MIC values
#' as.rsi(x = as.mic(2),
#'        mo = as.mo("S. pneumoniae"),
#'        ab = "AMX",
#'        guideline = "EUCAST")
#' as.rsi(x = as.mic(4),
#'        mo = as.mo("S. pneumoniae"),
#'        ab = "AMX",
#'        guideline = "EUCAST")
#'
#' plot(rsi_data)    # for percentages
#' barplot(rsi_data) # for frequencies
#' freq(rsi_data)    # frequency table with informative header
#'
#' # using dplyr's mutate
#' library(dplyr)
#' septic_patients %>%
#'   mutate_at(vars(PEN:RIF), as.rsi)
#'
#'
#' # fastest way to transform all columns with already valid AB results to class `rsi`:
#' septic_patients %>%
#'   mutate_if(is.rsi.eligible,
#'             as.rsi)
#'
#' # default threshold of `is.rsi.eligible` is 5%.
#' is.rsi.eligible(WHONET$`First name`) # fails, >80% is invalid
#' is.rsi.eligible(WHONET$`First name`, threshold = 0.99) # succeeds
as.rsi <- function(x, ...) {
  UseMethod("as.rsi")
}

#' @export
as.rsi.default <- function(x, ...) {
  if (is.rsi(x)) {
    x
  } else if (identical(levels(x), c("S", "I", "R"))) {
    structure(x, class = c('rsi', 'ordered', 'factor'))
  } else {
    # if (input_resembles_mic(x) > 0.5) {
    #   warning("`as.rsi` is intended to clean antimicrobial interpretations - not to interpret MIC values.", call. = FALSE)
    # }

    x <- x %>% unlist()
    x.bak <- x

    na_before <- x[is.na(x) | x == ''] %>% length()
    # remove all spaces
    x <- gsub(' +', '', x)
    # remove all MIC-like values: numbers, operators and periods
    x <- gsub('[0-9.,;:<=>]+', '', x)
    # remove everything between brackets, and 'high' and 'low'
    x <- gsub("([(].*[)])", "", x)
    x <- gsub("(high|low)", "", x, ignore.case = TRUE)
    # disallow more than 3 characters
    x[nchar(x) > 3] <- NA
    # set to capitals
    x <- toupper(x)
    # remove all invalid characters
    x <- gsub('[^RSI]+', '', x)
    # in cases of "S;S" keep S, but in case of "S;I" make it NA
    x <- gsub('^S+$', 'S', x)
    x <- gsub('^I+$', 'I', x)
    x <- gsub('^R+$', 'R', x)
    x[!x %in% c('S', 'I', 'R')] <- NA
    na_after <- x[is.na(x) | x == ''] %>% length()

    if (na_before != na_after) {
      list_missing <- x.bak[is.na(x) & !is.na(x.bak) & x.bak != ''] %>%
        unique() %>%
        sort()
      list_missing <- paste0('"', list_missing , '"', collapse = ", ")
      warning(na_after - na_before, ' results truncated (',
              round(((na_after - na_before) / length(x)) * 100),
              '%) that were invalid antimicrobial interpretations: ',
              list_missing, call. = FALSE)
    }

    structure(.Data = factor(x, levels = c("S", "I", "R"), ordered = TRUE),
              class =  c('rsi', 'ordered', 'factor'))
  }
}

input_resembles_mic <- function(x) {
  mic <- x %>%
    gsub("[^0-9.,]+", "", .) %>%
    unique()
  mic_valid <- suppressWarnings(as.mic(mic))
  result <- sum(!is.na(mic_valid)) / length(mic)
  if (is.na(result)) {
    0
  } else {
    result
  }
}

#' @rdname as.rsi
#' @importFrom dplyr case_when
#' @export
as.rsi.mic <- function(x, mo, ab, guideline = "EUCAST", ...) {
  exec_as.rsi(method = "mic",
              x = x,
              mo = mo,
              ab = ab,
              guideline = guideline)
}

#' @rdname as.rsi
#' @export
as.rsi.disk <- function(x, mo, ab, guideline = "EUCAST", ...) {
  exec_as.rsi(method = "disk",
              x = x,
              mo = mo,
              ab = ab,
              guideline = guideline)
}

exec_as.rsi <- function(method, x, mo, ab, guideline) {
  if (method == "mic") {
    x <- as.mic(x) # when as.rsi.mic is called directly
  } else if (method == "disk") {
    x <- as.disk(x) # when as.rsi.disk is called directly
  }

  mo <- as.mo(mo)
  ab <- as.ab(ab)

  mo_genus <- as.mo(mo_genus(mo))
  mo_family <- as.mo(mo_family(mo))
  mo_order <- as.mo(mo_order(mo))
  mo_becker <- as.mo(mo, Becker = TRUE)
  mo_lancefield <- as.mo(mo, Lancefield = TRUE)

  guideline <- toupper(guideline)
  if (guideline %in% c("CLSI", "EUCAST")) {
    guideline <- AMR::rsi_translation %>%
      filter(guideline %like% guideline) %>%
      pull(guideline) %>%
      sort() %>%
      rev() %>%
      .[1]
  }

  if (!guideline %in% AMR::rsi_translation$guideline) {
    stop(paste0("invalid guideline: '", guideline,
                "'.\nValid guidelines are: ", paste0("'", rev(sort(unique(AMR::rsi_translation$guideline))), "'", collapse = ", ")),
         call. = FALSE)
  }

  new_rsi <- rep(NA_character_, length(x))
  trans <- AMR::rsi_translation %>%
    filter(guideline == guideline) %>%
    mutate(lookup = paste(mo, ab))

  lookup_mo <- paste(mo, ab)
  lookup_genus <- paste(mo_genus, ab)
  lookup_family <- paste(mo_family, ab)
  lookup_order <- paste(mo_order, ab)
  lookup_becker <- paste(mo_becker, ab)
  lookup_lancefield <- paste(mo_lancefield, ab)

  for (i in 1:length(x)) {
    get_record <- trans %>%
      filter(lookup %in% c(lookup_mo[i],
                           lookup_genus[i],
                           lookup_family[i],
                           lookup_order[i],
                           lookup_becker[i],
                           lookup_lancefield[i])) %>%
      # be as specific as possible (i.e. prefer species over genus):
      arrange(desc(nchar(mo))) %>%
      .[1L,]

    if (NROW(get_record) > 0) {
      if (method == "mic") {
        new_rsi[i] <- case_when(is.na(get_record$S_mic) | is.na(get_record$R_mic) ~ NA_character_,
                                x[i] <= get_record$S_mic ~ "S",
                                x[i] >= get_record$R_mic ~ "R",
                                TRUE ~ "I")
      } else if (method == "disk") {
        new_rsi[i] <- case_when(is.na(get_record$S_disk) | is.na(get_record$R_disk) ~ NA_character_,
                                x[i] <= get_record$S_disk ~ "S",
                                x[i] >= get_record$R_disk ~ "R",
                                TRUE ~ "I")
      }

    }
  }
  structure(.Data = factor(new_rsi, levels = c("S", "I", "R"), ordered = TRUE),
            class =  c('rsi', 'ordered', 'factor'))
}

#' @rdname as.rsi
#' @importFrom crayon red blue
#' @export
as.rsi.data.frame <- function(x, col_mo = NULL, guideline = "EUCAST", ...) {
  tbl_ <- x

  ab_cols <- colnames(tbl_)[sapply(tbl_, function(x) is.mic(x) | is.disk(x))]
  if (length(ab_cols) == 0) {
    stop("No columns with MIC values or disk zones found in this data set. Use as.mic or as.disk to transform antibiotic columns.", call. = FALSE)
  }

  # try to find columns based on type
  # -- mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(tbl = tbl_, type = "mo")
  }
  if (is.null(col_mo)) {
    stop("`col_mo` must be set.", call. = FALSE)
  }

  # transform all MICs
  ab_cols <- colnames(tbl_)[sapply(tbl_, is.mic)]
  if (length(ab_cols) > 0) {
    for (i in 1:length(ab_cols)) {
      if (is.na(suppressWarnings(as.ab(ab_cols[i])))) {
        message(red(paste0("Unknown drug: `", bold(ab_cols[i]), "`. Rename this column to a drug name or code, and check the output with as.ab().")))
        next
      }
      message(blue(paste0("Interpreting column `", bold(ab_cols[i]), "` (", ab_name(ab_cols[i], tolower = TRUE), ")...")), appendLF = FALSE)
      tbl_[, ab_cols[i]] <- exec_as.rsi(method = "mic",
                                        x = tbl_ %>% pull(ab_cols[i]),
                                        mo = tbl_ %>% pull(col_mo),
                                        ab = as.ab(ab_cols[i]),
                                        guideline = guideline)
      message(blue(" OK."))
    }
  }
  # transform all disks
  ab_cols <- colnames(tbl_)[sapply(tbl_, is.disk)]
  if (length(ab_cols) > 0) {
    for (i in 1:length(ab_cols)) {
      if (is.na(suppressWarnings(as.ab(ab_cols[i])))) {
        message(red(paste0("Unknown drug: `", bold(ab_cols[i]), "`. Rename this column to a drug name or code, and check the output with as.ab().")))
        next
      }
      message(blue(paste0("Interpreting column `", bold(ab_cols[i]), "` (", ab_name(ab_cols[i], tolower = TRUE), ")...")), appendLF = FALSE)
      tbl_[, ab_cols[i]] <- exec_as.rsi(method = "disk",
                                        x = tbl_ %>% pull(ab_cols[i]),
                                        mo = tbl_ %>% pull(col_mo),
                                        ab = as.ab(ab_cols[i]),
                                        guideline = guideline)
      message(blue(" OK."))
    }
  }

  tbl_
}

#' @rdname as.rsi
#' @export
is.rsi <- function(x) {
  identical(class(x),
            c('rsi', 'ordered', 'factor'))
}

#' @rdname as.rsi
#' @export
is.rsi.eligible <- function(x, threshold = 0.05) {
  if (NCOL(x) > 1) {
    stop('`x` must be a one-dimensional vector.')
  }
  if (any(c("logical",
            "numeric",
            "integer",
            "mo",
            "Date",
            "POSIXct",
            "rsi",
            "raw",
            "hms")
          %in% class(x))) {
    # no transformation needed
    FALSE
  } else {
    x <- x[!is.na(x) & !is.null(x) & !identical(x, "")]
    if (length(x) == 0) {
      return(FALSE)
    }
    checked <- suppressWarnings(as.rsi(x))
    outcome <- sum(is.na(checked)) / length(x)
    outcome <= threshold
  }
}

#' @exportMethod print.rsi
#' @export
#' @importFrom dplyr %>%
#' @noRd
print.rsi <- function(x, ...) {
  cat("Class 'rsi'\n")
  print(as.character(x), quote = FALSE)
}

#' @exportMethod droplevels.rsi
#' @export
#' @noRd
droplevels.rsi <- function(x, exclude = if(anyNA(levels(x))) NULL else NA, ...) {
  x <- droplevels.factor(x, exclude = exclude, ...)
  class(x) <- c('rsi', 'ordered', 'factor')
  x
}

#' @exportMethod summary.rsi
#' @export
#' @noRd
summary.rsi <- function(object, ...) {
  x <- object
  c(
    "Class" = 'rsi',
    "<NA>" = sum(is.na(x)),
    "Sum S" = sum(x == "S", na.rm = TRUE),
    "Sum IR" = sum(x %in% c("I", "R"), na.rm = TRUE),
    "-Sum R" = sum(x == "R", na.rm = TRUE),
    "-Sum I" = sum(x == "I", na.rm = TRUE)
  )
}

#' @exportMethod plot.rsi
#' @export
#' @importFrom dplyr %>% group_by summarise filter mutate if_else n_distinct
#' @importFrom graphics plot text
#' @noRd
plot.rsi <- function(x, ...) {
  x_name <- deparse(substitute(x))

  suppressWarnings(
    data <- data.frame(x = x,
                       y = 1,
                       stringsAsFactors = TRUE) %>%
      group_by(x) %>%
      summarise(n = sum(y)) %>%
      filter(!is.na(x)) %>%
      mutate(s = round((n / sum(n)) * 100, 1))
  )
  if (!"S" %in% data$x) {
    data <- rbind(data, data.frame(x = "S", n = 0, s = 0))
  }
  if (!"I" %in% data$x) {
    data <- rbind(data, data.frame(x = "I", n = 0, s = 0))
  }
  if (!"R" %in% data$x) {
    data <- rbind(data, data.frame(x = "R", n = 0, s = 0))
  }

  data$x <- factor(data$x, levels = c('S', 'I', 'R'), ordered = TRUE)

  ymax <- if_else(max(data$s) > 95, 105, 100)

  plot(x = data$x,
       y = data$s,
       lwd = 2,
       col = c('green', 'orange', 'red'),
       ylim = c(0, ymax),
       ylab = 'Percentage',
       xlab = 'Antimicrobial Interpretation',
       main = paste('Susceptibility Analysis of', x_name),
       axes = FALSE,
       ...)
  # x axis
  axis(side = 1, at = 1:n_distinct(data$x), labels = levels(data$x), lwd = 0)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, 100, 5))

  text(x = data$x,
       y = data$s + 4,
       labels = paste0(data$s, '% (n = ', data$n, ')'))
}


#' @exportMethod barplot.rsi
#' @export
#' @importFrom dplyr %>% group_by summarise
#' @importFrom graphics barplot axis
#' @noRd
barplot.rsi <- function(height, ...) {
  x <- height
  x_name <- deparse(substitute(height))

  suppressWarnings(
    data <- data.frame(rsi = x, cnt = 1) %>%
      group_by(rsi) %>%
      summarise(cnt = sum(cnt)) %>%
      droplevels()
  )

  barplot(table(x),
          col = c('green3', 'orange2', 'red3'),
          xlab = 'Antimicrobial Interpretation',
          main = paste('Susceptibility Analysis of', x_name),
          ylab = 'Frequency',
          axes = FALSE,
          ...)
  # y axis, 0-100%
  axis(side = 2, at = seq(0, max(data$cnt) + max(data$cnt) * 1.1, by = 25))
}
