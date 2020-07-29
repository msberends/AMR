# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.github.io/AMR.    #
# ==================================================================== #

# functions from dplyr, will perhaps become poorman
distinct <- function(.data, ..., .keep_all = FALSE) {
  check_is_dataframe(.data)
  if ("grouped_data" %in% class(.data)) {
    distinct.grouped_data(.data, ..., .keep_all = .keep_all)
  } else {
    distinct.default(.data, ..., .keep_all = .keep_all)
  }
}
distinct.default <- function(.data, ..., .keep_all = FALSE) {
  names <- rownames(.data)
  rownames(.data) <- NULL
  if (length(deparse_dots(...)) == 0) {
    selected <- .data
  } else {
    selected <- select(.data, ...)
  }
  rows <- as.integer(rownames(unique(selected)))
  if (isTRUE(.keep_all)) {
    res <- .data[rows, , drop = FALSE]
  } else {
    res <- selected[rows, , drop = FALSE]
  }
  rownames(res) <- names[rows]
  res
}
distinct.grouped_data <- function(.data, ..., .keep_all = FALSE) {
  apply_grouped_function(.data, "distinct", ..., .keep_all = .keep_all)
}
filter_join_worker <- function(x, y, by = NULL, type = c("anti", "semi")) {
  type <- match.arg(type, choices = c("anti", "semi"), several.ok = FALSE)
  if (is.null(by)) {
    by <- intersect(names(x), names(y))
    join_message(by)
  }
  rows <- interaction(x[, by]) %in% interaction(y[, by])
  if (type == "anti") rows <- !rows
  res <- x[rows, , drop = FALSE]
  rownames(res) <- NULL
  res
}

# No export, no Rd
addin_insert_in <- function() {
  import_fn("insertText", "rstudioapi")(" %in% ")
}

# No export, no Rd
addin_insert_like <- function() {
  import_fn("insertText", "rstudioapi")(" %like% ")
}

check_dataset_integrity <- function() {
  # check if user overwrote our data sets in their global environment
  data_in_pkg <- data(package = "AMR", envir = asNamespace("AMR"))$results[, "Item", drop = TRUE]
  data_in_globalenv <- ls(envir = globalenv())
  overwritten <- data_in_pkg[data_in_pkg %in% data_in_globalenv]
  stop_if(length(overwritten) > 0,
          "the following data set is overwritten by your global environment and prevents the AMR package from working correctly:\n",
          paste0("'", overwritten, "'", collapse = ", "),
          ".\nPlease rename your object before using this function.", call = FALSE)
  # check if other packages did not overwrite our data sets
  tryCatch({
    check_microorganisms <- all(c("mo", "fullname", "kingdom", "phylum",
                                  "class", "order", "family", "genus", 
                                  "species", "subspecies", "rank",
                                  "species_id", "source", "ref", "prevalence") %in% colnames(microorganisms),
                                na.rm = TRUE)
    check_antibiotics <- all(c("ab", "atc", "cid", "name", "group", 
                               "atc_group1", "atc_group2", "abbreviations",
                               "synonyms", "oral_ddd", "oral_units", 
                               "iv_ddd", "iv_units", "loinc") %in% colnames(antibiotics),
                             na.rm = TRUE)
  }, error = function(e)
    stop_('please use the command \'library("AMR")\' before using this function, to load the required reference data.', call = FALSE)
  )
  invisible(TRUE)
}

search_type_in_df <- function(x, type) {
  # try to find columns based on type
  found <- NULL
  
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  colnames(x) <- trimws(colnames(x))
  
  # -- mo
  if (type == "mo") {
    if (any(sapply(x, is.mo))) {
      found <- sort(colnames(x)[sapply(x, is.mo)])[1]
    } else if ("mo" %in% colnames(x) &
               suppressWarnings(
                 all(x$mo %in% c(NA,
                                 microorganisms$mo,
                                 microorganisms.translation$mo_old)))) {
      found <- "mo"
    } else if (any(colnames(x) %like% "^(mo|microorganism|organism|bacteria|ba[ck]terie)s?$")) {
      found <- sort(colnames(x)[colnames(x) %like% "^(mo|microorganism|organism|bacteria|ba[ck]terie)s?$"])[1]
    } else if (any(colnames(x) %like% "^(microorganism|organism|bacteria|ba[ck]terie)")) {
      found <- sort(colnames(x)[colnames(x) %like% "^(microorganism|organism|bacteria|ba[ck]terie)"])[1]
    } else if (any(colnames(x) %like% "species")) {
      found <- sort(colnames(x)[colnames(x) %like% "species"])[1]
    }
    
  }
  # -- key antibiotics
  if (type == "keyantibiotics") {
    if (any(colnames(x) %like% "^key.*(ab|antibiotics)")) {
      found <- sort(colnames(x)[colnames(x) %like% "^key.*(ab|antibiotics)"])[1]
    }
  }
  # -- date
  if (type == "date") {
    if (any(colnames(x) %like% "^(specimen date|specimen_date|spec_date)")) {
      # WHONET support
      found <- sort(colnames(x)[colnames(x) %like% "^(specimen date|specimen_date|spec_date)"])[1]
      if (!any(class(pull(x, found)) %in% c("Date", "POSIXct"))) {
        stop(font_red(paste0("ERROR: Found column `", font_bold(found), "` to be used as input for `col_", type,
                             "`, but this column contains no valid dates. Transform its values to valid dates first.")),
             call. = FALSE)
      }
    } else if (any(sapply(x, function(x) inherits(x, c("Date", "POSIXct"))))) {
      found <- sort(colnames(x)[sapply(x, function(x) inherits(x, c("Date", "POSIXct")))])[1]
    }
  }
  # -- patient id
  if (type == "patient_id") {
    if (any(colnames(x) %like% "^(identification |patient|patid)")) {
      found <- sort(colnames(x)[colnames(x) %like% "^(identification |patient|patid)"])[1]
    }
  }
  # -- specimen
  if (type == "specimen") {
    if (any(colnames(x) %like% "(specimen type|spec_type)")) {
      found <- sort(colnames(x)[colnames(x) %like% "(specimen type|spec_type)"])[1]
    } else if (any(colnames(x) %like% "^(specimen)")) {
      found <- sort(colnames(x)[colnames(x) %like% "^(specimen)"])[1]
    }
  }
  # -- UTI (urinary tract infection)
  if (type == "uti") {
    if (any(colnames(x) == "uti")) {
      found <- colnames(x)[colnames(x) == "uti"][1]
    } else if (any(colnames(x) %like% "(urine|urinary)")) {
      found <- sort(colnames(x)[colnames(x) %like% "(urine|urinary)"])[1]
    }
    if (!is.null(found)) {
      # this column should contain logicals
      if (!is.logical(x[, found, drop = TRUE])) {
        message(font_red(paste0("NOTE: Column `", font_bold(found), "` found as input for `col_", type,
                                "`, but this column does not contain 'logical' values (TRUE/FALSE) and was ignored.")))
        found <- NULL
      }
    }
  }
  
  if (!is.null(found)) {
    msg <- paste0("NOTE: Using column `", font_bold(found), "` as input for `col_", type, "`.")
    if (type %in% c("keyantibiotics", "specimen")) {
      msg <- paste(msg, "Use", font_bold(paste0("col_", type), "= FALSE"), "to prevent this.")
    }
    message(font_blue(msg))
  }
  found
}

stop_ifnot_installed <- function(package) {
  # no "utils::installed.packages()" since it requires non-staged install since R 3.6.0
  # https://developer.r-project.org/Blog/public/2019/02/14/staged-install/index.html
  sapply(package, function(pkg)
    tryCatch(get(".packageName", envir = asNamespace(pkg)),
             error = function(e) {
               if (package == "rstudioapi") {
                 stop("This function only works in RStudio.", call. = FALSE)
               } else if (pkg != "base") {
                 stop("package '", pkg, "' required but not installed.",
                      "\nTry to install it with: install.packages(\"", pkg, "\")",
                      call. = FALSE)
               }
             }))
  return(invisible())
}

import_fn <- function(name, pkg) {
  stop_ifnot_installed(pkg)
  tryCatch(
    get(name, envir = asNamespace(pkg)),
    error = function(e) stop_("an error occurred in import_fn() while using this function", call = FALSE))
}

stop_ <- function(..., call = TRUE) {
  msg <- paste0(c(...), collapse = "")
  if (!isFALSE(call)) {
    if (isTRUE(call)) {
      call <- as.character(sys.call(-1)[1])
    } else {
      # so you can go back more than 1 call, as used in rsi_calc(), that now throws a reference to e.g. n_rsi()
      call <- as.character(sys.call(call)[1])
    }
    msg <- paste0("in ", call, "(): ", msg)
  }
  stop(msg, call. = FALSE)
}

stop_if <- function(expr, ..., call = TRUE) {
  if (isTRUE(expr)) {
    if (isTRUE(call)) {
      call <- -1
    }
    if (!isFALSE(call)) {
      # since we're calling stop_(), which is another call
      call <- call - 1
    }
    stop_(..., call = call)
  }
}

stop_ifnot <- function(expr, ..., call = TRUE) {
  if (!isTRUE(expr)) {
    if (isTRUE(call)) {
      call <- -1
    }
    if (!isFALSE(call)) {
      # since we're calling stop_(), which is another call
      call <- call - 1
    }
    stop_(..., call = call)
  }
}

"%or%" <- function(x, y) {
  if (is.null(x) | is.null(y)) {
    if (is.null(x)) {
      return(y)
    } else {
      return(x)
    }
  }
  ifelse(!is.na(x),
         x,
         ifelse(!is.na(y), y, NA))
}

class_integrity_check <- function(value, type, check_vector) {
  if (!all(value[!is.na(value)] %in% check_vector)) {
    warning(paste0("invalid ", type, ", NA generated"), call. = FALSE)
    value[!value %in% check_vector] <- NA
  }
  value
}

# transforms data set to data.frame with only ASCII values, to comply with CRAN policies
dataset_UTF8_to_ASCII <- function(df) {
  trans <- function(vect) {
    iconv(vect, from = "UTF-8", to = "ASCII//TRANSLIT")
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (i in seq_len(NCOL(df))) {
    col <- df[, i]
    if (is.list(col)) {
      col <- lapply(col, function(j) trans(j))
      df[, i] <- list(col)
    } else {
      if (is.factor(col)) {
        levels(col) <- trans(levels(col))
      } else if (is.character(col)) {
        col <- trans(col)
      } else {
        col
      }
      df[, i] <- col
    }
  }
  df
}

has_colour <- function() {
  # this is a base R version of crayon::has_color
  enabled <- getOption("crayon.enabled")
  if (!is.null(enabled)) {
    return(isTRUE(enabled))
  }
  rstudio_with_ansi_support <- function(x) {
    if (Sys.getenv("RSTUDIO", "") == "") {
      return(FALSE)
    }
    if ((cols <- Sys.getenv("RSTUDIO_CONSOLE_COLOR", "")) != "" && !is.na(as.numeric(cols))) {
      return(TRUE)
    }
    tryCatch(get("isAvailable", envir = asNamespace("rstudioapi"))(), error = function(e) return(FALSE)) &&
      tryCatch(get("hasFun", envir = asNamespace("rstudioapi"))("getConsoleHasColor"), error = function(e) return(FALSE))
  }
  if (rstudio_with_ansi_support() && sink.number() == 0) {
    return(TRUE)
  }
  if (!isatty(stdout())) {
    return(FALSE)
  }
  if (tolower(Sys.info()["sysname"]) == "windows") {
    if (Sys.getenv("ConEmuANSI") == "ON") {
      return(TRUE)
    }
    if (Sys.getenv("CMDER_ROOT") != "") {
      return(TRUE)
    }
    return(FALSE)
  }
  emacs_version <- function() {
    ver <- Sys.getenv("INSIDE_EMACS")
    if (ver == "") {
      return(NA_integer_)
    }
    ver <- gsub("'", "", ver)
    ver <- strsplit(ver, ",", fixed = TRUE)[[1]]
    ver <- strsplit(ver, ".", fixed = TRUE)[[1]]
    as.numeric(ver)
  }
  if ((Sys.getenv("EMACS") != "" || Sys.getenv("INSIDE_EMACS") != "") &&
      !is.na(emacs_version()[1]) && emacs_version()[1] >= 23) {
    return(TRUE)
  }
  if ("COLORTERM" %in% names(Sys.getenv())) {
    return(TRUE)
  }
  if (Sys.getenv("TERM") == "dumb") {
    return(FALSE)
  }
  grepl(pattern = "^screen|^xterm|^vt100|color|ansi|cygwin|linux", 
        x = Sys.getenv("TERM"),
        ignore.case = TRUE,
        perl = TRUE)
}

# the crayon colours
try_colour <- function(..., before, after, collapse = " ") {
  txt <- paste0(unlist(list(...)), collapse = collapse)
  if (isTRUE(has_colour())) {
    if (is.null(collapse)) {
      paste0(before, txt, after, collapse = NULL)
    } else {
      paste0(before, txt, after, collapse = "")
    }
  } else {
    txt
  }
}
font_black <- function(..., collapse = " ") {
  try_colour(..., before = "\033[38;5;232m", after = "\033[39m", collapse = collapse)
}
font_blue <- function(..., collapse = " ") {
  try_colour(..., before = "\033[34m", after = "\033[39m", collapse = collapse)
}
font_green <- function(..., collapse = " ") {
  try_colour(..., before = "\033[32m", after = "\033[39m", collapse = collapse)
}
font_magenta <- function(..., collapse = " ") {
  try_colour(..., before = "\033[35m", after = "\033[39m", collapse = collapse)
}
font_red <- function(..., collapse = " ") {
  try_colour(..., before = "\033[31m", after = "\033[39m", collapse = collapse)
}
font_silver <- function(..., collapse = " ") {
  try_colour(..., before = "\033[90m", after = "\033[39m", collapse = collapse)
}
font_white <- function(..., collapse = " ") {
  try_colour(..., before = "\033[37m", after = "\033[39m", collapse = collapse)
}
font_yellow <- function(..., collapse = " ") {
  try_colour(..., before = "\033[33m", after = "\033[39m", collapse = collapse)
}
font_subtle <- function(..., collapse = " ") {
  try_colour(..., before = "\033[38;5;246m", after = "\033[39m", collapse = collapse)
}
font_grey <- function(..., collapse = " ") {
  try_colour(..., before = "\033[38;5;249m", after = "\033[39m", collapse = collapse)
}
font_green_bg <- function(..., collapse = " ") {
  try_colour(..., before = "\033[42m", after = "\033[49m", collapse = collapse)
}
font_red_bg <- function(..., collapse = " ") {
  try_colour(..., before = "\033[41m", after = "\033[49m", collapse = collapse)
}
font_yellow_bg <- function(..., collapse = " ") {
  try_colour(..., before = "\033[43m", after = "\033[49m", collapse = collapse)
}
font_bold <- function(..., collapse = " ") {
  try_colour(..., before = "\033[1m", after = "\033[22m", collapse = collapse)
}
font_italic <- function(..., collapse = " ") {
  try_colour(..., before = "\033[3m", after = "\033[23m", collapse = collapse)
}
font_underline <- function(..., collapse = " ") {
  try_colour(..., before = "\033[4m", after = "\033[24m", collapse = collapse)
}
font_stripstyle <- function(x) {
  # from crayon:::ansi_regex
  gsub("(?:(?:\\x{001b}\\[)|\\x{009b})(?:(?:[0-9]{1,3})?(?:(?:;[0-9]{0,3})*)?[A-M|f-m])|\\x{001b}[A-M]", "", x, perl = TRUE)
}

progress_estimated <- function(n = 1, n_min = 0, ...) {
  if (!interactive() || n < n_min) {
    pb <- list()
    pb$tick <- function() {
      invisible()
    }
    pb$kill <- function() {
      invisible()
    }
    structure(pb, class = "txtProgressBar")
  } else if (n >= n_min) {
    pb <- utils::txtProgressBar(max = n, style = 3)
    pb$tick <- function() {
      pb$up(pb$getVal() + 1)
    }
    pb
  }
}

# works exactly like round(), but rounds `round2(44.55, 1)` to 44.6 instead of 44.5
# and adds decimal zeroes until `digits` is reached when force_zero = TRUE
round2 <- function(x, digits = 0, force_zero = TRUE) {
  x <- as.double(x)
  # https://stackoverflow.com/a/12688836/4575331
  val <- (trunc((abs(x) * 10 ^ digits) + 0.5) / 10 ^ digits) * sign(x)
  if (digits > 0 & force_zero == TRUE) {
    values_trans <- val[val != as.integer(val) & !is.na(val)]
    val[val != as.integer(val) & !is.na(val)] <- paste0(values_trans,
                                                        strrep("0", 
                                                               max(0, 
                                                                   digits - nchar(
                                                                     format(
                                                                       as.double(
                                                                         gsub(".*[.](.*)$", 
                                                                              "\\1",
                                                                              values_trans)),
                                                                       scientific = FALSE)))))
  }
  as.double(val)
}


# percentage from our other package: 'cleaner'
percentage <- function(x, digits = NULL, ...) {
  
  # getdecimalplaces() function
  getdecimalplaces <- function(x, minimum = 0, maximum = 3) {
    if (maximum < minimum) {
      maximum <- minimum
    }
    if (minimum > maximum) {
      minimum <- maximum
    }
    max_places <- max(unlist(lapply(strsplit(sub("0+$", "", 
                                                 as.character(x * 100)), ".", fixed = TRUE),
                                    function(y) ifelse(length(y) == 2, nchar(y[2]), 0))), na.rm = TRUE)
    max(min(max_places,
            maximum, na.rm = TRUE),
        minimum, na.rm = TRUE)
  }
  
  # format_percentage() function
  format_percentage <- function(x, digits = NULL, ...) {
    if (is.null(digits)) {
      digits <- getdecimalplaces(x)
    }
    
    # round right: percentage(0.4455) and format(as.percentage(0.4455), 1) should return "44.6%", not "44.5%"
    x_formatted <- format(round2(as.double(x), digits = digits + 2) * 100,
                          scientific = FALSE,
                          digits = digits,
                          nsmall = digits,
                          ...)
    x_formatted <- paste0(x_formatted, "%")
    x_formatted[!grepl(pattern = "^[0-9.,e-]+$", x = x)] <- NA_character_
    x_formatted
  }
  
  # the actual working part
  x <- as.double(x)
  if (is.null(digits)) {
    # max one digit if undefined
    digits <- getdecimalplaces(x, minimum = 0, maximum = 1)
  }
  format_percentage(structure(.Data = as.double(x),
                              class = c("percentage", "numeric")),
                    digits = digits, ...)
}

# prevent dependency on package 'backports'
# these functions were not available in previous versions of R (last checked: R 4.0.0)
# see here for the full list: https://github.com/r-lib/backports
strrep <- function(x, times) {
  x <- as.character(x)
  if (length(x) == 0L) 
    return(x)
  unlist(.mapply(function(x, times) {
    if (is.na(x) || is.na(times)) 
      return(NA_character_)
    if (times <= 0L) 
      return("")
    paste0(replicate(times, x), collapse = "")
  }, list(x = x, times = times), MoreArgs = list()), use.names = FALSE)
}
trimws <- function(x, which = c("both", "left", "right")) {
  which <- match.arg(which)
  mysub <- function(re, x) sub(re, "", x, perl = TRUE)
  if (which == "left") 
    return(mysub("^[ \t\r\n]+", x))
  if (which == "right") 
    return(mysub("[ \t\r\n]+$", x))
  mysub("[ \t\r\n]+$", mysub("^[ \t\r\n]+", x))
}
isFALSE <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x) && !x
}
deparse1 <- function(expr, collapse = " ", width.cutoff = 500L, ...) {
  paste(deparse(expr, width.cutoff, ...), collapse = collapse)
}
