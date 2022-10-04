# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# doi:10.18637/jss.v104.i03                                            #
#                                                                      #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# faster implementation of left_join than using merge() by poorman - we use match():
pm_left_join <- function(x, y, by = NULL, suffix = c(".x", ".y")) {
  if (is.null(by)) {
    by <- intersect(names(x), names(y))[1L]
    if (is.na(by)) {
      stop_("no common column found for pm_left_join()")
    }
    pm_join_message(by)
  } else if (!is.null(names(by))) {
    by <- unname(c(names(by), by))
  }
  if (length(by) == 1) {
    by <- rep(by, 2)
  }

  int_x <- colnames(x) %in% colnames(y) & colnames(x) != by[1]
  int_y <- colnames(y) %in% colnames(x) & colnames(y) != by[2]
  colnames(x)[int_x] <- paste0(colnames(x)[int_x], suffix[1L])
  colnames(y)[int_y] <- paste0(colnames(y)[int_y], suffix[2L])

  merged <- cbind(
    x,
    y[match(
      x[, by[1], drop = TRUE],
      y[, by[2], drop = TRUE]
    ),
    colnames(y)[!colnames(y) %in% colnames(x) & !colnames(y) == by[2]],
    drop = FALSE
    ]
  )

  rownames(merged) <- NULL
  merged
}

# support where() like tidyverse:
# adapted from https://github.com/nathaneastwood/poorman/blob/52eb6947e0b4430cd588976ed8820013eddf955f/R/where.R#L17-L32
where <- function(fn) {
  if (!is.function(fn)) {
    stop(pm_deparse_var(fn), " is not a valid predicate function.")
  }
  preds <- unlist(lapply(
    pm_select_env$.data,
    function(x, fn) {
      do.call("fn", list(x))
    },
    fn
  ))
  if (!is.logical(preds)) stop("`where()` must be used with functions that return `TRUE` or `FALSE`.")
  data_cols <- pm_select_env$get_colnames()
  cols <- data_cols[preds]
  which(data_cols %in% cols)
}

# copied and slightly rewritten from poorman under same license (2021-10-15)
quick_case_when <- function(...) {
  fs <- list(...)
  lapply(fs, function(x) {
    if (!inherits(x, "formula")) {
      stop("`case_when()` requires formula inputs.")
    }
  })
  n <- length(fs)
  if (n == 0L) {
    stop("No cases provided.")
  }

  validate_case_when_length <- function(query, value, fs) {
    lhs_lengths <- lengths(query)
    rhs_lengths <- lengths(value)
    all_lengths <- unique(c(lhs_lengths, rhs_lengths))
    if (length(all_lengths) <= 1L) {
      return(all_lengths[[1L]])
    }
    non_atomic_lengths <- all_lengths[all_lengths != 1L]
    len <- non_atomic_lengths[[1L]]
    if (length(non_atomic_lengths) == 1L) {
      return(len)
    }
    inconsistent_lengths <- non_atomic_lengths[-1L]
    lhs_problems <- lhs_lengths %in% inconsistent_lengths
    rhs_problems <- rhs_lengths %in% inconsistent_lengths
    problems <- lhs_problems | rhs_problems
    if (any(problems)) {
      stop("The following formulas must be length ", len, " or 1, not ",
        paste(inconsistent_lengths, collapse = ", "), ".\n    ",
        paste(fs[problems], collapse = "\n    "),
        call. = FALSE
      )
    }
  }

  replace_with <- function(x, i, val, arg_name) {
    if (is.null(val)) {
      return(x)
    }
    i[is.na(i)] <- FALSE
    if (length(val) == 1L) {
      x[i] <- val
    } else {
      x[i] <- val[i]
    }
    x
  }

  query <- vector("list", n)
  value <- vector("list", n)
  default_env <- parent.frame()
  for (i in seq_len(n)) {
    query[[i]] <- eval(fs[[i]][[2]], envir = default_env)
    value[[i]] <- eval(fs[[i]][[3]], envir = default_env)
    if (!is.logical(query[[i]])) {
      stop(fs[[i]][[2]], " does not return a `logical` vector.")
    }
  }
  m <- validate_case_when_length(query, value, fs)
  out <- value[[1]][rep(NA_integer_, m)]
  replaced <- rep(FALSE, m)
  for (i in seq_len(n)) {
    out <- replace_with(
      out, query[[i]] & !replaced, value[[i]],
      NULL
    )
    replaced <- replaced | (query[[i]] & !is.na(query[[i]]))
  }
  out
}

# No export, no Rd
addin_insert_in <- function() {
  import_fn("insertText", "rstudioapi")(" %in% ")
}

# No export, no Rd
addin_insert_like <- function() {
  # we want Shift + Ctrl/Cmd + L to iterate over %like%, %unlike%, %like_case%, and %unlike_case%

  getActiveDocumentContext <- import_fn("getActiveDocumentContext", "rstudioapi")
  insertText <- import_fn("insertText", "rstudioapi")
  modifyRange <- import_fn("modifyRange", "rstudioapi")
  document_range <- import_fn("document_range", "rstudioapi")
  document_position <- import_fn("document_position", "rstudioapi")

  context <- getActiveDocumentContext()
  current_row <- context$selection[[1]]$range$end[1]
  current_col <- context$selection[[1]]$range$end[2]
  current_row_txt <- context$contents[current_row]
  if (is.null(current_row) || current_row_txt %unlike% "%(un)?like") {
    insertText(" %like% ")
    return(invisible())
  }

  pos_preceded_by <- function(txt) {
    if (tryCatch(substr(current_row_txt, current_col - nchar(trimws(txt, which = "right")), current_col) == trimws(txt, which = "right"),
      error = function(e) FALSE
    )) {
      return(TRUE)
    }
    tryCatch(substr(current_row_txt, current_col - nchar(txt), current_col) %like% paste0("^", txt),
      error = function(e) FALSE
    )
  }
  replace_pos <- function(old, with) {
    modifyRange(document_range(
      document_position(current_row, current_col - nchar(old)),
      document_position(current_row, current_col)
    ),
    text = with,
    id = context$id
    )
  }

  if (pos_preceded_by(" %like% ")) {
    replace_pos(" %like% ", with = " %unlike% ")
  } else if (pos_preceded_by(" %unlike% ")) {
    replace_pos(" %unlike% ", with = " %like_case% ")
  } else if (pos_preceded_by(" %like_case% ")) {
    replace_pos(" %like_case% ", with = " %unlike_case% ")
  } else if (pos_preceded_by(" %unlike_case% ")) {
    replace_pos(" %unlike_case% ", with = " %like% ")
  } else {
    insertText(" %like% ")
  }
}

search_type_in_df <- function(x, type, info = TRUE) {
  meet_criteria(x, allow_class = "data.frame")
  meet_criteria(type, allow_class = "character", has_length = 1)

  # try to find columns based on type
  found <- NULL

  # remove attributes from other packages
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  colnames_formatted <- tolower(generalise_antibiotic_name(colnames(x)))

  # -- mo
  if (type == "mo") {
    if (any(vapply(FUN.VALUE = logical(1), x, is.mo))) {
      # take first <mo> column
      found <- colnames(x)[vapply(FUN.VALUE = logical(1), x, is.mo)]
    } else if ("mo" %in% colnames_formatted &&
      suppressWarnings(all(x$mo %in% c(NA, AMR::microorganisms$mo)))) {
      found <- "mo"
    } else if (any(colnames_formatted %like_case% "^(mo|microorganism|organism|bacteria|ba[ck]terie)s?$")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "^(mo|microorganism|organism|bacteria|ba[ck]terie)s?$"])
    } else if (any(colnames_formatted %like_case% "^(microorganism|organism|bacteria|ba[ck]terie)")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "^(microorganism|organism|bacteria|ba[ck]terie)"])
    } else if (any(colnames_formatted %like_case% "species")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "species"])
    }
  }
  # -- key antibiotics
  if (type %in% c("keyantibiotics", "keyantimicrobials")) {
    if (any(colnames_formatted %like_case% "^key.*(ab|antibiotics|antimicrobials)")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "^key.*(ab|antibiotics|antimicrobials)"])
    }
  }
  # -- date
  if (type == "date") {
    if (any(colnames_formatted %like_case% "^(specimen date|specimen_date|spec_date)")) {
      # WHONET support
      found <- sort(colnames(x)[colnames_formatted %like_case% "^(specimen date|specimen_date|spec_date)"])
      if (!inherits(pm_pull(x, found), c("Date", "POSIXct"))) {
        stop(font_red(paste0(
          "Found column '", font_bold(found), "' to be used as input for `col_", type,
          "`, but this column contains no valid dates. Transform its values to valid dates first."
        )),
        call. = FALSE
        )
      }
    } else if (any(vapply(FUN.VALUE = logical(1), x, function(x) inherits(x, c("Date", "POSIXct"))))) {
      # take first <Date> column
      found <- colnames(x)[vapply(FUN.VALUE = logical(1), x, function(x) inherits(x, c("Date", "POSIXct")))]
    }
  }
  # -- patient id
  if (type == "patient_id") {
    crit1 <- colnames_formatted %like_case% "^(patient|patid)"
    if (any(crit1)) {
      found <- colnames(x)[crit1]
    } else {
      crit2 <- colnames_formatted %like_case% "(identification |patient|pat.*id)"
      if (any(crit2)) {
        found <- colnames(x)[crit2]
      }
    }
  }
  # -- specimen
  if (type == "specimen") {
    if (any(colnames_formatted %like_case% "(specimen type|spec_type)")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "(specimen type|spec_type)"])
    } else if (any(colnames_formatted %like_case% "^(specimen)")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "^(specimen)"])
    }
  }
  # -- UTI (urinary tract infection)
  if (type == "uti") {
    if (any(colnames_formatted == "uti")) {
      found <- colnames(x)[colnames_formatted == "uti"]
    } else if (any(colnames_formatted %like_case% "(urine|urinary)")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "(urine|urinary)"])
    }
    if (!is.null(found)) {
      # this column should contain logicals
      if (!is.logical(x[, found, drop = TRUE])) {
        message_("Column '", font_bold(found), "' found as input for `col_", type,
          "`, but this column does not contain 'logical' values (TRUE/FALSE) and was ignored.",
          add_fn = font_red
        )
        found <- NULL
      }
    }
  }

  found <- found[1]

  if (!is.null(found) && info == TRUE) {
    if (message_not_thrown_before("search_in_type", type)) {
      msg <- paste0("Using column '", font_bold(found), "' as input for `col_", type, "`.")
      if (type %in% c("keyantibiotics", "keyantimicrobials", "specimen")) {
        msg <- paste(msg, "Use", font_bold(paste0("col_", type), "= FALSE"), "to prevent this.")
      }
      message_(msg)
    }
  }
  found
}

is_valid_regex <- function(x) {
  regex_at_all <- tryCatch(vapply(
    FUN.VALUE = logical(1),
    X = strsplit(x, "", fixed = TRUE),
    FUN = function(y) {
      any(y %in% c(
        "$", "(", ")", "*", "+", "-",
        ".", "?", "[", "]", "^", "{",
        "|", "}", "\\"
      ),
      na.rm = TRUE
      )
    },
    USE.NAMES = FALSE
  ),
  error = function(e) rep(TRUE, length(x))
  )
  regex_valid <- vapply(
    FUN.VALUE = logical(1),
    X = x,
    FUN = function(y) {
      !inherits(try(grepl(y, "", perl = TRUE), silent = TRUE), "try-error")
    },
    USE.NAMES = FALSE
  )
  regex_at_all & regex_valid
}

stop_ifnot_installed <- function(package) {
  # no "utils::installed.packages()" since it requires non-staged install since R 3.6.0
  # https://developer.r-project.org/Blog/public/2019/02/14/staged-install/index.html
  vapply(FUN.VALUE = character(1), package, function(pkg) {
    tryCatch(get(".packageName", envir = asNamespace(pkg)),
      error = function(e) {
        if (pkg == "rstudioapi") {
          stop("This function only works in RStudio when using R >= 3.2.", call. = FALSE)
        } else if (pkg != "base") {
          stop("This requires the '", pkg, "' package.",
            "\nTry to install it with: install.packages(\"", pkg, "\")",
            call. = FALSE
          )
        }
      }
    )
  })
  return(invisible())
}

pkg_is_available <- function(pkg, also_load = TRUE, min_version = NULL) {
  if (also_load == TRUE) {
    out <- suppressWarnings(require(pkg, character.only = TRUE, warn.conflicts = FALSE))
  } else {
    out <- requireNamespace(pkg, quietly = TRUE)
  }
  if (!is.null(min_version)) {
    out <- out && utils::packageVersion(pkg) >= min_version
  }
  isTRUE(out)
}

import_fn <- function(name, pkg, error_on_fail = TRUE) {
  if (isTRUE(error_on_fail)) {
    stop_ifnot_installed(pkg)
  }
  tryCatch(
    # don't use get() to avoid fetching non-API functions
    getExportedValue(name = name, ns = asNamespace(pkg)),
    error = function(e) {
      if (isTRUE(error_on_fail)) {
        stop_("function ", name, "() is not an exported object from package '", pkg,
          "'. Please create an issue at https://github.com/msberends/AMR/issues. Many thanks!",
          call = FALSE
        )
      } else {
        return(NULL)
      }
    }
  )
}

# this alternative wrapper to the message(), warning() and stop() functions:
# - wraps text to never break lines within words
# - ignores formatted text while wrapping
# - adds indentation dependent on the type of message (such as NOTE)
# - can add additional formatting functions like blue or bold text
word_wrap <- function(...,
                      add_fn = list(),
                      as_note = FALSE,
                      width = 0.95 * getOption("width"),
                      extra_indent = 0) {
  msg <- paste0(c(...), collapse = "")

  if (isTRUE(as_note)) {
    msg <- paste0(AMR_env$info_icon, " ", gsub("^note:? ?", "", msg, ignore.case = TRUE))
  }

  if (msg %like% "\n") {
    # run word_wraps() over every line here, bind them and return again
    return(paste0(vapply(
      FUN.VALUE = character(1),
      trimws(unlist(strsplit(msg, "\n", fixed = TRUE)), which = "right"),
      word_wrap,
      add_fn = add_fn,
      as_note = FALSE,
      width = width,
      extra_indent = extra_indent
    ),
    collapse = "\n"
    ))
  }

  # correct for operators (will add the space later on)
  ops <- "([,./><\\]\\[])"
  msg <- gsub(paste0(ops, " ", ops), "\\1\\2", msg, perl = TRUE)
  # we need to correct for already applied style, that adds text like "\033[31m\"
  msg_stripped <- font_stripstyle(msg)
  # where are the spaces now?
  msg_stripped_wrapped <- paste0(strwrap(msg_stripped,
    simplify = TRUE,
    width = width
  ),
  collapse = "\n"
  )
  msg_stripped_wrapped <- paste0(unlist(strsplit(msg_stripped_wrapped, "(\n|\\*\\|\\*)")),
    collapse = "\n"
  )
  msg_stripped_spaces <- which(unlist(strsplit(msg_stripped, "", fixed = TRUE)) == " ")
  msg_stripped_wrapped_spaces <- which(unlist(strsplit(msg_stripped_wrapped, "", fixed = TRUE)) != "\n")
  # so these are the indices of spaces that need to be replaced
  replace_spaces <- which(!msg_stripped_spaces %in% msg_stripped_wrapped_spaces)
  # put it together
  msg <- unlist(strsplit(msg, " ", fixed = TRUE))
  msg[replace_spaces] <- paste0(msg[replace_spaces], "\n")
  # add space around operators again
  msg <- gsub(paste0(ops, ops), "\\1 \\2", msg, perl = TRUE)
  msg <- paste0(msg, collapse = " ")
  msg <- gsub("\n ", "\n", msg, fixed = TRUE)

  if (msg_stripped %like% "\u2139 ") {
    indentation <- 2 + extra_indent
  } else if (msg_stripped %like% "^=> ") {
    indentation <- 3 + extra_indent
  } else {
    indentation <- 0 + extra_indent
  }
  msg <- gsub("\n", paste0("\n", strrep(" ", indentation)), msg, fixed = TRUE)
  # remove trailing empty characters
  msg <- gsub("(\n| )+$", "", msg)

  if (length(add_fn) > 0) {
    if (!is.list(add_fn)) {
      add_fn <- list(add_fn)
    }
    for (i in seq_len(length(add_fn))) {
      msg <- add_fn[[i]](msg)
    }
  }

  # format backticks
  msg <- gsub("(`.+?`)", font_grey_bg("\\1"), msg)

  # clean introduced whitespace between fullstops
  msg <- gsub("[.] +[.]", "..", msg)
  # remove extra space that was introduced (case: "Smith et al., 2022")
  msg <- gsub(". ,", ".,", msg, fixed = TRUE)

  msg
}

message_ <- function(...,
                     appendLF = TRUE,
                     add_fn = list(font_blue),
                     as_note = TRUE) {
  message(word_wrap(...,
    add_fn = add_fn,
    as_note = as_note
  ),
  appendLF = appendLF
  )
}

warning_ <- function(...,
                     add_fn = list(),
                     immediate = FALSE,
                     call = FALSE) {
  warning(word_wrap(...,
    add_fn = add_fn,
    as_note = FALSE
  ),
  immediate. = immediate,
  call. = call
  )
}

# this alternative to the stop() function:
# - adds the function name where the error was thrown
# - wraps text to never break lines within words
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
  msg <- word_wrap(msg, add_fn = list(), as_note = FALSE)
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
  if (isFALSE(expr)) {
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
  if (is.null(x) || is.null(y)) {
    if (is.null(x)) {
      return(y)
    } else {
      return(x)
    }
  }
  ifelse(is.na(x), y, x)
}

return_after_integrity_check <- function(value, type, check_vector) {
  if (!all(value[!is.na(value)] %in% check_vector)) {
    warning_(paste0("invalid ", type, ", NA generated"))
    value[!value %in% check_vector] <- NA
  }
  value
}

# transforms data set to a tibble with only ASCII values, to comply with CRAN policies
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
  import_fn("as_tibble", "tibble")(df)
}

documentation_date <- function(d) {
  paste0(trimws(format(d, "%e")), " ", month.name[as.integer(format(d, "%m"))], ", ", format(d, "%Y"))
}

format_included_data_number <- function(data) {
  if (is.data.frame(data)) {
    n <- nrow(data)
  } else {
    n <- length(unique(data))
  }
  if (n > 10000) {
    rounder <- -3 # round on thousands
  } else if (n > 1000) {
    rounder <- -2 # round on hundreds
  } else {
    rounder <- -1 # round on tens
  }
  paste0("~", format(round(n, rounder), decimal.mark = ".", big.mark = ","))
}

# for eucast_rules() and mdro(), creates markdown output with URLs and names
create_eucast_ab_documentation <- function() {
  x <- trimws(unique(toupper(unlist(strsplit(EUCAST_RULES_DF$then_change_these_antibiotics, ",", fixed = TRUE)))))
  ab <- character()
  for (val in x) {
    if (paste0("AB_", val) %in% ls(envir = asNamespace("AMR"))) {
      # antibiotic group names, as defined in data-raw/_pre_commit_hook.R, such as `CARBAPENEMS`
      val <- eval(parse(text = paste0("AB_", val)), envir = asNamespace("AMR"))
    } else if (val %in% AB_lookup$ab) {
      # separate drugs, such as `AMX`
      val <- as.ab(val)
    } else {
      val <- as.rsi(NA)
    }
    ab <- c(ab, val)
  }
  ab <- unique(ab)
  atcs <- ab_atc(ab, only_first = TRUE)
  # only keep ABx with an ATC code:
  ab <- ab[!is.na(atcs)]
  ab_names <- ab_name(ab, language = NULL, tolower = TRUE)
  ab <- ab[order(ab_names)]
  ab_names <- ab_names[order(ab_names)]
  atc_txt <- paste0("[", atcs[!is.na(atcs)], "](", ab_url(ab), ")")
  out <- paste0(ab_names, " (`", ab, "`, ", atc_txt, ")", collapse = ", ")
  substr(out, 1, 1) <- toupper(substr(out, 1, 1))
  out
}

vector_or <- function(v, quotes = TRUE, reverse = FALSE, sort = TRUE, initial_captital = FALSE, last_sep = " or ") {
  # makes unique and sorts, and this also removed NAs
  v <- unique(v)
  if (isTRUE(sort)) {
    v <- sort(v)
  }
  if (isTRUE(reverse)) {
    v <- rev(v)
  }
  if (isTRUE(quotes)) {
    quotes <- '"'
  } else if (isFALSE(quotes)) {
    quotes <- ""
  } else {
    quotes <- quotes[1L]
  }
  if (isTRUE(initial_captital)) {
    v[1] <- gsub("^([a-z])", "\\U\\1", v[1], perl = TRUE)
  }
  if (length(v) <= 1) {
    return(paste0(quotes, v, quotes))
  }
  if (identical(v, c("I", "R", "S"))) {
    # class <rsi> should be sorted like this
    v <- c("R", "S", "I")
  }
  # all commas except for last item, so will become '"val1", "val2", "val3" or "val4"'
  paste0(
    paste0(quotes, v[seq_len(length(v) - 1)], quotes, collapse = ", "),
    last_sep, paste0(quotes, v[length(v)], quotes)
  )
}

vector_and <- function(v, quotes = TRUE, reverse = FALSE, sort = TRUE, initial_captital = FALSE) {
  vector_or(
    v = v, quotes = quotes, reverse = reverse, sort = sort,
    initial_captital = initial_captital, last_sep = " and "
  )
}

format_class <- function(class, plural = FALSE) {
  class.bak <- class
  class[class == "numeric"] <- "number"
  class[class == "integer"] <- "whole number"
  if (all(c("numeric", "integer") %in% class.bak, na.rm = TRUE)) {
    class[class %in% c("number", "whole number")] <- "(whole) number"
  }
  class[class == "character"] <- "text string"
  class[class == "Date"] <- "date"
  class[class %in% c("POSIXt", "POSIXct", "POSIXlt")] <- "date/time"
  class[class != class.bak] <- paste0(
    ifelse(plural, "", "a "),
    class[class != class.bak],
    ifelse(plural, "s", "")
  )
  # exceptions
  class[class == "logical"] <- ifelse(plural, "a vector of `TRUE`/`FALSE`", "`TRUE` or `FALSE`")
  class[class == "data.frame"] <- "a data set"
  if ("list" %in% class) {
    class <- "a list"
  }
  if ("matrix" %in% class) {
    class <- "a matrix"
  }
  if ("custom_eucast_rules" %in% class) {
    class <- "input created with `custom_eucast_rules()`"
  }
  if (any(c("mo", "ab", "rsi") %in% class)) {
    class <- paste0("of class <", class[1L], ">")
  }
  class[class == class.bak] <- paste0("of class <", class[class == class.bak], ">")
  # output
  vector_or(class, quotes = FALSE, sort = FALSE)
}

# a check for every single argument in all functions
meet_criteria <- function(object,
                          allow_class = NULL,
                          has_length = NULL,
                          looks_like = NULL,
                          is_in = NULL,
                          is_positive = NULL,
                          is_positive_or_zero = NULL,
                          is_finite = NULL,
                          contains_column_class = NULL,
                          allow_NULL = FALSE,
                          allow_NA = FALSE,
                          ignore.case = FALSE,
                          .call_depth = 0) { # depth in calling

  obj_name <- deparse(substitute(object))
  call_depth <- -2 - abs(.call_depth)

  # if object is missing, or another error:
  tryCatch(invisible(object),
    error = function(e) AMR_env$meet_criteria_error_txt <- e$message
  )
  if (!is.null(AMR_env$meet_criteria_error_txt)) {
    error_txt <- AMR_env$meet_criteria_error_txt
    AMR_env$meet_criteria_error_txt <- NULL
    stop(error_txt, call. = FALSE) # don't use stop_() here, our pkg may not be loaded yet
  }
  AMR_env$meet_criteria_error_txt <- NULL

  if (is.null(object)) {
    stop_if(allow_NULL == FALSE, "argument `", obj_name, "` must not be NULL", call = call_depth)
    return(invisible())
  }
  if (is.null(dim(object)) && length(object) == 1 && suppressWarnings(is.na(object))) { # suppressWarnings for functions
    stop_if(allow_NA == FALSE, "argument `", obj_name, "` must not be NA", call = call_depth)
    return(invisible())
  }

  if (!is.null(allow_class)) {
    stop_ifnot(inherits(object, allow_class), "argument `", obj_name,
      "` must be ", format_class(allow_class, plural = isTRUE(has_length > 1)),
      ", i.e. not be ", format_class(class(object), plural = isTRUE(has_length > 1)),
      call = call_depth
    )
    # check data.frames for data
    if (inherits(object, "data.frame")) {
      stop_if(any(dim(object) == 0),
        "the data provided in argument `", obj_name,
        "` must contain rows and columns (current dimensions: ",
        paste(dim(object), collapse = "x"), ")",
        call = call_depth
      )
    }
  }
  if (!is.null(has_length)) {
    stop_ifnot(length(object) %in% has_length, "argument `", obj_name,
      "` must ", # ifelse(allow_NULL, "be NULL or must ", ""),
      "be of length ", vector_or(has_length, quotes = FALSE),
      ", not ", length(object),
      call = call_depth
    )
  }
  if (!is.null(looks_like)) {
    stop_ifnot(object %like% looks_like, "argument `", obj_name,
      "` must ", # ifelse(allow_NULL, "be NULL or must ", ""),
      "resemble the regular expression \"", looks_like, "\"",
      call = call_depth
    )
  }
  if (!is.null(is_in)) {
    if (ignore.case == TRUE) {
      object <- tolower(object)
      is_in <- tolower(is_in)
    }
    stop_ifnot(all(object %in% is_in, na.rm = TRUE), "argument `", obj_name, "` ",
      ifelse(!is.null(has_length) && length(has_length) == 1 && has_length == 1,
        "must be either ",
        "must only contain values "
      ),
      vector_or(is_in, quotes = !isTRUE(any(c("double", "numeric", "integer") %in% allow_class))),
      ifelse(allow_NA == TRUE, ", or NA", ""),
      call = call_depth
    )
  }
  if (isTRUE(is_positive)) {
    stop_if(is.numeric(object) && !all(object > 0, na.rm = TRUE), "argument `", obj_name,
      "` must ",
      ifelse(!is.null(has_length) && length(has_length) == 1 && has_length == 1,
        "be a number higher than zero",
        "all be numbers higher than zero"
      ),
      call = call_depth
    )
  }
  if (isTRUE(is_positive_or_zero)) {
    stop_if(is.numeric(object) && !all(object >= 0, na.rm = TRUE), "argument `", obj_name,
      "` must ",
      ifelse(!is.null(has_length) && length(has_length) == 1 && has_length == 1,
        "be zero or a positive number",
        "all be zero or numbers higher than zero"
      ),
      call = call_depth
    )
  }
  if (isTRUE(is_finite)) {
    stop_if(is.numeric(object) && !all(is.finite(object[!is.na(object)]), na.rm = TRUE), "argument `", obj_name,
      "` must ",
      ifelse(!is.null(has_length) && length(has_length) == 1 && has_length == 1,
        "be a finite number",
        "all be finite numbers"
      ),
      " (i.e. not be infinite)",
      call = call_depth
    )
  }
  if (!is.null(contains_column_class)) {
    stop_ifnot(any(vapply(
      FUN.VALUE = logical(1),
      object,
      function(col, columns_class = contains_column_class) {
        inherits(col, columns_class)
      }
    ), na.rm = TRUE),
    "the data provided in argument `", obj_name,
    "` must contain at least one column of class <", contains_column_class, ">. ",
    "See ?as.", contains_column_class, ".",
    call = call_depth
    )
  }
  return(invisible())
}

get_current_data <- function(arg_name, call) {
  valid_df <- function(x) {
    !is.null(x) && is.data.frame(x)
  }
  # try dplyr::cur_data_all() first to support dplyr groups
  # only useful for e.g. dplyr::filter(), dplyr::mutate() and dplyr::summarise()
  # not useful (throws error) with e.g. dplyr::select() - but that will be caught later in this function
  cur_data_all <- import_fn("cur_data_all", "dplyr", error_on_fail = FALSE)
  if (!is.null(cur_data_all)) {
    out <- tryCatch(cur_data_all(), error = function(e) NULL)
    if (valid_df(out)) {
      return(out)
    }
  }

  # try a manual (base R) method, by going over all underlying environments with sys.frames()
  for (env in sys.frames()) {
    if (!is.null(env$`.Generic`)) {
      # don't check `".Generic" %in% names(env)`, because in R < 3.2, `names(env)` is always NULL

      if (valid_df(env$`.data`)) {
        # an element `.data` will be in the environment when using `dplyr::select()`
        # (but not when using `dplyr::filter()`, `dplyr::mutate()` or `dplyr::summarise()`)
        return(env$`.data`)
      } else if (valid_df(env$xx)) {
        # an element `xx` will be in the environment for rows + cols, e.g. `example_isolates[c(1:3), carbapenems()]`
        return(env$xx)
      } else if (valid_df(env$x)) {
        # an element `x` will be in the environment for only cols, e.g. `example_isolates[, carbapenems()]`
        return(env$x)
      }
    }
  }

  # no data.frame found, so an error  must be returned:
  if (is.na(arg_name)) {
    if (isTRUE(is.numeric(call))) {
      fn <- as.character(sys.call(call + 1)[1])
      examples <- paste0(
        ", e.g.:\n",
        "  your_data %>% select(", fn, "())\n",
        "  your_data %>% select(column_a, column_b, ", fn, "())\n",
        "  your_data[, ", fn, "()]\n",
        '  your_data[, c("column_a", "column_b", ', fn, "())]"
      )
    } else {
      examples <- ""
    }
    stop_("this function must be used inside a `dplyr` verb or `data.frame` call",
      examples,
      call = call
    )
  } else {
    # mimic a base R error that the argument is missing
    stop_("argument `", arg_name, "` is missing with no default", call = call)
  }
}

get_current_column <- function() {
  # try dplyr::cur_columns() first
  cur_column <- import_fn("cur_column", "dplyr", error_on_fail = FALSE)
  if (!is.null(cur_column)) {
    out <- tryCatch(cur_column(), error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
  }

  # cur_column() doesn't always work (only allowed for certain conditions set by dplyr), but it's probably still possible:
  frms <- lapply(sys.frames(), function(env) {
    if (!is.null(env$i)) {
      if (!is.null(env$tibble_vars)) {
        # for mutate_if()
        env$tibble_vars[env$i]
      } else {
        # for mutate(across())
        df <- tryCatch(get_current_data(NA, 0), error = function(e) NULL)
        if (is.data.frame(df)) {
          colnames(df)[env$i]
        } else {
          env$i
        }
      }
    } else {
      NULL
    }
  })

  vars <- unlist(frms)
  if (length(vars) > 0) {
    vars[length(vars)]
  } else {
    # not found, so:
    NULL
  }
}

is_null_or_grouped_tbl <- function(x) {
  # class "grouped_df" might change at one point, so only set in one place; here.
  is.null(x) || inherits(x, "grouped_df")
}

unique_call_id <- function(entire_session = FALSE, match_fn = NULL) {
  if (entire_session == TRUE) {
    return(c(envir = "session", call = "session"))
  }

  # combination of environment ID (such as "0x7fed4ee8c848")
  # and relevant system call (where 'match_fn' is being called in)
  calls <- sys.calls()
  in_test <- any(as.character(calls[[1]]) %like_case% "run_test_dir|run_test_file|test_all|tinytest|test_package|testthat", na.rm = TRUE)
  if (!isTRUE(in_test)) {
    for (i in seq_len(length(calls))) {
      call_clean <- gsub("[^a-zA-Z0-9_().-]", "", as.character(calls[[i]]), perl = TRUE)
      if (any(call_clean %like% paste0(match_fn, "\\("), na.rm = TRUE)) {
        return(c(
          envir = gsub("<environment: (.*)>", "\\1", utils::capture.output(sys.frames()[[1]]), perl = TRUE),
          call = paste0(deparse(calls[[i]]), collapse = "")
        ))
      }
    }
  }
  c(
    envir = paste0(sample(c(0:9, letters[1:6]), size = 32, replace = TRUE), collapse = ""),
    call = paste0(sample(c(0:9, letters[1:6]), size = 32, replace = TRUE), collapse = "")
  )
}

#' @noRd
#' @param fn name of the function as a character
#' @param ... character elements to be pasted together as a 'salt'
#' @param entire_session show message once per session
message_not_thrown_before <- function(fn, ..., entire_session = FALSE) {
  # this is to prevent that messages/notes will be printed for every dplyr group or more than once per session
  # e.g. this would show a msg 4 times: example_isolates %>% group_by(ward) %>% filter(mo_is_gram_negative())
  salt <- gsub("[^a-zA-Z0-9|_-]", "?", substr(paste(c(...), sep = "|", collapse = "|"), 1, 512), perl = TRUE)
  not_thrown_before <- is.null(AMR_env[[paste0("thrown_msg.", fn, ".", salt)]]) ||
    !identical(
      AMR_env[[paste0("thrown_msg.", fn, ".", salt)]],
      unique_call_id(
        entire_session = entire_session,
        match_fn = fn
      )
    )
  if (isTRUE(not_thrown_before)) {
    # message was not thrown before - remember this so on the next run it will return FALSE:
    assign(
      x = paste0("thrown_msg.", fn, ".", salt),
      value = unique_call_id(entire_session = entire_session, match_fn = fn),
      envir = AMR_env
    )
  }
  not_thrown_before
}

has_colour <- function() {
  # this is a base R version of crayon::has_color, but disables colours on emacs

  if (Sys.getenv("EMACS") != "" || Sys.getenv("INSIDE_EMACS") != "") {
    # disable on emacs, which only supports 8 colours
    return(FALSE)
  }
  enabled <- getOption("crayon.enabled")
  if (!is.null(enabled)) {
    return(isTRUE(enabled))
  }
  rstudio_with_ansi_support <- function(x) {
    if (Sys.getenv("RSTUDIO", "") == "") {
      return(FALSE)
    }
    if ((cols <- Sys.getenv("RSTUDIO_CONSOLE_COLOR", "")) != "" && !is.na(as.double(cols))) {
      return(TRUE)
    }
    tryCatch(getExportedValue("isAvailable", ns = asNamespace("rstudioapi"))(), error = function(e) {
      return(FALSE)
    }) &&
      tryCatch(getExportedValue("hasFun", ns = asNamespace("rstudioapi"))("getConsoleHasColor"), error = function(e) {
        return(FALSE)
      })
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
  if ("COLORTERM" %in% names(Sys.getenv())) {
    return(TRUE)
  }
  if (Sys.getenv("TERM") == "dumb") {
    return(FALSE)
  }
  grepl(
    pattern = "^screen|^xterm|^vt100|color|ansi|cygwin|linux",
    x = Sys.getenv("TERM"),
    ignore.case = TRUE,
    perl = TRUE
  )
}

# set colours if console has_colour()
try_colour <- function(..., before, after, collapse = " ") {
  if (length(c(...)) == 0) {
    return(character(0))
  }
  txt <- paste0(c(...), collapse = collapse)
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
  before <- "\033[38;5;232m"
  after <- "\033[39m"
  theme_info <- import_fn("getThemeInfo", "rstudioapi", error_on_fail = FALSE)
  if (!is.null(theme_info) && isTRUE(theme_info()$dark)) {
    # white
    before <- "\033[37m"
    after <- "\033[39m"
  }
  try_colour(..., before = before, after = after, collapse = collapse)
}
font_white <- function(..., collapse = " ") {
  before <- "\033[37m"
  after <- "\033[39m"
  theme_info <- import_fn("getThemeInfo", "rstudioapi", error_on_fail = FALSE)
  if (!is.null(theme_info) && isTRUE(theme_info()$dark)) {
    # black
    before <- "\033[38;5;232m"
    after <- "\033[39m"
  }
  try_colour(..., before = before, after = after, collapse = collapse)
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
font_yellow <- function(..., collapse = " ") {
  try_colour(..., before = "\033[33m", after = "\033[39m", collapse = collapse)
}
font_subtle <- function(..., collapse = " ") {
  try_colour(..., before = "\033[38;5;246m", after = "\033[39m", collapse = collapse)
}
font_grey <- function(..., collapse = " ") {
  try_colour(..., before = "\033[38;5;249m", after = "\033[39m", collapse = collapse)
}
font_grey_bg <- function(..., collapse = " ") {
  if (tryCatch(import_fn("getThemeInfo", "rstudioapi", error_on_fail = FALSE)()$dark, error = function(e) FALSE)) {
    # similar to HTML #444444
    try_colour(..., before = "\033[48;5;238m", after = "\033[49m", collapse = collapse)
  } else {
    # similar to HTML #f0f0f0
    try_colour(..., before = "\033[48;5;255m", after = "\033[49m", collapse = collapse)
  }
}
font_red_bg <- function(..., collapse = " ") {
  # this is #ed553b (picked to be colourblind-safe with other RSI colours)
  try_colour(font_black(..., collapse = collapse), before = "\033[48;5;203m", after = "\033[49m", collapse = collapse)
}
font_orange_bg <- function(..., collapse = " ") {
  # this is #f6d55c (picked to be colourblind-safe with other RSI colours)
  try_colour(font_black(..., collapse = collapse), before = "\033[48;5;222m", after = "\033[49m", collapse = collapse)
}
font_yellow_bg <- function(..., collapse = " ") {
  try_colour(font_black(..., collapse = collapse), before = "\033[48;5;228m", after = "\033[49m", collapse = collapse)
}
font_green_bg <- function(..., collapse = " ") {
  # this is #3caea3 (picked to be colourblind-safe with other RSI colours)
  try_colour(font_black(..., collapse = collapse), before = "\033[48;5;79m", after = "\033[49m", collapse = collapse)
}
font_purple_bg <- function(..., collapse = " ") {
  try_colour(font_black(..., collapse = collapse), before = "\033[48;5;89m", after = "\033[49m", collapse = collapse)
}
font_rose_bg <- function(..., collapse = " ") {
  try_colour(font_black(..., collapse = collapse), before = "\033[48;5;217m", after = "\033[49m", collapse = collapse)
}
font_na <- function(..., collapse = " ") {
  font_red(..., collapse = collapse)
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

progress_ticker <- function(n = 1, n_min = 0, print = TRUE, ...) {
  if (print == FALSE || n < n_min) {
    pb <- list()
    pb$tick <- function() {
      invisible()
    }
    pb$kill <- function() {
      invisible()
    }
    set_clean_class(pb, new_class = "txtProgressBar")
  } else if (n >= n_min) {
    # rely on the progress package if it is available - it has a more verbose output
    progress_bar <- import_fn("progress_bar", "progress", error_on_fail = FALSE)
    if (!is.null(progress_bar)) {
      # so we use progress::progress_bar
      # a close() method was also added, see below this function
      pb <- progress_bar$new(
        format = "[:bar] :percent (:current/:total)",
        total = n
      )
    } else {
      pb <- utils::txtProgressBar(max = n, style = 3)
      pb$tick <- function() {
        pb$up(pb$getVal() + 1)
      }
    }
    pb
  }
}

#' @method close progress_bar
#' @export
#' @noRd
close.progress_bar <- function(con, ...) {
  con$terminate()
}

set_clean_class <- function(x, new_class) {
  # return the object with only the new class and no additional attributes where possible
  if (is.null(x)) {
    x <- NA_character_
  }
  if (is.factor(x)) {
    # keep only levels and remove all other attributes
    lvls <- levels(x)
    attributes(x) <- NULL
    levels(x) <- lvls
  } else if (!is.list(x) && !is.function(x)) {
    attributes(x) <- NULL
  }
  class(x) <- new_class
  x
}

formatted_filesize <- function(...) {
  size_kb <- file.size(...) / 1024
  if (size_kb < 1) {
    paste(round(size_kb, 1), "kB")
  } else if (size_kb < 100) {
    paste(round(size_kb, 0), "kB")
  } else {
    paste(round(size_kb / 1024, 1), "MB")
  }
}

create_pillar_column <- function(x, ...) {
  new_pillar_shaft_simple <- import_fn("new_pillar_shaft_simple", "pillar")
  new_pillar_shaft_simple(x, ...)
}

as_original_data_class <- function(df, old_class = NULL) {
  if ("tbl_df" %in% old_class && pkg_is_available("tibble", also_load = FALSE)) {
    fn <- import_fn("as_tibble", "tibble")
  } else if ("tbl_ts" %in% old_class && pkg_is_available("tsibble", also_load = FALSE)) {
    fn <- import_fn("as_tsibble", "tsibble")
  } else if ("data.table" %in% old_class && pkg_is_available("data.table", also_load = FALSE)) {
    fn <- import_fn("as.data.table", "data.table")
  } else if ("tabyl" %in% old_class && pkg_is_available("janitor", also_load = FALSE)) {
    fn <- import_fn("as_tabyl", "janitor")
  } else {
    fn <- base::as.data.frame
  }
  fn(df)
}

# works exactly like round(), but rounds `round2(44.55, 1)` to 44.6 instead of 44.5
# and adds decimal zeroes until `digits` is reached when force_zero = TRUE
round2 <- function(x, digits = 1, force_zero = TRUE) {
  x <- as.double(x)
  # https://stackoverflow.com/a/12688836/4575331
  val <- (trunc((abs(x) * 10^digits) + 0.5) / 10^digits) * sign(x)
  if (digits > 0 && force_zero == TRUE) {
    values_trans <- val[val != as.integer(val) & !is.na(val)]
    val[val != as.integer(val) & !is.na(val)] <- paste0(
      values_trans,
      strrep(
        "0",
        max(
          0,
          digits - nchar(
            format(
              as.double(
                gsub(
                  ".*[.](.*)$",
                  "\\1",
                  values_trans
                )
              ),
              scientific = FALSE
            )
          )
        )
      )
    )
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
    max_places <- max(unlist(lapply(
      strsplit(sub(
        "0+$", "",
        as.character(x * 100)
      ), ".", fixed = TRUE),
      function(y) ifelse(length(y) == 2, nchar(y[2]), 0)
    )), na.rm = TRUE)
    max(min(max_places,
      maximum,
      na.rm = TRUE
    ),
    minimum,
    na.rm = TRUE
    )
  }

  # format_percentage() function
  format_percentage <- function(x, digits = NULL, ...) {
    if (is.null(digits)) {
      digits <- getdecimalplaces(x)
    }
    if (is.null(digits) || is.na(digits) || !is.numeric(digits)) {
      digits <- 2
    }

    # round right: percentage(0.4455) and format(as.percentage(0.4455), 1) should return "44.6%", not "44.5%"
    x_formatted <- format(round2(as.double(x), digits = digits + 2) * 100,
      scientific = FALSE,
      digits = max(1, digits),
      nsmall = digits,
      ...
    )
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
  format_percentage(structure(
    .Data = as.double(x),
    class = c("percentage", "numeric")
  ),
  digits = digits, ...
  )
}

time_start_tracking <- function() {
  AMR_env$time_start <- round(as.double(Sys.time()) * 1000)
}

time_track <- function(name = NULL) {
  paste("(until now:", trimws(round(as.double(Sys.time()) * 1000) - AMR_env$time_start), "ms)")
}

trimws2 <- function(..., whitespace = "[\u0009\u000A\u000B\u000C\u000D\u0020\u0085\u00A0\u1680\u180E\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200A\u200B\u200C\u200D\u2028\u2029\u202F\u205F\u2060\u3000\uFEFF]") {
  # this is even faster than trimws() itself which sets " \t\n\r".
  trimws(..., whitespace = whitespace)
}


# Faster data.table implementations ----

match <- function(x, ...) {
  if (isTRUE(AMR_env$has_data.table) && is.character(x)) {
    # data.table::chmatch() is 35% faster than base::match() for character
    getExportedValue(name = "chmatch", ns = asNamespace("data.table"))(x, ...)
  } else {
    base::match(x, ...)
  }
}
`%in%` <- function(x, table) {
  if (isTRUE(AMR_env$has_data.table) && is.character(x) && is.character(table)) {
    # data.table::`%chin%`() is 20-50% faster than base::`%in%`() for character
    getExportedValue(name = "%chin%", ns = asNamespace("data.table"))(x, table)
  } else {
    base::`%in%`(x, table)
  }
}

# nolint start

# Register S3 methods ----
# copied from vctrs::s3_register by their permission:
# https://github.com/r-lib/vctrs/blob/05968ce8e669f73213e3e894b5f4424af4f46316/R/register-s3.R
s3_register <- function(generic, class, method = NULL) {
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)
  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]
  caller <- parent.frame()
  get_method_env <- function() {
    top <- topenv(caller)
    if (isNamespace(top)) {
      asNamespace(environmentName(top))
    } else {
      caller
    }
  }
  get_method <- function(method, env) {
    if (is.null(method)) {
      get(paste0(generic, ".", class), envir = get_method_env())
    } else {
      method
    }
  }
  method_fn <- get_method(method)
  stopifnot(is.function(method_fn))
  setHook(packageEvent(package, "onLoad"), function(...) {
    ns <- asNamespace(package)
    method_fn <- get_method(method)
    registerS3method(generic, class, method_fn, envir = ns)
  })
  if (!isNamespaceLoaded(package)) {
    return(invisible())
  }
  envir <- asNamespace(package)
  if (exists(generic, envir)) {
    registerS3method(generic, class, method_fn, envir = envir)
  }
  invisible()
}


# Support old R versions ----
# these functions were not available in previous versions of R
# see here for the full list: https://github.com/r-lib/backports
if (getRversion() < "3.1.0") {
  # R-3.0 does not contain these functions, set them here to prevent installation failure
  # (required for extension of the <mic> class)
  cospi <- function(...) 1
  sinpi <- function(...) 1
  tanpi <- function(...) 1
}

if (getRversion() < "3.2.0") {
  anyNA <- function(x, recursive = FALSE) {
    if (isTRUE(recursive) && (is.list(x) || is.pairlist(x))) {
      return(any(rapply(x, anyNA, how = "unlist", recursive = FALSE)))
    }
    any(is.na(x))
  }
  dir.exists <- function(paths) {
    x <- base::file.info(paths)$isdir
    !is.na(x) & x
  }
  file.size <- function(...) {
    file.info(...)$size
  }
  file.mtime <- function(...) {
    file.info(...)$mtime
  }
  isNamespaceLoaded <- function(pkg) {
    pkg %in% loadedNamespaces()
  }
  lengths <- function(x, use.names = TRUE) {
    vapply(x, length, FUN.VALUE = NA_integer_, USE.NAMES = use.names)
  }
}

if (getRversion() < "3.3.0") {
  strrep <- function(x, times) {
    x <- as.character(x)
    if (length(x) == 0L) {
      return(x)
    }
    unlist(.mapply(function(x, times) {
      if (is.na(x) || is.na(times)) {
        return(NA_character_)
      }
      if (times <= 0L) {
        return("")
      }
      paste0(replicate(times, x), collapse = "")
    }, list(x = x, times = times), MoreArgs = list()), use.names = FALSE)
  }
}

if (getRversion() < "3.5.0") {
  isFALSE <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x) && !x
  }
  # trims() was introduced in 3.3.0, but its argument `whitespace` only in 3.5.0
  trimws <- function(x, which = c("both", "left", "right"), whitespace = "[ \t\r\n]") {
    which <- match.arg(which)
    mysub <- function(re, x) sub(re, "", x, perl = TRUE)
    switch(which,
      left = mysub(paste0("^", whitespace, "+"), x),
      right = mysub(paste0(whitespace, "+$"), x),
      both = mysub(paste0(whitespace, "+$"), mysub(paste0("^", whitespace, "+"), x))
    )
  }
}

if (getRversion() < "3.6.0") {
  str2lang <- function(s) {
    stopifnot(length(s) == 1L)
    ex <- parse(text = s, keep.source = FALSE)
    stopifnot(length(ex) == 1L)
    ex[[1L]]
  }
}

if (getRversion() < "4.0.0") {
  deparse1 <- function(expr, collapse = " ", width.cutoff = 500L, ...) {
    paste(deparse(expr, width.cutoff, ...), collapse = collapse)
  }
}

# nolint end
