# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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
# how to conduct AMR data analysis: https://amr-for-r.org              #
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
    y[
      match(
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

# support where() like tidyverse (this function will also be used when running `antibiogram()`):
where <- function(fn) {
  # based on https://github.com/nathaneastwood/poorman/blob/52eb6947e0b4430cd588976ed8820013eddf955f/R/where.R#L17-L32
  if (!is.function(fn)) {
    stop_("`", deparse(substitute(fn)), "()` is not a valid predicate function.")
  }
  df <- pm_select_env$.data
  cols <- pm_select_env$get_colnames()
  if (is.null(df)) {
    df <- get_current_data("where", call = FALSE)
    cols <- colnames(df)
  }
  preds <- unlist(lapply(
    df,
    function(x, fn) {
      do.call("fn", list(x))
    },
    fn
  ))
  if (!is.logical(preds)) stop_("`where()` must be used with functions that return `TRUE` or `FALSE`.")
  data_cols <- cols
  cols <- data_cols[preds]
  which(data_cols %in% cols)
}

# copied and slightly rewritten from {poorman} under permissive license (2021-10-15)
# https://github.com/nathaneastwood/poorman, MIT licensed, Nathan Eastwood, 2020
case_when_AMR <- function(...) {
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

rbind_AMR <- function(...) {
  # this is just rbind(), but with the functionality of dplyr::bind_rows(),
  # to allow differences in available columns
  l <- list(...)
  l_names <- unique(unlist(lapply(l, names)))
  l_new <- lapply(l, function(df) {
    rownames(df) <- NULL
    for (col in l_names[!l_names %in% colnames(df)]) {
      # create the new column, could also be length 0
      df[, col] <- rep(NA, NROW(df))
    }
    df
  })
  do.call(rbind, l_new)
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
    modifyRange(
      document_range(
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

search_type_in_df <- function(x, type, info = TRUE, add_col_prefix = TRUE) {
  meet_criteria(x, allow_class = "data.frame")
  meet_criteria(type, allow_class = "character", has_length = 1)

  # try to find columns based on type
  found <- NULL

  # remove attributes from other packages
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  colnames_formatted <- tolower(generalise_antibiotic_name(colnames(x)))

  # -- mo
  if (type == "mo") {
    add_MO_lookup_to_AMR_env()

    if (any(vapply(FUN.VALUE = logical(1), x, is.mo))) {
      # take first 'mo' column
      found <- colnames(x)[vapply(FUN.VALUE = logical(1), x, is.mo)]
    } else if ("mo" %in% colnames_formatted &&
      suppressWarnings(all(x$mo %in% c(NA, AMR_env$MO_lookup$mo)))) {
      found <- "mo"
    } else if (any(colnames_formatted %like_case% "^(mo|microorganism|organism|bacteria|ba[ck]terie)s?$")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "^(mo|microorganism|organism|bacteria|ba[ck]terie)s?$"])
    } else if (any(colnames_formatted %like_case% "^(microorganism|organism|bacteria|ba[ck]terie)")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "^(microorganism|organism|bacteria|ba[ck]terie)"])
    } else if (any(colnames_formatted %like_case% "species")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "species"])
    }
  }
  # -- key antimicrobials
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
        stop(
          font_red(paste0(
            "Found column '", font_bold(found), "' to be used as input for `", ifelse(add_col_prefix, "col_", ""), type,
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
  # -- host (animals)
  if (type == "host") {
    if (any(colnames_formatted %like_case% "^(host|animal)")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "^(host|animal)"])
    } else if (any(colnames_formatted %like_case% "((^|[^A-Za-z])host($|[^A-Za-z])|animal)")) {
      found <- sort(colnames(x)[colnames_formatted %like_case% "((^|[^A-Za-z])host($|[^A-Za-z])|animal)"])
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
        message_("Column '", font_bold(found), "' found as input for `", ifelse(add_col_prefix, "col_", ""), type,
          "`, but this column does not contain 'logical' values (TRUE/FALSE) and was ignored.",
          add_fn = font_red
        )
        found <- NULL
      }
    }
  }

  found <- found[1]

  if (!is.null(found) && isTRUE(info)) {
    if (message_not_thrown_before("search_in_type", type)) {
      msg <- paste0("Using column '", font_bold(found), "' as input for `", ifelse(add_col_prefix, "col_", ""), type, "`.")
      if (type %in% c("keyantibiotics", "keyantimicrobials", "specimen")) {
        msg <- paste(msg, "Use", font_bold(paste0(ifelse(add_col_prefix, "col_", ""), type), "= FALSE"), "to prevent this.")
      }
      message_(msg)
    }
  }
  found
}

is_valid_regex <- function(x) {
  regex_at_all <- tryCatch(
    vapply(
      FUN.VALUE = logical(1),
      X = strsplit(x, "", fixed = TRUE),
      FUN = function(y) {
        any(
          y %in% c(
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
  installed <- vapply(FUN.VALUE = logical(1), package, requireNamespace, quietly = TRUE)
  if (any(!installed) && any(package == "rstudioapi")) {
    stop("This function only works in RStudio when using R >= 3.2.", call. = FALSE)
  } else if (any(!installed)) {
    stop("This requires the ", vector_and(package[!installed]), " package.",
      "\nTry to install with install.packages().",
      call. = FALSE
    )
  } else {
    return(invisible())
  }
}

pkg_is_available <- function(pkg, also_load = FALSE, min_version = NULL) {
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
        stop_("function `", name, "()` is not an exported object from package '", pkg,
          "'. Please create an issue at ", font_url("https://github.com/msberends/AMR/issues"), ". Many thanks!",
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
    return(paste0(
      vapply(
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
  msg_stripped <- gsub("(.*)?\\033\\]8;;.*\\a(.*?)\\033\\]8;;\\a(.*)", "\\1\\2\\3", msg, perl = TRUE) # for font_url()
  msg_stripped <- font_stripstyle(msg_stripped)
  # where are the spaces now?
  msg_stripped_wrapped <- paste0(
    strwrap(msg_stripped,
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
  if (pkg_is_available("cli") &&
    tryCatch(isTRUE(getExportedValue("ansi_has_hyperlink_support", ns = asNamespace("cli"))()), error = function(e) FALSE) &&
    tryCatch(getExportedValue("isAvailable", ns = asNamespace("rstudioapi"))(), error = function(e) {
      return(FALSE)
    }) &&
    tryCatch(getExportedValue("versionInfo", ns = asNamespace("rstudioapi"))()$version > "2023.6.0.0", error = function(e) {
      return(FALSE)
    })) {
    # we are in a recent version of RStudio, so do something nice: add links to our help pages in the console.
    parts <- strsplit(msg, "`", fixed = TRUE)[[1]]
    cmds <- parts %in% paste0(ls(envir = asNamespace("AMR")), "()")
    # functions with a dot are not allowed: https://github.com/rstudio/rstudio/issues/11273#issuecomment-1156193252
    # lead them to the help page of our package
    parts[cmds & parts %like% "[.]"] <- font_url(
      url = paste0("ide:help:AMR::", gsub("()", "", parts[cmds & parts %like% "[.]"], fixed = TRUE)),
      txt = parts[cmds & parts %like% "[.]"]
    )
    # otherwise, give a 'click to run' popup
    parts[cmds & parts %unlike% "[.]"] <- font_url(
      url = paste0("ide:run:AMR::", parts[cmds & parts %unlike% "[.]"]),
      txt = parts[cmds & parts %unlike% "[.]"]
    )
    # datasets should give help page as well
    parts[parts %in% c("antimicrobials", "microorganisms", "microorganisms.codes", "microorganisms.groups")] <- font_url(
      url = paste0("ide:help:AMR::", gsub("()", "", parts[parts %in% c("antimicrobials", "microorganisms", "microorganisms.codes", "microorganisms.groups")], fixed = TRUE)),
      txt = parts[parts %in% c("antimicrobials", "microorganisms", "microorganisms.codes", "microorganisms.groups")]
    )
    # text starting with `?` must also lead to the help page
    parts[parts %like% "^[?].+"] <- font_url(
      url = paste0("ide:help:AMR::", gsub("?", "", parts[parts %like% "^[?].+"], fixed = TRUE)),
      txt = parts[parts %like% "^[?].+"]
    )
    msg <- paste0(parts, collapse = "`")
  }
  msg <- gsub("`(.+?)`", font_grey_bg("\\1"), msg)

  # clean introduced whitespace in between fullstops
  msg <- gsub("[.] +[.]", "..", msg)
  # remove extra space that was introduced (e.g. "Smith et al. , 2022")
  msg <- gsub(". ,", ".,", msg, fixed = TRUE)
  msg <- gsub("[ ,", "[,", msg, fixed = TRUE)
  msg <- gsub("/ /", "//", msg, fixed = TRUE)

  msg
}

message_ <- function(...,
                     appendLF = TRUE,
                     add_fn = list(font_blue),
                     as_note = TRUE) {
  message(
    word_wrap(...,
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
  warning(
    trimws2(word_wrap(...,
      add_fn = add_fn,
      as_note = FALSE
    )),
    immediate. = immediate,
    call. = call
  )
}

# this alternative to the stop() function:
# - adds the function name where the error was thrown
# - wraps text to never break lines within words
stop_ <- function(..., call = TRUE) {
  msg <- paste0(c(...), collapse = "")
  msg_call <- ""
  if (!isFALSE(call)) {
    if (isTRUE(call)) {
      call <- as.character(sys.call(-1)[1])
    } else {
      # so you can go back more than 1 call, as used in sir_calc(), that now throws a reference to e.g. n_sir()
      call <- as.character(sys.call(call)[1])
    }
    msg_call <- paste0("in ", call, "():")
  }
  msg <- trimws2(word_wrap(msg, add_fn = list(), as_note = FALSE))
  if (!is.null(AMR_env$cli_abort) && length(unlist(strsplit(msg, "\n", fixed = TRUE))) <= 1) {
    if (is.character(call)) {
      call <- as.call(str2lang(paste0(call, "()")))
    } else {
      call <- NULL
    }
    AMR_env$cli_abort(msg, call = call)
  } else {
    stop(paste(msg_call, msg), call. = FALSE)
  }
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
  day <- as.integer(format(d, "%e"))
  suffix <- rep("th", length(day))
  suffix[day %in% c(1, 21, 31)] <- "st"
  suffix[day %in% c(2, 22)] <- "nd"
  suffix[day %in% c(3, 23)] <- "rd"
  paste0(month.name[as.integer(format(d, "%m"))], " ", day, suffix, ", ", format(d, "%Y"))
}

format_included_data_number <- function(data) {
  if (is.numeric(data) && length(data) == 1) {
    n <- data
  } else if (is.data.frame(data)) {
    n <- nrow(data)
  } else {
    n <- length(unique(data))
  }
  if (n > 10000) {
    rounder <- -3 # round on thousands
  } else if (n > 1000) {
    rounder <- -2 # round on hundreds
  } else if (n < 50) {
    # do not round
    rounder <- 0
  } else {
    rounder <- -1 # round on tens
  }
  paste0(ifelse(rounder == 0, "", "~"), format(round(n, rounder), decimal.mark = ".", big.mark = " "))
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
    # class 'sir' should be sorted like this
    v <- c("S", "I", "R")
  }
  if (identical(v, c("I", "NI", "R", "S", "SDD"))) {
    # class 'sir' should be sorted like this
    v <- c("S", "SDD", "I", "R", "NI")
  }
  # oxford comma
  if (last_sep %in% c(" or ", " and ") && length(v) > 2) {
    last_sep <- paste0(",", last_sep)
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
  if (any(c("mo", "ab", "sir") %in% class)) {
    class <- paste0("of class '", class[1L], "'")
  }
  class[class == class.bak] <- paste0("of class '", class[class == class.bak], "'")
  # output
  vector_or(class, quotes = FALSE, sort = FALSE)
}

# a check for every single argument in all functions
meet_criteria <- function(object, # can be literally `list(...)` for `allow_arguments_from`
                          allow_class = NULL,
                          has_length = NULL,
                          looks_like = NULL,
                          is_in = NULL,
                          is_positive = NULL,
                          is_positive_or_zero = NULL,
                          is_finite = NULL,
                          allow_NULL = FALSE,
                          allow_NA = FALSE,
                          ignore.case = FALSE,
                          allow_arguments_from = NULL, # 1 function, or a list of functions
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

  if (identical(class(object), "list") && !"list" %in% allow_class) {
    # coming from Python, possibly - turn lists (not data.frame) to the underlying data type
    object <- unlist(object)
  }

  if (!is.null(allow_class) && !(suppressWarnings(all(is.na(object))) && allow_NA == TRUE)) {
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
    is_in.bak <- is_in
    if ("logical" %in% allow_class) {
      is_in <- is_in[!is_in %in% c("TRUE", "FALSE")]
    }
    or_values <- vector_or(is_in, quotes = !isTRUE(any(c("numeric", "integer") %in% allow_class)))
    if ("logical" %in% allow_class) {
      or_values <- paste0(or_values, ", or TRUE or FALSE")
    }
    stop_ifnot(all(object %in% is_in.bak, na.rm = TRUE), "argument `", obj_name, "` ",
      ifelse(!is.null(has_length) && length(has_length) == 1 && has_length == 1,
        "must be either ",
        "must only contain values "
      ),
      or_values,
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
  if (!is.null(allow_arguments_from) && !is.null(names(object))) {
    args_given <- names(object)
    if (is.function(allow_arguments_from)) {
      allow_arguments_from <- list(allow_arguments_from)
    }
    args_allowed <- sort(unique(unlist(lapply(allow_arguments_from, function(x) names(formals(x))))))
    args_allowed <- args_allowed[args_allowed != "..."]
    disallowed <- args_given[!args_given %in% args_allowed]
    stop_if(length(disallowed) > 0,
      ifelse(length(disallowed) == 1,
        paste("the argument", vector_and(disallowed), "is"),
        paste("the arguments", vector_and(disallowed), "are")
      ),
      " not valid. Valid arguments are: ",
      vector_and(args_allowed), ".",
      call = call_depth
    )
  }
  return(invisible())
}

ascertain_sir_classes <- function(x, obj_name) {
  sirs <- vapply(FUN.VALUE = logical(1), x, is.sir)
  if (!any(sirs, na.rm = TRUE)) {
    warning_(
      "the data provided in argument `", obj_name,
      "` should contain at least one column of class 'sir'. Eligible SIR column were now guessed. ",
      "See `?as.sir`.",
      immediate = TRUE
    )
    sirs_eligible <- is_sir_eligible(x)
    for (col in colnames(x)[sirs_eligible]) {
      x[[col]] <- as.sir(x[[col]])
    }
  }
  x
}

get_current_data <- function(arg_name, call) {
  # This function enables AMR selectors (e.g., AMR::carbapenems()) to work seamlessly across different environments, including dplyr, base R, data.table, and tidymodels.
  # It identifies and extracts the appropriate data frame from the current execution context.
  valid_df <- function(x) {
    !is.null(x) && is.data.frame(x)
  }

  frms <- sys.frames()

  # check dplyr environments to support dplyr groups
  with_mask <- vapply(FUN.VALUE = logical(1), frms, function(e) !is.null(e$mask))
  for (env in frms[which(with_mask)]) {
    if (is.function(env$mask$current_rows) && (valid_df(env$data) || valid_df(env$`.data`))) {
      # an element `.data` or `data` (containing all data) and `mask` (containing functions) will be in the environment when using dplyr verbs
      # we use their mask$current_rows() below to get the group rows, since dplyr::cur_data_all() is deprecated and will be removed in the future
      # e.g. for `example_isolates %>% group_by(ward) %>% mutate(first = first_isolate(.))`
      if (valid_df(env$data)) {
        # support for dplyr 1.1.x
        df <- env$data
      } else {
        # support for dplyr 1.0.x
        df <- env$`.data`
      }
      rows <- tryCatch(env$mask$current_rows(), error = function(e) seq_len(NROW(df)))
      return(df[rows, , drop = FALSE])
    }
  }

  # now go over all underlying environments looking for other dplyr, tidymodels, data.table and base R selection environments
  with_generic <- vapply(FUN.VALUE = logical(1), frms, function(e) !is.null(e$`.Generic`))
  for (env in frms[which(with_generic)]) {
    if (valid_df(env$`.data`)) {
      # an element `.data` will be in the environment when using dplyr::select()
      return(env$`.data`)
    } else if (valid_df(env$training)) {
      # an element `training` will be in the environment when using some tidymodels functions such as `prep()`
      return(env$training)
    } else if (valid_df(env$data)) {
      # an element `data` will be in the environment when using older dplyr versions, or some tidymodels functions such as `fit()`
      return(env$data)
    } else if (valid_df(env$xx)) {
      # an element `xx` will be in the environment for rows + cols in base R, e.g. `example_isolates[c(1:3), carbapenems()]`
      return(env$xx)
    } else if (valid_df(env$x)) {
      # an element `x` will be in the environment for only cols in base R, e.g. `example_isolates[, carbapenems()]`
      # this element will also be present in data.table environments where there's a .Generic available
      return(env$x)
    }
  }

  # now a special case for dplyr's 'scoped' variants
  with_tbl <- vapply(FUN.VALUE = logical(1), frms, function(e) valid_df(e$`.tbl`))
  for (env in frms[which(with_tbl)]) {
    if (!is.null(names(env)) && all(c(".tbl", ".vars", ".cols") %in% names(env), na.rm = TRUE)) {
      # an element `.tbl` will be in the environment when using scoped dplyr variants, with or without `dplyr::vars()`
      # e.g. `dplyr::summarise_at(carbapenems(), ...)` or `dplyr::mutate_at(vars(carbapenems()), ...)`
      return(env$`.tbl`)
    }
  }

  # now check if it was run with eval(), which has arguments `expr`, `envir`, and `enclos`
  from_eval_parse <- vapply(FUN.VALUE = logical(1), frms, function(e) all(c("expr", "envir", "enclos") %in% names(e)))
  for (env in frms[which(from_eval_parse)]) {
    if (valid_df(env$envir)) {
      # the element `envir` could contain the data in case of
      # e.g. `eval(parse(text = "any(cephalosporins_3rd() == 'R')"), envir = example_isolates)`
      # this is also used by run_custom_mdro_guideline() to support antimicrobial selectors in the part before `~`
      return(env$envir)
    }
  }

  # no data.frame found, so an error  must be returned:
  if (is.na(arg_name)) {
    if (isTRUE(is.numeric(call))) {
      fn <- as.character(sys.call(call + 1)[1])
      examples <- paste0(
        ", e.g.:\n",
        " ", AMR_env$bullet_icon, " your_data %>% select(", fn, "())\n",
        " ", AMR_env$bullet_icon, " your_data %>% select(column_a, column_b, ", fn, "())\n",
        " ", AMR_env$bullet_icon, " your_data %>% filter(any(", fn, "() == \"R\"))\n",
        " ", AMR_env$bullet_icon, " your_data[, ", fn, "()]\n",
        " ", AMR_env$bullet_icon, " your_data[, c(\"column_a\", \"column_b\", ", fn, "())]"
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
  out <- tryCatch(cur_column(), error = function(e) NULL)
  if (!is.null(out)) {
    return(out)
  }

  # cur_column() doesn't always work (only allowed for certain conditions set by dplyr), but it's probably still possible:
  frms <- lapply(sys.frames(), function(env) {
    if (tryCatch(!is.null(env$i), error = function(e) FALSE)) {
      if (!is.null(env$tibble_vars)) {
        # for mutate_if()
        # TODO remove later, was part of older dplyr versions (at least not in dplyr 1.1.4)
        env$tibble_vars[env$i]
      } else {
        # for mutate(across())
        if (!is.null(env$data) && is.data.frame(env$data)) {
          df <- env$data
        } else {
          df <- tryCatch(get_current_data(NA, 0), error = function(e) NULL)
        }
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
  # class "grouped_data" is from {poorman}, see aa_helper_pm_functions.R
  # class "grouped_df" is from {dplyr} and might change at one point, so only set in one place; here.
  is.null(x) || inherits(x, "grouped_data") || inherits(x, "grouped_df")
}

get_group_names <- function(x) {
  if ("pm_groups" %in% names(attributes(x))) {
    pm_get_groups(x)
  } else if (!is.null(x) && is_null_or_grouped_tbl(x)) {
    grps <- colnames(attributes(x)$groups)
    grps[!grps %in% c(".group_id", ".rows")]
  } else {
    character(0)
  }
}

format_custom_query_rule <- function(query, colours = has_colour()) {
  # this is used by custom EUCAST and custom MDRO rules

  # font_black() is a bit expensive so do it once:
  txt <- font_black("{text}")
  query <- gsub(" & ", sub("{text}", font_bold(" and "), txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" | ", sub("{text}", " or ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" + ", sub("{text}", " plus ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" - ", sub("{text}", " minus ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" / ", sub("{text}", " divided by ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" * ", sub("{text}", " times ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" == ", sub("{text}", " is ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" > ", sub("{text}", " is higher than ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" < ", sub("{text}", " is lower than ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" >= ", sub("{text}", " is higher than or equal to ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" <= ", sub("{text}", " is lower than or equal to ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" ^ ", sub("{text}", " to the power of ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" %in% ", sub("{text}", " is one of ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub(" %like% ", sub("{text}", " resembles ", txt, fixed = TRUE), query, fixed = TRUE)
  query <- gsub("any\\((.*)\\)$", paste0(font_black("any of "), "\\1"), query)
  query <- gsub("all\\((.*)\\)$", paste0(font_black("all of "), "\\1"), query)
  if (colours == TRUE) {
    query <- gsub("[\"']R[\"']", font_rose_bg(" R "), query)
    query <- gsub("[\"']SDD[\"']", font_orange_bg(" SDD "), query)
    query <- gsub("[\"']S[\"']", font_green_bg(" S "), query)
    query <- gsub("[\"']NI[\"']", font_grey_bg(font_black(" NI ")), query)
    query <- gsub("[\"']I[\"']", font_orange_bg(" I "), query)
  }
  # replace the black colour 'stops' with blue colour 'starts'
  query <- gsub("\033[39m", "\033[34m", as.character(query), fixed = TRUE)
  # start with blue
  query <- paste0("\033[34m", query)
  if (colours == FALSE) {
    query <- font_stripstyle(query)
  }
  query
}

unique_call_id <- function(entire_session = FALSE, match_fn = NULL) {
  if (entire_session == TRUE) {
    return(c(envir = "session", call = "session"))
  }

  # combination of environment ID (such as "0x7fed4ee8c848")
  # and relevant system call (where 'match_fn' is being called in)
  calls <- sys.calls()
  in_test <- any(as.character(calls[[1]]) %like_case% "run_test_dir|run_test_file|test_all|tinytest|test_package|testthat", na.rm = TRUE)
  if (!isTRUE(in_test) && !is.null(match_fn)) {
    for (i in seq_len(length(calls))) {
      call_clean <- gsub("[^a-zA-Z0-9_().-]", "", as.character(calls[[i]]), perl = TRUE)
      if (match_fn %in% call_clean || any(call_clean %like% paste0(match_fn, "\\("), na.rm = TRUE)) {
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
#' @param fn Name of the function as a character.
#' @param ... Character elements to be pasted together as a 'salt'.
#' @param entire_session Show message once per session.
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

reset_all_thrown_messages <- function() {
  rm(
    list = grep("^thrown_msg", ls(envir = AMR_env), value = TRUE),
    envir = AMR_env
  )
}

has_colour <- function() {
  if (is.null(AMR_env$supports_colour)) {
    if (Sys.getenv("EMACS") != "" || Sys.getenv("INSIDE_EMACS") != "") {
      # disable on emacs, which only supports 8 colours
      AMR_env$supports_colour <- FALSE
    } else {
      has_color <- import_fn("has_color", "crayon", error_on_fail = FALSE)
      AMR_env$supports_colour <- !is.null(has_color) && isTRUE(has_color())
    }
  }
  # always FALSE for GitHub documents (`index.Rmd` and `README.Rmd`)
  isTRUE(AMR_env$supports_colour) && !identical(getOption("rmarkdown.output.format"), "github_document")
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
is_dark <- function() {
  if (is.null(AMR_env$is_dark_theme)) {
    AMR_env$is_dark_theme <- !has_colour() || tryCatch(isTRUE(getExportedValue("getThemeInfo", ns = asNamespace("rstudioapi"))()$dark), error = function(e) FALSE)
  }
  isTRUE(AMR_env$is_dark_theme)
}
font_black <- function(..., collapse = " ", adapt = TRUE) {
  before <- "\033[38;5;232m"
  after <- "\033[39m"
  if (isTRUE(adapt) && is_dark()) {
    # white
    before <- "\033[37m"
    after <- "\033[39m"
  }
  try_colour(..., before = before, after = after, collapse = collapse)
}
font_white <- function(..., collapse = " ", adapt = TRUE) {
  before <- "\033[37m"
  after <- "\033[39m"
  if (isTRUE(adapt) && is_dark()) {
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
  if (is_dark()) {
    # similar to HTML #444444
    try_colour(..., before = "\033[48;5;238m", after = "\033[49m", collapse = collapse)
  } else {
    # similar to HTML #f0f0f0
    try_colour(..., before = "\033[48;5;255m", after = "\033[49m", collapse = collapse)
  }
}
font_red_bg <- function(..., collapse = " ") {
  # this is #ed553b (picked to be colourblind-safe with other SIR colours)
  try_colour(font_black(..., collapse = collapse, adapt = FALSE), before = "\033[48;5;203m", after = "\033[49m", collapse = collapse)
}
font_orange_bg <- function(..., collapse = " ") {
  # this is #f6d55c (picked to be colourblind-safe with other SIR colours)
  try_colour(font_black(..., collapse = collapse, adapt = FALSE), before = "\033[48;5;222m", after = "\033[49m", collapse = collapse)
}
font_yellow_bg <- function(..., collapse = " ") {
  try_colour(font_black(..., collapse = collapse, adapt = FALSE), before = "\033[48;5;228m", after = "\033[49m", collapse = collapse)
}
font_green_bg <- function(..., collapse = " ") {
  # this is #3caea3 (picked to be colourblind-safe with other SIR colours)
  try_colour(font_black(..., collapse = collapse, adapt = FALSE), before = "\033[48;5;79m", after = "\033[49m", collapse = collapse)
}
font_purple_bg <- function(..., collapse = " ") {
  try_colour(font_black(..., collapse = collapse, adapt = FALSE), before = "\033[48;5;89m", after = "\033[49m", collapse = collapse)
}
font_rose_bg <- function(..., collapse = " ") {
  if (is_dark()) {
    # this is #ed553b (picked to be colourblind-safe with other SIR colours)
    try_colour(font_black(..., collapse = collapse, adapt = FALSE), before = "\033[48;5;203m", after = "\033[49m", collapse = collapse)
  } else {
    # also colourblind-safe but softer
    try_colour(font_black(..., collapse = collapse, adapt = FALSE), before = "\033[48;5;217m", after = "\033[49m", collapse = collapse)
  }
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
font_url <- function(url, txt = url) {
  if (tryCatch(isTRUE(getExportedValue("ansi_has_hyperlink_support", ns = asNamespace("cli"))()), error = function(e) FALSE)) {
    paste0("\033]8;;", url, "\a", txt, "\033]8;;\a")
  } else {
    url
  }
}
font_stripstyle <- function(x) {
  # remove URLs
  x <- gsub("\033]8;;(.*?)\a.*?\033]8;;\a", "\\1", x)
  # from crayon:::ansi_regex
  x <- gsub("(?:(?:\\x{001b}\\[)|\\x{009b})(?:(?:[0-9]{1,3})?(?:(?:;[0-9]{0,3})*)?[A-M|f-m])|\\x{001b}[A-M]", "", x, perl = TRUE)
  x
}

progress_ticker <- function(n = 1, n_min = 0, print = TRUE, clear = TRUE, title = "", only_bar_percent = FALSE, ...) {
  if (print == FALSE || n < n_min) {
    # create fake/empty object
    pb <- list()
    pb$tick <- function() {
      invisible()
    }
    pb$kill <- function() {
      invisible()
    }
    set_clean_class(pb, new_class = "txtProgressBar")
  } else if (n >= n_min) {
    title <- trimws2(title)
    if (title != "") {
      title <- paste0(title, " ")
    }
    progress_bar <- import_fn("progress_bar", "progress", error_on_fail = FALSE)
    if (!is.null(progress_bar)) {
      # so we use progress::progress_bar
      # a close()-method was also added, see below for that
      pb <- progress_bar$new(
        show_after = 0,
        format = paste0(
          title,
          ifelse(only_bar_percent == TRUE, "[:bar] :percent", "[:bar] :percent (:current/:total,:eta)")
        ),
        clear = clear,
        total = n
      )
    } else {
      # use base R's txtProgressBar
      cat(title, "\n", sep = "")
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
  # for progress::progress_bar$new()
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

as_original_data_class <- function(df, old_class = NULL, extra_class = NULL) {
  if ("tbl_df" %in% old_class && pkg_is_available("tibble")) {
    # this will then also remove groups
    fn <- import_fn("as_tibble", "tibble")
  } else if ("data.table" %in% old_class && pkg_is_available("data.table")) {
    fn <- import_fn("as.data.table", "data.table")
  } else {
    fn <- function(x) base::as.data.frame(df, stringsAsFactors = FALSE)
  }
  out <- fn(df)
  # don't keep row names
  rownames(out) <- NULL
  # add additional class if needed
  if (!is.null(extra_class)) {
    class(out) <- c(extra_class, class(out))
  }
  out
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
    max(
      min(max_places,
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
  format_percentage(
    structure(
      .Data = as.double(x),
      class = c("percentage", "numeric")
    ),
    digits = digits, ...
  )
}

add_intrinsic_resistance_to_AMR_env <- function() {
  # for mo_is_intrinsic_resistant() - saves a lot of time when executed on this vector
  if (is.null(AMR_env$intrinsic_resistant)) {
    AMR_env$intrinsic_resistant <- paste(AMR::intrinsic_resistant$mo, AMR::intrinsic_resistant$ab)
  }
}

add_MO_lookup_to_AMR_env <- function() {
  # for all MO functions, saves a lot of time on package load and in package size
  if (is.null(AMR_env$MO_lookup)) {
    MO_lookup <- AMR::microorganisms

    MO_lookup$kingdom_index <- NA_real_
    MO_lookup[which(MO_lookup$kingdom == "Bacteria" | as.character(MO_lookup$mo) == "UNKNOWN"), "kingdom_index"] <- 1
    MO_lookup[which(MO_lookup$kingdom == "Fungi"), "kingdom_index"] <- 1.25
    MO_lookup[which(MO_lookup$kingdom == "Protozoa"), "kingdom_index"] <- 1.5
    MO_lookup[which(MO_lookup$kingdom == "Chromista"), "kingdom_index"] <- 1.75
    MO_lookup[which(MO_lookup$kingdom == "Archaea"), "kingdom_index"] <- 2
    # all the rest
    MO_lookup[which(is.na(MO_lookup$kingdom_index)), "kingdom_index"] <- 3

    # the fullname lowercase, important for the internal algorithms in as.mo()
    MO_lookup$fullname_lower <- tolower(trimws2(paste(
      MO_lookup$genus,
      MO_lookup$species,
      MO_lookup$subspecies
    )))
    ind <- MO_lookup$genus == "" | grepl("^[(]unknown ", MO_lookup$fullname, perl = TRUE)
    MO_lookup[ind, "fullname_lower"] <- tolower(MO_lookup[ind, "fullname", drop = TRUE])
    MO_lookup$fullname_lower <- trimws2(gsub("[^.a-z0-9/ \\-]+", "", MO_lookup$fullname_lower, perl = TRUE))
    # special for Salmonella - they have cities as subspecies but not the species (enterica) in the fullname:
    MO_lookup$fullname_lower[which(MO_lookup$subspecies %like_case% "^[A-Z]")] <- gsub(" enterica ", " ", MO_lookup$fullname_lower[which(MO_lookup$subspecies %like_case% "^[A-Z]")], fixed = TRUE)

    MO_lookup$genus_lower <- tolower(MO_lookup$genus)

    MO_lookup$full_first <- substr(MO_lookup$fullname_lower, 1, 1)
    MO_lookup$species_first <- tolower(substr(MO_lookup$species, 1, 1)) # tolower for groups (Streptococcus, Salmonella)
    MO_lookup$subspecies_first <- tolower(substr(MO_lookup$subspecies, 1, 1)) # tolower for Salmonella serovars
    AMR_env$MO_lookup <- MO_lookup
  }
}

trimws2 <- function(..., whitespace = "[\u0009\u000A\u000B\u000C\u000D\u0020\u0085\u00A0\u1680\u180E\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200A\u200B\u200C\u200D\u2028\u2029\u202F\u205F\u2060\u3000\uFEFF]") {
  # this is even faster than trimws() itself which sets "[ \t\r\n]".
  trimws(..., whitespace = whitespace)
}

totitle <- function(x) {
  gsub("^(.)", "\\U\\1", x, perl = TRUE)
}

readRDS_AMR <- function(file, refhook = NULL) {
  # this is readRDS with remote file support
  con <- file(file)
  on.exit(close(con))
  readRDS(con, refhook = refhook)
}

get_n_cores <- function(max_cores = Inf) {
  if (pkg_is_available("parallelly", min_version = "0.8.0", also_load = FALSE)) {
    available_cores <- import_fn("availableCores", "parallelly")
    n_cores <- min(available_cores(), na.rm = TRUE)
  } else {
    # `parallel` is part of base R since 2.14.0, but detectCores() is not very precise on exotic systems like Docker and quota-set Linux environments
    n_cores <- parallel::detectCores()[1]
    if (is.na(n_cores)) {
      n_cores <- 1
    }
  }
  max_cores <- floor(max_cores)
  if (max_cores == 0) {
    n_cores <- 1
  } else if (max_cores < 0) {
    n_cores <- max(1, n_cores - abs(max_cores))
  } else if (max_cores > 0) {
    n_cores <- min(n_cores, max_cores)
  }
  n_cores
}

# Faster data.table implementations ----

match <- function(x, table, ...) {
  if (!is.null(AMR_env$chmatch) && inherits(x, "character") && inherits(table, "character")) {
    # data.table::chmatch() is much faster than base::match() for character
    tryCatch(AMR_env$chmatch(x, table, ...), error = function(e) base::match(x, table, ...))
  } else {
    base::match(x, table, ...)
  }
}
`%in%` <- function(x, table) {
  if (!is.null(AMR_env$chin) && inherits(x, "character") && inherits(table, "character")) {
    # data.table::`%chin%`() is much faster than base::`%in%`() for character
    tryCatch(AMR_env$chin(x, table), error = function(e) base::`%in%`(x, table))
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
  # (required for extension of the 'mic' class)
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
}

if (getRversion() < "3.6.0") {
  str2lang <- function(s) {
    stopifnot(length(s) == 1L)
    ex <- parse(text = s, keep.source = FALSE)
    stopifnot(length(ex) == 1L)
    ex[[1L]]
  }
  # trims() was introduced in 3.3.0, but its argument `whitespace` only in 3.6.0
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

if (getRversion() < "4.0.0") {
  deparse1 <- function(expr, collapse = " ", width.cutoff = 500L, ...) {
    paste(deparse(expr, width.cutoff, ...), collapse = collapse)
  }
}

# nolint end
