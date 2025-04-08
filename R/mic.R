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
# how to conduct AMR data analysis: https://amr-for-r.org/             #
# ==================================================================== #

# these are allowed MIC values and will become factor levels
VALID_MIC_LEVELS <- c(
  as.double(paste0("0.000", c(1:9))),
  as.double(paste0("0.00", c(1:99, 1953125, 390625, 78125))),
  as.double(paste0("0.0", c(1:99, 125, 128, 156, 165, 256, 512, 625, 3125, 15625))),
  as.double(paste0("0.", c(1:99, 125, 128, 256, 512))),
  1:9, 1.5,
  c(10:98)[9:98 %% 2 == TRUE],
  2^c(7:12), 192 * c(1:5), 80 * c(2:12)
)
VALID_MIC_LEVELS <- trimws(gsub("[.]?0+$", "", format(unique(sort(VALID_MIC_LEVELS)), scientific = FALSE), perl = TRUE))
operators <- c("<", "<=", "", ">=", ">")
VALID_MIC_LEVELS <- c(t(vapply(
  FUN.VALUE = character(length(VALID_MIC_LEVELS)),
  c("<", "<=", "", ">=", ">"),
  paste0,
  VALID_MIC_LEVELS
)))
COMMON_MIC_VALUES <- c(
  0.0001, 0.0002, 0.0005,
  0.001, 0.002, 0.004, 0.008,
  0.016, 0.032, 0.064,
  0.125, 0.25, 0.5,
  1, 2, 4, 8,
  16, 32, 64,
  128, 256, 512,
  1024, 2048, 4096
)

#' Transform Input to Minimum Inhibitory Concentrations (MIC)
#'
#' This transforms vectors to a new class [`mic`], which treats the input as decimal numbers, while maintaining operators (such as ">=") and only allowing valid MIC values known to the field of (medical) microbiology.
#' @rdname as.mic
#' @param x A [character] or [numeric] vector.
#' @param na.rm A [logical] indicating whether missing values should be removed.
#' @param keep_operators A [character] specifying how to handle operators (such as `>` and `<=`) in the input. Accepts one of three values: `"all"` (or `TRUE`) to keep all operators, `"none"` (or `FALSE`) to remove all operators, or `"edges"` to keep operators only at both ends of the range.
#' @param ... Arguments passed on to methods.
#' @details To interpret MIC values as SIR values, use [as.sir()] on MIC values. It supports guidelines from EUCAST (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`) and CLSI (`r min(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(clinical_breakpoints, guideline %like% "CLSI")$guideline)))`).
#'
#' This class for MIC values is a quite a special data type: formally it is an ordered [factor] with valid MIC values as [factor] levels (to make sure only valid MIC values are retained), but for any mathematical operation it acts as decimal numbers:
#'
#' ```
#' x <- random_mic(10)
#' x
#' #> Class 'mic'
#' #>  [1] 16     1      8      8      64     >=128  0.0625 32     32     16
#'
#' is.factor(x)
#' #> [1] TRUE
#'
#' x[1] * 2
#' #> [1] 32
#'
#' median(x)
#' #> [1] 26
#' ```
#'
#' This makes it possible to maintain operators that often come with MIC values, such ">=" and "<=", even when filtering using [numeric] values in data analysis, e.g.:
#'
#' ```
#' x[x > 4]
#' #> Class 'mic'
#' #> [1] 16    8     8     64    >=128 32    32    16
#'
#' df <- data.frame(x, hospital = "A")
#' subset(df, x > 4) # or with dplyr: df %>% filter(x > 4)
#' #>        x hospital
#' #> 1     16        A
#' #> 5     64        A
#' #> 6  >=128        A
#' #> 8     32        A
#' #> 9     32        A
#' #> 10    16        A
#' ```
#'
#' All so-called [group generic functions][groupGeneric()] are implemented for the MIC class (such as `!`, `!=`, `<`, `>=`, [exp()], [log2()]). Some mathematical functions are also implemented (such as [quantile()], [median()], [fivenum()]). Since [sd()] and [var()] are non-generic functions, these could not be extended. Use [mad()] as an alternative, or use e.g. `sd(as.numeric(x))` where `x` is your vector of MIC values.
#'
#' Using [as.double()] or [as.numeric()] on MIC values will remove the operators and return a numeric vector. Do **not** use [as.integer()] on MIC values as by the \R convention on [factor]s, it will return the index of the factor levels (which is often useless for regular users).
#'
#' The function [is.mic()] detects if the input contains class `mic`. If the input is a [data.frame] or [list], it iterates over all columns/items and returns a [logical] vector.
#'
#' Use [droplevels()] to drop unused levels. At default, it will return a plain factor. Use `droplevels(..., as.mic = TRUE)` to maintain the `mic` class.
#'
#' With [rescale_mic()], existing MIC ranges can be limited to a defined range of MIC values. This can be useful to better compare MIC distributions.
#'
#' For `ggplot2`, use one of the [`scale_*_mic()`][scale_x_mic()] functions to plot MIC values. They allows custom MIC ranges and to plot intermediate log2 levels for missing MIC values.
#' @return Ordered [factor] with additional class [`mic`], that in mathematical operations acts as a [numeric] vector. Bear in mind that the outcome of any mathematical operation on MICs will return a [numeric] value.
#' @aliases mic
#' @export
#' @seealso [as.sir()]
#' @examples
#' mic_data <- as.mic(c(">=32", "1.0", "1", "1.00", 8, "<=0.128", "8", "16", "16"))
#' mic_data
#' is.mic(mic_data)
#'
#' # this can also coerce combined MIC/SIR values:
#' as.mic("<=0.002; S")
#'
#' # mathematical processing treats MICs as numeric values
#' fivenum(mic_data)
#' quantile(mic_data)
#' all(mic_data < 512)
#'
#' # rescale MICs using rescale_mic()
#' rescale_mic(mic_data, mic_range = c(4, 16))
#'
#' # interpret MIC values
#' as.sir(
#'   x = as.mic(2),
#'   mo = as.mo("Streptococcus pneumoniae"),
#'   ab = "AMX",
#'   guideline = "EUCAST"
#' )
#' as.sir(
#'   x = as.mic(c(0.01, 2, 4, 8)),
#'   mo = as.mo("Streptococcus pneumoniae"),
#'   ab = "AMX",
#'   guideline = "EUCAST"
#' )
#'
#' # plot MIC values, see ?plot
#' plot(mic_data)
#' plot(mic_data, mo = "E. coli", ab = "cipro")
#'
#' if (require("ggplot2")) {
#'   autoplot(mic_data, mo = "E. coli", ab = "cipro")
#' }
#' if (require("ggplot2")) {
#'   autoplot(mic_data, mo = "E. coli", ab = "cipro", language = "nl") # Dutch
#' }
as.mic <- function(x, na.rm = FALSE, keep_operators = "all") {
  meet_criteria(x, allow_NA = TRUE)
  meet_criteria(na.rm, allow_class = "logical", has_length = 1)
  meet_criteria(keep_operators, allow_class = c("character", "logical"), is_in = c("all", "none", "edges", FALSE, TRUE), has_length = 1)
  if (isTRUE(keep_operators)) {
    keep_operators <- "all"
  } else if (isFALSE(keep_operators)) {
    keep_operators <- "none"
  }

  if (is.mic(x) && (keep_operators == "all" || !any(x %like% "[>=<]", na.rm = TRUE))) {
    if (!identical(levels(x), VALID_MIC_LEVELS)) {
      # might be from an older AMR version - just update MIC factor levels
      x <- set_clean_class(factor(as.character(x), levels = VALID_MIC_LEVELS, ordered = TRUE),
        new_class = c("mic", "ordered", "factor")
      )
    }
    return(x)
  }

  x.bak <- NULL
  if (is.numeric(x)) {
    x.bak <- format(x, scientific = FALSE)
    # MICs never have more than 9 decimals, so:
    x <- format(round(x, 9), scientific = FALSE)
  } else {
    x <- as.character(unlist(x))
  }
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  }
  x <- trimws2(x)
  x[x == ""] <- NA
  if (is.null(x.bak)) {
    x.bak <- x
  }
  # remove NAs on beforehand to not count them
  x.bak <- gsub("(NA)+", "", x.bak)
  # and trim
  x.bak <- trimws2(x.bak)

  # comma to period
  x <- gsub(",", ".", x, fixed = TRUE)
  # transform Unicode for >= and <=
  x <- gsub("\u2264", "<=", x, fixed = TRUE)
  x <- gsub("\u2265", ">=", x, fixed = TRUE)
  if (any(x %like% "[0-9]/.*[0-9]", na.rm = TRUE)) {
    warning_("Some MICs were combined values, only the first values are kept")
    x[x %like% "[0-9]/.*[0-9]"] <- gsub("/.*", "", x[x %like% "[0-9]/.*[0-9]"])
  }
  # remove other invalid characters
  x <- gsub("[^a-zA-Z0-9.><= -]+", "", x, perl = TRUE)
  # transform => to >= and =< to <=
  x <- gsub("=<", "<=", x, fixed = TRUE)
  x <- gsub("=>", ">=", x, fixed = TRUE)
  # Remove leading == and =
  x <- gsub("^=+", "", x)
  # retrieve signs and remove them from input
  x_signs <- trimws(gsub("[^>=<]", "", x))
  x <- trimws(gsub("[>=<]", "", x))
  # dots without a leading zero must start with 0
  x <- gsub("([^0-9]|^)[.]", "\\10.", x, perl = TRUE)
  # values like "<=0.2560.512" should be 0.512
  x <- gsub(".*[.].*[.]", "0.", x, perl = TRUE)
  # remove ending .0
  x <- gsub("[.]+0$", "", x, perl = TRUE)
  # remove all after last digit
  x <- gsub("[^0-9]+$", "", x, perl = TRUE)
  # keep only one zero before dot
  x <- gsub("^0+[.]", "0.", x, perl = TRUE)
  # starting 00 is probably 0.0 if there's no dot yet
  x[x %unlike% "[.]"] <- gsub("^00", "0.0", x[!x %like% "[.]"])
  # remove last zeroes
  x <- gsub("([.].?)0+$", "\\1", x, perl = TRUE)
  x <- gsub("(.*[.])0+$", "\\10", x, perl = TRUE)
  # remove ending .0 again
  x[x %like% "[.]"] <- gsub("0+$", "", x[x %like% "[.]"])
  # never end with dot
  x <- gsub("[.]$", "", x, perl = TRUE)
  # remove scientific notation
  x[x %like% "[0-9]e[-]?[0-9]"] <- trimws(format(suppressWarnings(as.double(x[x %like% "[0-9]e[-]?[0-9]"])), scientific = FALSE))
  # add signs again
  x <- paste0(x_signs, x)
  # remove NAs introduced by format()
  x <- gsub("(NA)+", "", x)
  # trim it
  x <- trimws2(x)

  ## previously unempty values now empty - should return a warning later on
  x[x.bak != "" & x == ""] <- "invalid"

  na_before <- x[is.na(x) | x == ""] %pm>% length()
  x[!as.character(x) %in% VALID_MIC_LEVELS] <- NA
  na_after <- x[is.na(x) | x == ""] %pm>% length()

  if (na_before != na_after) {
    list_missing <- x.bak[is.na(x) & !is.na(x.bak) & x.bak != ""] %pm>%
      unique() %pm>%
      sort() %pm>%
      vector_and(quotes = TRUE)
    cur_col <- get_current_column()
    warning_("in `as.mic()`: ", na_after - na_before, " result",
      ifelse(na_after - na_before > 1, "s", ""),
      ifelse(is.null(cur_col), "", paste0(" in index '", cur_col, "'")),
      " truncated (",
      round(((na_after - na_before) / length(x)) * 100),
      "%) that were invalid MICs: ",
      list_missing,
      call = FALSE
    )
  }

  if (keep_operators == "none" && !all(is.na(x))) {
    x <- gsub("[>=<]", "", x)
  } else if (keep_operators == "edges" && !all(is.na(x))) {
    dbls <- as.double(gsub("[>=<]", "", x))
    x[dbls == min(dbls, na.rm = TRUE)] <- paste0("<=", min(dbls, na.rm = TRUE))
    x[dbls == max(dbls, na.rm = TRUE)] <- paste0(">=", max(dbls, na.rm = TRUE))
    keep <- x[dbls == max(dbls, na.rm = TRUE) | dbls == min(dbls, na.rm = TRUE)]
    x[!x %in% keep] <- gsub("[>=<]", "", x[!x %in% keep])
  }

  set_clean_class(factor(x, levels = VALID_MIC_LEVELS, ordered = TRUE),
    new_class = c("mic", "ordered", "factor")
  )
}

#' @rdname as.mic
#' @export
is.mic <- function(x) {
  if (identical(typeof(x), "list")) {
    unname(vapply(FUN.VALUE = logical(1), x, is.mic))
  } else {
    isTRUE(inherits(x, "mic"))
  }
}

#' @rdname as.mic
#' @details `NA_mic_` is a missing value of the new `mic` class, analogous to e.g. base \R's [`NA_character_`][base::NA].
#' @format NULL
#' @export
NA_mic_ <- set_clean_class(factor(NA, levels = VALID_MIC_LEVELS, ordered = TRUE),
  new_class = c("mic", "ordered", "factor")
)

#' @rdname as.mic
#' @param mic_range A manual range to rescale the MIC values, e.g., `mic_range = c(0.001, 32)`. Use `NA` to prevent rescaling on one side, e.g., `mic_range = c(NA, 32)`.
#' @export
rescale_mic <- function(x, mic_range, keep_operators = "edges", as.mic = TRUE) {
  meet_criteria(mic_range, allow_class = c("numeric", "integer", "logical", "mic"), has_length = 2, allow_NA = TRUE, allow_NULL = TRUE)
  if (is.numeric(mic_range)) {
    mic_range <- trimws(format(mic_range, scientific = FALSE))
    mic_range <- gsub("[.]0+$", "", mic_range)
    mic_range[mic_range == "NA"] <- NA_character_
  } else if (is.mic(mic_range)) {
    mic_range <- as.character(mic_range)
  }
  stop_ifnot(
    all(mic_range %in% c(VALID_MIC_LEVELS, NA)),
    "Values in `mic_range` must be valid MIC values. ",
    "The allowed range is ", format(as.double(as.mic(VALID_MIC_LEVELS)[1]), scientific = FALSE), " to ", format(as.double(as.mic(VALID_MIC_LEVELS)[length(VALID_MIC_LEVELS)]), scientific = FALSE), ". ",
    "Unvalid: ", vector_and(mic_range[!mic_range %in% c(VALID_MIC_LEVELS, NA)], quotes = FALSE), "."
  )

  x <- as.mic(x)
  if (is.null(mic_range)) {
    mic_range <- c(NA, NA)
  }
  mic_range <- as.mic(mic_range)

  min_mic <- mic_range[1]
  max_mic <- mic_range[2]
  if (!is.na(min_mic)) {
    x[x < min_mic] <- min_mic
  }
  if (!is.na(max_mic)) {
    x[x > max_mic] <- max_mic
  }

  x <- as.mic(x, keep_operators = ifelse(keep_operators == "edges", "none", keep_operators))

  if (isTRUE(as.mic)) {
    if (keep_operators == "edges" && length(unique(x)) > 1) {
      x[x == min(x, na.rm = TRUE)] <- paste0("<=", x[x == min(x, na.rm = TRUE)])
      x[x == max(x, na.rm = TRUE)] <- paste0(">=", x[x == max(x, na.rm = TRUE)])
    }
    return(x)
  }

  # create a manual factor with levels only within desired range
  expanded <- plotrange_as_table(x,
    expand = TRUE,
    keep_operators = ifelse(keep_operators == "edges", "none", keep_operators),
    mic_range = mic_range
  )
  if (keep_operators == "edges") {
    names(expanded)[1] <- paste0("<=", names(expanded)[1])
    names(expanded)[length(expanded)] <- paste0(">=", names(expanded)[length(expanded)])
  }
  # MICs contain all MIC levels, so strip this to only existing levels and their intermediate values
  out <- factor(names(expanded),
    levels = names(expanded),
    ordered = TRUE
  )
  # and only keep the ones in the data
  if (keep_operators == "edges") {
    out <- out[match(x, as.double(as.mic(out, keep_operators = "all")))]
  } else {
    out <- out[match(x, out)]
  }
  out
}

#' @rdname as.mic
#' @details Use [mic_p50()] and [mic_p90()] to get the 50th and 90th percentile of MIC values. They return 'normal' [numeric] values.
#' @export
mic_p50 <- function(x, na.rm = FALSE, ...) {
  x <- as.mic(x)
  as.double(stats::quantile(x, probs = 0.5, na.rm = na.rm))
}

#' @rdname as.mic
#' @export
mic_p90 <- function(x, na.rm = FALSE, ...) {
  x <- as.mic(x)
  as.double(stats::quantile(x, probs = 0.9, na.rm = na.rm))
}

#' @method as.double mic
#' @export
#' @noRd
as.double.mic <- function(x, ...) {
  as.double(gsub("[<=>]+", "", as.character(x), perl = TRUE))
}

#' @method as.numeric mic
#' @export
#' @noRd
as.numeric.mic <- function(x, ...) {
  as.numeric(gsub("[<=>]+", "", as.character(x), perl = TRUE))
}

#' @rdname as.mic
#' @method droplevels mic
#' @param as.mic A [logical] to indicate whether the `mic` class should be kept - the default is `TRUE` for [rescale_mic()] and `FALSE` for [droplevels()]. When setting this to `FALSE` in [rescale_mic()], the output will have factor levels that acknowledge `mic_range`.
#' @export
droplevels.mic <- function(x, as.mic = FALSE, ...) {
  x <- as.mic(x) # make sure that currently implemented MIC levels are used
  x <- droplevels.factor(x, ...)
  if (as.mic == TRUE) {
    class(x) <- c("mic", "ordered", "factor")
  }
  x
}

all_valid_mics <- function(x) {
  if (!inherits(x, c("mic", "character", "factor", "numeric", "integer"))) {
    return(FALSE)
  }
  x_mic <- tryCatch(suppressWarnings(as.mic(x[!is.na(x)])),
    error = function(e) NA
  )
  !any(is.na(x_mic)) && !all(is.na(x))
}

# will be exported using s3_register() in R/zzz.R
pillar_shaft.mic <- function(x, ...) {
  if (!identical(levels(x), VALID_MIC_LEVELS) && message_not_thrown_before("pillar_shaft.mic")) {
    warning_(AMR_env$sup_1_icon, " These columns contain an outdated or altered structure - convert with `as.mic()` to update",
      call = FALSE
    )
  }
  crude_numbers <- as.double(x)
  operators <- gsub("[^<=>]+", "", as.character(x))
  operators[!is.na(operators) & operators != ""] <- font_silver(operators[!is.na(operators) & operators != ""], collapse = NULL)
  out <- trimws(paste0(operators, trimws(format(crude_numbers))))
  out[is.na(x)] <- font_na(NA)
  # make trailing zeroes less visible
  out[out %like% "[.]"] <- gsub("([.]?0+)$", font_silver("\\1"), out[out %like% "[.]"], perl = TRUE)
  create_pillar_column(out, align = "right", width = max(nchar(font_stripstyle(out))))
}

# will be exported using s3_register() in R/zzz.R
type_sum.mic <- function(x, ...) {
  if (!identical(levels(x), VALID_MIC_LEVELS)) {
    paste0("mic", AMR_env$sup_1_icon)
  } else {
    "mic"
  }
}

#' @method print mic
#' @export
#' @noRd
print.mic <- function(x, ...) {
  cat("Class 'mic'")
  if (!identical(levels(x), VALID_MIC_LEVELS)) {
    cat(font_red(" with an outdated or altered structure - convert with `as.mic()` to update"))
  }
  cat("\n")
  print(as.character(x), quote = FALSE)
  att <- attributes(x)
  if ("na.action" %in% names(att)) {
    cat(font_silver(paste0("(NA ", class(att$na.action), ": ", paste0(att$na.action, collapse = ", "), ")\n")))
  }
}

#' @method summary mic
#' @export
#' @noRd
summary.mic <- function(object, ...) {
  summary(as.double(object), ...)
}

#' @method as.matrix mic
#' @export
#' @noRd
as.matrix.mic <- function(x, ...) {
  as.matrix(as.double(x), ...)
}
#' @method as.vector mic
#' @export
#' @noRd
as.vector.mic <- function(x, mode = "numneric", ...) {
  y <- NextMethod()
  y <- as.mic(y)
  calls <- unlist(lapply(sys.calls(), as.character))
  if (any(calls %in% c("rbind", "cbind")) && message_not_thrown_before("as.vector.mic")) {
    warning_("Functions `rbind()` and `cbind()` cannot preserve the structure of MIC values. Use dplyr's `bind_rows()` or `bind_cols()` instead.", call = FALSE)
  }
  y
}
#' @method as.list mic
#' @export
#' @noRd
as.list.mic <- function(x, ...) {
  lapply(as.list(as.character(x), ...), as.mic)
}
#' @method as.data.frame mic
#' @export
#' @noRd
as.data.frame.mic <- function(x, ...) {
  as.data.frame.vector(as.mic(x), ...)
}

#' @method [ mic
#' @export
#' @noRd
"[.mic" <- function(x, ...) {
  y <- NextMethod()
  as.mic(y)
}
#' @method [[ mic
#' @export
#' @noRd
"[[.mic" <- function(x, ...) {
  y <- NextMethod()
  as.mic(y)
}
#' @method [<- mic
#' @export
#' @noRd
"[<-.mic" <- function(i, j, ..., value) {
  value <- as.mic(value)
  y <- NextMethod()
  as.mic(y)
}
#' @method [[<- mic
#' @export
#' @noRd
"[[<-.mic" <- function(i, j, ..., value) {
  value <- as.mic(value)
  y <- NextMethod()
  as.mic(y)
}
#' @method c mic
#' @export
#' @noRd
c.mic <- function(...) {
  as.mic(unlist(lapply(list(...), as.character)))
}

#' @method unique mic
#' @export
#' @noRd
unique.mic <- function(x, incomparables = FALSE, ...) {
  y <- NextMethod()
  as.mic(y)
}

#' @method rep mic
#' @export
#' @noRd
rep.mic <- function(x, ...) {
  y <- NextMethod()
  as.mic(y)
}

#' @method sort mic
#' @export
#' @noRd
sort.mic <- function(x, decreasing = FALSE, ...) {
  x <- as.mic(x) # make sure that currently implemented MIC levels are used
  dbl <- as.double(x)
  # make sure that e.g. '<0.001' comes before '0.001', and '>0.001' comes after
  dbl[as.character(x) %like% "<[0-9]"] <- dbl[as.character(x) %like% "<[0-9]"] - 0.000002
  dbl[as.character(x) %like% "<="] <- dbl[as.character(x) %like% "<="] - 0.000001
  dbl[as.character(x) %like% ">="] <- dbl[as.character(x) %like% ">="] + 0.000001
  dbl[as.character(x) %like% ">[0-9]"] <- dbl[as.character(x) %like% ">[0-9]"] + 0.000002
  if (decreasing == TRUE) {
    x[order(-dbl)]
  } else {
    x[order(dbl)]
  }
}

#' @method hist mic
#' @importFrom graphics hist
#' @export
#' @noRd
hist.mic <- function(x, ...) {
  warning_("in `hist()`: use `plot()` or ggplot2's `autoplot()` for optimal plotting of MIC values")
  hist(log2(x))
}

# will be exported using s3_register() in R/zzz.R
get_skimmers.mic <- function(column) {
  column <- as.mic(column) # make sure that currently implemented MIC levels are used
  skimr::sfl(
    skim_type = "mic",
    p0 = ~ stats::quantile(., probs = 0, na.rm = TRUE, names = FALSE),
    p25 = ~ stats::quantile(., probs = 0.25, na.rm = TRUE, names = FALSE),
    p50 = ~ stats::quantile(., probs = 0.5, na.rm = TRUE, names = FALSE),
    p75 = ~ stats::quantile(., probs = 0.75, na.rm = TRUE, names = FALSE),
    p100 = ~ stats::quantile(., probs = 1, na.rm = TRUE, names = FALSE),
    hist = ~ skimr::inline_hist(log2(stats::na.omit(.)), 5)
  )
}

# Miscellaneous mathematical functions ------------------------------------

#' @method mean mic
#' @export
#' @noRd
mean.mic <- function(x, trim = 0, na.rm = FALSE, ...) {
  mean(as.double(x), trim = trim, na.rm = na.rm, ...)
}

#' @method median mic
#' @importFrom stats median
#' @export
#' @noRd
median.mic <- function(x, na.rm = FALSE, ...) {
  median(as.double(x), na.rm = na.rm, ...)
}

#' @method quantile mic
#' @importFrom stats quantile
#' @export
#' @noRd
quantile.mic <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE,
                         names = TRUE, type = 7, ...) {
  quantile(as.double(x), probs = probs, na.rm = na.rm, names = names, type = type, ...)
}

# Math (see ?groupGeneric) ------------------------------------------------

#' @export
Math.mic <- function(x, ...) {
  x <- as.double(x)
  # set class to numeric, because otherwise NextMethod will be factor (since mic is a factor)
  .Class <- class(x)
  NextMethod(.Generic)
}

# Ops (see ?groupGeneric) -------------------------------------------------

#' @export
Ops.mic <- function(e1, e2) {
  e1_chr <- as.character(e1)
  e2_chr <- character(0)
  e1 <- as.double(e1)
  if (!missing(e2)) {
    # when .Generic is `!`, e2 is missing
    e2_chr <- as.character(e2)
    e2 <- as.double(e2)
  }
  if (as.character(.Generic) %in% c("<", "<=", "==", "!=", ">", ">=")) {
    # make sure that <0.002 is lower than 0.002
    # and that >32 is higher than 32, but equal to >=32
    e1[e1_chr %like% "<" & e1_chr %unlike% "="] <- e1[e1_chr %like% "<" & e1_chr %unlike% "="] - 0.000001
    e1[e1_chr %like% ">" & e1_chr %unlike% "="] <- e1[e1_chr %like% ">" & e1_chr %unlike% "="] + 0.000001
    e2[e2_chr %like% "<" & e2_chr %unlike% "="] <- e2[e2_chr %like% "<" & e2_chr %unlike% "="] - 0.000001
    e2[e2_chr %like% ">" & e2_chr %unlike% "="] <- e2[e2_chr %like% ">" & e2_chr %unlike% "="] + 0.000001
  }
  # set .Class to numeric, because otherwise NextMethod will be factor (since mic is a factor)
  .Class <- class(e1)
  NextMethod(.Generic)
}

# Complex (see ?groupGeneric) ---------------------------------------------

#' @export
Complex.mic <- function(z) {
  z <- as.double(z)
  # set class to numeric, because otherwise NextMethod will be factor (since mic is a factor)
  .Class <- class(z)
  NextMethod(.Generic)
}

# Summary (see ?groupGeneric) ---------------------------------------------

#' @export
Summary.mic <- function(..., na.rm = FALSE) {
  # NextMethod() cannot be called from an anonymous function (`...`), so we get() the generic directly:
  fn <- get(.Generic, envir = .GenericCallEnv)
  fn(as.double(c(...)),
    na.rm = na.rm
  )
}
