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

# ------------------------------------------------
# THIS FILE WAS CREATED AUTOMATICALLY!
# Source file: data-raw/reproduction_of_poorman.R
# ------------------------------------------------

# Poorman: a package to replace all dplyr functions with base R so we can lose dependency on dplyr.
# These functions were downloaded from https://github.com/nathaneastwood/poorman,
# from this commit: https://github.com/nathaneastwood/poorman/tree/7d76d77f8f7bc663bf30fb5a161abb49801afa17
#
# All code below was released under MIT license, that permits 'free of charge, to any person obtaining a 
# copy of the software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
# is furnished to do so', given that a copyright notice is given in the software.
#
# Copyright notice as found on https://github.com/nathaneastwood/poorman/blob/master/LICENSE on 2 May 2020:
# YEAR: 2020
# COPYRIGHT HOLDER: Nathan Eastwood

arrange <- function(.data, ...) {
  check_is_dataframe(.data)
  if ("grouped_data" %in% class(.data)) {
    arrange.grouped_data(.data, ...)
  } else {
    arrange.default(.data, ...)
  }
}

arrange.default <- function(.data, ...) {
  rows <- eval.parent(substitute(with(.data, order(...))))
  .data[rows, , drop = FALSE]
}

arrange.grouped_data <- function(.data, ...) {
  apply_grouped_function(.data, "arrange", ...)
}
between <- function(x, left, right) {
  if (!is.null(attr(x, "class")) && !inherits(x, c("Date", "POSIXct"))) {
    warning("`between()` called on numeric vector with S3 class")
  }
  if (!is.double(x)) x <- as.numeric(x)
  x >= as.numeric(left) & x <= as.numeric(right)
}
count <- function(x, ..., wt = NULL, sort = FALSE, name = NULL) {
  groups <- get_groups(x)
  if (!missing(...)) x <- group_by(x, ..., .add = TRUE)
  wt <- deparse_var(wt)
  res <- do.call(tally, list(x, wt, sort, name))
  if (length(groups) > 0L) res <- do.call(group_by, list(res, as.name(groups)))
  res
}

tally <- function(x, wt = NULL, sort = FALSE, name = NULL) {
  name <- check_name(x, name)
  wt <- deparse_var(wt)
  res <- do.call(summarise, set_names(list(x, as.name(tally_n(x, wt))), c(".data", name)))
  res <- ungroup(res)
  if (isTRUE(sort)) res <- do.call(arrange, list(res, call("desc", as.name(name))))
  rownames(res) <- NULL
  res
}

add_count <- function(x, ..., wt = NULL, sort = FALSE, name = NULL) {
  name <- check_name(x, name)
  row_names <- rownames(x)
  wt <- deparse_var(wt)
  if (!missing(...)) x <- group_by(x, ..., .add = TRUE)
  res <- do.call(add_tally, list(x, wt, sort, name))
  res[row_names, ]
}

add_tally <- function(x, wt = NULL, sort = FALSE, name = NULL) {
  wt <- deparse_var(wt)
  n <- tally_n(x, wt)
  name <- check_name(x, name)
  res <- do.call(mutate, set_names(list(x, as.name(n)), c(".data", name)))

  if (isTRUE(sort)) {
    do.call(arrange, list(res, call("desc", as.name(name))))
  } else {
    res
  }
}

tally_n <- function(x, wt) {
  if (is.null(wt) && "n" %in% colnames(x)) {
    message("Using `n` as weighting variable")
    wt <- "n"
  }
  context$.data <- x
  on.exit(rm(list = ".data", envir = context))
  if (is.null(wt)) {
    "n()"
  } else {
    paste0("sum(", wt, ", na.rm = TRUE)")
  }
}

check_name <- function(df, name) {
  if (is.null(name)) {
    if ("n" %in% colnames(df)) {
      stop(
        "Column 'n' is already present in output\n",
        "* Use `name = \"new_name\"` to pick a new name"
      )
    }
    return("n")
  }

  if (!is.character(name) || length(name) != 1) {
    stop("`name` must be a single string")
  }

  name
}
desc <- function(x) -xtfrm(x)
select_env <- new.env()

peek_vars <- function() {
  get(".col_names", envir = select_env)
}

context <- new.env()

n <- function() {
  do.call(nrow, list(quote(.data)), envir = context)
}
filter <- function(.data, ...) {
  check_is_dataframe(.data)
  if ("grouped_data" %in% class(.data)) {
    filter.grouped_data(.data, ...)
  } else {
    filter.default(.data, ...)
  }
}

filter.default <- function(.data, ...) {
  conditions <- paste(deparse_dots(...), collapse = " & ")
  context$.data <- .data
  on.exit(rm(.data, envir = context))
  .data[do.call(with, list(.data, str2lang(unname(conditions)))), ]
}

filter.grouped_data <- function(.data, ...) {
  rows <- rownames(.data)
  res <- apply_grouped_function(.data, "filter", ...)
  res[rows[rows %in% rownames(res)], ]
}
group_by <- function(.data, ..., .add = FALSE) {
  check_is_dataframe(.data)
  pre_groups <- get_groups(.data)
  groups <- deparse_dots(...)
  if (isTRUE(.add)) groups <- unique(c(pre_groups, groups))
  unknown <- !(groups %in% colnames(.data))
  if (any(unknown)) stop("Invalid groups: ", groups[unknown])
  structure(.data, class = c("grouped_data", class(.data)), groups = groups)
}

ungroup <- function(x, ...) {
  check_is_dataframe(x)
  rm_groups <- deparse_dots(...)
  groups <- attr(x, "groups")
  if (length(rm_groups) == 0L) rm_groups <- groups
  attr(x, "groups") <- groups[!(groups %in% rm_groups)]
  if (length(attr(x, "groups")) == 0L) {
    attr(x, "groups") <- NULL
    class(x) <- class(x)[!(class(x) %in% "grouped_data")]
  }
  x
}

get_groups <- function(x) {
  attr(x, "groups", exact = TRUE)
}

has_groups <- function(x) {
  groups <- get_groups(x)
  if (is.null(groups)) FALSE else TRUE
}

set_groups <- function(x, groups) {
  attr(x, "groups") <- groups
  x
}

apply_grouped_function <- function(.data, fn, ...) {
  groups <- get_groups(.data)
  grouped <- split_into_groups(.data, groups)
  res <- do.call(rbind, unname(lapply(grouped, fn, ...)))
  if (any(groups %in% colnames(res))) {
    class(res) <- c("grouped_data", class(res))
    attr(res, "groups") <- groups[groups %in% colnames(res)]
  }
  res
}

split_into_groups <- function(.data, groups) {
  class(.data) <- "data.frame"
  group_factors <- lapply(groups, function(x, .data) as.factor(.data[, x]), .data)
  res <- split(x = .data, f = group_factors)
  res
}

print.grouped_data <- function(x, ..., digits = NULL, quote = FALSE, right = TRUE, row.names = TRUE, max = NULL) {
  class(x) <- "data.frame"
  print(x, ..., digits = digits, quote = quote, right = right, row.names = row.names, max = max)
  cat("\nGroups: ", paste(attr(x, "groups", exact = TRUE), collapse = ", "), "\n\n")
}
if_else <- function(condition, true, false, missing = NULL) {
  if (!is.logical(condition)) stop("`condition` must be a logical vector.")
  cls_true <- class(true)
  cls_false <- class(false)
  cls_missing <- class(missing)
  if (!identical(cls_true, cls_false)) {
    stop("The class of `true` <", class(true), "> is not the same as the class of `false` <", class(false), ">")
  }
  if (!is.null(missing) && !identical(cls_true, cls_missing)) {
    stop("`missing` must be a ", cls_true, " vector, not a ", cls_missing, " vector.")
  }
  res <- ifelse(condition, true, false)
  if (!is.null(missing)) res[is.na(res)] <- missing
  attributes(res) <- attributes(true)
  res
}

inner_join <- function(x, y, by = NULL, suffix = c(".x", ".y")) {
  join_worker(x = x, y = y, by = by, suffix = suffix, sort = FALSE)
}

left_join <- function(x, y, by = NULL, suffix = c(".x", ".y")) {
  join_worker(x = x, y = y, by = by, suffix = suffix, all.x = TRUE)
}

right_join <- function(x, y, by = NULL, suffix = c(".x", ".y")) {
  join_worker(x = x, y = y, by = by, suffix = suffix, all.y = TRUE)
}

full_join <- function(x, y, by = NULL, suffix = c(".x", ".y")) {
  join_worker(x = x, y = y, by = by, suffix = suffix, all = TRUE)
}

join_worker <- function(x, y, by = NULL, suffix = c(".x", ".y"), ...) {
  x[, ".join_id"] <- seq_len(nrow(x))
  if (is.null(by)) {
    by <- intersect(names(x), names(y))
    join_message(by)
    merged <- merge(x = x, y = y, by = by, suffixes = suffix, ...)[, union(names(x), names(y))]
  } else if (is.null(names(by))) {
    merged <- merge(x = x, y = y, by = by, suffixes = suffix, ...)
  } else {
    merged <- merge(x = x, y = y, by.x = names(by), by.y = by, suffixes = suffix, ...)
  }
  merged <- merged[order(merged[, ".join_id"]), colnames(merged) != ".join_id"]
  rownames(merged) <- NULL
  merged
}

join_message <- function(by) {
  if (length(by) > 1L) {
    message("Joining, by = c(\"", paste0(by, collapse = "\", \""), "\")\n", sep = "")
  } else {
    message("Joining, by = \"", by, "\"\n", sep = "")
  }
}

anti_join <- function(x, y, by = NULL) {
  filter_join_worker(x, y, by, type = "anti")
}

semi_join <- function(x, y, by = NULL) {
  filter_join_worker(x, y, by, type = "semi")
}

# filter_join_worker <- function(x, y, by = NULL, type = c("anti", "semi")) {
#   type <- match.arg(type, choices = c("anti", "semi"), several.ok = FALSE)
#   if (is.null(by)) {
#     by <- intersect(names(x), names(y))
#     join_message(by)
#   }
#   rows <- interaction(x[, by]) %in% interaction(y[, by])
#   if (type == "anti") rows <- !rows
#   res <- x[rows, ]
#   rownames(res) <- NULL
#   res
# }
lag <- function (x, n = 1L, default = NA) {
  if (inherits(x, "ts")) stop("`x` must be a vector, not a `ts` object, do you want `stats::lag()`?")
  if (length(n) != 1L || !is.numeric(n) || n < 0L) stop("`n` must be a nonnegative integer scalar")
  if (n == 0L) return(x)
  tryCatch(
    storage.mode(default) <- typeof(x),
    warning = function(w) {
      stop("Cannot convert `default` <", typeof(default), "> to `x` <", typeof(x), ">")
    }
  )
  xlen <- length(x)
  n <- pmin(n, xlen)
  res <- c(rep(default, n), x[seq_len(xlen - n)])
  attributes(res) <- attributes(x)
  res
}

lead <- function (x, n = 1L, default = NA) {
  if (length(n) != 1L || !is.numeric(n) || n < 0L) stop("n must be a nonnegative integer scalar")
  if (n == 0L) return(x)
  tryCatch(
    storage.mode(default) <- typeof(x),
    warning = function(w) {
      stop("Cannot convert `default` <", typeof(default), "> to `x` <", typeof(x), ">")
    }
  )
  xlen <- length(x)
  n <- pmin(n, xlen)
  res <- c(x[-seq_len(n)], rep(default, n))
  attributes(res) <- attributes(x)
  res
}
mutate <- function(.data, ...) {
  check_is_dataframe(.data)
  if ("grouped_data" %in% class(.data)) {
    mutate.grouped_data(.data, ...)
  } else {
    mutate.default(.data, ...)
  }
}

mutate.default <- function(.data, ...) {
  conditions <- deparse_dots(...)
  cond_names <- names(conditions)
  unnamed <- which(nchar(cond_names) == 0L)
  if (is.null(cond_names)) {
    names(conditions) <- conditions
  } else if (length(unnamed) > 0L) {
    names(conditions)[unnamed] <- conditions[unnamed]
  }
  not_matched <- names(conditions)[!names(conditions) %in% names(.data)]
  .data[, not_matched] <- NA
  context$.data <- .data
  on.exit(rm(.data, envir = context))
  for (i in seq_along(conditions)) {
    .data[, names(conditions)[i]] <- do.call(with, list(.data, str2lang(unname(conditions)[i])))
  }
  .data
}

mutate.grouped_data <- function(.data, ...) {
  rows <- rownames(.data)
  res <- apply_grouped_function(.data, "mutate", ...)
  res[rows, ]
}
n_distinct <- function(..., na.rm = FALSE) {
  res <- c(...)
  if (is.list(res)) return(nrow(unique(as.data.frame(res, stringsAsFactors = FALSE))))
  if (isTRUE(na.rm)) res <- res[!is.na(res)]
  length(unique(res))
}
`%>%` <- function(lhs, rhs) {
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  eval(as.call(c(rhs[[1L]], lhs, as.list(rhs[-1L]))), envir = parent.frame())
}
pull <- function(.data, var = -1) {
  var_deparse <- deparse_var(var)
  col_names <- colnames(.data)
  if (!(var_deparse %in% col_names) & grepl("^[[:digit:]]+L|[[:digit:]]", var_deparse)) {
    var <- as.integer(gsub("L", "", var_deparse))
    var <- if_else(var < 1L, rev(col_names)[abs(var)], col_names[var])
  } else if (var_deparse %in% col_names) {
    var <- var_deparse
  }
  .data[, var]
}
relocate <- function(.data, ..., .before = NULL, .after = NULL) {
  check_is_dataframe(.data)
  data_names <- colnames(.data)
  col_pos <- select_positions(.data, ...)

  .before <- deparse_var(.before)
  .after <- deparse_var(.after)
  has_before <- !is.null(.before)
  has_after <- !is.null(.after)

  if (has_before && has_after) {
    stop("You must supply only one of `.before` and `.after`")
  } else if (has_before) {
    where <- min(match(.before, data_names))
    col_pos <- c(setdiff(col_pos, where), where)
  } else if (has_after) {
    where <- max(match(.after, data_names))
    col_pos <- c(where, setdiff(col_pos, where))
  } else {
    where <- 1L
    col_pos <- union(col_pos, where)
  }
  lhs <- setdiff(seq(1L, where - 1L), col_pos)
  rhs <- setdiff(seq(where + 1L, ncol(.data)), col_pos)
  col_pos <- unique(c(lhs, col_pos, rhs))
  col_pos <- col_pos[col_pos <= length(data_names)]

  res <- .data[col_pos]
  if (has_groups(.data)) res <- set_groups(res, get_groups(.data))
  res
}
rename <- function(.data, ...) {
  check_is_dataframe(.data)
  new_names <- names(deparse_dots(...))
  if (length(new_names) == 0L) {
    warning("You didn't give any new names")
    return(.data)
  }
  col_pos <- select_positions(.data, ...)
  old_names <- colnames(.data)[col_pos]
  new_names_zero <- nchar(new_names) == 0L
  if (any(new_names_zero)) {
    warning("You didn't provide new names for: ", paste0("`", old_names[new_names_zero], collapse = ", "), "`")
    new_names[new_names_zero] <- old_names[new_names_zero]
  }
  colnames(.data)[col_pos] <- new_names
  .data
}
rownames_to_column <- function(.data, var = "rowname") {
  check_is_dataframe(.data)
  col_names <- colnames(.data)
  if (var %in% col_names) stop("Column `", var, "` already exists in `.data`")
  .data[, var] <- rownames(.data)
  rownames(.data) <- NULL
  .data[, c(var, setdiff(col_names, var))]
}

select <- function(.data, ...) {
  map <- names(deparse_dots(...))
  col_pos <- select_positions(.data, ..., group_pos = TRUE)
  res <- .data[, col_pos, drop = FALSE]
  to_map <- nchar(map) > 0L
  colnames(res)[to_map] <- map[to_map]
  if (has_groups(.data)) res <- set_groups(res, get_groups(.data))
  res
}
starts_with <- function(match, ignore.case = TRUE, vars = peek_vars()) {
  grep(pattern = paste0("^", paste0(match, collapse = "|^")), x = vars, ignore.case = ignore.case)
}

ends_with <- function(match, ignore.case = TRUE, vars = peek_vars()) {
  grep(pattern = paste0(paste0(match, collapse = "$|"), "$"), x = vars, ignore.case = ignore.case)
}

contains <- function(match, ignore.case = TRUE, vars = peek_vars()) {
  matches <- lapply(
    match,
    function(x) {
      if (isTRUE(ignore.case)) {
        match_u <- toupper(x)
        match_l <- tolower(x)
        pos_u <- grep(pattern = match_u, x = toupper(vars), fixed = TRUE)
        pos_l <- grep(pattern = match_l, x = tolower(vars), fixed = TRUE)
        unique(c(pos_l, pos_u))
      } else {
        grep(pattern = x, x = vars, fixed = TRUE)
      }
    }
  )
  unique(matches)
}

matches <- function(match, ignore.case = TRUE, perl = FALSE, vars = peek_vars()) {
  grep(pattern = match, x = vars, ignore.case = ignore.case, perl = perl)
}

num_range <- function(prefix, range, width = NULL, vars = peek_vars()) {
  if (!is.null(width)) {
    range <- sprintf(paste0("%0", width, "d"), range)
  }
  find <- paste0(prefix, range)
  if (any(duplicated(vars))) {
    stop("Column names must be unique")
  } else {
    x <- match(find, vars)
    x[!is.na(x)]
  }
}

all_of <- function(x, vars = peek_vars()) {
  x_ <- !x %in% vars
  if (any(x_)) {
    which_x_ <- which(x_)
    if (length(which_x_) == 1L) {
      stop("The column ", x[which_x_], " does not exist.")
    } else {
      stop("The columns ", paste(x[which_x_], collapse = ", "), " do not exist.")
    }
  } else {
    which(vars %in% x)
  }
}

any_of <- function(x, vars = peek_vars()) {
  which(vars %in% x)
}

everything <- function(vars = peek_vars()) {
  seq_along(vars)
}

last_col <- function(offset = 0L, vars = peek_vars()) {
  if (!is_wholenumber(offset)) stop("`offset` must be an integer")
  n <- length(vars)
  if (offset && n <= offset) {
    stop("`offset` must be smaller than the number of `vars`")
  } else if (n == 0) {
    stop("Can't select last column when `vars` is empty")
  } else {
    n - offset
  }
}
select_positions <- function(.data, ..., group_pos = FALSE) {
  cols <- eval(substitute(alist(...)))
  data_names <- colnames(.data)
  select_env$.col_names <- data_names
  on.exit(rm(list = ".col_names", envir = select_env))
  exec_env <- parent.frame(2L)
  pos <- unlist(lapply(cols, eval_expr, exec_env = exec_env))
  if (isTRUE(group_pos)) {
    groups <- get_groups(.data)
    missing_groups <- !(groups %in% cols)
    if (any(missing_groups)) {
      message("Adding missing grouping variables: `", paste(groups[missing_groups], collapse = "`, `"), "`")
      pos <- c(match(groups[missing_groups], data_names), pos)
    }
  }
  unique(pos)
}

eval_expr <- function(x, exec_env) {
  type <- typeof(x)
  switch(
    type,
    "integer" = x,
    "double" = as.integer(x),
    "character" = select_char(x),
    "symbol" = select_symbol(x, exec_env = exec_env),
    "language" = eval_call(x),
    stop("Expressions of type <", typeof(x), "> cannot be evaluated for use when subsetting.")
  )
}

select_char <- function(expr) {
  pos <- match(expr, select_env$.col_names)
  if (is.na(pos)) stop("Column `", expr, "` does not exist")
  pos
}

select_symbol <- function(expr, exec_env) {
  res <- try(select_char(as.character(expr)), silent = TRUE)
  if (inherits(res, "try-error")) {
    res <- tryCatch(
      select_char(eval(expr, envir = exec_env)),
      error = function(e) stop("Column ", expr, " does not exist.")
    )
  }
  res
}

eval_call <- function(x) {
  type <- as.character(x[[1]])
  switch(
    type,
    `:` = select_seq(x),
    `!` = select_negate(x),
    `-` = select_minus(x),
    `c` = select_c(x),
    `(` = select_bracket(x),
    select_context(x)
  )
}

select_seq <- function(expr) {
  x <- eval_expr(expr[[2]])
  y <- eval_expr(expr[[3]])
  x:y
}

select_negate <- function(expr) {
  x <- if (is_negated_colon(expr)) {
    expr <- call(":", expr[[2]][[2]], expr[[2]][[3]][[2]])
    eval_expr(expr)
  } else {
    eval_expr(expr[[2]])
  }
  x * -1L
}

is_negated_colon <- function(expr) {
  expr[[1]] == "!" && length(expr[[2]]) > 1L && expr[[2]][[1]] == ":" && expr[[2]][[3]][[1]] == "!"
}

select_minus <- function(expr) {
  x <- eval_expr(expr[[2]])
  x * -1L
}

select_c <- function(expr) {
  lst_expr <- as.list(expr)
  lst_expr[[1]] <- NULL
  unlist(lapply(lst_expr, eval_expr))
}

select_bracket <- function(expr) {
  eval_expr(expr[[2]])
}

select_context <- function(expr) {
  eval(expr, envir = context$.data)
}
slice <- function(.data, ...) {
  check_is_dataframe(.data)
  if ("grouped_data" %in% class(.data)) {
    slice.grouped_data(.data, ...)
  } else {
    slice.default(.data, ...)
  }
}

slice.default <- function(.data, ...) {
  rows <- c(...)
  stopifnot(is.numeric(rows) | is.integer(rows))
  if (all(rows > 0L)) rows <- rows[rows <= nrow(.data)]
  .data[rows, ]
}

slice.grouped_data <- function(.data, ...) {
  apply_grouped_function(.data, "slice", ...)
}
summarise <- function(.data, ...) {
  check_is_dataframe(.data)
  if ("grouped_data" %in% class(.data)) {
    summarise.grouped_data(.data, ...)
  } else {
    summarise.default(.data, ...)
  }
}

summarise.default <- function(.data, ...) {
  fns <- vapply(substitute(...()), deparse, NA_character_)
  context$.data <- .data
  on.exit(rm(.data, envir = context))
  if (has_groups(.data)) {
    group <- unique(.data[, get_groups(.data), drop = FALSE])
    if (nrow(group) == 0L) return(NULL)
  }
  res <- lapply(fns, function(x) do.call(with, list(.data, str2lang(x))))
  res <- as.data.frame(res)
  fn_names <- names(fns)
  colnames(res) <- if (is.null(fn_names)) fns else fn_names
  if (has_groups(.data)) res <- cbind(group, res)
  res
}

summarise.grouped_data <- function(.data, ...) {
  groups <- get_groups(.data)
  res <- apply_grouped_function(.data, "summarise", ...)
  res <- res[do.call(order, lapply(groups, function(x) res[, x])), ]
  rownames(res) <- NULL
  res
}

summarize <- summarise
summarize.default <- summarise.default
summarize.grouped_data <- summarise.grouped_data
transmute <- function(.data, ...) {
  check_is_dataframe(.data)
  if ("grouped_data" %in% class(.data)) {
    transmute.grouped_data(.data, ...)
  } else {
    transmute.default(.data, ...)
  }
}

transmute.default <- function(.data, ...) {
  conditions <- deparse_dots(...)
  mutated <- mutate(.data, ...)
  mutated[, names(conditions), drop = FALSE]
}

transmute.grouped_data <- function(.data, ...) {
  rows <- rownames(.data)
  res <- apply_grouped_function(.data, "transmute", ...)
  res[rows, ]
}
deparse_dots <- function(...) {
  vapply(substitute(...()), deparse, NA_character_)
}

deparse_var <- function(var) {
  sub_var <- eval(substitute(substitute(var)), parent.frame())
  if (is.symbol(sub_var)) var <- as.character(sub_var)
  var
}

check_is_dataframe <- function(.data) {
  parent_fn <- all.names(sys.call(-1L), max.names = 1L)
  if (!is.data.frame(.data)) stop(parent_fn, " must be given a data.frame")
  invisible()
}

is_wholenumber <- function(x) {
  x %% 1L == 0L
}

set_names <- function(object = nm, nm) {
  names(object) <- nm
  object
}

cume_dist <- function(x) {
  rank(x, ties.method = "max", na.last = "keep") / sum(!is.na(x))
}

dense_rank <- function(x) {
  match(x, sort(unique(x)))
}

min_rank <- function(x) {
  rank(x, ties.method = "min", na.last = "keep")
}

ntile <- function (x = row_number(), n) {
  if (!missing(x)) x <- row_number(x)
  len <- length(x) - sum(is.na(x))
  n <- as.integer(floor(n))
  if (len == 0L) {
    rep(NA_integer_, length(x))
  } else {
    n_larger <- as.integer(len %% n)
    n_smaller <- as.integer(n - n_larger)
    size <- len / n
    larger_size <- as.integer(ceiling(size))
    smaller_size <- as.integer(floor(size))
    larger_threshold <- larger_size * n_larger
    bins <- if_else(
      x <= larger_threshold,
      (x + (larger_size - 1L)) / larger_size,
      (x + (-larger_threshold + smaller_size - 1L)) / smaller_size + n_larger
    )
    as.integer(floor(bins))
  }
}

percent_rank <- function(x) {
  (min_rank(x) - 1) / (sum(!is.na(x)) - 1)
}

row_number <- function(x) {
  if (missing(x)) seq_len(n()) else rank(x, ties.method = "first", na.last = "keep")
}
