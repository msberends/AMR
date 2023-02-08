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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

# ------------------------------------------------
# THIS FILE WAS CREATED AUTOMATICALLY!
# Source file: data-raw/reproduction_of_poorman.R
# ------------------------------------------------

# {poorman}: a package to replace all dplyr functions with base R so we can lose dependency on {dplyr}.
# These functions were downloaded from https://github.com/nathaneastwood/poorman,
# from this commit: https://github.com/nathaneastwood/poorman/tree/3cc0a9920b1eb559dd166f548561244189586b3a.
#
# All functions are prefixed with 'pm_' to make it obvious that they are {dplyr} substitutes.
#
# All code below was released under MIT license, that permits 'free of charge, to any person obtaining a
# copy of the software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
# is furnished to do so', given that a copyright notice is given in the software.
#
# Copyright notice on 8 February 2023, the day this code was downloaded, as found on
# https://github.com/nathaneastwood/poorman/blob/3cc0a9920b1eb559dd166f548561244189586b3a/LICENSE:
# YEAR: 2020
# COPYRIGHT HOLDER: Nathan Eastwood

pm_across <- function(.cols = everything(), .fns = NULL, ..., .names = NULL) {
  setup <- pm_setup_across(substitute(.cols), .fns, .names)
  if (length(setup$names) == 1 && grepl("\\{\\.col\\}|\\{\\.fn\\}", setup$names)) {
    ref <- setup$names
    id <- 1
    fn_names <- unique(names(setup$funs))
    for (i in seq_along(setup$cols)) {
      .col <- setup$cols[i]
      for (j in seq_along(fn_names)) {
        .fn <- fn_names[j]
        setup$names[id] <- gluestick(ref)
        id <- id + 1
      }
    }
  }
  cols <- setup$cols
  n_cols <- length(cols)
  if (n_cols == 0L) return(data.frame())
  funs <- setup$funs
  data <- pm_context$get_columns(cols)
  names <- setup$names
  if (is.null(funs)) {
    data <- data.frame(data)
    if (is.null(names)) {
      return(data)
    } else {
      return(setNames(data, names))
    }
  }
  n_fns <- length(funs)
  res <- vector(mode = "list", length = n_fns * n_cols)
  k <- 1L
  for (i in seq_len(n_cols)) {
    pm_context$cur_column <- cols[[i]]
    col <- data[[i]]
    for (j in seq_len(n_fns)) {
      res[[k]] <- funs[[j]](col, ...)
      k <- k + 1L
    }
  }
  if (is.null(names(res))) names(res) <- names
  as.data.frame(res)
}
pm_if_any <- function(.cols, .fns = NULL, ..., .names = NULL) {
  df <- do.call(across, list(.cols = substitute(.cols), .fns = .fns, ..., .names = .names))
  if (nrow(df) == 0L) return(FALSE)
  check_if_types(df)
  Reduce(`|`, df)
}
pm_if_all <- function(.cols, .fns = NULL, ..., .names = NULL) {
  df <- do.call(across, list(.cols = substitute(.cols), .fns = .fns, ..., .names = .names))
  if (nrow(df) == 0L) return(FALSE)
  check_if_types(df)
  Reduce(`&`, df)
}
pm_check_if_types <- function(df) {
  types <- vapply(df, class, NA_character_)
  not_logical <- types != "logical"
  if (any(not_logical)) {
    stop(
      "Cannot convert the following columns to <logical>:\n    ",
      paste0(colnames(df)[not_logical], " <", types, "> ", collapse = "\n    ")
    )
  }
}
pm_setup_across <- function(.cols, .fns, .names) {
  cols <- pm_eval_select_pos(.data = pm_context$.data, .cols, .pm_group_pos = FALSE)
  cols <- pm_context$get_colnames()[cols]
  if (pm_context$is_grouped()) cols <- setdiff(cols, pm_group_vars(pm_context$.data))
  funs <- if (is.null(.fns)) NULL else if (!is.list(.fns)) list(.fns) else .fns
  if (is.null(funs)) return(list(cols = cols, funs = funs, names = .names))
  f_nms <- names(funs)
  if (is.null(f_nms) && !is.null(.fns)) names(funs) <- seq_along(funs)
  if (any(nchar(f_nms) == 0L)) {
    miss <- which(nchar(f_nms) == 0L)
    names(funs)[miss] <- miss
    f_nms <- names(funs)
  }
  funs <- lapply(funs, as_function)
  names <- if (!is.null(.names)) {
    .names
  } else {
    if (length(funs) == 1L && is.null(f_nms)) {
      cols
    } else {
      nms <- do.call(paste, c(rev(expand.grid(names(funs), cols)), sep = "_"))
      if (length(nms) == 0L) nms <- NULL
      nms
    }
  }
  list(cols = cols, funs = funs, names = names)
}
pm_arrange <- function(.data, ...) {
  pm_arrange.data.frame(.data, ...)
}
pm_arrange.data.frame <- function(.data, ..., .by_group = FALSE) {
  dots <- pm_dotdotdot(...)
  is_grouped <- pm_has_groups(.data)
  if (isTRUE(.by_group)) dots <- c(groups(.data), dots)
  rows <- pm_arrange_rows(.data = .data, dots)
  row_number <- attr(.data, "row.names") 
  out <- .data[rows, , drop = FALSE]
  if (is.numeric(row_number)) {
    row.names(out) <- row_number
  }
  if (is_grouped) {
    attr(out, "groups") <- pm_calculate_groups(out, pm_group_vars(out))
  }
  out
}
pm_arrange_rows <- function(.data, dots) {
  if (length(dots) == 0L) return(seq_len(nrow(.data)))
  for (i in seq_along(dots)) {
    tmp <- deparse(dots[[i]])
    if (startsWith(tmp, "desc(")) {
      tmp <- gsub("^desc\\(", "-", tmp)
      tmp <- gsub("\\)$", "", tmp)
    }
    dots[[i]] <- parse(text = tmp, keep.source = FALSE)[[1]]
  }
  used <- unname(do.call(c, lapply(dots, pm_find_used)))
  used <- used[used %in% colnames(.data)]
  for (i in seq_along(dots)) {
    if (is.character(.data[[used[[i]]]])) {
      .data[[used[[i]]]] <- factor(.data[[used[[i]]]])
    }
    if (is.factor(.data[[used[[i]]]]) &&
        (startsWith(deparse(dots[[i]]), "desc(") ||
         startsWith(deparse(dots[[i]]), "-"))) {
      dots[[i]] <- bquote(-xtfrm(.(as.name(used[[i]]))))
    }
  }
  data <- do.call(pm_transmute, c(list(.data = pm_ungroup(.data)), dots))
  do.call(order, c(data, list(decreasing = FALSE, na.last = TRUE)))
}
pm_bind_cols <- function(...) {
  lsts <- list(...)
  lsts <- squash(lsts)
  lsts <- Filter(Negate(is.null), lsts)
  if (length(lsts) == 0L) return(data.frame())
  lapply(lsts, function(x) is_df_or_vector(x))
  lsts <- do.call(cbind, lsts)
  if (!is.data.frame(lsts)) lsts <- as.data.frame(lsts)
  lsts
}
pm_bind_rows <- function(..., .id = NULL) {
  lsts <- list(...)
  lsts <- flatten(lsts)
  lsts <- Filter(Negate(is.null), lsts)
  lapply(lsts, function(x) is_df_or_vector(x))
  lapply(lsts, function(x) if (is.atomic(x) && !is_named(x)) stop("Vectors must be named."))
  if (!missing(.id)) {
    lsts <- lapply(seq_along(lsts), function(i) {
      nms <- names(lsts)
      id_df <- data.frame(id = if (is.null(nms)) as.character(i) else nms[i], stringsAsFactors = FALSE)
      colnames(id_df) <- .id
      cbind(id_df, lsts[[i]])
    })
  }
  nms <- unique(unlist(lapply(lsts, names)))
  lsts <- lapply(
    lsts,
    function(x) {
      if (!is.data.frame(x)) x <- data.frame(as.list(x), stringsAsFactors = FALSE)
      for (i in nms[!nms %in% names(x)]) x[[i]] <- NA
      x
    }
  )
  names(lsts) <- NULL
  do.call(rbind, lsts)
}
pm_case_when <- function(...) {
  fs <- list(...)
  lapply(fs, function(x) if (!inherits(x, "formula")) stop("`case_when()` requires formula inputs."))
  n <- length(fs)
  if (n == 0L) stop("No cases provided.")
  query <- vector("list", n)
  value <- vector("list", n)
  default_env <- parent.frame()
  for (i in seq_len(n)) {
    query[[i]] <- eval(fs[[i]][[2]], envir = default_env)
    value[[i]] <- eval(fs[[i]][[3]], envir = default_env)
    if (!is.logical(query[[i]])) stop(fs[[i]][[2]], " does not return a `logical` vector.")
  }
  m <- validate_case_when_length(query, value, fs)
  out <- value[[1]][rep(NA_integer_, m)]
  replaced <- rep(FALSE, m)
  for (i in seq_len(n)) {
    out <- replace_with(out, query[[i]] & !replaced, value[[i]], NULL)
    replaced <- replaced | (query[[i]] & !is.na(query[[i]]))
  }
  out
}
pm_validate_case_when_length <- function(query, value, fs) {
  lhs_lengths <- lengths(query)
  rhs_lengths <- lengths(value)
  all_lengths <- unique(c(lhs_lengths, rhs_lengths))
  if (length(all_lengths) <= 1L) return(all_lengths[[1L]])
  non_atomic_lengths <- all_lengths[all_lengths != 1L]
  len <- non_atomic_lengths[[1L]]
  if (length(non_atomic_lengths) == 1L) return(len)
  inconsistent_lengths <- non_atomic_lengths[-1L]
  lhs_problems <- lhs_lengths %in% inconsistent_lengths
  rhs_problems <- rhs_lengths %in% inconsistent_lengths
  problems <- lhs_problems | rhs_problems
  if (any(problems)) {
    stop(
      "The following formulas must be length ", len, " or 1, not ",
      paste(inconsistent_lengths, collapse = ", "), ".\n    ",
      paste(fs[problems], collapse = "\n    ")
    )
  }
}
pm_context <- new.env()
pm_context$setup <- function(.data) pm_context$.data <- .data
pm_context$get_data <- function() pm_context$.data
pm_context$get_columns <- function(cols) pm_context$.data[, cols, drop = FALSE]
pm_context$cur_column <- NULL
pm_context$get_nrow <- function() nrow(pm_context$.data)
pm_context$get_colnames <- function() colnames(pm_context$.data)
pm_context$is_grouped <- function() pm_has_groups(pm_context$.data)
pm_context$as_env <- function() {
  if (any(pm_is_nested(pm_context$.data))) {
    lapply(as.list(pm_context$.data), function(x) if (is.data.frame(x[[1]])) x[[1]] else x)
  } else {
    pm_context$.data
  }
}
pm_context$pm_group_env <- NULL
pm_context$clean <- function() {
  rm(list = c(".data"), envir = pm_context)
  if (!is.null(pm_context$cur_column)) rm(list = c("cur_column"), envir = pm_context)
}
pm_n <- function() {
  check_pm_context("`n()`", pm_context$.data)
  pm_context$get_nrow()
}
pm_cur_data <- function() {
  check_pm_context("`cur_data()`", pm_context$.data)
  data <- pm_context$get_data()
  data[, !(colnames(data) %in% pm_group_vars(data)), drop = FALSE]
}
pm_cur_data_all <- function() {
  check_pm_context("`cur_data_all()`", pm_context$.data)
  pm_ungroup(pm_context$get_data())
}
pm_cur_group <- function() {
  check_pm_context("`cur_group()`", pm_context$.data)
  data <- pm_context$get_data()
  res <- data[1L, pm_group_vars(data), drop = FALSE]
  rownames(res) <- NULL
  res
}
pm_cur_pm_group_id <- function() {
  check_pm_context("`cur_pm_group_id()`", pm_context$.data)
  data <- pm_context$get_data()
  res <- data[1L, pm_group_vars(data), drop = FALSE]
  details <- get_pm_group_details(data)
  details[, ".pm_group_id"] <- seq_len(nrow(details))
  res <- suppressMessages(semi_join(details, res))
  res[, ".pm_group_id"]
}
pm_cur_pm_group_rows <- function() {
  check_pm_context("`cur_pm_group_rows()`", pm_context$.data)
  data <- pm_context$get_data()
  res <- data[1L, pm_group_vars(data), drop = FALSE]
  res <- suppressMessages(semi_join(get_pm_group_details(data), res))
  unlist(res[, ".rows"])
}
pm_cur_column <- function() {
  check_pm_context("`cur_column()`", pm_context$cur_column, "`across`")
  pm_context$cur_column
}
pm_check_pm_context <- function(fn, pm_context, name = NULL) {
  if (is.null(pm_context)) {
    stop(fn, " must only be used inside ", if (is.null(name)) "poorman verbs" else name)
  }
}
pm_count <- function(x, ..., wt = NULL, sort = FALSE, name = NULL) {
  groups <- pm_group_vars(x)
  if (!missing(...)) x <- pm_group_by(x, ..., .add = TRUE)
  wt <- pm_deparse_var(wt)
  res <- do.call(tally, list(x, wt, sort, name))
  if (length(groups) > 0L) res <- do.call(pm_group_by, list(res, as.name(groups)))
  res
}
pm_tally <- function(x, wt = NULL, sort = FALSE, name = NULL) {
  name <- check_name(x, name)
  wt <- pm_deparse_var(wt)
  res <- do.call(pm_summarise, setNames(list(x, tally_n(x, wt)), c(".data", name)))
  res <- pm_ungroup(res)
  if (isTRUE(sort)) res <- do.call(pm_arrange, list(res, call("desc", as.name(name))))
  rownames(res) <- NULL
  res
}
pm_add_count <- function(x, ..., wt = NULL, sort = FALSE, name = NULL) {
  name <- check_name(x, name)
  row_names <- rownames(x)
  wt <- pm_deparse_var(wt)
  if (!missing(...)) x <- pm_group_by(x, ..., .add = TRUE)
  res <- do.call(add_tally, list(x, wt, sort, name))
  res[row_names, ]
}
pm_add_tally <- function(x, wt = NULL, sort = FALSE, name = NULL) {
  wt <- pm_deparse_var(wt)
  n <- tally_n(x, wt)
  name <- check_name(x, name)
  res <- do.call(pm_mutate, setNames(list(x, n), c(".data", name)))
  if (isTRUE(sort)) {
    do.call(pm_arrange, list(res, call("desc", as.name(name))))
  } else {
    res
  }
}
pm_tally_n <- function(x, wt) {
  if (is.null(wt) && "n" %in% colnames(x)) {
    message("Using `n` as weighting variable")
    wt <- "n"
  }
  pm_context$setup(.data = x)
  on.exit(pm_context$clean(), add = TRUE)
  if (is.null(wt)) {
    call("n")
  } else {
    call("sum", as.name(wt), na.rm = TRUE)
  }
}
pm_check_name <- function(df, name) {
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
pm_desc <- function(x) -xtfrm(x)
pm_distinct <- function(.data, ..., .keep_all = FALSE) {
  if ("grouped_df" %in% class(.data)) pm_distinct.grouped_df(.data, ..., .keep_all = FALSE) else pm_distinct.data.frame(.data, ..., .keep_all = FALSE)
}
pm_distinct.data.frame <- function(.data, ..., .keep_all = FALSE) {
  if (ncol(.data) == 0L) return(.data[1, ])
  cols <- pm_dotdotdot(...)
  col_names <- names(cols)
  col_len <- length(cols)
  if (is.null(col_names) && col_len > 0L) names(cols) <- cols
  if (col_len == 0L) {
    res <- .data
  } else {
    mut <- pm_mutate_df(.data, ...)
    res <- mut$data
    col_names <- names(cols)
    res <- if (!is.null(col_names)) {
      zero_names <- nchar(col_names) == 0L
      if (any(zero_names)) {
        names(cols)[zero_names] <- cols[zero_names]
        col_names <- names(cols)
      }
      suppressMessages(select(res, col_names))
    } else {
      suppressMessages(select(res, cols))
    }
  }
  res <- unique(res)
  if (isTRUE(.keep_all)) {
    res <- cbind(res, .data[rownames(res), setdiff(colnames(.data), colnames(res)), drop = FALSE])
  }
  common_cols <- c(intersect(colnames(.data), colnames(res)), setdiff(col_names, colnames(.data)))
  if (is.numeric(attr(res, "row.names"))) {
    row.names(res) <- seq_len(nrow(res))
  }
  if (length(common_cols) > 0L) res[, common_cols, drop = FALSE] else res
}
pm_distinct.grouped_df <- function(.data, ..., .keep_all = FALSE) {
  pm_apply_grouped_function("pm_distinct", .data, drop = TRUE, ..., .keep_all = .keep_all)
}
pm_dotdotdot <- function(..., .impute_names = FALSE) {
  dots <- eval(substitute(alist(...)))
  if (isTRUE(.impute_names)) {
    pm_deparse_dots <- lapply(dots, deparse)
    names_dots <- names(dots)
    unnamed <- if (is.null(names_dots)) rep(TRUE, length(dots)) else nchar(names_dots) == 0L
    names(dots)[unnamed] <- pm_deparse_dots[unnamed]
  }
  dots
}
pm_deparse_dots <- function(...) {
  vapply(substitute(...()), deparse, NA_character_)
}
pm_deparse_var <- function(var, frame = if (is.null(pm_eval_env$env)) parent.frame() else pm_eval_env$env) {
  sub_var <- eval(substitute(substitute(var)), frame)
  if (is.symbol(sub_var)) var <- as.character(sub_var)
  var
}
pm_eval_env <- new.env()
pm_filter <- function(.data, ..., .preserve = FALSE) {
  if ("grouped_df" %in% class(.data)) pm_filter.grouped_df(.data, ..., .preserve = FALSE) else pm_filter.data.frame(.data, ..., .preserve = FALSE)
}
pm_filter.data.frame <- function(.data, ..., .preserve = FALSE) {
  conditions <- pm_dotdotdot(...)
  if (length(conditions) == 0L) return(.data)
  check_filter(conditions)
  cond_class <- vapply(conditions, typeof, NA_character_)
  cond_class <- cond_class[!cond_class %in% c("language", "logical")]
  if (length(cond_class) > 0L) stop("Conditions must be logical vectors")
  pm_context$setup(.data)
  on.exit(pm_context$clean(), add = TRUE)
  pm_eval_env$env <- parent.frame()
  on.exit(rm(list = "env", envir = pm_eval_env), add = TRUE)
  rows <- lapply(
    conditions,
    function(cond, frame) eval(cond, pm_context$.data, frame),
    frame = pm_eval_env$env
  )
  rows <- Reduce("&", rows)
  .data[rows & !is.na(rows), ]
}
pm_filter.grouped_df <- function(.data, ..., .preserve = FALSE) {
  rows <- rownames(.data)
  res <- pm_apply_grouped_function("pm_filter", .data, drop = TRUE, ...)
  res <- res[rows[rows %in% rownames(res)], ]
  groups <- pm_group_vars(.data)
  pre_filtered_groups <- pm_group_data(.data)
  post_filtered_groups <- pm_calculate_groups(res, groups)
  if (!(!.preserve && isTRUE(attr(pre_filtered_groups, ".drop")))) {
    filtered_groups <- anti_join(pre_filtered_groups, post_filtered_groups, by = groups)
    filtered_groups <- filtered_groups[, groups, drop = FALSE]
    filtered_groups[[".rows"]] <- rep(list(integer()), length.out = nrow(filtered_groups))
    post_filtered_groups <- bind_rows(post_filtered_groups, filtered_groups)
    ordered <- do.call(pm_arrange_rows, list(post_filtered_groups, pm_as_symbols(groups)))
    post_filtered_groups <- post_filtered_groups[ordered, ]
  }
  attr(res, "groups") <- post_filtered_groups
  res
}
pm_check_filter <- function(conditions) {
  named <- have_name(conditions)
  for (i in which(named)) {
    if (!is.logical(conditions[[i]])) {
      stop(
        sprintf("Problem with `pm_filter()` input `..%s`.\n", i),
        sprintf("Input `..%s` is named.\n", i),
        "This usually means that you've used `=` instead of `==`.\n",
        sprintf("Did you mean `%s == %s`?", names(conditions)[[i]], conditions[[i]])
      )
    }
  }
}
pm_group_by <- function(.data, ..., .add = FALSE, .drop = pm_group_by_drop_default(.data)) {
  pm_group_by.data.frame(.data, ..., .add = FALSE, .drop = pm_group_by_drop_default(.data))
}
pm_group_by.data.frame <- function(.data, ..., .add = FALSE, .drop = pm_group_by_drop_default(.data)) {
  vars <- pm_dotdotdot(..., .impute_names = TRUE)
  if (all(vapply(vars, is.null, FALSE))) {
    res <- pm_groups_set(.data, NULL)
    class(res) <- class(res)[!(class(res) %in% "grouped_df")]
    return(res)
  }
  new_cols <- pm_add_group_columns(.data, vars)
  res <- new_cols$data
  groups <- new_cols$groups
  if (isTRUE(.add)) groups <- union(pm_group_vars(.data), groups)
  unknown <- !(groups %in% colnames(res))
  if (any(unknown)) stop("Invalid groups: ", groups[unknown])
  if (length(groups) > 0L) {
    res <- pm_groups_set(res, groups, .drop)
    class(res) <- union("grouped_df", class(res))
  }
  res
}
pm_group_by_drop_default <- function(.tbl) {
  if ("grouped_df" %in% class(.tbl)) pm_group_by_drop_default.grouped_df(.tbl) else pm_group_by_drop_default.data.frame(.tbl)
}
pm_group_by_drop_default.data.frame <- function(.tbl) {
  TRUE
}
pm_group_by_drop_default.grouped_df <- function(.tbl) {
  tryCatch({
    !identical(attr(pm_group_data(.tbl), ".drop"), FALSE)
  }, error = function(e) {
    TRUE
  })
}
pm_add_group_columns <- function(.data, vars) {
  vars <- vars[!vapply(vars, is.null, FALSE)]
  types <- do.call(c, lapply(vars, typeof))
  test <- any(types == "language")
  needs_mutate <- if (test) unname(which(types == "language")) else NULL
  if (!is.null(needs_mutate)) {
    .data <- do.call(pm_mutate, c(list(.data = pm_ungroup(.data)), vars[needs_mutate]))
  }
  list(data = .data, groups = names(vars))
}
pm_group_data <- function(.data) {
  if ("grouped_df" %in% class(.data)) pm_group_data.grouped_df(.data) else pm_group_data.data.frame(.data)
}
pm_group_data.data.frame <- function(.data) {
  structure(list(.rows = list(seq_len(nrow(.data)))), class = "data.frame", row.names = c(NA, -1L))
}
pm_group_data.grouped_df <- function(.data) {
  attr(.data, "groups")
}
pm_group_rows <- function(.data) {
  pm_group_data(.data)[[".rows"]]
}
pm_group_indices <- function(.data) {
  if (!pm_has_groups(.data)) return(rep(1L, nrow(.data)))
  groups <- pm_group_vars(.data)
  res <- unique(.data[, groups, drop = FALSE])
  res <- res[do.call(order, lapply(groups, function(x) res[, x])), , drop = FALSE]
  class(res) <- "data.frame"
  nrow_data <- nrow(.data)
  rows <- rep(NA, nrow_data)
  for (i in seq_len(nrow_data)) {
    rows[i] <- which(interaction(res[, groups]) %in% interaction(.data[i, groups]))
  }
  rows
}
pm_group_vars <- function(x) {
  groups <- attr(x, "groups", exact = TRUE)
  if (is.null(groups)) character(0) else colnames(groups)[!colnames(groups) %in% c(".pm_group_id", ".rows")]
}
pm_groups <- function(x) {
  pm_as_symbols(pm_group_vars(x))
}
pm_group_size <- function(x) {
  lengths(pm_group_rows(x))
}
pm_n_groups <- function(x) {
  nrow(pm_group_data(x))
}
pm_group_split <- function(.data, ..., .keep = TRUE) {
  dots_len <- length(pm_dotdotdot(...)) > 0L
  if (pm_has_groups(.data) && isTRUE(dots_len)) {
    warning("... is ignored in pm_group_split(<grouped_df>), please use pm_group_by(..., .add = TRUE) %pm>% pm_group_split()")
  }
  if (!pm_has_groups(.data) && isTRUE(dots_len)) {
    .data <- pm_group_by(.data, ...)
  }
  if (!pm_has_groups(.data) && !isTRUE(dots_len)) {
    return(list(.data))
  }
  pm_context$setup(.data)
  on.exit(pm_context$clean(), add = TRUE)
  groups <- pm_group_vars(.data)
  attr(pm_context$.data, "groups") <- NULL
  res <- pm_split_into_groups(pm_context$.data, groups)
  names(res) <- NULL
  if (!isTRUE(.keep)) {
    res <- lapply(res, function(x) x[, !colnames(x) %in% groups])
  }
  any_empty <- unlist(lapply(res, function(x) !(nrow(x) == 0L)))
  res[any_empty]
}
pm_group_keys <- function(.data) {
  groups <- pm_group_vars(.data)
  pm_context$setup(.data)
  res <- pm_context$get_columns(pm_context$get_colnames() %in% groups)
  res <- res[!duplicated(res), , drop = FALSE]
  if (nrow(res) == 0L) return(res)
  class(res) <- "data.frame"
  res <- res[do.call(order, lapply(groups, function(x) res[, x])), , drop = FALSE]
  rownames(res) <- NULL
  res
}
pm_split_into_groups <- function(.data, groups, drop = FALSE, ...) {
  class(.data) <- "data.frame"
  pm_group_factors <- lapply(groups, function(x, .data) as.factor(.data[, x]), .data)
  split(x = .data, f = pm_group_factors, drop = drop, ...)
}
pm_groups_set <- function(x, groups, drop = pm_group_by_drop_default(x)) {
  attr(x, "groups") <- if (is.null(groups) || length(groups) == 0L) {
    NULL
  } else {
    pm_calculate_groups(x, groups, drop)
  }
  x
}
pm_get_pm_group_details <- function(x) {
  groups <- attr(x, "groups", exact = TRUE)
  if (is.null(groups)) character(0) else groups
}
pm_has_groups <- function(x) {
  groups <- pm_group_vars(x)
  if (length(groups) == 0L) FALSE else TRUE
}
pm_apply_grouped_function <- function(fn, .data, drop = FALSE, ...) {
  groups <- pm_group_vars(.data)
  grouped <- pm_split_into_groups(.data, groups, drop)
  res <- do.call(rbind, unname(lapply(grouped, fn, ...)))
  if (any(groups %in% colnames(res))) {
    class(res) <- c("grouped_df", class(res))
    res <- pm_groups_set(res, groups[groups %in% colnames(res)])
  }
  res
}
pm_calculate_groups <- function(data, groups, drop = pm_group_by_drop_default(data)) {
  data <- pm_ungroup(data)
  unknown <- setdiff(groups, colnames(data))
  if (length(unknown) > 0L) {
    stop(sprintf("`groups` missing from `data`: %s.", paste0(groups, collapse = ", ")))
  }
  unique_groups <- unique(data[, groups, drop = FALSE])
  is_factor <- do.call(c, lapply(unique_groups, function(x) is.factor(x)))
  n_comb <- nrow(unique_groups)
  rows <- rep(list(NA), n_comb)
  data_groups <- interaction(data[, groups, drop = TRUE])
  for (i in seq_len(n_comb)) {
    rows[[i]] <- which(data_groups %in% interaction(unique_groups[i, groups]))
  }
  if (!isTRUE(drop) && any(is_factor)) {
    na_lvls <- do.call(
      expand.grid,
      lapply(unique_groups, function(x) if (is.factor(x)) levels(x)[!(levels(x) %in% x)] else NA)
    )
    unique_groups <- rbind(unique_groups, na_lvls)
    for (i in seq_len(nrow(na_lvls))) {
      rows[[length(rows) + 1]] <- integer(0)
    }
  }
  unique_groups[[".rows"]] <- rows
  unique_groups <- unique_groups[do.call(order, lapply(groups, function(x) unique_groups[, x])), , drop = FALSE]
  rownames(unique_groups) <- NULL
  structure(unique_groups, .drop = drop)
}
pm_is.grouped_df <- function(x) {
  inherits(x, "grouped_df")
}
pm_if_else <- function(condition, true, false, missing = NULL) {
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
pm_anti_join <- function(x, y, by = NULL) {
  pm_filter_join_worker(x, y, by, type = "anti")
}
pm_semi_join <- function(x, y, by = NULL) {
  pm_filter_join_worker(x, y, by, type = "semi")
}
pm_filter_join_worker <- function(x, y, by = NULL, type = c("anti", "semi")) {
  type <- match.arg(type, choices = c("anti", "semi"), several.ok = FALSE)
  if (is.null(by)) {
    by <- intersect(names(x), names(y))
    join_message(by)
  }
  rows <- interaction(x[, by]) %in% interaction(y[, by])
  if (type == "anti") rows <- !rows
  res <- x[rows, , drop = FALSE]
  rownames(res) <- NULL
  reconstruct_attrs(res, x)
}
pm_inner_join <- function(x, y, by = NULL, suffix = c(".x", ".y"), ..., na_matches = c("na", "never")) {
  join_worker(x = x, y = y, by = by, suffix = suffix, sort = FALSE, ..., keep = FALSE, na_matches = na_matches)
}
pm_left_join <- function(x, y, by = NULL, suffix = c(".x", ".y"), ..., keep = FALSE, na_matches = c("na", "never")) {
  join_worker(x = x, y = y, by = by, suffix = suffix, all.x = TRUE, ..., keep = keep, na_matches = na_matches)
}
pm_right_join <- function(x, y, by = NULL, suffix = c(".x", ".y"), ..., keep = FALSE, na_matches = c("na", "never")) {
  join_worker(x = x, y = y, by = by, suffix = suffix, all.y = TRUE, ..., keep = keep, na_matches = na_matches)
}
pm_full_join <- function(x, y, by = NULL, suffix = c(".x", ".y"), ..., keep = FALSE, na_matches = c("na", "never")) {
  join_worker(x = x, y = y, by = by, suffix = suffix, all = TRUE, ..., keep = keep, na_matches = na_matches)
}
pm_join_worker <- function(x, y, by = NULL, suffix = c(".x", ".y"), keep = FALSE, na_matches = c("na", "never"), ...) {
  na_matches <- match.arg(arg = na_matches, choices = c("na", "never"), several.ok = FALSE)
  incomparables <- if (na_matches == "never") NA else NULL
  x[, ".join_id"] <- seq_len(nrow(x))
  merged <- if (is.null(by)) {
    by <- intersect(names(x), names(y))
    join_message(by)
    merge(
      x = x, y = y, by = by, suffixes = suffix, incomparables = incomparables, ...
    )[, union(names(x), names(y)), drop = FALSE]
  } else if (is.null(names(by))) {
    merge(x = x, y = y, by = by, suffixes = suffix, incomparables = incomparables, ...)
  } else {
    merge(x = x, y = y, by.x = names(by), by.y = by, suffixes = suffix, incomparables = incomparables, ...)
  }
  merged <- merged[order(merged[, ".join_id"]), colnames(merged) != ".join_id", drop = FALSE]
  if (isTRUE(keep)) {
    keep_pos <- match(by, names(merged))
    x_by <- paste0(by, suffix[1L])
    colnames(merged)[keep_pos] <- x_by
    merged[, paste0(by, suffix[2L])] <- merged[, x_by]
  }
  rownames(merged) <- NULL
  reconstruct_attrs(merged, x)
}
pm_join_message <- function(by) {
  if (length(by) > 1L) {
    message("Joining, by = c(\"", paste0(by, collapse = "\", \""), "\")\n", sep = "")
  } else {
    message("Joining, by = \"", by, "\"\n", sep = "")
  }
}
pm_as_function <- function(x, env = parent.frame()) {
  if (is.function(x)) return(x)
  if (is_formula(x)) {
    if (length(x) > 2) stop("Can't convert a two-sided formula to a function")
    env <- attr(x, ".Environment", exact = TRUE)
    rhs <- as.list(x)[[2]]
    return(as.function(list(... = substitute(), .x = quote(..1), .y = quote(..2), . = quote(..1), rhs), envir = env))
  }
  if (is_string(x)) return(get(x, envir = env, mode = "function"))
  stop("Can't convert an object of class ", class(x), " to a function.")
}
pm_is_formula <- function(x) {
  inherits(x, "formula")
}
pm_is_string <- function(x) {
  is.character(x) && length(x) == 1L
}
pm_is_wholenumber <- function(x) {
  x %% 1L == 0L
}
pm_names_are_invalid <- function(x) {
  x == "" | is.na(x)
}
pm_is_named <- function(x) {
  nms <- names(x)
  if (is.null(nms)) return(FALSE)
  if (any(names_are_invalid(nms))) return(FALSE)
  TRUE
}
pm_have_name <- function(x) {
  nms <- names(x)
  if (is.null(nms)) rep(FALSE, length(x)) else !names_are_invalid(nms)
}
pm_is_empty_list <- function(x) {
  inherits(x, "list") && length(x) == 0L
}
pm_as_symbols <- function(x) {
  lapply(x, as.symbol)
}
pm_is_df_or_vector <- function(x) {
  res <- is.data.frame(x) || is.atomic(x)
  if (!isTRUE(res)) stop("You must pass vector(s) and/or data.frame(s).")
  TRUE
}
pm_lag <- function(x, n = 1L, default = NA) {
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
pm_lead <- function(x, n = 1L, default = NA) {
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
pm_lst <- function(...) {
  fn_call <- match.call()
  list_to_eval <- as.list(fn_call)[-1]
  out <- vector(mode = "list", length = length(list_to_eval))
  names(out) <- names(list_to_eval)
  exprs <- lapply(substitute(list(...)), deparse)[-1]
  for (element in seq_along(list_to_eval)) {
    value <- list_to_eval[[element]]
    if (is.language(value)) {
      value <- eval(
        value,
        envir = if (length(out) == 0) {
          list_to_eval
        } else {
          drop_dup_list(out[1:(element - 1)])
        }
      )
    }
    if (is.null(value)) {
      out[element] <- list(NULL)
    } else {
      out[[element]] <- value
    }
    invalid_name <- is.null(names(out)[element]) ||
      is.na(names(out)[element]) ||
      names(out)[element] == ""
    if (invalid_name) {
      if (exprs[[element]] != "NULL" || (exprs[[element]] == "NULL" && is.null(out[[element]]))) {
        names(out)[element] <- exprs[[element]]
      }
    }
  }
  out
}
pm_drop_dup_list <- function(x) {
  list_names <- names(x)
  if (identical(list_names, unique(list_names))) return(x)
  count <- table(list_names)
  dupes <- names(count[count > 1])
  uniques <- names(count[count == 1])
  to_drop <- do.call(c, lapply(
    dupes,
    function(x) {
      matches <- which(list_names == x)
      matches[-length(matches)]
    }
  ))
  x[uniques] <- Filter(Negate(is.null), x[uniques])
  return(x[-to_drop])
}
pm_mutate <- function(.data, ...) {
  if ("grouped_df" %in% class(.data)) pm_mutate.grouped_df(.data, ...) else pm_mutate.data.frame(.data, ...)
}
pm_mutate.data.frame <- function(
  .data,
  ...,
  .keep = c("all", "used", "unused", "none"),
  .before = NULL,
  .after = NULL
) {
  keep <- match.arg(arg = .keep, choices = c("all", "used", "unused", "none"), several.ok = FALSE)
  res <- pm_mutate_df(.data = .data, ...)
  data <- res$data
  new_cols <- res$new_cols
  .before <- substitute(.before)
  .after <- substitute(.after)
  if (!is.null(.before) || !is.null(.after)) {
    new <- setdiff(new_cols, names(.data))
    data <- do.call(pm_relocate, c(list(.data = data), new, .before = .before, .after = .after))
  }
  if (keep == "all") {
    data
  } else if (keep == "unused") {
    unused <- setdiff(colnames(.data), res$used_cols)
    keep <- intersect(colnames(data), c(pm_group_vars(.data), unused, new_cols))
    select(.data = data, keep)
  } else if (keep == "used") {
    keep <- intersect(colnames(data), c(pm_group_vars(.data), res$used_cols, new_cols))
    select(.data = data, keep)
  } else if (keep == "none") {
    keep <- c(setdiff(pm_group_vars(.data), new_cols), intersect(new_cols, colnames(data)))
    select(.data = data, keep)
  }
}
pm_mutate.grouped_df <- function(.data, ...) {
  pm_context$pm_group_env <- parent.frame(n = 1)
  on.exit(rm(list = c("pm_group_env"), envir = pm_context), add = TRUE)
  rows <- rownames(.data)
  res <- pm_apply_grouped_function("pm_mutate", .data, drop = TRUE, ...)
  res[rows, , drop = FALSE]
}
pm_mutate_df <- function(.data, ...) {
  conditions <- pm_dotdotdot(..., .impute_names = TRUE)
  cond_nms <- names(pm_dotdotdot(..., .impute_names = FALSE))
  if (length(conditions) == 0L) {
    return(list(
      data = .data,
      used_cols = NULL,
      new_cols = NULL
    ))
  }
  used <- unname(do.call(c, lapply(conditions, pm_find_used)))
  used <- used[used %in% colnames(.data)]
  pm_context$setup(.data)
  on.exit(pm_context$clean(), add = TRUE)
  for (i in seq_along(conditions)) {
    not_named <- (is.null(cond_nms) || cond_nms[i] == "")
    res <- eval(
      conditions[[i]],
      envir = pm_context$as_env(),
      enclos = if (!is.null(pm_context$pm_group_env)) pm_context$pm_group_env else parent.frame(n = 2)
    )
    res_nms <- names(res)
    if (is.data.frame(res)) {
      if (not_named) {
        pm_context$.data[, res_nms] <- res
      } else {
        pm_context$.data[[cond_nms[i]]] <- res
      }
    } else if (is.atomic(res)) {
      cond_nms[i] <- names(conditions)[[i]]
      pm_context$.data[[cond_nms[i]]] <- res
    } else {
      if (is.null(res_nms)) names(res) <- names(conditions)[[i]]
      pm_context$.data[[names(res)]] <- res
    }
  }
  list(
    data = pm_context$.data,
    used_cols = used,
    new_cols = cond_nms
  )
}
pm_find_used <- function(expr) {
  if (is.symbol(expr)) {
    as.character(expr)
  } else {
    unique(unlist(lapply(expr[-1], pm_find_used)))
  }
}
pm_n_distinct <- function(..., na.rm = FALSE) {
  res <- do.call(cbind, list(...))
  if (isTRUE(na.rm)) res <- res[!is.na(res), , drop = FALSE]
  nrow(unique(res))
}
pm_nth <- function(x, n, order_by = NULL, default = pm_default_missing(x)) {
  if (length(n) != 1 || !is.numeric(n)) stop("`n` must be a single integer.")
  n <- trunc(n)
  if (n == 0 || n > length(x) || n < -length(x)) return(default)
  if (n < 0) n <- length(x) + n + 1
  if (is.null(order_by)) x[[n]] else x[[order(order_by)[[n]]]]
}
pm_first <- function(x, order_by = NULL, default = pm_default_missing(x)) {
  nth(x, 1L, order_by = order_by, default = default)
}
pm_last <- function(x, order_by = NULL, default = pm_default_missing(x)) {
  nth(x, -1L, order_by = order_by, default = default)
}
pm_default_missing <- function(x) {
  pm_default_missing.data.frame(x)
}
pm_default_missing.data.frame <- function(x) {
  if (!is.object(x) && is.list(x)) NULL else x[NA_real_]
}
pm_default_missing.data.frame <- function(x) {
  rep(NA, nrow(x))
}
`%pm>%` <- function(lhs, rhs) {
  rhs_call <- pm_insert_dot(substitute(rhs))
  eval(rhs_call, envir = list(`.` = lhs), enclos = parent.frame())
}
pm_insert_dot <- function(expr) {
  if (is.symbol(expr) || expr[[1]] == quote(`(`)) {
    expr <- as.call(c(expr, quote(`.`)))
  } else if (length(expr) == 1) {
    expr <- as.call(c(expr[[1]], quote(`.`)))
  } else if (
    expr[[1]] != quote(`{`) &&
    !any(vapply(expr[-1], identical, quote(`.`), FUN.VALUE = logical(1))) &&
    !any(vapply(expr[-1], identical, quote(`!!!.`), FUN.VALUE = logical(1)))
  ) {
    expr <- as.call(c(expr[[1]], quote(`.`), as.list(expr[-1])))
  }
  expr
}
pm_pivot_longer <- function(
  data,
  cols,
  names_to = "name",
  names_prefix = NULL,
  names_sep = NULL,
  names_pattern = NULL,
  values_to = "value",
  values_drop_na = FALSE,
  ...
) {
  if (missing(cols)) {
    stop("`cols` must select at least one column.")
  }
  cols <- names(pm_eval_select_pos(data, substitute(cols)))
  if (any(names_to %in% setdiff(names(data), cols))) {
    stop(
      paste0(
        "Some values of the columns specified in 'names_to' are already present
        as column names. Either use another value in `names_to` or pm_rename the
        following columns: ",
        paste(names_to[which(names_to %in% setdiff(names(data), cols))], sep = ", ")
      ),
      call. = FALSE)
  }
  if (length(cols) == 0L) {
    stop("No columns found for reshaping data.", call. = FALSE)
  }
  data[["_Row"]] <- as.numeric(rownames(data))
  names_to_2 <- paste(names_to, collapse = "_")
  long <- stats::reshape(
    as.data.frame(data, stringsAsFactors = FALSE),
    varying = cols,
    idvar = "_Row",
    v.names = values_to,
    timevar = names_to_2,
    direction = "long"
  )
  long <- long[do.call(order, long[, c("_Row", names_to_2)]), ]
  long[["_Row"]] <- NULL
  long[[names_to_2]] <- cols[long[[names_to_2]]]
  if (length(names_to) > 1) {
    if (is.null(names_pattern)) {
      for (i in seq_along(names_to)) {
        new_vals <- unlist(lapply(
          strsplit(unique(long[[names_to_2]]), names_sep, fixed = TRUE),
          function(x) x[i]
        ))
        long[[names_to[i]]] <- new_vals
      }
    } else {
      tmp <- regmatches(
        unique(long[[names_to_2]]),
        regexec(names_pattern, unique(long[[names_to_2]]))
      )
      tmp <- as.data.frame(do.call(rbind, tmp), stringsAsFactors = FALSE)
      names(tmp) <- c(names_to_2, names_to)
      long <- cbind(long, tmp[match(long[[names_to_2]], tmp[[names_to_2]]), -1])
    }
    long[[names_to_2]] <- NULL
  }
  long <- pm_relocate(.data = long, "value", .after = -1)
  if (!is.null(names_prefix)) {
    if (length(names_to) > 1) {
      stop("`names_prefix` only works when `names_to` is of length 1.", call. = FALSE)
    }
    long[[names_to]] <- gsub(paste0("^", names_prefix), "", long[[names_to]])
  }
  if (values_drop_na) {
    long <- long[!is.na(long[, values_to]), ]
  }
  rownames(long) <- NULL
  attributes(long)$reshapeLong <- NULL
  long
}
pm_pivot_wider <- function(
  data,
  id_cols = NULL,
  values_from = "Value",
  names_from = "Name",
  names_sep = "_",
  names_prefix = "",
  names_glue = NULL,
  values_fill = NULL,
  ...
) {
  old_names <- names(data)
  names_from <- names(pm_eval_select_pos(data, substitute(names_from)))
  values_from <- names(pm_eval_select_pos(data, substitute(values_from)))
  variable_attr <- lapply(data, attributes)
  if (is.null(id_cols)) {
    row_index <- do.call(
      paste, 
      c(data[, !names(data) %in% c(values_from, names_from), drop = FALSE], sep = "_")
    )
    if (length(row_index) == 0) row_index <- rep("", nrow(data))
    data[["_Rows"]] <- row_index
    id_cols <- "_Rows"
  }
  current_colnames <- colnames(data)
  current_colnames <- current_colnames[current_colnames != "_Rows"]
  if (is.null(names_glue)) {
    future_colnames <- unique(do.call(paste, c(data[, names_from, drop = FALSE], sep = names_sep)))
  } else {
    vars <- regmatches(names_glue, gregexpr("\\{\\K[^{}]+(?=\\})", names_glue, perl = TRUE))[[1]]
    tmp_data <- unique(data[, vars])
    future_colnames <- unique(apply(tmp_data, 1, function(x) {
      tmp_vars <- list()
      for (i in seq_along(vars)) {
        tmp_vars[[i]] <- x[vars[i]]
      }
      tmp_colname <- gsub("\\{\\K[^{}]+(?=\\})", "", names_glue, perl = TRUE)
      tmp_colname <- gsub("\\{\\}", "%s", tmp_colname)
      do.call(sprintf, c(fmt = tmp_colname, tmp_vars))
    }))
  }
  if (any(future_colnames %in% current_colnames)) {
    stop(
      paste0(
        "Some values of the columns specified in 'names_from' are already present
        as column names. Either use `name_prefix` or pm_rename the following columns: ",
        paste(current_colnames[which(current_colnames %in% future_colnames)], sep = ", ")
      ),
      call. = FALSE
    )
  }
  data$new_time <- do.call(paste, c(data[, names_from, drop = FALSE], sep = "_"))
  data[, names_from] <- NULL
  wide <- stats::reshape(
    as.data.frame(data, stringsAsFactors = FALSE),
    v.names = values_from,
    idvar = id_cols,
    timevar = "new_time",
    sep = names_sep,
    direction = "wide"
  )
  if ("_Rows" %in% names(wide)) wide[["_Rows"]] <- NULL
  rownames(wide) <- NULL
  if (length(values_from) == 1) {
    to_rename <- which(startsWith(names(wide), paste0(values_from, names_sep)))
    names(wide)[to_rename] <- future_colnames
  }
  if (length(values_from) > 1) {
    for (i in values_from) {
      tmp1 <- wide[, which(!startsWith(names(wide), i))]
      tmp2 <- wide[, which(startsWith(names(wide), i))]
      wide <- cbind(tmp1, tmp2)
    }
  }
  new_cols <- setdiff(names(wide), old_names)
  names(wide)[which(names(wide) %in% new_cols)] <- paste0(names_prefix, new_cols)
  if (!is.null(values_fill)) {
    if (length(values_fill) == 1) {
      if (is.numeric(wide[[new_cols[1]]])) {
        if (!is.numeric(values_fill)) {
          stop(paste0("`values_fill` must be of type numeric."), call. = FALSE)
        } else {
          for (i in new_cols) {
            wide[[i]] <- replace_na(wide[[i]], replace = values_fill)
          }
        }
      } else if (is.character(wide[[new_cols[1]]])) {
        if (!is.character(values_fill)) {
          stop(paste0("`values_fill` must be of type character."), call. = FALSE)
        } else {
          for (i in new_cols) {
            wide[[i]] <- replace_na(wide[[i]], replace = values_fill)
          }
        }
      } else if (is.factor(wide[[new_cols[1]]])) {
        if (!is.factor(values_fill)) {
          stop(paste0("`values_fill` must be of type factor."), call. = FALSE)
        } else {
          for (i in new_cols) {
            wide[[i]] <- replace_na(wide[[i]], replace = values_fill)
          }
        }
      }
    } else {
      stop("`values_fill` must be of length 1.", call. = FALSE)
    }
  }
  attributes(wide)$reshapeWide <- NULL
  for (i in colnames(wide)) {
    attributes(wide[[i]]) <- variable_attr[[i]]
  }
  wide
}
pm_pull <- function(.data, var = -1) {
  var_list <- as.list(seq_along(.data))
  names(var_list) <- names(.data)
  .var <- eval(substitute(var), var_list)
  if (.var < 0L) .var <- length(var_list) + .var + 1L
  .data[[.var]]
}
pm_relocate <- function(.data, ..., .before = NULL, .after = NULL) {
  pm_relocate.data.frame(.data, ..., .before = NULL, .after = NULL)
}
pm_relocate.data.frame <- function(.data, ..., .before = NULL, .after = NULL) {
  data_names <- colnames(.data)
  col_pos <- pm_select_positions(.data, ...)
  if (!missing(.before) && !is.null(.before)) .before <- colnames(.data)[pm_eval_select_pos(.data, substitute(.before))]
  if (!missing(.after) && !is.null(.after)) .after <- colnames(.data)[pm_eval_select_pos(.data, substitute(.after))]
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
  if (pm_has_groups(.data)) res <- pm_groups_set(res, pm_group_vars(.data))
  res
}
pm_rename <- function(.data, ...) {
  pm_rename.data.frame(.data, ...)
}
pm_rename.data.frame <- function(.data, ...) {
  new_names <- names(pm_dotdotdot(...))
  if (length(new_names) == 0L) {
    warning("You didn't give any new names")
    return(.data)
  }
  col_pos <- pm_select_positions(.data, ...)
  old_names <- colnames(.data)[col_pos]
  new_names_zero <- nchar(new_names) == 0L
  if (any(new_names_zero)) {
    warning("You didn't provide new names for: ", paste0("`", old_names[new_names_zero], collapse = ", "), "`")
    new_names[new_names_zero] <- old_names[new_names_zero]
  }
  colnames(.data)[col_pos] <- new_names
  .data
}
pm_rename_with <- function(.data, .fn, .cols = everything(), ...) {
  pm_rename_with.data.frame(.data, .fn, .cols = everything(), ...)
}
pm_rename_with.data.frame <- function(.data, .fn, .cols = everything(), ...) {
  if (!is.function(.fn)) stop("`", .fn, "` is not a valid function")
  grouped <- is.grouped_df(.data)
  if (grouped) grp_pos <- which(colnames(.data) %in% pm_group_vars(.data))
  col_pos <- pm_eval_select_pos(.data = .data, .pm_group_pos = TRUE, .cols = substitute(.cols))
  cols <- colnames(.data)[col_pos]
  new_cols <- .fn(cols, ...)
  if (any(duplicated(new_cols))) {
    stop("New names must be unique however `", deparse(substitute(.fn)), "` returns duplicate column names")
  }
  colnames(.data)[col_pos] <- new_cols
  if (grouped) .data <- pm_groups_set(.data, colnames(.data)[grp_pos])
  .data
}
pm_starts_with <- function(match, ignore.case = TRUE, vars = peek_vars()) {
  grep(pattern = paste0("^", paste0(match, collapse = "|^")), x = vars, ignore.case = ignore.case)
}
pm_ends_with <- function(match, ignore.case = TRUE, vars = peek_vars()) {
  grep(pattern = paste0(paste0(match, collapse = "$|"), "$"), x = vars, ignore.case = ignore.case)
}
pm_contains <- function(match, ignore.case = TRUE, vars = peek_vars()) {
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
  unique(unlist(matches))
}
pm_matches <- function(match, ignore.case = TRUE, perl = FALSE, vars = peek_vars()) {
  grep(pattern = match, x = vars, ignore.case = ignore.case, perl = perl)
}
pm_num_range <- function(prefix, range, width = NULL, vars = peek_vars()) {
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
pm_all_of <- function(x, vars = peek_vars()) {
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
pm_any_of <- function(x, vars = peek_vars()) {
  which(vars %in% x)
}
pm_everything <- function(vars = peek_vars()) {
  seq_along(vars)
}
pm_last_col <- function(offset = 0L, vars = peek_vars()) {
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
pm_peek_vars <- function() {
  pm_select_env$get_colnames()
}
pm_select_positions <- function(.data, ..., .pm_group_pos = FALSE) {
  cols <- pm_dotdotdot(...)
  cols <- cols[!vapply(cols, is.null, FALSE)]
  if (length(cols) == 0L) return(integer(0))
  pm_select_env$setup(.data = .data, calling_frame = parent.frame(2L))
  on.exit(pm_select_env$clean(), add = TRUE)
  data_names <- pm_select_env$get_colnames()
  pos <- unlist(lapply(cols, pm_eval_expr))
  if (length(pos) > 0) pos <- if (pos[1] >= 0) pos[pos >= 0] else pos[pos < 0]
  col_len <- pm_select_env$get_ncol()
  if (any(pos > col_len)) {
    oor <- pos[which(pos > col_len)]
    oor_len <- length(oor)
    stop(
      "Location", if (oor_len > 1) "s " else " ", collapse_to_sentence(oor),
      if (oor_len > 1) " don't " else " doesn't ", "exist. There are only ", col_len, " columns."
    )
  }
  if (isTRUE(.pm_group_pos)) {
    groups <- pm_group_vars(.data)
    missing_groups <- !(groups %in% cols)
    if (any(missing_groups)) {
      sel_missing <- groups[missing_groups]
      readd <- match(sel_missing, data_names)
      readd <- readd[!(readd %in% pos)]
      if (length(readd) > 0L) {
        message("Adding missing grouping variables: `", paste(sel_missing, collapse = "`, `"), "`")
        if (length(names(cols)) > 0L) names(readd) <- data_names[readd]
        pos <- c(readd, pos)
      }
    }
  }
  if (length(data_names[pos]) != 0L) {
    nm_pos <- names(pos)
    if (any(nm_pos == "")) {
      names(pos)[which(nm_pos == "")] <- data_names[pos[which(nm_pos == "")]]
    }
    if (is.null(nm_pos)) {
      names(pos) <- data_names[abs(pos)]
    }
  }
  uniques <- pos[!duplicated(pos)]
  res_nms <- data_names[uniques]
  res <- match(res_nms, data_names)
  if (length(res) != 0L) {
    res <- if (length(setdiff(names(uniques), data_names)) > 0L) {
      if (all(uniques > 0L)) structure(res, .Names = names(uniques)) else structure(res, .Names = res_nms)
    } else {
      structure(res, .Names = res_nms)
    }
  }
  res
}
pm_eval_expr <- function(x) {
  type <- typeof(x)
  switch(
    type,
    "integer" = x,
    "double" = as.integer(x),
    "character" = pm_select_char(x),
    "symbol" = pm_select_symbol(x),
    "language" = pm_eval_call(x),
    stop("Expressions of type <", typeof(x), "> cannot be evaluated for use when subsetting.")
  )
}
pm_select_char <- function(expr) {
  pos <- match(expr, pm_select_env$get_colnames())
  if (any(is.na(pos))) stop("The following columns do not exist:\n    ", paste(expr, collapse = "\n    "))
  pos
}
pm_select_symbol <- function(expr) {
  expr_name <- as.character(expr)
  if (grepl("^is\\.", expr_name) && is.function(expr)) {
    stop(
      "Predicate functions must be wrapped in `where()`.\n\n",
      sprintf("  data %%pm>%% select(where(%s))", expr_name)
    )
  }
  res <- try(pm_select_char(as.character(expr)), silent = TRUE)
  if (inherits(res, "try-error")) {
    res <- tryCatch(
      unlist(lapply(eval(expr, envir = pm_select_env$calling_frame), pm_eval_expr)),
      error = function(e) stop("Column ", expr, " does not exist.")
    )
  }
  res
}
pm_eval_call <- function(x) {
  type <- as.character(x[[1]])
  if (length(type) > 1L) {
    type <- "pm_context"
  }
  switch(
    type,
    `:` = pm_select_seq(x),
    `!` = pm_select_negate(x),
    `-` = pm_select_minus(x),
    `c` = pm_select_c(x),
    `(` = pm_select_bracket(x),
    `&` = pm_select_and(x),
    pm_select_context(x)
  )
}
pm_select_and <- function(expr) {
  exprs <- as.list(expr)[-1]
  res <- do.call(c, lapply(exprs, pm_eval_expr))
  if (all(res > 0) || all(res < 0)) return(unique(res))
  res <- res[!(duplicated(abs(res)) | duplicated(abs(res), fromLast = TRUE))]
  res[res > 0]
}
pm_select_seq <- function(expr) {
  x <- pm_eval_expr(expr[[2]])
  y <- pm_eval_expr(expr[[3]])
  x:y
}
pm_select_negate <- function(expr) {
  x <- if (is_negated_colon(expr)) {
    expr <- call(":", expr[[2]][[2]], expr[[2]][[3]][[2]])
    pm_eval_expr(expr)
  } else {
    pm_eval_expr(expr[[2]])
  }
  x * -1L
}
pm_is_negated_colon <- function(expr) {
  expr[[1]] == "!" && length(expr[[2]]) > 1L && expr[[2]][[1]] == ":" && expr[[2]][[3]][[1]] == "!"
}
pm_select_minus <- function(expr) {
  x <- pm_eval_expr(expr[[2]])
  x * -1L
}
pm_select_c <- function(expr) {
  lst_expr <- as.list(expr)
  lst_expr[[1]] <- NULL
  unlist(lapply(lst_expr, pm_eval_expr))
}
pm_select_bracket <- function(expr) {
  pm_eval_expr(expr[[2]])
}
pm_select_context <- function(expr) {
  eval(expr, envir = pm_select_env$.data)
}
pm_select_env <- new.env()
pm_select_env$setup <- function(.data, calling_frame) {
  pm_select_env$.data <- .data
  pm_select_env$calling_frame <- calling_frame
}
pm_select_env$clean <- function() {
  rm(list = c(".data", "calling_frame"), envir = pm_select_env)
}
pm_select_env$get_colnames <- function() colnames(pm_select_env$.data)
pm_select_env$get_nrow <- function() nrow(pm_select_env$.data)
pm_select_env$get_ncol <- function() ncol(pm_select_env$.data)
pm_eval_select_pos <- function(.data, .cols, .pm_group_pos = FALSE) {
  do.call(pm_select_positions, list(.data = .data, .cols, .pm_group_pos = .pm_group_pos))
}
pm_select <- function(.data, ...) {
  col_pos <- pm_select_positions(.data, ..., .pm_group_pos = TRUE)
  res <- .data[, col_pos, drop = FALSE]
  if (length(names(res)) != 0) colnames(res) <- names(col_pos)
  if (pm_has_groups(.data)) res <- pm_groups_set(res, pm_group_vars(.data))
  res
}
pm_summarise <- function(.data, ..., .groups = NULL) {
  if ("grouped_df" %in% class(.data)) pm_summarise.grouped_df(.data, ..., .groups = NULL) else pm_summarise.data.frame(.data, ..., .groups = NULL)
}
pm_summarise.data.frame <- function(.data, ..., .groups = NULL) {
  fns <- pm_dotdotdot(...)
  pm_context$setup(.data)
  on.exit(pm_context$clean(), add = TRUE)
  groups_exist <- pm_context$is_grouped()
  if (groups_exist) {
    group <- unique(pm_context$get_columns(pm_group_vars(pm_context$.data)))
  }
  if (is_empty_list(fns)) {
    if (groups_exist) return(group) else return(data.frame())
  }
  res <- vector(mode = "list", length = length(fns))
  pm_eval_env <- c(as.list(pm_context$.data), vector(mode = "list", length = length(fns)))
  new_pos <- seq(length(pm_context$.data) + 1L, length(pm_eval_env), 1L)
  for (i in seq_along(fns)) {
    pm_eval_env[[new_pos[i]]] <- do.call(with, list(pm_eval_env, fns[[i]]))
    nms <- if (!is_named(pm_eval_env[[new_pos[i]]])) {
      if (!is.null(names(fns)[[i]])) names(fns)[[i]] else deparse(fns[[i]])
    } else {
      NULL
    }
    if (!is.null(nms)) names(pm_eval_env)[[new_pos[i]]] <- nms
    res[[i]] <- build_data_frame(pm_eval_env[[new_pos[i]]], nms = nms)
  }
  res <- do.call(cbind, res)
  if (groups_exist) res <- cbind(group, res, row.names = NULL)
  res
}
pm_summarise.grouped_df <- function(.data, ..., .groups = NULL) {
  if (!is.null(.groups)) {
    .groups <- match.arg(arg = .groups, choices = c("drop", "drop_last", "keep"), several.ok = FALSE)
  }
  groups <- pm_group_vars(.data)
  res <- pm_apply_grouped_function("pm_summarise", .data, drop = TRUE, ...)
  res <- res[pm_arrange_rows(res, pm_as_symbols(groups)), , drop = FALSE]
  verbose <- pm_summarise_verbose(.groups)
  if (is.null(.groups)) {
    all_one <- as.data.frame(table(res[, groups]))
    all_one <- all_one[all_one$Freq != 0, ]
    .groups <- if (all(all_one$Freq == 1)) "drop_last" else "keep"
  }
  if (.groups == "drop_last") {
    n <- length(groups)
    if (n > 1) {
      if (verbose) pm_summarise_inform(groups[-n])
      res <- pm_groups_set(res, groups[-n], pm_group_by_drop_default(.data))
    }
  } else if (.groups == "keep") {
    if (verbose) pm_summarise_inform(groups)
    res <- pm_groups_set(res, groups, pm_group_by_drop_default(.data))
  } else if (.groups == "drop") {
    attr(res, "groups") <- NULL
  }
  rownames(res) <- NULL
  res
}
pm_summarise_inform <- function(new_groups) {
  message(sprintf(
    "`pm_summarise()` has grouped output by %s. You can override using the `.groups` argument.",
    paste0("'", new_groups, "'", collapse = ", ")
  ))
}
pm_summarise_verbose <- function(.groups) {
  is.null(.groups) &&
    !identical(getOption("poorman.summarise.inform"), FALSE)
}
pm_transmute <- function(.data, ...) {
  if ("grouped_df" %in% class(.data)) pm_transmute.grouped_df(.data, ...) else pm_transmute.data.frame(.data, ...)
}
pm_transmute.data.frame <- function(.data, ...) {
  pm_mutate(.data, ..., .keep = "none")
}
pm_transmute.grouped_df <- function(.data, ...) {
  rows <- rownames(.data)
  res <- pm_apply_grouped_function("pm_transmute", .data, drop = TRUE, ...)
  res[rows, ]
}
pm_ungroup <- function(x, ...) {
  if ("grouped_df" %in% class(x)) pm_ungroup.grouped_df(x, ...) else pm_ungroup.data.frame(x, ...)
}
pm_ungroup.data.frame <- function(x, ...) {
  rm_groups <- pm_deparse_dots(...)
  groups <- pm_group_vars(x)
  if (length(rm_groups) == 0L) rm_groups <- groups
  x <- pm_groups_set(x, groups[!(groups %in% rm_groups)])
  if (length(attr(x, "groups")) == 0L) {
    attr(x, "groups") <- NULL
    class(x) <- class(x)[!(class(x) %in% "grouped_df")]
  }
  x
}
pm_ungroup.grouped_df <- function(x, ...) {
  pm_ungroup.data.frame(...)
}
pm_check_is_dataframe <- function(.data) {
  parent_fn <- all.names(sys.call(-1L), max.names = 1L)
  if (!is.data.frame(.data)) stop(parent_fn, " must be given a data.frame")
  invisible()
}
pm_seq2 <- function(from, to) {
  if (length(from) != 1) stop("`from` must be length one")
  if (length(to) != 1) stop("`to` must be length one")
  if (from > to) integer() else seq.int(from, to)
}
pm_collapse_to_sentence <- function(x) {
  len_x <- length(x)
  if (len_x == 0L) {
    stop("Length of `x` is 0")
  } else if (len_x == 1L) {
    as.character(x)
  } else if (len_x == 2L) {
    paste(x, collapse = " and ")
  } else {
    paste(paste(x[1:(len_x - 1)], collapse = ", "), x[len_x], sep = " and ")
  }
}
pm_build_data_frame <- function(x, nms = NULL) {
  res <- if (is.atomic(x)) {
    data.frame(x)
  } else if (is.list(x) && !is.data.frame(x)) {
    structure(list(x = x), class = "data.frame", row.names = c(NA, -1L))
  } else if (is.data.frame(x)) {
    x
  }
  if (!is.null(nms)) colnames(res) <- nms
  res
}
pm_is_nested <- function(lst) vapply(lst, function(x) inherits(x[1L], "list"), FALSE)
pm_squash <- function(lst) {
  do.call(c, lapply(lst, function(x) if (is.list(x) && !is.data.frame(x)) squash(x) else list(x)))
}
pm_flatten <- function(lst) {
  nested <- pm_is_nested(lst)
  res <- c(lst[!nested], unlist(lst[nested], recursive = FALSE))
  if (sum(nested)) Recall(res) else return(res)
}
pm_where <- function(fn) {
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
