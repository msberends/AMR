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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Determine Clinical or Epidemic Episodes
#'
#' These functions determine which items in a vector can be considered (the start of) a new episode. This can be used to determine clinical episodes for any epidemiological analysis. The [get_episode()] function returns the index number of the episode per group, while the [is_new_episode()] function returns `TRUE` for every new [get_episode()] index. Both absolute and relative episode determination are supported.
#' @param x vector of dates (class `Date` or `POSIXt`), will be sorted internally to determine episodes
#' @param episode_days episode length in days to specify the time period after which a new episode begins, can also be less than a day or `Inf`, see *Details*
#' @param case_free_days (inter-epidemic) interval length in days after which a new episode will start, can also be less than a day or `Inf`, see *Details*
#' @param ... ignored, only in place to allow future extensions
#' @details Episodes can be determined in two ways: absolute and relative.
#'
#' 1. Absolute
#'
#'    This method uses `episode_days` to define an episode length in days, after which a new episode will start. A common use case in AMR data analysis is microbial epidemiology: episodes of *S. aureus* bacteraemia in ICU patients for example. The episode length could then be 30 days, so that new *S. aureus* isolates after an ICU episode of 30 days will be considered a different (or new) episode.
#'
#'    Thus, this method counts **since the start of the previous episode**.
#'
#' 2. Relative
#'
#'    This method uses `case_free_days` to quantify the duration of case-free days (the inter-epidemic interval), after which a new episode will start. A common use case is infectious disease epidemiology: episodes of norovirus outbreaks in a hospital for example. The case-free period could then be 14 days, so that new norovirus cases after that time will be considered a different (or new) episode.
#'
#'    Thus, this methods counts **since the last case in the previous episode**.
#'
#' In a table:
#'
#' |    Date    | Using `episode_days = 7` | Using `case_free_days = 7` |
#' |:----------:|:------------------------:|:--------------------------:|
#' | 2023-01-01 |             1            |              1             |
#' | 2023-01-02 |             1            |              1             |
#' | 2023-01-05 |             1            |              1             |
#' | 2023-01-08 |             2**          |              1             |
#' | 2023-02-21 |             3            |              2***          |
#' | 2023-02-22 |             3            |              2             |
#' | 2023-02-23 |             3            |              2             |
#' | 2023-02-24 |             3            |              2             |
#' | 2023-03-01 |             4            |              2             |
#'
#' ** This marks the start of a new episode, because 8 January 2023 is more than 7 days since the start of the previous episode (1 January 2023). \cr
#' *** This marks the start of a new episode, because 21 January 2023 is more than 7 days since the last case in the previous episode (8 January 2023).
#'
#' Either `episode_days` or `case_free_days` must be provided in the function.
#'
#' ### Difference between `get_episode()` and `is_new_episode()`
#'
#' The [get_episode()] function returns the index number of the episode, so all cases/patients/isolates in the first episode will have the number 1, all cases/patients/isolates in the second episode will have the number 2, etc.
#'
#' The [is_new_episode()] function on the other hand, returns `TRUE` for every new [get_episode()] index.
#'
#' To specify, when setting `episode_days = 365` (using method 1 as explained above), this is how the two functions differ:
#'
#' | patient   | date       | `get_episode()` | `is_new_episode()` |
#' |:---------:|:----------:|:---------------:|:------------------:|
#' | A         | 2019-01-01 |               1 | TRUE               |
#' | A         | 2019-03-01 |               1 | FALSE              |
#' | A         | 2021-01-01 |               2 | TRUE               |
#' | B         | 2008-01-01 |               1 | TRUE               |
#' | B         | 2008-01-01 |               1 | FALSE              |
#' | C         | 2020-01-01 |               1 | TRUE               |
#'
#' ### Other
#'
#' The [first_isolate()] function is a wrapper around the [is_new_episode()] function, but is more efficient for data sets containing microorganism codes or names and allows for different isolate selection methods.
#'
#' The `dplyr` package is not required for these functions to work, but these episode functions do support [variable grouping][dplyr::group_by()] and work conveniently inside `dplyr` verbs such as [`filter()`][dplyr::filter()], [`mutate()`][dplyr::mutate()] and [`summarise()`][dplyr::summarise()].
#' @return
#' * [get_episode()]: an [integer] vector
#' * [is_new_episode()]: a [logical] vector
#' @seealso [first_isolate()]
#' @rdname get_episode
#' @export
#' @examples
#' # difference between absolute and relative determination of episodes:
#' x <- data.frame(dates = as.Date(c(
#'   "2021-01-01",
#'   "2021-01-02",
#'   "2021-01-05",
#'   "2021-01-08",
#'   "2021-02-21",
#'   "2021-02-22",
#'   "2021-02-23",
#'   "2021-02-24",
#'   "2021-03-01",
#'   "2021-03-01"
#' )))
#' x$absolute <- get_episode(x$dates, episode_days = 7)
#' x$relative <- get_episode(x$dates, case_free_days = 7)
#' x
#'
#'
#' # `example_isolates` is a data set available in the AMR package.
#' # See ?example_isolates
#' df <- example_isolates[sample(seq_len(2000), size = 100), ]
#'
#' get_episode(df$date, episode_days = 60) # indices
#' is_new_episode(df$date, episode_days = 60) # TRUE/FALSE
#'
#' # filter on results from the third 60-day episode only, using base R
#' df[which(get_episode(df$date, 60) == 3), ]
#'
#' # the functions also work for less than a day, e.g. to include one per hour:
#' get_episode(
#'   c(
#'     Sys.time(),
#'     Sys.time() + 60 * 60
#'   ),
#'   episode_days = 1 / 24
#' )
#'
#' \donttest{
#' if (require("dplyr")) {
#'   # is_new_episode() can also be used in dplyr verbs to determine patient
#'   # episodes based on any (combination of) grouping variables:
#'   df %>%
#'     mutate(condition = sample(
#'       x = c("A", "B", "C"),
#'       size = 100,
#'       replace = TRUE
#'     )) %>%
#'     group_by(patient, condition) %>%
#'     mutate(new_episode = is_new_episode(date, 365)) %>%
#'     select(patient, date, condition, new_episode) %>%
#'     arrange(patient, condition, date)
#' }
#'
#' if (require("dplyr")) {
#'   df %>%
#'     group_by(ward, patient) %>%
#'     transmute(date,
#'       patient,
#'       new_index = get_episode(date, 60),
#'       new_logical = is_new_episode(date, 60)
#'     ) %>%
#'     arrange(patient, ward, date)
#' }
#'
#' if (require("dplyr")) {
#'   df %>%
#'     group_by(ward) %>%
#'     summarise(
#'       n_patients = n_distinct(patient),
#'       n_episodes_365 = sum(is_new_episode(date, episode_days = 365)),
#'       n_episodes_60 = sum(is_new_episode(date, episode_days = 60)),
#'       n_episodes_30 = sum(is_new_episode(date, episode_days = 30))
#'     )
#' }
#'
#' # grouping on patients and microorganisms leads to the same
#' # results as first_isolate() when using 'episode-based':
#' if (require("dplyr")) {
#'   x <- df %>%
#'     filter_first_isolate(
#'       include_unknown = TRUE,
#'       method = "episode-based"
#'     )
#'
#'   y <- df %>%
#'     group_by(patient, mo) %>%
#'     filter(is_new_episode(date, 365)) %>%
#'     ungroup()
#'
#'   identical(x, y)
#' }
#'
#' # but is_new_episode() has a lot more flexibility than first_isolate(),
#' # since you can now group on anything that seems relevant:
#' if (require("dplyr")) {
#'   df %>%
#'     group_by(patient, mo, ward) %>%
#'     mutate(flag_episode = is_new_episode(date, 365)) %>%
#'     select(group_vars(.), flag_episode)
#' }
#' }
get_episode <- function(x, episode_days = NULL, case_free_days = NULL, ...) {
  meet_criteria(x, allow_class = c("Date", "POSIXt"), allow_NA = TRUE)
  meet_criteria(episode_days, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = FALSE, allow_NULL = TRUE)
  meet_criteria(case_free_days, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = FALSE, allow_NULL = TRUE)
  as.integer(exec_episode(x, episode_days, case_free_days, ...))
}

#' @rdname get_episode
#' @export
is_new_episode <- function(x, episode_days = NULL, case_free_days = NULL, ...) {
  meet_criteria(x, allow_class = c("Date", "POSIXt"), allow_NA = TRUE)
  meet_criteria(episode_days, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = FALSE, allow_NULL = TRUE)
  meet_criteria(case_free_days, allow_class = c("numeric", "integer"), has_length = 1, is_positive = TRUE, is_finite = FALSE, allow_NULL = TRUE)
  !duplicated(exec_episode(x, episode_days, case_free_days, ...))
}

exec_episode <- function(x, episode_days, case_free_days, ...) {
  stop_ifnot(is.null(episode_days) || is.null(case_free_days),
    "either argument `episode_days` or argument `case_free_days` must be set.",
    call = -2
  )

  # running as.double() on a POSIXct object will return its number of seconds since 1970-01-01
  x <- as.double(as.POSIXct(x)) # as.POSIXct() required for Date classes

  # since x is now in seconds, get seconds from episode_days as well
  episode_seconds <- episode_days * 60 * 60 * 24
  case_free_seconds <- case_free_days * 60 * 60 * 24

  if (length(x) == 1) { # this will also match 1 NA, which is fine
    return(1)
  } else if (length(x) == 2 && all(!is.na(x))) {
    if ((length(episode_seconds) > 0 && (max(x) - min(x)) >= episode_seconds) ||
      (length(case_free_seconds) > 0 && (max(x) - min(x)) >= case_free_seconds)) {
      if (x[1] <= x[2]) {
        return(c(1, 2))
      } else {
        return(c(2, 1))
      }
    } else {
      return(c(1, 1))
    }
  }

  run_episodes <- function(x, episode_seconds, case_free) {
    NAs <- which(is.na(x))
    x[NAs] <- 0

    indices <- integer(length = length(x))
    start <- x[1]
    ind <- 1
    indices[ind] <- 1
    for (i in 2:length(x)) {
      if ((length(episode_seconds) > 0 && (x[i] - start) >= episode_seconds) ||
        (length(case_free_seconds) > 0 && (x[i] - x[i - 1]) >= case_free_seconds)) {
        ind <- ind + 1
        start <- x[i]
      }
      indices[i] <- ind
    }
    indices[NAs] <- NA
    indices
  }

  ord <- order(x)
  out <- run_episodes(x[ord], episode_seconds, case_free_seconds)[order(ord)]
  out[is.na(x) & ord != 1] <- NA # every NA expect for the first must remain NA
  out
}
