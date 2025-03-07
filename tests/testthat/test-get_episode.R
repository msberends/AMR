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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

test_that("get_episode works", {
  x <- data.frame(dates = as.Date(c("2021-01-01", "2021-01-02", "2021-01-05", "2021-01-08", "2021-02-21", "2021-02-22", "2021-02-23", "2021-02-24", "2021-03-01", "2021-03-01")))
  x$absolute <- get_episode(x$dates, episode_days = 7)
  x$relative <- get_episode(x$dates, case_free_days = 7)
  expect_equal(x$absolute, c(1, 1, 1, 2, 3, 3, 3, 3, 4, 4))
  expect_equal(x$relative, c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2))
  expect_equal(get_episode(as.Date(c("2022-01-01", "2020-01-01")), 365), c(2, 1))
  expect_equal(get_episode(as.Date(c("2020-01-01", "2022-01-01")), 365), c(1, 2))

  test_df <- rbind(
    data.frame(
      date = as.Date(c("2015-01-01", "2015-10-01", "2016-02-04", "2016-12-31", "2017-01-01", "2017-02-01", "2017-02-05", "2020-01-01")),
      patient_id = "A"
    ),
    data.frame(
      date = as.Date(c("2015-01-01", "2016-02-01", "2016-12-31", "2017-01-01", "2017-02-03")),
      patient_id = "B"
    )
  )

  expect_equal(
    get_episode(test_df$date, 365),
    c(1, 1, 2, 2, 2, 3, 3, 4, 1, 2, 2, 2, 3)
  )
  expect_equal(
    get_episode(test_df$date[which(test_df$patient_id == "A")], 365),
    c(1, 1, 2, 2, 2, 2, 3, 4)
  )
  expect_equal(
    get_episode(test_df$date[which(test_df$patient_id == "B")], 365),
    c(1, 2, 2, 2, 3)
  )

  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    expect_identical(
      test_df %>% group_by(patient_id) %>% mutate(f = is_new_episode(date, 365)) %>% pull(f),
      c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE)
    )

    suppressMessages(
      x <- example_isolates %>%
        mutate(out = first_isolate(., include_unknown = TRUE, method = "episode-based", info = FALSE))
    )
    y <- example_isolates %>%
      group_by(patient, mo) %>%
      mutate(out = is_new_episode(date, 365))

    expect_identical(which(x$out), which(y$out))
  }
})
