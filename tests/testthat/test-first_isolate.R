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

test_that("test-first_isolate.R", {
  skip_on_cran()

  # all four methods
  expect_equal(
    sum(first_isolate(x = example_isolates, method = "isolate-based", info = TRUE), na.rm = TRUE),
    1984
  )
  expect_equal(
    sum(first_isolate(x = example_isolates, method = "patient-based", info = TRUE), na.rm = TRUE),
    1265
  )
  expect_equal(
    sum(first_isolate(x = example_isolates, method = "episode-based", info = TRUE), na.rm = TRUE),
    1300
  )
  expect_equal(
    sum(first_isolate(x = example_isolates, method = "phenotype-based", info = TRUE), na.rm = TRUE),
    1387
  )

  # for phenotype determination
  expect_equal(
    AMR:::duplicated_antibiogram("SSSS", points_threshold = 2, ignore_I = TRUE, type = "points"),
    FALSE
  )
  expect_equal(
    AMR:::duplicated_antibiogram(c("RRR", "SSS"),
      points_threshold = 2, ignore_I = TRUE, type = "points"
    ),
    c(FALSE, FALSE)
  )
  expect_equal(
    AMR:::duplicated_antibiogram(c("RRR", "RRR", "SSS"),
      points_threshold = 2, ignore_I = TRUE, type = "points"
    ),
    c(FALSE, TRUE, FALSE)
  )
  expect_equal(
    AMR:::duplicated_antibiogram(c("RRR", "RSS", "SSS", "RSS", "RRR", "RRR", "SSS", "RSS", "RSR", "RRR"),
      points_threshold = 2, ignore_I = TRUE, type = "points"
    ),
    c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  )

  # Phenotype-based, using key antimicrobials
  expect_equal(
    sum(first_isolate(
      x = example_isolates,
      method = "phenotype-based",
      type = "keyantimicrobials",
      antifungal = NULL, info = TRUE
    ), na.rm = TRUE),
    1383
  )
  expect_equal(
    sum(first_isolate(
      x = example_isolates,
      method = "phenotype-based",
      type = "keyantimicrobials",
      antifungal = NULL, info = TRUE, ignore_I = FALSE
    ), na.rm = TRUE),
    1397
  )


  # first non-ICU isolates
  expect_true(
    sum(
      first_isolate(example_isolates,
        col_mo = "mo",
        col_date = "date",
        col_patient_id = "patient",
        col_icu = example_isolates$ward == "ICU",
        info = TRUE,
        icu_exclude = TRUE
      ),
      na.rm = TRUE
    ) < 950
  )

  # set 1500 random observations to be of specimen type 'Urine'
  random_rows <- sample(x = 1:2000, size = 1500, replace = FALSE)
  x <- example_isolates
  x$specimen <- "Other"
  x[random_rows, "specimen"] <- "Urine"
  expect_true(
    sum(first_isolate(
      x = x,
      col_date = "date",
      col_patient_id = "patient",
      col_mo = "mo",
      col_specimen = "specimen",
      filter_specimen = "Urine",
      info = TRUE
    ), na.rm = TRUE) < 1400
  )
  # same, but now exclude ICU
  expect_true(
    sum(first_isolate(
      x = x,
      col_date = "date",
      col_patient_id = "patient",
      col_mo = "mo",
      col_specimen = "specimen",
      filter_specimen = "Urine",
      col_icu = x$ward == "ICU",
      icu_exclude = TRUE,
      info = TRUE
    ), na.rm = TRUE) < 1000
  )

  # "No isolates found"
  test_iso <- example_isolates
  test_iso$specimen <- "test"
  expect_message(first_isolate(test_iso,
    "date",
    "patient",
    col_mo = "mo",
    col_specimen = "specimen",
    filter_specimen = "something_unexisting",
    info = TRUE
  ))

  # printing of exclusion message
  expect_message(first_isolate(example_isolates,
    col_date = "date",
    col_mo = "mo",
    col_patient_id = "patient",
    col_testcode = "gender",
    testcodes_exclude = "M",
    info = TRUE
  ))

  # errors
  expect_error(first_isolate("date", "patient", col_mo = "mo"))
  expect_error(first_isolate(example_isolates,
    col_date = "non-existing col",
    col_mo = "mo"
  ))

  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    # if mo is not an mo class, result should be the same
    expect_identical(
      example_isolates %>%
        mutate(mo = as.character(mo)) %>%
        first_isolate(
          col_date = "date",
          col_mo = "mo",
          col_patient_id = "patient",
          info = FALSE
        ),
      example_isolates %>%
        first_isolate(
          col_date = "date",
          col_mo = "mo",
          col_patient_id = "patient",
          info = FALSE
        )
    )

    # support for WHONET
    expect_message(example_isolates %>%
      select(-patient) %>%
      mutate(
        `First name` = "test",
        `Last name` = "test",
        Sex = "Female"
      ) %>%
      first_isolate(info = TRUE))

    # groups
    x <- example_isolates %>%
      group_by(ward) %>%
      mutate(first = first_isolate())
    y <- example_isolates %>%
      group_by(ward) %>%
      mutate(first = first_isolate(.))
    expect_identical(x, y)
  }

  # missing dates should be no problem
  df <- example_isolates
  df[1:100, "date"] <- NA
  expect_equal(
    sum(
      first_isolate(
        x = df,
        col_date = "date",
        col_patient_id = "patient",
        col_mo = "mo",
        info = TRUE
      ),
      na.rm = TRUE
    ),
    1390
  )

  # unknown MOs
  test_unknown <- example_isolates
  test_unknown$mo <- ifelse(test_unknown$mo == "B_ESCHR_COLI", "UNKNOWN", test_unknown$mo)
  expect_equal(
    sum(first_isolate(test_unknown, include_unknown = FALSE)),
    1116
  )
  expect_equal(
    sum(first_isolate(test_unknown, include_unknown = TRUE)),
    1599
  )

  test_unknown$mo <- ifelse(test_unknown$mo == "UNKNOWN", NA, test_unknown$mo)
  expect_equal(
    sum(first_isolate(test_unknown)),
    1116
  )

  # empty sir results
  expect_equal(
    sum(first_isolate(example_isolates, include_untested_sir = FALSE)),
    1374
  )

  # shortcuts
  expect_identical(
    filter_first_isolate(example_isolates),
    subset(example_isolates, first_isolate(example_isolates))
  )


  # notice that all mo's are distinct, so all are TRUE
  expect_true(all(first_isolate(AMR:::pm_distinct(example_isolates, mo, .keep_all = TRUE), info = TRUE) == TRUE))

  # only one isolate, so return fast
  expect_true(first_isolate(data.frame(mo = "Escherichia coli", date = Sys.Date(), patient = "patient"), info = TRUE))
})
