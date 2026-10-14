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
# Center Groningen in the Netherlands, in collaboration with many      #
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

test_that("test-mo_source.R", {
  # set_mo_source() and the interactive branches of get_mo_source() only work in interactive sessions, see .github/scripts/run_tests.sh
  skip_if(!interactive())
  skip_on_cran()

  dir <- tempfile("mo_source_test")
  dir.create(dir)
  src_a <- file.path(dir, "source_a.csv")
  src_b <- file.path(dir, "source_b.csv")
  dest_a <- file.path(dir, "dest_a.rds")
  dest_b <- file.path(dir, "dest_b.rds")
  write.csv(data.frame(lab_code = "TESTCODE_A", mo = "B_ESCHR_COLI"), src_a, row.names = FALSE)
  write.csv(data.frame(lab_code = "TESTCODE_B", mo = "B_STPHY_AURS"), src_b, row.names = FALSE)
  # existing destination files are updated without asking for permission, so these tests never prompt
  saveRDS(NULL, dest_a)
  saveRDS(NULL, dest_b)
  amr_env <- AMR:::AMR_env
  mo_source_bak <- amr_env$mo_source
  on.exit(amr_env$mo_source <- mo_source_bak, add = TRUE)

  suppressMessages(set_mo_source(src_b, destination = dest_b))
  suppressMessages(set_mo_source(src_a, destination = dest_a))
  expect_identical(get_mo_source(destination = dest_a)$lab_code, "TESTCODE_A")
  # another destination returns its own mo_source, not the one already loaded
  expect_identical(get_mo_source(destination = dest_b)$lab_code, "TESTCODE_B")
  expect_identical(as.character(as.mo("TESTCODE_B", reference_df = get_mo_source(destination = dest_b), info = FALSE)), "B_STPHY_AURS")

  # a destination that does not exist removes no other mo_source file, nor the loaded mo_source or the cache
  suppressMessages(get_mo_source(destination = dest_a))
  suppressWarnings(as.mo("Klebsiela pneumonia", info = FALSE))
  # as.mo() uses the mo_source at the default destination, which can exist on the system that runs these tests
  loaded_destination <- attributes(amr_env$mo_source)$mo_source_destination
  n_cached <- NROW(amr_env$mo_previously_coerced)
  expect_null(get_mo_source(destination = file.path(dir, "does_not_exist.rds")))
  expect_true(file.exists(dest_a))
  expect_true(file.exists(dest_b))
  expect_identical(attributes(amr_env$mo_source)$mo_source_destination, loaded_destination)
  expect_identical(NROW(amr_env$mo_previously_coerced), n_cached)

  # an updated original file updates the mo_source file at the same destination
  write.csv(data.frame(lab_code = "TESTCODE_A2", mo = "B_ESCHR_COLI"), src_a, row.names = FALSE)
  Sys.setFileTime(src_a, Sys.time() + 60)
  expect_identical(suppressMessages(get_mo_source(destination = dest_a))$lab_code, "TESTCODE_A2")
  expect_true(file.exists(dest_b))

  # a deleted mo_source file removes the reference to it
  unlink(dest_a)
  expect_null(get_mo_source(destination = dest_a))
  expect_null(amr_env$mo_source)
  expect_true(file.exists(dest_b))

  # a moved or deleted original file removes the mo_source file, as documented, without an error
  suppressMessages(get_mo_source(destination = dest_b))
  unlink(src_b)
  expect_null(suppressMessages(get_mo_source(destination = dest_b)))
  expect_false(file.exists(dest_b))

  unlink(dir, recursive = TRUE)
})
