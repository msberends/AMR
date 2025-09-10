# ==================================================================== #
# TITLE:                                                               #
# Antidiskrobial Resistance (AMR) Analysis                              #
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

test_that("test-disk.R", {
  skip_on_cran()

  expect_true(as.disk(8) == as.disk("8"))
  expect_true(is.disk(as.disk(8)))

  expect_equal(suppressWarnings(as.logical(as.disk("INVALID VALUE"))), NA)

  # all levels should be valid disks
  x <- as.disk(c(20, 40))
  expect_inherits(x[1], "disk")
  expect_inherits(x[[1]], "disk")
  expect_inherits(c(x[1], x[9]), "disk")
  expect_inherits(unique(x[1], x[9]), "disk")
  # expect_warning(as.disk("INVALID VALUE"))
  x[2] <- 32
  expect_inherits(x, "disk")

  pdf(NULL) # prevent Rplots.pdf being created
  expect_silent(barplot(as.disk(c(10, 20, 40))))
  expect_silent(plot(as.disk(c(10, 20, 40))))
  expect_silent(plot(as.disk(c(10, 20, 40)), expand = FALSE))
  expect_silent(plot(as.disk(c(10, 20, 40)), mo = "Escherichia coli", ab = "cipr"))
  if (AMR:::pkg_is_available("ggplot2")) {
    expect_inherits(ggplot2::autoplot(as.disk(c(10, 20, 40))), "gg")
    expect_inherits(ggplot2::autoplot(as.disk(c(10, 20, 40)), expand = FALSE), "gg")
    expect_inherits(ggplot2::autoplot(as.disk(c(10, 20, 40)), mo = "Escherichia coli", ab = "cipr"), "gg")
  }
  expect_output(print(as.disk(12)))

  if (AMR:::pkg_is_available("tibble")) {
    expect_output(print(tibble::tibble(d = as.disk(12))))
  }

  # skimr
  if (AMR:::pkg_is_available("skimr", min_version = "2.0.0", also_load = TRUE)) {
    expect_named(
      skim(random_disk(100)),
      c("skim_type", "skim_variable", "n_missing", "complete_rate", "disk.p0", "disk.p25", "disk.p50", "disk.p75", "disk.p100", "disk.hist")
    )
  }
})
