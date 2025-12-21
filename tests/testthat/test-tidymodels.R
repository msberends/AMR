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

test_that("tidymodels.R", {
  skip_on_cran()

  if (AMR:::pkg_is_available("recipes", also_load = TRUE) && AMR:::pkg_is_available("dplyr", also_load = TRUE)) {
    # SIR
    df <- tibble(
      sir1 = as.sir(c("S", "I", "R", "S", "R")),
      sir2 = as.sir(c("I", "R", "S", "R", "I")),
      not_sir = c("S", "R", "R", "S", "I")
    )
    rec <- recipe(~., data = df) %>% step_sir_numeric(all_sir())
    prepped <- prep(rec)
    baked <- bake(prepped, new_data = df)
    expect_inherits(baked$sir1, "numeric")
    expect_inherits(baked$sir2, "numeric")
    expect_equal(baked$not_sir, as.factor(df$not_sir))

    # MIC
    df <- tibble(
      mic_col1 = as.mic(c("<=0.002", "0.002", "0.004", "0.016", "32")),
      mic_col2 = as.mic(c("0.5", "1", "2", "4", "8")),
      non_mic = c(1, 2, 3, 4, 5)
    )
    rec <- recipe(~., data = df) %>% step_mic_log2(all_mic())
    prepped <- prep(rec)
    baked <- bake(prepped, new_data = df)
    expect_inherits(baked$mic_col1, "numeric")
    expect_inherits(baked$mic_col2, "numeric")
    expect_equal(baked$non_mic, df$non_mic)
    expect_equal(baked$mic_col2, log2(as.numeric(df$mic_col2)))

    # disk
    df <- tibble(
      disk_col = as.disk(c(21, 22, 23, 24, 25)),
      non_disk = c(21, 22, 23, 24, 25)
    )
    rec <- recipe(~., data = df) %>% step_rm(!all_disk())
    prepped <- prep(rec)
    baked <- bake(prepped, new_data = df)
    expect_inherits(baked$disk_col, "disk")

    # steps check
    df <- tibble(x = as.mic(c("1", "2", "4")))
    rec <- recipe(~x, data = df) %>% step_mic_log2(all_mic())
    prepped <- prep(rec)
    tidy_df <- tidy(prepped, number = 1)
    expect_equal(unname(tidy_df$terms), "x")

    df <- tibble(x = as.sir(c("S", "I", "R")))
    rec <- recipe(~x, data = df) %>% step_sir_numeric(all_sir())
    prepped <- prep(rec)
    tidy_df <- tidy(prepped, number = 1)
    expect_equal(unname(tidy_df$terms), "x")
  }
})
