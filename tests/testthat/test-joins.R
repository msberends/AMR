context("joins.R")

test_that("joins work", {
  unjoined <- septic_patients
  inner <- septic_patients %>% inner_join_microorganisms()
  left <- septic_patients %>% left_join_microorganisms()
  suppressWarnings(right <- septic_patients %>% right_join_microorganisms())
  suppressWarnings(full <- septic_patients %>% full_join_microorganisms())
  
  expect_true(ncol(unjoined) < ncol(inner))
  expect_true(nrow(unjoined) == nrow(inner))
  
  expect_true(ncol(unjoined) < ncol(left))
  expect_true(nrow(unjoined) == nrow(left))
  
  expect_true(nrow(unjoined) < nrow(right))
  expect_true(nrow(unjoined) < nrow(full))
})
