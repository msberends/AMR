context("like.R")

test_that("`like` works", {
  expect_true(suppressWarnings("test" %like% c("^t", "^s")))
  expect_true("test" %like% "test")
  expect_true("test" %like% "TEST")
  expect_true(as.factor("test") %like% "TEST")
  expect_identical(factor(c("Test case", "Something different", "Yet another thing")) %like% c("case", "diff", "yet"),
                   c(TRUE, TRUE, TRUE))
})
