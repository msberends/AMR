context("abname.R")

test_that("abname works", {
  expect_equal(abname("AMOX"), "Amoxicillin")
  expect_equal(abname(c("AMOX", "GENT")), c("Amoxicillin", "Gentamicin"))
  expect_equal(abname(c("AMOX+GENT")), "Amoxicillin + gentamicin")
  expect_equal(abname("AMOX", from = 'umcg'), "Amoxicillin")
  expect_equal(abname("amox", from = 'molis', tolower = TRUE), "amoxicillin")
  expect_equal(abname("J01CA04", from = 'atc'), "Amoxicillin")
  expect_equal(abname(c("amox", "J01CA04", "Trimox", "dispermox", "Amoxil")),
                      rep("Amoxicillin", 5))
  expect_equal(abname("AMOX", to = 'atc'), "J01CA04")

  expect_error(abname("AMOX", to = c(1:3)))
  expect_error(abname("AMOX", to = "test"))
  expect_warning(abname("TEST
       "))
  expect_warning(abname("AMOX or GENT"))
})
