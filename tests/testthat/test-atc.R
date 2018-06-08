context("atc.R")

test_that("atc_property works", {
  expect_equal(tolower(atc_property("J01CA04", property = "Name")), "amoxicillin")
  expect_equivalent(atc_property("J01CA04", "DDD"), 1)
})

test_that("abname works", {
  expect_equal(abname("AMOX"), "Amoxicillin")
  expect_equal(abname(c("AMOX", "GENT")), c("Amoxicillin", "Gentamicin"))
  expect_equal(abname(c("AMOX+GENT")), "Amoxicillin + gentamicin")
  expect_equal(abname("AMOX", from = 'umcg'), "Amoxicillin")
  expect_equal(abname("amox", from = 'molis'), "Amoxicillin")
  expect_equal(abname("J01CA04", from = 'atc'), "Amoxicillin")
})

test_that("guess_atc works", {
  expect_equal(guess_atc(c("J01FA01",
                           "Erythromycin",
                           "eryt",
                           "ERYT",
                           "ERY",
                           "Erythrocin",
                           "Eryzole",
                           "Pediamycin")),
               rep("J01FA01", 8))

})
