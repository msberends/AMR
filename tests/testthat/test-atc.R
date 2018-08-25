context("atc.R")

test_that("atc_property works", {
  skip_on_cran() # relies on internet connection of server, don't test
  skip_on_appveyor() # security error on AppVeyor

  if (!is.null(curl::nslookup("www.whocc.no", error = FALSE))) {
    expect_equal(tolower(atc_property("J01CA04", property = "Name")), "amoxicillin")
    expect_equal(atc_property("J01CA04", property = "unit"), "g")
    expect_equal(atc_property("J01CA04", property = "DDD"),
                 atc_ddd("J01CA04"))

    expect_identical(atc_property("J01CA04", property = "Groups"),
                     atc_groups("J01CA04"))

    expect_warning(atc_property("ABCDEFG", property = "DDD"))

    expect_error(atc_property("J01CA04", property = c(1:5)))
    expect_error(atc_property("J01CA04", property = "test"))
    expect_error(atc_property("J01CA04", property = "test", administration = c(1:5)))
  }
})

test_that("guess_atc works", {
  expect_equal(as.character(guess_atc(c("J01FA01",
                           "Erythromycin",
                           "eryt",
                           "ERYT",
                           "ERY",
                           "Erythrocin",
                           "Eryzole",
                           "Pediamycin"))),
               rep("J01FA01", 8))

  expect_identical(class(as.atc("amox")), "atc")


})

test_that("atc.property works", {
  expect_equal(atc.certe("J01CA04"), "amox")
  expect_equal(atc.umcg("J01CA04"), "AMOX")
  expect_equal(atc.official("J01CA04"), "Amoxicillin")
  expect_equal(atc.official_nl("J01CA04"), "Amoxicilline")
  expect_equal(atc.trivial_nl("J01CA04"), "Amoxicilline")
})
