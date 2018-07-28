context("atc.R")

test_that("atc_property works", {
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

test_that("abname works", {
  expect_equal(abname("AMOX"), "Amoxicillin")
  expect_equal(abname(c("AMOX", "GENT")), c("Amoxicillin", "Gentamicin"))
  expect_equal(abname(c("AMOX+GENT")), "Amoxicillin + gentamicin")
  expect_equal(abname("AMOX", from = 'umcg'), "Amoxicillin")
  expect_equal(abname("amox", from = 'molis', tolower = TRUE), "amoxicillin")
  expect_equal(abname("J01CA04", from = 'atc'), "Amoxicillin")
  expect_equal(abname("AMOX", to = 'atc'), "J01CA04")
  expect_equal(abname("AMOX en GENT"), "Amoxicillin + gentamicin")
  expect_error(abname("AMOX", to = c(1:3)))
  expect_error(abname("AMOX", to = "test"))
  expect_warning(abname("TEST
       "))
  expect_warning(abname("AMOX or GENT"))
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
