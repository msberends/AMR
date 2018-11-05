context("get_locale.R")

test_that("get_locale works", {
  expect_identical(mo_genus("B_GRAMP", language = "pt"),
                   "(Gram positivos desconhecidos)")

  expect_identical(mo_fullname("CoNS", "en"), "Coagulase Negative Staphylococcus (CoNS)")
  expect_identical(mo_fullname("CoNS", "de"), "Koagulase-negative Staphylococcus (KNS)")
  expect_identical(mo_fullname("CoNS", "nl"), "Coagulase-negatieve Staphylococcus (CNS)")
  expect_identical(mo_fullname("CoNS", "es"), "Staphylococcus coagulasa negativo (CoNS)")
  expect_identical(mo_fullname("CoNS", "it"), "Staphylococcus negativo coagulasi (CoNS)")
  # expect_identical(mo_fullname("CoNS", "fr"), "Staphylococcus \u00e0 coagulase n\u00e9gative (CoNS)")
  expect_identical(mo_fullname("CoNS", "pt"), "Staphylococcus coagulase negativo (CoNS)")

})
