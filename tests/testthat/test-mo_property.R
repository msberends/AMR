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

test_that("test-mo_property.R", {
  skip_on_cran()

  expect_equal(mo_kingdom("Escherichia coli"), "Bacteria")
  expect_equal(mo_kingdom("Escherichia coli"), mo_domain("Escherichia coli"))
  expect_equal(mo_phylum("Escherichia coli"), "Pseudomonadota")
  expect_equal(mo_class("Escherichia coli"), "Gammaproteobacteria")
  expect_equal(mo_order("Escherichia coli"), "Enterobacterales")
  expect_equal(mo_family("Escherichia coli"), "Enterobacteriaceae")
  expect_equal(mo_fullname("Escherichia coli"), "Escherichia coli")
  expect_equal(mo_genus("Escherichia coli"), "Escherichia")
  expect_equal(mo_name("Escherichia coli"), "Escherichia coli")
  expect_equal(mo_shortname("Escherichia coli"), "E. coli")
  expect_equal(mo_shortname("Escherichia"), "Escherichia")
  expect_equal(mo_shortname("Staphylococcus aureus"), "S. aureus")
  expect_equal(mo_shortname("Staphylococcus aureus", Becker = TRUE), "S. aureus")
  expect_equal(mo_shortname("Staphylococcus aureus", Becker = "all", language = "en"), "CoPS")
  expect_equal(mo_shortname("Streptococcus agalactiae"), "S. agalactiae")
  expect_equal(mo_shortname("Streptococcus agalactiae", Lancefield = TRUE), "GBS")

  # check gram stain determination, to prevent we lag after a taxonomic renaming
  current_grampos_phyla <- c(
    "Actinomycetota", # since 2021, old name was Actinobacteria
    "Chloroflexota", # since 2021, old name was Chloroflexi
    "Bacillota", # since 2021, old name was Firmicutes
    "Mycoplasmatota" # since 2021, old name was Tenericutes
  )
  expect_true(all(current_grampos_phyla %in% microorganisms$phylum, na.rm = TRUE))
  current_grampos_classes <- c(
    "",
    "Acidimicrobiia",
    "Actinomycetes",
    "Anaerolineae",
    "Ardenticatenia",
    "Bacilli",
    "Caldilineae",
    "Chloroflexia",
    "Clostridia",
    "Coriobacteriia",
    "Culicoidibacteria",
    "Dehalococcoidia",
    "Erysipelotrichia",
    "Ktedonobacteria",
    "Limnochordia",
    "Limnocylindria",
    "Mollicutes",
    "Negativicutes",
    "Nitriliruptoria",
    "Rubrobacteria",
    "Tepidiformia",
    "Thermoflexia",
    "Thermoleophilia",
    "Thermolithobacteria"
  )
  expect_identical(
    sort(unique(microorganisms[which(microorganisms$phylum %in% current_grampos_phyla), "class", drop = TRUE])),
    current_grampos_classes
  )

  expect_equal(mo_species("Escherichia coli"), "coli")
  expect_equal(mo_subspecies("Escherichia coli"), "")
  expect_equal(mo_type("Escherichia coli", language = "en"), "Bacteria")
  expect_equal(mo_gramstain("Escherichia coli", language = "en"), "Gram-negative")
  expect_inherits(mo_taxonomy("Escherichia coli"), "list")
  expect_equal(names(mo_taxonomy("Escherichia coli")), c(
    "kingdom", "phylum", "class", "order",
    "family", "genus", "species", "subspecies"
  ))
  expect_equal(mo_synonyms("Escherichia coli"), NULL)
  expect_true(length(mo_synonyms("Candida albicans")) > 1)
  expect_inherits(mo_synonyms(c("Candida albicans", "Escherichia coli")), "list")
  expect_equal(names(mo_info("Escherichia coli")), c(
    "mo", "rank",
    "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies",
    "status", "synonyms", "gramstain", "oxygen_tolerance",
    "url", "ref", "snomed", "lpsn", "mycobank", "gbif", "group_members"
  ))
  expect_inherits(mo_info(c("Escherichia coli", "Staphylococcus aureus")), "list")
  expect_true(length(mo_group_members("B_HACEK")) > 1)
  expect_inherits(mo_group_members(c("Candida albicans", "Escherichia coli")), "list")

  expect_identical(
    mo_oxygen_tolerance(c("Klebsiella pneumoniae", "Clostridioides difficile")),
    c("facultative anaerobe", "anaerobe")
  )

  expect_equal(
    as.character(table(mo_pathogenicity(example_isolates$mo))),
    c("1911", "72", "1", "16")
  )

  expect_equal(mo_ref("Escherichia coli"), "Castellani et al., 1919")
  expect_equal(mo_authors("Escherichia coli"), "Castellani et al.")
  expect_equal(mo_year("Escherichia coli"), 1919)

  expect_true(mo_url("Amoeba dysenteriae") %like% "gbif.org")
  expect_true(mo_url("Candida albicans") %like% "mycobank.org")
  expect_true(mo_url("Escherichia coli") %like% "lpsn.dsmz.de")

  # test integrity of getting back full names
  expect_identical(
    microorganisms$fullname[microorganisms$fullname %unlike% "(Fungi|{)"],
    suppressWarnings(mo_fullname(microorganisms$fullname[microorganisms$fullname %unlike% "(Fungi|{)"], language = "en", keep_synonyms = TRUE))
  )

  # check languages
  expect_equal(mo_type("Escherichia coli", language = "de"), "Bakterien")
  expect_equal(mo_gramstain("Escherichia coli", language = "nl"), "Gram-negatief")

  gr <- mo_gramstain("Escherichia coli", language = NULL)
  for (l in AMR:::LANGUAGES_SUPPORTED[-1]) {
    expect_false(mo_gramstain("Escherichia coli", language = l) == gr, info = paste("Gram-stain in language", l))
  }

  # test languages
  expect_error(mo_gramstain("Escherichia coli", language = "UNKNOWN"))
  fullnames <- microorganisms$fullname[which(microorganisms$fullname %unlike% "unknown|coagulase|Fungi|[(]class[)]|[{]")]
  to_dutch <- suppressWarnings(mo_name(fullnames, language = "nl", keep_synonyms = TRUE))
  back_to_english <- suppressWarnings(mo_name(to_dutch, language = NULL, keep_synonyms = TRUE))
  diffs <- paste0('"', fullnames[fullnames != back_to_english], '"', collapse = ", ")
  expect_identical(fullnames, back_to_english, info = diffs) # gigantic test - will run ALL names


  # manual property function
  expect_error(mo_property("Escherichia coli", property = c("genus", "fullname")))
  expect_error(mo_property("Escherichia coli", property = "UNKNOWN"))
  expect_identical(
    mo_property("Escherichia coli", property = "fullname"),
    mo_fullname("Escherichia coli")
  )
  expect_identical(
    mo_property("Escherichia coli", property = "genus"),
    mo_genus("Escherichia coli")
  )
  expect_identical(
    mo_property("Escherichia coli", property = "species"),
    mo_species("Escherichia coli")
  )
  expect_identical(
    mo_property("Escherichia coli", property = "lpsn"),
    mo_lpsn("Escherichia coli")
  )
  expect_identical(
    mo_property("Escherichia coli", property = "gbif"),
    mo_gbif("Escherichia coli")
  )
  expect_identical(
    mo_property("Absidia abundans", property = "mycobank"),
    mo_mycobank("Absidia abundans")
  )

  expect_true("Escherichia blattae" %in% mo_synonyms("Shimwellia blattae"))
  expect_true(is.list(mo_synonyms(rep("Shimwellia blattae", 2))))
  expect_identical(
    mo_current(c("Escherichia blattae", "Escherichia coli")),
    c("Shimwellia blattae", "Escherichia coli")
  )

  expect_identical(mo_ref("Chlamydia psittaci"), "Garcia-Lopez et al., 2019")
  expect_identical(mo_ref("Chlamydophila psittaci", keep_synonyms = TRUE), "Everett et al., 1999")

  expect_true(112283007 %in% mo_snomed("Escherichia coli")[[1]])
  # outcome of mo_fullname must always return the fullname from the data set
  x <- data.frame(
    mo = microorganisms$mo,
    # fullname from the original data:
    f1 = microorganisms$fullname,
    # newly created fullname based on MO code:
    f2 = mo_fullname(microorganisms$mo, language = "en", keep_synonyms = TRUE),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(subset(x, f1 != f2)), 0)
  # is gram pos/neg (also return FALSE for all non-bacteria)
  expect_equal(
    mo_is_gram_negative(c("Escherichia coli", "Staphylococcus aureus", "Candida albicans")),
    c(TRUE, FALSE, FALSE)
  )
  expect_equal(
    mo_is_gram_positive(c("Escherichia coli", "Staphylococcus aureus", "Candida albicans")),
    c(FALSE, TRUE, FALSE)
  )
  expect_equal(
    mo_is_yeast(c("Candida", "Trichophyton", "Klebsiella")),
    c(TRUE, FALSE, FALSE)
  )
  # is intrinsic resistant
  expect_equal(
    mo_is_intrinsic_resistant(
      c("Escherichia coli", "Staphylococcus aureus", "Candida albicans"),
      "vanco"
    ),
    c(TRUE, FALSE, FALSE)
  )
  # with reference data
  expect_equal(
    mo_name("test", reference_df = data.frame(col1 = "test", mo = "B_ESCHR_COLI")),
    "Escherichia coli"
  )
  if (AMR:::pkg_is_available("dplyr", min_version = "1.0.0", also_load = TRUE)) {
    expect_equal(example_isolates %>% filter(mo_is_gram_negative()) %>% nrow(),
      730,
      tolerance = 0.5
    )
    expect_equal(example_isolates %>% filter(mo_is_gram_positive()) %>% nrow(),
      1238,
      tolerance = 0.5
    )
    expect_equal(example_isolates %>% filter(mo_is_intrinsic_resistant(ab = "Vancomycin")) %>% nrow(),
      710,
      tolerance = 0.5
    )
  }
})
