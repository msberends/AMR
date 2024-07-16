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
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #


# Traditional antibiogram ----------------------------------------------

ab1 <- antibiogram(example_isolates,
                   antibiotics = c(aminoglycosides(), carbapenems()))

ab2 <- antibiogram(example_isolates,
                   antibiotics = aminoglycosides(),
                   ab_transform = "atc",
                   mo_transform = "gramstain")

ab3 <- antibiogram(example_isolates,
                   antibiotics = carbapenems(),
                   ab_transform = "name",
                   mo_transform = "name")

expect_inherits(ab1, "antibiogram")
expect_inherits(ab2, "antibiogram")
expect_inherits(ab3, "antibiogram")
expect_equal(colnames(ab1), c("Pathogen (N min-max)", "AMK", "GEN", "IPM", "KAN", "MEM", "TOB"))
expect_equal(colnames(ab2), c("Pathogen (N min-max)", "J01GB01", "J01GB03", "J01GB04", "J01GB06"))
expect_equal(colnames(ab3), c("Pathogen (N min-max)", "Imipenem", "Meropenem"))
expect_equal(ab3$Meropenem, c(52, NA, 100, 100, NA))

# Combined antibiogram -------------------------------------------------

# combined antibiotics yield higher empiric coverage
ab4 <- antibiogram(example_isolates,
                   antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"),
                   mo_transform = "gramstain")

ab5 <- antibiogram(example_isolates,
                   antibiotics = c("TZP", "TZP+TOB"),
                   mo_transform = "gramstain",
                   ab_transform = "name",
                   sep = " & ",
                   add_total_n = FALSE)

expect_inherits(ab4, "antibiogram")
expect_inherits(ab5, "antibiogram")
expect_equal(colnames(ab4), c("Pathogen (N min-max)", "TZP", "TZP + GEN", "TZP + TOB"))
expect_equal(colnames(ab5), c("Pathogen", "Piperacillin/tazobactam", "Piperacillin/tazobactam & Tobramycin"))

# Syndromic antibiogram ------------------------------------------------

# the data set could contain a filter for e.g. respiratory specimens
ab6 <- antibiogram(example_isolates,
                   antibiotics = c(aminoglycosides(), carbapenems()),
                   syndromic_group = "ward")

# with a custom language, though this will be determined automatically
# (i.e., this table will be in Dutch on Dutch systems)
ex1 <- example_isolates[which(mo_genus() == "Escherichia"), ]
ab7 <- antibiogram(ex1,
                   antibiotics = aminoglycosides(),
                   ab_transform = "name",
                   syndromic_group = ifelse(ex1$ward == "ICU",
                                            "IC", "Geen IC"),
                   language = "nl")

expect_inherits(ab6, "antibiogram")
expect_inherits(ab7, "antibiogram")
expect_equal(colnames(ab6), c("Syndromic Group", "Pathogen (N min-max)", "AMK", "GEN", "IPM", "KAN", "MEM", "TOB"))
expect_equal(colnames(ab7), c("Syndroomgroep", "Pathogeen (N min-max)", "Amikacine", "Gentamicine", "Tobramycine"))

# Weighted-incidence syndromic combination antibiogram (WISCA) ---------

# the data set could contain a filter for e.g. respiratory specimens
ab8 <- antibiogram(example_isolates,
                   antibiotics = c("AMC", "AMC+CIP", "TZP", "TZP+TOB"),
                   mo_transform = "gramstain",
                   minimum = 10, # this should be >= 30, but now just as example
                   syndromic_group = ifelse(example_isolates$age >= 65 &
                                              example_isolates$gender == "M",
                                            "WISCA Group 1", "WISCA Group 2"))

expect_inherits(ab8, "antibiogram")
expect_equal(colnames(ab8), c("Syndromic Group", "Pathogen (N min-max)", "AMC", "AMC + CIP", "TZP", "TZP + TOB"))

# Generate plots with ggplot2 or base R --------------------------------

pdf(NULL) # prevent Rplots.pdf being created

expect_silent(plot(ab1))
expect_silent(plot(ab2))
expect_silent(plot(ab3))
expect_silent(plot(ab4))
expect_silent(plot(ab5))
expect_silent(plot(ab6))
expect_silent(plot(ab7))
expect_silent(plot(ab8))

if (AMR:::pkg_is_available("ggplot2")) {
  expect_inherits(ggplot2::autoplot(ab1), "gg")
  expect_inherits(ggplot2::autoplot(ab2), "gg")
  expect_inherits(ggplot2::autoplot(ab3), "gg")
  expect_inherits(ggplot2::autoplot(ab4), "gg")
  expect_inherits(ggplot2::autoplot(ab5), "gg")
  expect_inherits(ggplot2::autoplot(ab6), "gg")
  expect_inherits(ggplot2::autoplot(ab7), "gg")
  expect_inherits(ggplot2::autoplot(ab8), "gg")
}
