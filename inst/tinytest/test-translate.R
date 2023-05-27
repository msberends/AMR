# ==================================================================== #
# TITLE                                                                #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# CITE AS                                                              #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
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

expect_identical(mo_genus("B_GRAMP", language = "pt"), "(gênero desconhecido)")

expect_identical(mo_fullname("CoNS", "cs"), "Koaguláza-negativní stafylokok (KNS)")
expect_identical(mo_fullname("CoNS", "da"), "Koagulase-negative stafylokokker (KNS)")
expect_identical(mo_fullname("CoNS", "de"), "Koagulase-negative Staphylococcus (KNS)")
expect_identical(mo_fullname("CoNS", "el"), "Σταφυλόκοκκος με αρνητική πηκτικότητα (CoNS)")
expect_identical(mo_fullname("CoNS", "en"), "Coagulase-negative Staphylococcus (CoNS)")
expect_identical(mo_fullname("CoNS", "es"), "Staphylococcus coagulasa negativo (SCN)")
expect_identical(mo_fullname("CoNS", "fi"), "Koagulaasinegatiivinen stafylokokki (KNS)")
expect_identical(mo_fullname("CoNS", "fr"), "Staphylococcus à coagulase négative (CoNS)")
expect_identical(mo_fullname("CoNS", "it"), "Staphylococcus negativo coagulasi (CoNS)")
expect_identical(mo_fullname("CoNS", "ja"), "コアグラーゼ陰性ブドウ球菌 (グラム陰性)")
expect_identical(mo_fullname("CoNS", "nl"), "Coagulase-negatieve Staphylococcus (CNS)")
expect_identical(mo_fullname("CoNS", "no"), "Koagulase-negative stafylokokker (KNS)")
expect_identical(mo_fullname("CoNS", "pl"), "Staphylococcus koagulazoujemny (CoNS)")
expect_identical(mo_fullname("CoNS", "pt"), "Staphylococcus coagulase negativo (CoNS)")
expect_identical(mo_fullname("CoNS", "ro"), "Stafilococ coagulazo-negativ (SCN)")
expect_identical(mo_fullname("CoNS", "ru"), "Коагулазоотрицательный стафилококк (КОС)")
expect_identical(mo_fullname("CoNS", "sv"), "Koagulasnegativa stafylokocker (KNS)")
expect_identical(mo_fullname("CoNS", "tr"), "Koagülaz-negatif Stafilokok (KNS)")
expect_identical(mo_fullname("CoNS", "uk"), "Коагулазонегативний стафілокок (КНС)")
expect_identical(mo_fullname("CoNS", "zh"), "凝固酶阴性葡萄球菌 (CoNS)")

expect_error(mo_fullname("CoNS", "aa"))
