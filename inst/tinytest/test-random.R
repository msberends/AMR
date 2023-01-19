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
# doi:10.18637/jss.v104.i03                                            #
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

expect_inherits(random_mic(100), "mic")
expect_inherits(random_mic(100, mo = "Klebsiella pneumoniae"), "mic")
expect_inherits(random_mic(100, mo = "Klebsiella pneumoniae", ab = "meropenem"), "mic")
expect_inherits(random_mic(100, ab = "meropenem"), "mic")
# no normal factors of 2
expect_inherits(random_mic(100, "Haemophilus influenzae", "ceftaroline"), "mic")
expect_inherits(random_disk(100), "disk")
expect_inherits(random_disk(100, mo = "Klebsiella pneumoniae"), "disk")
expect_inherits(random_disk(100, mo = "Klebsiella pneumoniae", ab = "meropenem"), "disk")
expect_inherits(random_disk(100, ab = "meropenem"), "disk")
expect_inherits(random_sir(100), "sir")
