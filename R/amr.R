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
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       #
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

#' The `AMR` Package
#'
#' @description
#' Welcome to the `AMR` package.
#'
#' `AMR` is a free, open-source and independent \R package to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial data and properties, by using evidence-based methods. Our aim is to provide a standard for clean and reproducible antimicrobial resistance data analysis, that can therefore empower epidemiological analyses to continuously enable surveillance and treatment evaluation in any setting.
#'
#' This work was published in the Journal of Statistical Software (Volume 104(3); \doi{10.18637/jss.v104.i03}) and formed the basis of two PhD theses (\doi{10.33612/diss.177417131} and \doi{10.33612/diss.192486375}).
#'
#' After installing this package, \R knows `r format_included_data_number(microorganisms)` distinct microbial species and all `r format_included_data_number(rbind(antibiotics[, "atc", drop = FALSE], antivirals[, "atc", drop = FALSE]))` antibiotic, antimycotic and antiviral drugs by name and code (including ATC, EARS-NET, LOINC and SNOMED CT), and knows all about valid R/SI and MIC values. It supports any data format, including WHONET/EARS-Net data.
#'
#' This package is fully independent of any other \R package and works on Windows, macOS and Linux with all versions of \R since R-3.0.0 (April 2013). It was designed to work in any setting, including those with very limited resources. It was created for both routine data analysis and academic research at the Faculty of Medical Sciences of the University of Groningen, in collaboration with non-profit organisations Certe Medical Diagnostics and Advice and University Medical Center Groningen. This \R package is actively maintained and free software; you can freely use and distribute it for both personal and commercial (but not patent) purposes under the terms of the GNU General Public License version 2.0 (GPL-2), as published by the Free Software Foundation.
#'
#' This package can be used for:
#' - Reference for the taxonomy of microorganisms, since the package contains all microbial (sub)species from the List of Prokaryotic names with Standing in Nomenclature (LPSN) and the Global Biodiversity Information Facility (GBIF)
#' - Interpreting raw MIC and disk diffusion values, based on the latest CLSI or EUCAST guidelines
#' - Retrieving antimicrobial drug names, doses and forms of administration from clinical health care records
#' - Determining first isolates to be used for AMR data analysis
#' - Calculating antimicrobial resistance
#' - Determining multi-drug resistance (MDR) / multi-drug resistant organisms (MDRO)
#' - Calculating (empirical) susceptibility of both mono therapy and combination therapies
#' - Predicting future antimicrobial resistance using regression models
#' - Getting properties for any microorganism (such as Gram stain, species, genus or family)
#' - Getting properties for any antibiotic (such as name, code of EARS-Net/ATC/LOINC/PubChem, defined daily dose or trade name)
#' - Plotting antimicrobial resistance
#' - Applying EUCAST expert rules
#' - Getting SNOMED codes of a microorganism, or getting properties of a microorganism based on a SNOMED code
#' - Getting LOINC codes of an antibiotic, or getting properties of an antibiotic based on a LOINC code
#' - Machine reading the EUCAST and CLSI guidelines from 2011-2020 to translate MIC values and disk diffusion diameters to R/SI
#' - Principal component analysis for AMR
#'
#' @section Reference Data Publicly Available:
#' All data sets in this `AMR` package (about microorganisms, antibiotics, R/SI interpretation, EUCAST rules, etc.) are publicly and freely available for download in the following formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, SAS, and Stata. We also provide tab-separated plain text files that are machine-readable and suitable for input in any software program, such as laboratory information systems. Please visit [our website for the download links](https://msberends.github.io/AMR/articles/datasets.html). The actual files are of course available on [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw).
#' @source
#' To cite AMR in publications use:
#'
#' Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C (2022). "AMR: An R Package for Working with Antimicrobial Resistance Data." _Journal of Statistical Software_, *104*(3), 1-31. \doi{10.18637/jss.v104.i03}.
#'
#' A BibTeX entry for LaTeX users is:
#'
#' \preformatted{
#' `r format(citation("AMR"), style = "bib")`
#' }
#' @name AMR
#' @keywords internal
#' @rdname AMR
"_PACKAGE"
