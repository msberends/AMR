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

#' The `AMR` Package
#'
#' @description
#' Welcome to the `AMR` package.
#'
#' The `AMR` package is a peer-reviewed, [free and open-source](https://amr-for-r.org/#copyright) R package with [zero dependencies](https://en.wikipedia.org/wiki/Dependency_hell) to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial data and properties, by using evidence-based methods. **Our aim is to provide a standard** for clean and reproducible AMR data analysis, that can therefore empower epidemiological analyses to continuously enable surveillance and treatment evaluation in any setting. We are a team of [many different researchers](https://amr-for-r.org/authors.html) from around the globe to make this a successful and durable project!
#'
#' This work was published in the Journal of Statistical Software (Volume 104(3); \doi{10.18637/jss.v104.i03}) and formed the basis of two PhD theses (\doi{10.33612/diss.177417131} and \doi{10.33612/diss.192486375}).
#'
#' After installing this package, R knows [**`r AMR:::format_included_data_number(AMR::microorganisms)` distinct microbial species**](https://amr-for-r.org/reference/microorganisms.html) (updated June 2024) and all [**`r AMR:::format_included_data_number(NROW(AMR::antimicrobials) + NROW(AMR::antivirals))` antimicrobial and antiviral drugs**](https://amr-for-r.org/reference/antimicrobials.html) by name and code (including ATC, EARS-Net, ASIARS-Net, PubChem, LOINC and SNOMED CT), and knows all about valid SIR and MIC values. The integral clinical breakpoint guidelines from CLSI `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, grepl("CLSI", guideline))$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, grepl("CLSI", guideline))$guideline)))` and EUCAST `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, grepl("EUCAST", guideline))$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, grepl("EUCAST", guideline))$guideline)))` are included, even with epidemiological cut-off (ECOFF) values. It supports and can read any data format, including WHONET data. This package works on Windows, macOS and Linux with all versions of R since R-3.0 (April 2013). **It was designed to work in any setting, including those with very limited resources**. It was created for both routine data analysis and academic research at the Faculty of Medical Sciences of the [University of Groningen](https://www.rug.nl) and the [University Medical Center Groningen](https://www.umcg.nl).
#'
#' The `AMR` package is available in `r vector_and(vapply(FUN.VALUE = character(1), LANGUAGES_SUPPORTED_NAMES, function(x) x$exonym), quotes = FALSE, sort = FALSE)`. Antimicrobial drug (group) names and colloquial microorganism names are provided in these languages.
#' @section Download Our Reference Data:
#' All reference data sets in the AMR package - including information on microorganisms, antimicrobials, and clinical breakpoints - are freely available for download in multiple formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, and Stata.
#'
#' For maximum compatibility, we also provide machine-readable, tab-separated plain text files suitable for use in any software, including laboratory information systems.
#'
#' Visit [our website for direct download links](https://amr-for-r.org/articles/datasets.html), or explore the actual files in [our GitHub repository](https://github.com/msberends/AMR/tree/main/data-raw/datasets).
#' @source
#' To cite AMR in publications use:
#'
#' Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C (2022). "AMR: An R Package for Working with Antimicrobial Resistance Data." _Journal of Statistical Software_, *104*(3), 1-31. \doi{10.18637/jss.v104.i03}
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
