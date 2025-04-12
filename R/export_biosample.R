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

#' Export Data Set as NCBI BioSample Antibiogram
#'
#'
#' @param x A data set.
#' @param filename A character string specifying the file name.
#' @param type A character string specifying the type of data set, either "pathogen MIC" or "beta-lactamase MIC", see <https://www.ncbi.nlm.nih.gov/biosample/docs/>.
#' @keywords internal
export_ncbi_biosample <- function(x,
                                  filename = paste0("biosample_", format(Sys.time(), "%Y-%m-%d-%H%M%S"), ".xlsx"),
                                  type = "pathogen MIC",
                                  columns = where(is.mic),
                                  save_as_xlsx = TRUE) {
  meet_criteria(x, allow_class = "data.frame") # also checks dimensions to be >0
  meet_criteria(filename, allow_class = "character", has_length = 1)
  meet_criteria(type, allow_class = "character", has_length = 1, is_in = c("pathogen MIC", "beta-lactamase MIC"))
  meet_criteria(save_as_xlsx, allow_class = "logical", has_length = 1)

  out <- x %pm>%
    pm_select(columns)
  stop_if(NROW(out) == 0, "No columns found.")

  if (isTRUE(save_as_xlsx)) {
    export <- import_fn("write.xlsx", pkg = "openxlsx", error_on_fail = TRUE)
    export(out, file = filename, overwrite = TRUE, asTable = FALSE)
  } else {
    out
  }
}
