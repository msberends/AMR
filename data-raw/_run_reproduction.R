# Wrapper to run clinical breakpoints reproduction non-interactively
# Set UTF-8 locale so gsub() can handle Unicode patterns
Sys.setlocale("LC_CTYPE", "C.utf8")
Sys.setlocale("LC_ALL", "C.utf8")

# Overrides View() to just print a summary instead
View <- function(x, title = NULL) {
  if (is.data.frame(x) || is.matrix(x)) {
    cat("=== View() called:", if (!is.null(title)) title else deparse(substitute(x)), "===\n")
    cat("Dimensions:", nrow(x), "rows x", ncol(x), "cols\n")
    print(head(x, 10))
    cat("...\n\n")
  } else {
    print(x)
  }
  invisible(x)
}

setwd("/home/user/AMR")

cat("=== Step 1: Running reproduction_of_microorganisms.groups.R ===\n")
source("data-raw/_reproduction_scripts/reproduction_of_microorganisms.groups.R")

cat("\n=== Step 2: Running reproduction_of_clinical_breakpoints.R ===\n")
source("data-raw/_reproduction_scripts/reproduction_of_clinical_breakpoints.R")

cat("\n=== Done! ===\n")
