# Export Data Set as NCBI BioSample Antibiogram

Export Data Set as NCBI BioSample Antibiogram

## Usage

``` r
export_ncbi_biosample(x, filename = paste0("biosample_", format(Sys.time(),
  "%Y-%m-%d-%H%M%S"), ".xlsx"), type = "pathogen MIC",
  columns = where(is.mic), save_as_xlsx = TRUE)
```

## Arguments

- x:

  A data set.

- filename:

  A character string specifying the file name.

- type:

  A character string specifying the type of data set, either "pathogen
  MIC" or "beta-lactamase MIC", see
  <https://www.ncbi.nlm.nih.gov/biosample/docs/>.
