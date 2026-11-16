# Data Set with Unclean Data

A data set containing 3 000 microbial isolates that are not cleaned up
and consequently not ready for AMR data analysis. This data set can be
used for practice.

## Usage

``` r
example_isolates_unclean
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 3
000 observations and 8 variables:

- `patient_id`  
  ID of the patient

- `date`  
  date of receipt at the laboratory

- `hospital`  
  ID of the hospital, from A to C

- `bacteria`  
  info about microorganism that can be transformed with
  [`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md),
  see also
  [microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)

- `AMX:GEN`  
  4 different antibiotics that have to be transformed with
  [`as.sir()`](https://amr-for-r.org/release/2.1.1/reference/as.sir.md)

## Details

Like all data sets in this package, this data set is publicly available
for download in the following formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the
download links](https://msberends.github.io/AMR/articles/datasets.html).
The actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## Examples

``` r
example_isolates_unclean
```
