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
  [`as.mo()`](https://amr-for-r.org/release/3.0.0/reference/as.mo.md),
  see also
  [microorganisms](https://amr-for-r.org/release/3.0.0/reference/microorganisms.md)

- `AMX:GEN`  
  4 different antimicrobials that have to be transformed with
  [`as.sir()`](https://amr-for-r.org/release/3.0.0/reference/as.sir.md)

## Download Our Reference Data

All reference data sets in the AMR package - including information on
microorganisms, antimicrobials, and clinical breakpoints - are freely
available for download in multiple formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, and Stata.

For maximum compatibility, we also provide machine-readable,
tab-separated plain text files suitable for use in any software,
including laboratory information systems.

Visit [our website for direct download
links](https://amr-for-r.org/articles/datasets.html), or explore the
actual files in [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw/datasets).

## Examples

``` r
example_isolates_unclean
```
