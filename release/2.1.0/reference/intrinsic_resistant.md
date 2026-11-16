# Data Set with Bacterial Intrinsic Resistance

Data set containing defined intrinsic resistance by EUCAST of all
bug-drug combinations.

## Usage

``` r
intrinsic_resistant
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 134
634 observations and 2 variables:

- `mo`  
  Microorganism ID

- `ab`  
  Antibiotic ID

## Details

This data set is based on ['EUCAST Expert Rules' and 'EUCAST Intrinsic
Resistance and Unusual Phenotypes'
v3.3](https://www.eucast.org/expert_rules_and_expected_phenotypes)
(2021).

### Direct download

Like all data sets in this package, this data set is publicly available
for download in the following formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the
download links](https://msberends.github.io/AMR/articles/datasets.html).
The actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

They **allow for machine reading EUCAST and CLSI guidelines**, which is
almost impossible with the MS Excel and PDF files distributed by EUCAST
and CLSI.

## Examples

``` r
intrinsic_resistant
```
