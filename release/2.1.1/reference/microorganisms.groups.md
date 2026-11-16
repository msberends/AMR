# Data Set with 521 Microorganisms In Species Groups

A data set containing species groups and microbiological complexes,
which are used in [the clinical breakpoints
table](https://amr-for-r.org/release/2.1.1/reference/clinical_breakpoints.md).

## Usage

``` r
microorganisms.groups
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 521
observations and 4 variables:

- `mo_group`  
  ID of the species group / microbiological complex

- `mo`  
  ID of the microorganism belonging in the species group /
  microbiological complex

- `mo_group_name`  
  Name of the species group / microbiological complex, as retrieved with
  [`mo_name()`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md)

- `mo_name`  
  Name of the microorganism belonging in the species group /
  microbiological complex, as retrieved with
  [`mo_name()`](https://amr-for-r.org/release/2.1.1/reference/mo_property.md)

## Details

Like all data sets in this package, this data set is publicly available
for download in the following formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the
download links](https://msberends.github.io/AMR/articles/datasets.html).
The actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## See also

[`as.mo()`](https://amr-for-r.org/release/2.1.1/reference/as.mo.md)
[microorganisms](https://amr-for-r.org/release/2.1.1/reference/microorganisms.md)

## Examples

``` r
microorganisms.groups

# these are all species in the Bacteroides fragilis group, as per WHONET:
microorganisms.groups[microorganisms.groups$mo_group == "B_BCTRD_FRGL-C", ]
```
