# Data Set with 466 Microorganisms In Species Groups

A data set containing species groups and microbiological complexes,
which are used in [the clinical breakpoints
table](https://amr-for-r.org/release/3.1.0/reference/clinical_breakpoints.md).

## Usage

``` r
microorganisms.groups
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 466
observations and 4 variables:

- `mo_group`  
  ID of the species group / microbiological complex

- `mo`  
  ID of the microorganism belonging in the species group /
  microbiological complex

- `mo_group_name`  
  Name of the species group / microbiological complex, as retrieved with
  [`mo_name()`](https://amr-for-r.org/release/3.1.0/reference/mo_property.md)

- `mo_name`  
  Name of the microorganism belonging in the species group /
  microbiological complex, as retrieved with
  [`mo_name()`](https://amr-for-r.org/release/3.1.0/reference/mo_property.md)

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

## See also

[`as.mo()`](https://amr-for-r.org/release/3.1.0/reference/as.mo.md)
[microorganisms](https://amr-for-r.org/release/3.1.0/reference/microorganisms.md)

## Examples

``` r
microorganisms.groups

# these are all species in the Bacteroides fragilis group, as per WHONET:
microorganisms.groups[microorganisms.groups$mo_group == "B_BCTRD_FRGL-C", ]
```
