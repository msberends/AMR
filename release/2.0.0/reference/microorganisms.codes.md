# Data Set with 5 910 Common Microorganism Codes

A data set containing commonly used codes for microorganisms, from
laboratory systems and WHONET. Define your own with
[`set_mo_source()`](https://amr-for-r.org/release/2.0.0/reference/mo_source.md).
They will all be searched when using
[`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md) and
consequently all the
[`mo_*`](https://amr-for-r.org/release/2.0.0/reference/mo_property.md)
functions.

## Usage

``` r
microorganisms.codes
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 5
910 observations and 2 variables:

- `code`  
  Commonly used code of a microorganism

- `mo`  
  ID of the microorganism in the
  [microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
  data set

## Details

Like all data sets in this package, this data set is publicly available
for download in the following formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, SAS, and Stata. Please visit [our website for the
download links](https://msberends.github.io/AMR/articles/datasets.html).
The actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## See also

[`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md)
[microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)

## Examples

``` r
microorganisms.codes
```
