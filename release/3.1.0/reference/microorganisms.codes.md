# Data Set with 6 293 Common Microorganism Codes

A data set containing commonly used codes for microorganisms, from
laboratory systems and [WHONET](https://whonet.org). Define your own
with
[`set_mo_source()`](https://amr-for-r.org/release/3.1.0/reference/mo_source.md).
They will all be searched when using
[`as.mo()`](https://amr-for-r.org/release/3.1.0/reference/as.mo.md) and
consequently all the
[`mo_*`](https://amr-for-r.org/release/3.1.0/reference/mo_property.md)
functions.

## Usage

``` r
microorganisms.codes
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 6
293 observations and 3 variables:

- `code`  
  Commonly used code of a microorganism. ***This is a unique
  identifier.***

- `mo`  
  ID of the microorganism in the
  [microorganisms](https://amr-for-r.org/release/3.1.0/reference/microorganisms.md)
  data set. This is the record of the name that the code was named
  after, also if that name is outdated, so that the code and the
  microorganism ID remain relatable: e.g., the WHONET code `EAE` refers
  to *Enterobacter aerogenes* (`B_ENTRBC_AERG`).
  [`as.mo()`](https://amr-for-r.org/release/3.1.0/reference/as.mo.md)
  returns the current name, here *Klebsiella aerogenes*. Codes for
  organisms that are not in the
  [microorganisms](https://amr-for-r.org/release/3.1.0/reference/microorganisms.md)
  data set, such as viruses, refer to `UNKNOWN`.

- `source`  
  Source of the code, with the name it refers to, e.g.
  `"WHONET (Enterobacter aerogenes)"`, or `NA` if the source is unknown

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
microorganisms.codes

# 'ECO' or 'eco' is the WHONET code for E. coli:
microorganisms.codes[microorganisms.codes$code == "ECO", ]

# and therefore, 'eco' will be understood as E. coli in this package:
mo_info("eco")

# works for all AMR functions:
mo_is_intrinsic_resistant("eco", ab = "vancomycin")
```
