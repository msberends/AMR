# Data Set Denoting Bacterial Intrinsic Resistance

Data set containing 'EUCAST Expected Resistant Phenotypes' of *all*
bug-drug combinations between the
[microorganisms](https://amr-for-r.org/reference/microorganisms.md) and
[antimicrobials](https://amr-for-r.org/reference/antimicrobials.md) data
sets.

## Usage

``` r
intrinsic_resistant
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 271
905 observations and 2 variables:

- `mo`  
  Microorganism ID which occurs in
  [`microorganisms$mo`](https://amr-for-r.org/reference/microorganisms.md).
  Names can be retrieved using
  [`mo_name()`](https://amr-for-r.org/reference/mo_property.md).

- `ab`  
  Antimicrobial ID which occurs in
  [`antimicrobials$ab`](https://amr-for-r.org/reference/antimicrobials.md).
  Names can be retrieved using
  [`ab_name()`](https://amr-for-r.org/reference/ab_property.md).

## Details

This data set is currently based on ['EUCAST Expected Resistant
Phenotypes'
v1.2](https://www.eucast.org/expert_rules_and_expected_phenotypes)
(2023).

This data set is internally used by:

- [`not_intrinsic_resistant()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  (an [antimicrobial
  selector](https://amr-for-r.org/reference/antimicrobial_selectors.md))

- [`mo_is_intrinsic_resistant()`](https://amr-for-r.org/reference/mo_property.md)

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
intrinsic_resistant
#> # A tibble: 271,905 × 2
#>    mo          ab  
#>    <mo>        <ab>
#>  1 B_GRAMP     ATM 
#>  2 B_GRAMP     COL 
#>  3 B_GRAMP     NAL 
#>  4 B_GRAMP     PLB 
#>  5 B_GRAMP     TEM 
#>  6 B_ANAER-POS ATM 
#>  7 B_ANAER-POS COL 
#>  8 B_ANAER-POS NAL 
#>  9 B_ANAER-POS PLB 
#> 10 B_ANAER-POS TEM 
#> # ℹ 271,895 more rows
```
