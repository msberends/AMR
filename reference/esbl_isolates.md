# Data Set with 500 ESBL Isolates

A data set containing 500 microbial isolates with MIC values of common
antibiotics and a binary `esbl` column for extended-spectrum
beta-lactamase (ESBL) production. This data set contains randomised
fictitious data but reflects reality and can be used to practise
AMR-related machine learning, e.g., classification modelling with
[tidymodels](https://amr-for-r.org/articles/AMR_with_tidymodels.html).

## Usage

``` r
esbl_isolates
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 500
observations and 19 variables:

- `esbl`  
  Logical indicator if the isolate is ESBL-producing

- `genus`  
  Genus of the microorganism

- `AMC:COL`  
  MIC values for 17 antimicrobial agents, transformed to class
  [`mic`](https://amr-for-r.org/reference/as.mic.md) (see
  [`as.mic()`](https://amr-for-r.org/reference/as.mic.md))

## Details

See our [tidymodels
integration](https://amr-for-r.org/reference/amr-tidymodels.md) for an
example using this data set.

## Examples

``` r
esbl_isolates
#> # A tibble: 500 × 19
#>    esbl  genus   AMC   AMP   TZP   CXM   FOX   CTX   CAZ   GEN   TOB   TMP   SXT
#>    <lgl> <chr> <mic> <mic> <mic> <mic> <mic> <mic> <mic> <mic> <mic> <mic> <mic>
#>  1 FALSE Esch…    32    32     4    64    64  8.00  8.00     1     1  16.0    20
#>  2 FALSE Esch…    32    32     4    64    64  4.00  8.00     1     1  16.0   320
#>  3 FALSE Esch…     4     2    64     8     4  8.00  0.12    16    16   0.5    20
#>  4 FALSE Kleb…    32    32    16    64    64  8.00  8.00     1     1   0.5    20
#>  5 FALSE Esch…    32    32     4     4     4  0.25  2.00     1     1  16.0   320
#>  6 FALSE Citr…    32    32    16    64    64 64.00 32.00     1     1   0.5    20
#>  7 FALSE Morg…    32    32     4    64    64 16.00  2.00     1     1   0.5    20
#>  8 FALSE Prot…    16    32     4     1     4  8.00  0.12     1     1  16.0   320
#>  9 FALSE Ente…    32    32     8    64    64 32.00  4.00     1     1   0.5    20
#> 10 FALSE Citr…    32    32    32    64    64  8.00 64.00     1     1  16.0   320
#> # ℹ 490 more rows
#> # ℹ 6 more variables: NIT <mic>, FOS <mic>, CIP <mic>, IPM <mic>, MEM <mic>,
#> #   COL <mic>
```
