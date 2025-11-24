# Data Set with 2 000 Example Isolates

A data set containing 2 000 microbial isolates with their full
antibiograms. This data set contains randomised fictitious data, but
reflects reality and can be used to practise AMR data analysis. For
examples, please read [the tutorial on our
website](https://amr-for-r.org/articles/AMR.html).

## Usage

``` r
example_isolates
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 2
000 observations and 46 variables:

- `date`  
  Date of receipt at the laboratory

- `patient`  
  ID of the patient

- `age`  
  Age of the patient

- `gender`  
  Gender of the patient, either "F" or "M"

- `ward`  
  Ward type where the patient was admitted, either "Clinical", "ICU", or
  "Outpatient"

- `mo`  
  ID of microorganism created with
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md), see also the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set

- `PEN:RIF`  
  40 different antimicrobials with class
  [`sir`](https://amr-for-r.org/reference/as.sir.md) (see
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)); these column
  names occur in the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set and can be translated with
  [`set_ab_names()`](https://amr-for-r.org/reference/ab_property.md) or
  [`ab_name()`](https://amr-for-r.org/reference/ab_property.md)

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
example_isolates
#> # A tibble: 2,000 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-03 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  3 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  4 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  5 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  6 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  7 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  8 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  9 2002-01-16 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#> 10 2002-01-17 858515     79 F      ICU      B_STPHY_EPDR   R     NA    S     NA 
#> # ℹ 1,990 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …
```
