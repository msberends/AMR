# Data Set with 534 Microorganisms In Species Groups

A data set containing species groups and microbiological complexes,
which are used in [the clinical breakpoints
table](https://amr-for-r.org/reference/clinical_breakpoints.md).

## Usage

``` r
microorganisms.groups
```

## Format

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with 534
observations and 4 variables:

- `mo_group`  
  ID of the species group / microbiological complex

- `mo`  
  ID of the microorganism belonging in the species group /
  microbiological complex

- `mo_group_name`  
  Name of the species group / microbiological complex, as retrieved with
  [`mo_name()`](https://amr-for-r.org/reference/mo_property.md)

- `mo_name`  
  Name of the microorganism belonging in the species group /
  microbiological complex, as retrieved with
  [`mo_name()`](https://amr-for-r.org/reference/mo_property.md)

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

[`as.mo()`](https://amr-for-r.org/reference/as.mo.md)
[microorganisms](https://amr-for-r.org/reference/microorganisms.md)

## Examples

``` r
microorganisms.groups
#> # A tibble: 534 × 4
#>    mo_group       mo           mo_group_name                   mo_name          
#>    <mo>           <mo>         <chr>                           <chr>            
#>  1 B_ACNTB_BMNN-C B_ACNTB_BMNN Acinetobacter baumannii complex Acinetobacter ba…
#>  2 B_ACNTB_BMNN-C B_ACNTB_CLCC Acinetobacter baumannii complex Acinetobacter ca…
#>  3 B_ACNTB_BMNN-C B_ACNTB_LCTC Acinetobacter baumannii complex Acinetobacter di…
#>  4 B_ACNTB_BMNN-C B_ACNTB_NSCM Acinetobacter baumannii complex Acinetobacter no…
#>  5 B_ACNTB_BMNN-C B_ACNTB_PITT Acinetobacter baumannii complex Acinetobacter pi…
#>  6 B_ACNTB_BMNN-C B_ACNTB_SFRT Acinetobacter baumannii complex Acinetobacter se…
#>  7 B_BCTRD_FRGL-C B_BCTRD_FRGL Bacteroides fragilis complex    Bacteroides frag…
#>  8 B_BCTRD_FRGL-C B_BCTRD_OVTS Bacteroides fragilis complex    Bacteroides ovat…
#>  9 B_BCTRD_FRGL-C B_BCTRD_THTT Bacteroides fragilis complex    Bacteroides thet…
#> 10 B_BCTRD_FRGL-C B_PHCCL_VLGT Bacteroides fragilis complex    Bacteroides vulg…
#> # ℹ 524 more rows

# these are all species in the Bacteroides fragilis group, as per WHONET:
microorganisms.groups[microorganisms.groups$mo_group == "B_BCTRD_FRGL-C", ]
#> # A tibble: 5 × 4
#>   mo_group       mo           mo_group_name                mo_name              
#>   <mo>           <mo>         <chr>                        <chr>                
#> 1 B_BCTRD_FRGL-C B_BCTRD_FRGL Bacteroides fragilis complex Bacteroides fragilis 
#> 2 B_BCTRD_FRGL-C B_BCTRD_OVTS Bacteroides fragilis complex Bacteroides ovatus   
#> 3 B_BCTRD_FRGL-C B_BCTRD_THTT Bacteroides fragilis complex Bacteroides thetaiot…
#> 4 B_BCTRD_FRGL-C B_PHCCL_VLGT Bacteroides fragilis complex Bacteroides vulgatus 
#> 5 B_BCTRD_FRGL-C B_PRBCT_DSTS Bacteroides fragilis complex Parabacteroides dist…
```
