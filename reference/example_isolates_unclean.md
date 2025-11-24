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
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md), see also
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)

- `AMX:GEN`  
  4 different antimicrobials that have to be transformed with
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)

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
#> # A tibble: 3,000 × 8
#>    patient_id hospital date       bacteria      AMX   AMC   CIP   GEN  
#>    <chr>      <chr>    <date>     <chr>         <chr> <chr> <chr> <chr>
#>  1 J3         A        2012-11-21 E. coli       R     I     S     S    
#>  2 R7         A        2018-04-03 K. pneumoniae R     I     S     S    
#>  3 P3         A        2014-09-19 E. coli       R     S     S     S    
#>  4 P10        A        2015-12-10 E. coli       S     I     S     S    
#>  5 B7         A        2015-03-02 E. coli       S     S     S     S    
#>  6 W3         A        2018-03-31 S. aureus     R     S     R     S    
#>  7 J8         A        2016-06-14 E. coli       R     S     S     S    
#>  8 M3         A        2015-10-25 E. coli       R     S     S     S    
#>  9 J3         A        2019-06-19 E. coli       S     S     S     S    
#> 10 G6         A        2015-04-27 S. aureus     S     S     S     S    
#> # ℹ 2,990 more rows
```
