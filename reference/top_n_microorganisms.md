# Filter Top *n* Microorganisms

This function filters a data set to include only the top *n*
microorganisms based on a specified property, such as taxonomic family
or genus. For example, it can filter a data set to the top 3 species, or
to any species in the top 5 genera, or to the top 3 species in each of
the top 5 genera.

## Usage

``` r
top_n_microorganisms(x, n, property = "species", n_for_each = NULL,
  col_mo = NULL, ...)
```

## Arguments

- x:

  A data frame containing microbial data.

- n:

  An integer specifying the maximum number of unique values of the
  `property` to include in the output.

- property:

  A character string indicating the microorganism property to use for
  filtering. Must be one of the column names of the
  [microorganisms](https://amr-for-r.org/reference/microorganisms.md)
  data set: "mo", "fullname", "status", "kingdom", "phylum", "class",
  "order", "family", "genus", "species", "subspecies", "rank", "ref",
  "oxygen_tolerance", "source", "lpsn", "lpsn_parent",
  "lpsn_renamed_to", "mycobank", "mycobank_parent",
  "mycobank_renamed_to", "gbif", "gbif_parent", "gbif_renamed_to",
  "prevalence", or "snomed". If `NULL`, the raw values from `col_mo`
  will be used without transformation. When using `"species"` (default)
  or `"subpecies"`, the genus will be added to make sure each
  (sub)species still belongs to the right genus.

- n_for_each:

  An optional integer specifying the maximum number of rows to retain
  for each value of the selected property. If `NULL`, all rows within
  the top *n* groups will be included.

- col_mo:

  A character string indicating the column in `x` that contains
  microorganism names or codes. Defaults to the first column of class
  [`mo`](https://amr-for-r.org/reference/as.mo.md). Values will be
  coerced using [`as.mo()`](https://amr-for-r.org/reference/as.mo.md).

- ...:

  Additional arguments passed on to
  [`mo_property()`](https://amr-for-r.org/reference/mo_property.md) when
  `property` is not `NULL`.

## Details

This function is useful for preprocessing data before creating
[antibiograms](https://amr-for-r.org/reference/antibiogram.md) or other
analyses that require focused subsets of microbial data. For example, it
can filter a data set to only include isolates from the top 10 species.

## See also

[`mo_property()`](https://amr-for-r.org/reference/mo_property.md),
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md),
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md)

## Examples

``` r
# filter to the top 3 species:
top_n_microorganisms(example_isolates,
  n = 3
)
#> # A tibble: 1,015 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-03 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  3 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  4 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  5 2002-01-19 738003     71 M      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  6 2002-01-19 738003     71 M      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  7 2002-02-03 481442     76 M      ICU      B_STPHY_CONS   R     NA    S     NA 
#>  8 2002-02-14 067927     45 F      ICU      B_STPHY_CONS   R     NA    R     NA 
#>  9 2002-02-14 067927     45 F      ICU      B_STPHY_CONS   S     NA    S     NA 
#> 10 2002-02-21 A56499     64 M      Clinical B_STPHY_CONS   S     NA    S     NA 
#> # ℹ 1,005 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …

# filter to any species in the top 5 genera:
top_n_microorganisms(example_isolates,
  n = 5, property = "genus"
)
#> # A tibble: 1,742 × 46
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
#> # ℹ 1,732 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …

# filter to the top 3 species in each of the top 5 genera:
top_n_microorganisms(example_isolates,
  n = 5, property = "genus", n_for_each = 3
)
#> # A tibble: 1,497 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-02-21 4FC193     69 M      Clinical B_ENTRC_FACM   NA    NA    NA    NA 
#>  2 2002-04-08 130252     78 M      ICU      B_ENTRC_FCLS   NA    NA    NA    NA 
#>  3 2002-06-23 798871     82 M      Clinical B_ENTRC_FCLS   NA    NA    NA    NA 
#>  4 2002-06-23 798871     82 M      Clinical B_ENTRC_FCLS   NA    NA    NA    NA 
#>  5 2003-04-20 6BC362     62 M      ICU      B_ENTRC        NA    NA    NA    NA 
#>  6 2003-04-21 6BC362     62 M      ICU      B_ENTRC        NA    NA    NA    NA 
#>  7 2003-08-13 F35553     52 M      ICU      B_ENTRC_FCLS   NA    NA    NA    NA 
#>  8 2003-09-05 F35553     52 M      ICU      B_ENTRC        NA    NA    NA    NA 
#>  9 2003-09-05 F35553     52 M      ICU      B_ENTRC_FCLS   NA    NA    NA    NA 
#> 10 2003-09-28 1B0933     80 M      Clinical B_ENTRC        NA    NA    NA    NA 
#> # ℹ 1,487 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …
```
