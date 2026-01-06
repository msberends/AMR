# Determine Bug-Drug Combinations

Determine antimicrobial resistance (AMR) of all bug-drug combinations in
your data set where at least 30 (default) isolates are available per
species. Use [`format()`](https://rdrr.io/r/base/format.html) on the
result to prettify it to a publishable/printable format, see *Examples*.

## Usage

``` r
bug_drug_combinations(x, col_mo = NULL, FUN = mo_shortname,
  include_n_rows = FALSE, ...)

# S3 method for class 'bug_drug_combinations'
format(x, translate_ab = "name (ab, atc)",
  language = get_AMR_locale(), minimum = 30, combine_SI = TRUE,
  add_ab_group = TRUE, remove_intrinsic_resistant = FALSE,
  decimal.mark = getOption("OutDec"), big.mark = ifelse(decimal.mark ==
  ",", ".", ","), ...)
```

## Arguments

- x:

  A data set with antimicrobials columns, such as `amox`, `AMX` and
  `AMC`.

- col_mo:

  Column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)) - the default
  is the first column of class
  [`mo`](https://amr-for-r.org/reference/as.mo.md). Values will be
  coerced using [`as.mo()`](https://amr-for-r.org/reference/as.mo.md).

- FUN:

  The function to call on the `mo` column to transform the microorganism
  codes - the default is
  [`mo_shortname()`](https://amr-for-r.org/reference/mo_property.md).

- include_n_rows:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate if the
  total number of rows must be included in the output.

- ...:

  Arguments passed on to `FUN`.

- translate_ab:

  A [character](https://rdrr.io/r/base/character.html) of length 1
  containing column names of the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set.

- language:

  Language of the returned text - the default is the current system
  language (see
  [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md))
  and can also be set with the package option
  [`AMR_locale`](https://amr-for-r.org/reference/AMR-options.md). Use
  `language = NULL` or `language = ""` to prevent translation.

- minimum:

  The minimum allowed number of available (tested) isolates. Any isolate
  count lower than `minimum` will return `NA` with a warning. The
  default number of `30` isolates is advised by the Clinical and
  Laboratory Standards Institute (CLSI) as best practice, see *Source*.

- combine_SI:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  values S, SDD, and I should be summed, so resistance will be based on
  only R - the default is `TRUE`.

- add_ab_group:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate where the
  group of the antimicrobials must be included as a first column.

- remove_intrinsic_resistant:

  [logical](https://rdrr.io/r/base/logical.html) to indicate that rows
  and columns with 100% resistance for all tested antimicrobials must be
  removed from the table.

- decimal.mark:

  the character to be used to indicate the numeric decimal point.

- big.mark:

  character; if not empty used as mark between every `big.interval`
  decimals *before* (hence `big`) the decimal point.

## Value

The function `bug_drug_combinations()` returns a
[data.frame](https://rdrr.io/r/base/data.frame.html) with columns "mo",
"ab", "S", "SDD", "I", "R", and "total".

## Details

The function [`format()`](https://rdrr.io/r/base/format.html) calculates
the resistance per bug-drug combination and returns a table ready for
reporting/publishing. Use `combine_SI = TRUE` (default) to test R vs.
S+I and `combine_SI = FALSE` to test R+I vs. S. This table can also
directly be used in R Markdown / Quarto without the need for e.g.
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html).

## Examples

``` r
# example_isolates is a data set available in the AMR package.
# run ?example_isolates for more info.
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

# \donttest{
x <- bug_drug_combinations(example_isolates)
head(x)
#> # A tibble: 6 × 8
#>   mo                ab        S   SDD     I     R    NI total
#>   <chr>             <chr> <int> <int> <int> <int> <int> <int>
#> 1 (unknown species) AMC      15     0     0     0     0    15
#> 2 (unknown species) AMK       0     0     0     0     0     0
#> 3 (unknown species) AMP      15     0     0     1     0    16
#> 4 (unknown species) AMX      15     0     0     1     0    16
#> 5 (unknown species) AZM       3     0     0     3     0     6
#> 6 (unknown species) CAZ       0     0     0     0     0     0
#> Use 'format()' on this result to get a publishable/printable format.
format(x, translate_ab = "name (atc)")
#> # A tibble: 39 × 12
#>    Group     Drug  CoNS  `E. coli` `E. faecalis` `K. pneumoniae` `P. aeruginosa`
#>    <chr>     <chr> <chr> <chr>     <chr>         <chr>           <chr>          
#>  1 "Aminogl… Amik… "100… "  0.0% … "100.0% (39/… ""              ""             
#>  2 ""        Gent… " 13… "  2.0% … "100.0% (39/… " 10.3% (6/58)" "  0.0% (0/30)"
#>  3 ""        Kana… "100… ""        "100.0% (39/… ""              "100.0% (30/30…
#>  4 ""        Tobr… " 78… "  2.6% … "100.0% (39/… " 10.3% (6/58)" "  0.0% (0/30)"
#>  5 "Aminope… Amox… " 93… " 50.0% … ""            "100.0% (58/58… "100.0% (30/30…
#>  6 ""        Amox… " 42… " 13.1% … ""            " 10.3% (6/58)" "100.0% (30/30…
#>  7 ""        Ampi… " 93… " 50.0% … ""            "100.0% (58/58… "100.0% (30/30…
#>  8 "Carbape… Imip… " 47… "  0.0% … "  0.0% (0/3… "  0.0% (0/51)" ""             
#>  9 ""        Mero… " 47… "  0.0% … ""            "  0.0% (0/53)" ""             
#> 10 "Cephalo… Cefa… " 47… "  2.4% … "100.0% (39/… ""              "100.0% (30/30…
#> # ℹ 29 more rows
#> # ℹ 5 more variables: `P. mirabilis` <chr>, `S. aureus` <chr>,
#> #   `S. epidermidis` <chr>, `S. hominis` <chr>, `S. pneumoniae` <chr>

# Use FUN to change to transformation of microorganism codes
bug_drug_combinations(example_isolates,
  FUN = mo_gramstain
)
#> # A tibble: 80 × 8
#>    mo            ab        S   SDD     I     R    NI total
#>    <chr>         <chr> <int> <int> <int> <int> <int> <int>
#>  1 Gram-negative AMC     463     0    89   174     0   726
#>  2 Gram-negative AMK     251     0     0     5     0   256
#>  3 Gram-negative AMP     226     0     0   405     0   631
#>  4 Gram-negative AMX     226     0     0   405     0   631
#>  5 Gram-negative AZM       1     0     2   696     0   699
#>  6 Gram-negative CAZ     607     0     0    27     0   634
#>  7 Gram-negative CHL       1     0     0    30     0    31
#>  8 Gram-negative CIP     610     0    11    63     0   684
#>  9 Gram-negative CLI      18     0     1   709     0   728
#> 10 Gram-negative COL     309     0     0    78     0   387
#> # ℹ 70 more rows
#> Use 'format()' on this result to get a publishable/printable format.

bug_drug_combinations(example_isolates,
  FUN = function(x) {
    ifelse(x == as.mo("Escherichia coli"),
      "E. coli",
      "Others"
    )
  }
)
#> # A tibble: 80 × 8
#>    mo      ab        S   SDD     I     R    NI total
#>    <chr>   <chr> <int> <int> <int> <int> <int> <int>
#>  1 E. coli AMC     332     0    74    61     0   467
#>  2 E. coli AMK     171     0     0     0     0   171
#>  3 E. coli AMP     196     0     0   196     0   392
#>  4 E. coli AMX     196     0     0   196     0   392
#>  5 E. coli AZM       0     0     0   467     0   467
#>  6 E. coli CAZ     449     0     0    11     0   460
#>  7 E. coli CHL       0     0     0     0     0     0
#>  8 E. coli CIP     398     0     1    57     0   456
#>  9 E. coli CLI       0     0     0   467     0   467
#> 10 E. coli COL     240     0     0     0     0   240
#> # ℹ 70 more rows
#> Use 'format()' on this result to get a publishable/printable format.
# }
```
