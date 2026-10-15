# Determine Bug-Drug Combinations

Determine antimicrobial resistance (AMR) of all bug-drug combinations in
your data set where at least 30 (default) isolates are available per
species. Use [`format()`](https://rdrr.io/r/base/format.html) on the
result to prettify it to a publishable/printable format, see *Examples*.

## Usage

``` r
bug_drug_combinations(
  x,
  col_mo = NULL,
  FUN = mo_shortname,
  include_n_rows = FALSE,
  guideline = getOption("AMR_guideline", "EUCAST"),
  ...
)

# S3 method for class 'bug_drug_combinations'
format(
  x,
  translate_ab = "name (ab, atc)",
  language = get_AMR_locale(),
  minimum = 30,
  combine_SI,
  guideline = attr(x, "guideline"),
  add_ab_group = TRUE,
  remove_intrinsic_resistant = FALSE,
  decimal.mark = getOption("OutDec"),
  big.mark = ifelse(decimal.mark == ",", ".", ","),
  ...
)
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

- guideline:

  A guideline name, `"EUCAST"` or `"CLSI"` (or any version thereof, e.g.
  `"CLSI 2025"`), to determine whether the 'I' category is considered
  susceptible (EUCAST) or resistant (CLSI) when calculating resistance
  with [`format()`](https://rdrr.io/r/base/format.html). The default is
  the package option
  [`AMR_guideline`](https://amr-for-r.org/reference/AMR-options.md), or
  `"EUCAST"` if it is not set. In `bug_drug_combinations()`, this can
  also be `NULL` or `FALSE`, after which
  [`format()`](https://rdrr.io/r/base/format.html) requires its own
  `guideline` to be set. In
  [`format()`](https://rdrr.io/r/base/format.html), the default is the
  guideline used in `bug_drug_combinations()`, which can be overridden
  by setting it. See section *Interpretation of SIR* below.

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

  *(deprecated in favour of `guideline`)* A
  [logical](https://rdrr.io/r/base/logical.html). If set, it takes
  precedence over `guideline`: `TRUE` bases resistance on R only, and
  `FALSE` on SDD, I, and R.

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
"ab", "S", "SDD", "I", "R", "NI", "WT", "NWT", "NS", and "total".

## Details

The function [`format()`](https://rdrr.io/r/base/format.html) calculates
the resistance per bug-drug combination and returns a table ready for
reporting/publishing. Following the `guideline`, resistance is based on
R (EUCAST) or on I and R (CLSI), see section *Interpretation of SIR*
below. This table can also directly be used in R Markdown / Quarto
without the need for e.g.
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html).

## Interpretation of SIR

In 2019, the European Committee on Antimicrobial Susceptibility Testing
(EUCAST) has decided to change the definitions of susceptibility testing
categories S, I, and R
(<https://www.eucast.org/bacteria/clinical-breakpoints-and-interpretation/definition-of-s-i-and-r/>).

In aggregated analyses, such as proportions, counts, antibiograms, and
plots, the 'I' category must be considered either susceptible or
resistant. In epidemiological AMR data analysis, this follows the
guideline: EUCAST considers 'I' susceptible, while CLSI considers 'I'
resistant. The dose-dependent susceptibility of CLSI ('SDD') is always
considered susceptible, as it equals what EUCAST defines as 'I'. All
functions with a `guideline` argument for this purpose use these
categories:

|             |               |               |
|-------------|---------------|---------------|
| `guideline` | Susceptible   | Resistant     |
| `"EUCAST"`  | S, SDD, I, WT | R, NWT, NS    |
| `"CLSI"`    | S, SDD, WT    | I, R, NWT, NS |

The values 'WT' (wild type) and 'NWT' (non-wild type) result from
epidemiological cut-off values (ECOFFs) instead of clinical breakpoints.
Since they are no clinical categories, 'WT' is always considered
susceptible and 'NWT' always resistant, regardless of the guideline. The
value 'NS' (non-susceptible) is always considered resistant, and 'NI'
(non-interpretable) is always considered a missing result.

The `guideline` argument can be:

- **Left unset**, to use the package option
  [`AMR_guideline`](https://amr-for-r.org/reference/AMR-options.md), or
  `"EUCAST"` if this option is not set. In the latter case, a message is
  shown once per session (but not in R Markdown or Quarto documents);
  set the option to silence it.

- **`"EUCAST"`**, or any EUCAST version such as `"EUCAST 2025"`, to
  consider 'I' susceptible, e.g., to return `"SI"` vs. `"R"` in
  [`proportion_df()`](https://amr-for-r.org/reference/proportion.md).

- **`"CLSI"`**, or any CLSI version such as `"CLSI 2025"`, to consider
  'I' resistant, e.g., to return `"S"` vs. `"IR"` in
  [`proportion_df()`](https://amr-for-r.org/reference/proportion.md).

- **`NULL` or `FALSE`**, to keep S, I, and R separate. This is only
  possible in functions that return results per category:
  [`count_df()`](https://amr-for-r.org/reference/count.md),
  [`proportion_df()`](https://amr-for-r.org/reference/proportion.md),
  [`sir_df()`](https://amr-for-r.org/reference/proportion.md),
  [`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`geom_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`labels_sir_count()`](https://amr-for-r.org/reference/ggplot_sir.md),
  and `bug_drug_combinations()`. Functions that return a single
  proportion or count, such as
  [`susceptibility()`](https://amr-for-r.org/reference/proportion.md),
  [`count_resistant()`](https://amr-for-r.org/reference/count.md),
  [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md),
  [`wisca()`](https://amr-for-r.org/reference/antibiogram.md), and
  [`format()`](https://rdrr.io/r/base/format.html) on the result of
  `bug_drug_combinations()`, require 'I' to be considered either
  susceptible or resistant, and return an error.

Use [`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
and [`resistance()`](https://amr-for-r.org/reference/proportion.md) (or
[`count_susceptible()`](https://amr-for-r.org/reference/count.md) and
[`count_resistant()`](https://amr-for-r.org/reference/count.md)) to
follow the guideline, or e.g.
[`proportion_SI()`](https://amr-for-r.org/reference/proportion.md) and
[`proportion_R()`](https://amr-for-r.org/reference/proportion.md) to be
explicit regardless of the guideline.

The `combine_SI` argument of the functions that return results per
category,
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md), and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md) is
deprecated in favour of `guideline`. If set, it is followed literally
and takes precedence over `guideline`.

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
#> # A tibble: 6 × 11
#>   mo                ab        S   SDD     I     R    NI    WT   NWT    NS total
#>   <chr>             <chr> <int> <int> <int> <int> <int> <int> <int> <int> <int>
#> 1 (unknown species) AMC      15     0     0     0     0     0     0     0    15
#> 2 (unknown species) AMK       0     0     0     0     0     0     0     0     0
#> 3 (unknown species) AMP      15     0     0     1     0     0     0     0    16
#> 4 (unknown species) AMX      15     0     0     1     0     0     0     0    16
#> 5 (unknown species) AZM       3     0     0     3     0     0     0     0     6
#> 6 (unknown species) CAZ       0     0     0     0     0     0     0     0     0
#> Use 'format()' on this result to get a publishable/printable format.
format(x, translate_ab = "name (atc)")
#> ℹ `bug_drug_combinations()` assumes the EUCAST guideline and thus considers the
#>   'I' category susceptible when calculating resistance in `format()` (R vs.
#>   S+SDD+I). Use `guideline = "CLSI"` for I+R vs. S+SDD. Set the `guideline`
#>   argument, or set the `AMR_guideline` option to change the default (see
#>   `?AMR-options`). This message is shown once per session.
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
#> # A tibble: 80 × 11
#>    mo            ab        S   SDD     I     R    NI    WT   NWT    NS total
#>    <chr>         <chr> <int> <int> <int> <int> <int> <int> <int> <int> <int>
#>  1 Gram-negative AMC     463     0    89   174     0     0     0     0   726
#>  2 Gram-negative AMK     251     0     0     5     0     0     0     0   256
#>  3 Gram-negative AMP     226     0     0   405     0     0     0     0   631
#>  4 Gram-negative AMX     226     0     0   405     0     0     0     0   631
#>  5 Gram-negative AZM       1     0     2   696     0     0     0     0   699
#>  6 Gram-negative CAZ     607     0     0    27     0     0     0     0   634
#>  7 Gram-negative CHL       1     0     0    30     0     0     0     0    31
#>  8 Gram-negative CIP     610     0    11    63     0     0     0     0   684
#>  9 Gram-negative CLI      18     0     1   709     0     0     0     0   728
#> 10 Gram-negative COL     309     0     0    78     0     0     0     0   387
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
#> # A tibble: 80 × 11
#>    mo      ab        S   SDD     I     R    NI    WT   NWT    NS total
#>    <chr>   <chr> <int> <int> <int> <int> <int> <int> <int> <int> <int>
#>  1 E. coli AMC     332     0    74    61     0     0     0     0   467
#>  2 E. coli AMK     171     0     0     0     0     0     0     0   171
#>  3 E. coli AMP     196     0     0   196     0     0     0     0   392
#>  4 E. coli AMX     196     0     0   196     0     0     0     0   392
#>  5 E. coli AZM       0     0     0   467     0     0     0     0   467
#>  6 E. coli CAZ     449     0     0    11     0     0     0     0   460
#>  7 E. coli CHL       0     0     0     0     0     0     0     0     0
#>  8 E. coli CIP     398     0     1    57     0     0     0     0   456
#>  9 E. coli CLI       0     0     0   467     0     0     0     0   467
#> 10 E. coli COL     240     0     0     0     0     0     0     0   240
#> # ℹ 70 more rows
#> Use 'format()' on this result to get a publishable/printable format.
# }
```
