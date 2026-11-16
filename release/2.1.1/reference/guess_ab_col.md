# Guess Antibiotic Column

This tries to find a column name in a data set based on information from
the
[antibiotics](https://amr-for-r.org/release/2.1.1/reference/antibiotics.md)
data set. Also supports WHONET abbreviations.

## Usage

``` r
guess_ab_col(
  x = NULL,
  search_string = NULL,
  verbose = FALSE,
  only_sir_columns = FALSE
)
```

## Arguments

- x:

  a [data.frame](https://rdrr.io/r/base/data.frame.html)

- search_string:

  a text to search `x` for, will be checked with
  [`as.ab()`](https://amr-for-r.org/release/2.1.1/reference/as.ab.md) if
  this value is not a column in `x`

- verbose:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  additional info should be printed

- only_sir_columns:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only antibiotic columns must be detected that were transformed to
  class `sir` (see
  [`as.sir()`](https://amr-for-r.org/release/2.1.1/reference/as.sir.md))
  on beforehand (default is `FALSE`)

## Value

A column name of `x`, or `NULL` when no result is found.

## Details

You can look for an antibiotic (trade) name or abbreviation and it will
search `x` and the
[antibiotics](https://amr-for-r.org/release/2.1.1/reference/antibiotics.md)
data set for any column containing a name or code of that antibiotic.

## Examples

``` r
df <- data.frame(
  amox = "S",
  tetr = "R"
)

guess_ab_col(df, "amoxicillin")
guess_ab_col(df, "J01AA07") # ATC code of tetracycline

guess_ab_col(df, "J01AA07", verbose = TRUE)
# NOTE: Using column 'tetr' as input for J01AA07 (tetracycline).

# WHONET codes
df <- data.frame(
  AMP_ND10 = "R",
  AMC_ED20 = "S"
)
guess_ab_col(df, "ampicillin")
guess_ab_col(df, "J01CR02")
guess_ab_col(df, as.ab("augmentin"))
```
