# Guess Antibiotic Column

This tries to find a column name in a data set based on information from
the
[antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
data set. Also supports WHONET abbreviations.

## Usage

``` r
guess_ab_col(x = NULL, search_string = NULL, verbose = FALSE,
  only_sir_columns = FALSE)
```

## Arguments

- x:

  A [data.frame](https://rdrr.io/r/base/data.frame.html).

- search_string:

  A text to search `x` for, will be checked with
  [`as.ab()`](https://amr-for-r.org/release/3.0.1/reference/as.ab.md) if
  this value is not a column in `x`.

- verbose:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  additional info should be printed.

- only_sir_columns:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only antimicrobial columns must be included that were transformed to
  class [sir](https://amr-for-r.org/release/3.0.1/reference/as.sir.md)
  on beforehand. Defaults to `FALSE` if no columns of `x` have a class
  [sir](https://amr-for-r.org/release/3.0.1/reference/as.sir.md).

## Value

A column name of `x`, or `NULL` when no result is found.

## Details

You can look for an antibiotic (trade) name or abbreviation and it will
search `x` and the
[antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
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

# WHONET codes
df <- data.frame(
  AMP_ND10 = "R",
  AMC_ED20 = "S"
)
guess_ab_col(df, "ampicillin")
guess_ab_col(df, "J01CR02")
guess_ab_col(df, "augmentin")
```
