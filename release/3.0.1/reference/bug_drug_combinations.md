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
  [`as.mo()`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md)) -
  the default is the first column of class
  [`mo`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md).

- FUN:

  The function to call on the `mo` column to transform the microorganism
  codes - the default is
  [`mo_shortname()`](https://amr-for-r.org/release/3.0.1/reference/mo_property.md).

- include_n_rows:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate if the
  total number of rows must be included in the output.

- ...:

  Arguments passed on to `FUN`.

- translate_ab:

  A [character](https://rdrr.io/r/base/character.html) of length 1
  containing column names of the
  [antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
  data set.

- language:

  Language of the returned text - the default is the current system
  language (see
  [`get_AMR_locale()`](https://amr-for-r.org/release/3.0.1/reference/translate.md))
  and can also be set with the package option
  [`AMR_locale`](https://amr-for-r.org/release/3.0.1/reference/AMR-options.md).
  Use `language = NULL` or `language = ""` to prevent translation.

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

# \donttest{
x <- bug_drug_combinations(example_isolates)
head(x)
format(x, translate_ab = "name (atc)")

# Use FUN to change to transformation of microorganism codes
bug_drug_combinations(example_isolates,
  FUN = mo_gramstain
)

bug_drug_combinations(example_isolates,
  FUN = function(x) {
    ifelse(x == as.mo("Escherichia coli"),
      "E. coli",
      "Others"
    )
  }
)
# }
```
