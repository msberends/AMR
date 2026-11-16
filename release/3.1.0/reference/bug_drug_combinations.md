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
  [`as.mo()`](https://amr-for-r.org/release/3.1.0/reference/as.mo.md)) -
  the default is the first column of class
  [`mo`](https://amr-for-r.org/release/3.1.0/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/3.1.0/reference/as.mo.md).

- FUN:

  The function to call on the `mo` column to transform the microorganism
  codes - the default is
  [`mo_shortname()`](https://amr-for-r.org/release/3.1.0/reference/mo_property.md).

- include_n_rows:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate if the
  total number of rows must be included in the output.

- guideline:

  A guideline name, `"EUCAST"` or `"CLSI"` (or any version thereof, e.g.
  `"CLSI 2025"`), to determine whether the 'I' category is considered
  susceptible (EUCAST) or resistant (CLSI) when calculating resistance
  with [`format()`](https://rdrr.io/r/base/format.html). The default is
  the package option
  [`AMR_guideline`](https://amr-for-r.org/release/3.1.0/reference/AMR-options.md),
  or `"EUCAST"` if it is not set. In `bug_drug_combinations()`, this can
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
  [antimicrobials](https://amr-for-r.org/release/3.1.0/reference/antimicrobials.md)
  data set.

- language:

  Language of the returned text - the default is the current system
  language (see
  [`get_AMR_locale()`](https://amr-for-r.org/release/3.1.0/reference/translate.md))
  and can also be set with the package option
  [`AMR_locale`](https://amr-for-r.org/release/3.1.0/reference/AMR-options.md).
  Use `language = NULL` or `language = ""` to prevent translation.

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
  [`AMR_guideline`](https://amr-for-r.org/release/3.1.0/reference/AMR-options.md),
  or `"EUCAST"` if this option is not set. In the latter case, a message
  is shown once per session (but not in R Markdown or Quarto documents);
  set the option to silence it.

- **`"EUCAST"`**, or any EUCAST version such as `"EUCAST 2025"`, to
  consider 'I' susceptible, e.g., to return `"SI"` vs. `"R"` in
  [`proportion_df()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md).

- **`"CLSI"`**, or any CLSI version such as `"CLSI 2025"`, to consider
  'I' resistant, e.g., to return `"S"` vs. `"IR"` in
  [`proportion_df()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md).

- **`NULL` or `FALSE`**, to keep S, I, and R separate. This is only
  possible in functions that return results per category:
  [`count_df()`](https://amr-for-r.org/release/3.1.0/reference/count.md),
  [`proportion_df()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md),
  [`sir_df()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md),
  [`ggplot_sir()`](https://amr-for-r.org/release/3.1.0/reference/ggplot_sir.md),
  [`geom_sir()`](https://amr-for-r.org/release/3.1.0/reference/ggplot_sir.md),
  [`labels_sir_count()`](https://amr-for-r.org/release/3.1.0/reference/ggplot_sir.md),
  and `bug_drug_combinations()`. Functions that return a single
  proportion or count, such as
  [`susceptibility()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md),
  [`count_resistant()`](https://amr-for-r.org/release/3.1.0/reference/count.md),
  [`antibiogram()`](https://amr-for-r.org/release/3.1.0/reference/antibiogram.md),
  [`wisca()`](https://amr-for-r.org/release/3.1.0/reference/antibiogram.md),
  and [`format()`](https://rdrr.io/r/base/format.html) on the result of
  `bug_drug_combinations()`, require 'I' to be considered either
  susceptible or resistant, and return an error.

Use
[`susceptibility()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md)
and
[`resistance()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md)
(or
[`count_susceptible()`](https://amr-for-r.org/release/3.1.0/reference/count.md)
and
[`count_resistant()`](https://amr-for-r.org/release/3.1.0/reference/count.md))
to follow the guideline, or e.g.
[`proportion_SI()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md)
and
[`proportion_R()`](https://amr-for-r.org/release/3.1.0/reference/proportion.md)
to be explicit regardless of the guideline.

The `combine_SI` argument of the functions that return results per
category,
[`antibiogram()`](https://amr-for-r.org/release/3.1.0/reference/antibiogram.md),
and
[`wisca()`](https://amr-for-r.org/release/3.1.0/reference/antibiogram.md)
is deprecated in favour of `guideline`. If set, it is followed literally
and takes precedence over `guideline`.

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
