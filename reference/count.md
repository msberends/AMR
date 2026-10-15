# Count Available Isolates

These functions can be used to count resistant/susceptible microbial
isolates. All functions support quasiquotation with pipes, can be used
in [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)
from the `dplyr` package and also support grouped variables, see
*Examples*.

`count_resistant()` should be used to count resistant isolates,
`count_susceptible()` should be used to count susceptible isolates.

## Usage

``` r
count_resistant(
  ...,
  only_all_tested = FALSE,
  guideline = getOption("AMR_guideline", "EUCAST")
)

count_susceptible(
  ...,
  only_all_tested = FALSE,
  guideline = getOption("AMR_guideline", "EUCAST")
)

count_S(..., only_all_tested = FALSE)

count_SI(..., only_all_tested = FALSE)

count_I(..., only_all_tested = FALSE)

count_IR(..., only_all_tested = FALSE)

count_R(..., only_all_tested = FALSE)

count_all(..., only_all_tested = FALSE)

n_sir(..., only_all_tested = FALSE)

count_df(
  data,
  translate_ab = "name",
  language = get_AMR_locale(),
  combine_SI,
  guideline = getOption("AMR_guideline", "EUCAST")
)
```

## Arguments

- ...:

  One or more vectors (or columns) with antibiotic interpretations. They
  will be transformed internally with
  [`as.sir()`](https://amr-for-r.org/reference/as.sir.md) if needed.

- only_all_tested:

  (for combination therapies, i.e. using more than one variable for
  `...`): a [logical](https://rdrr.io/r/base/logical.html) to indicate
  that isolates must be tested for all antimicrobials, see section
  *Combination Therapy* below. Use
  [`sir_combine()`](https://amr-for-r.org/reference/sir_combine.md) to
  create a new variable for a combination therapy.

- guideline:

  A guideline name, `"EUCAST"` or `"CLSI"` (or any version thereof, e.g.
  `"CLSI 2025"`), to determine whether the 'I' category is considered
  susceptible (EUCAST) or resistant (CLSI). The default is the package
  option
  [`AMR_guideline`](https://amr-for-r.org/reference/AMR-options.md), or
  `"EUCAST"` if it is not set. For functions that return results per
  category, i.e., `count_df()`,
  [`proportion_df()`](https://amr-for-r.org/reference/proportion.md),
  [`sir_df()`](https://amr-for-r.org/reference/proportion.md),
  [`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`geom_sir()`](https://amr-for-r.org/reference/ggplot_sir.md), and
  [`labels_sir_count()`](https://amr-for-r.org/reference/ggplot_sir.md),
  this can also be `NULL` or `FALSE` to keep S, I, and R separate. See
  section *Interpretation of SIR* below.

- data:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) containing
  columns with class [`sir`](https://amr-for-r.org/reference/as.sir.md)
  (see [`as.sir()`](https://amr-for-r.org/reference/as.sir.md)).

- translate_ab:

  A column name of the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set to translate the antibiotic abbreviations to, using
  [`ab_property()`](https://amr-for-r.org/reference/ab_property.md).

- language:

  Language of the returned text - the default is the current system
  language (see
  [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md))
  and can also be set with the package option
  [`AMR_locale`](https://amr-for-r.org/reference/AMR-options.md). Use
  `language = NULL` or `language = ""` to prevent translation.

- combine_SI:

  *(deprecated in favour of `guideline`)* A
  [logical](https://rdrr.io/r/base/logical.html). If set, it takes
  precedence over `guideline`: `TRUE` merges S, SDD, and I into `"SI"`
  vs. `"R"`, and `FALSE` keeps S, I, and R separate.

## Value

An [integer](https://rdrr.io/r/base/integer.html)

## Details

These functions are meant to count isolates. Use the
[`resistance()`](https://amr-for-r.org/reference/proportion.md)/[`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
functions to calculate microbial resistance/susceptibility.

The function `n_sir()` is an alias of `count_all()`. They can be used to
count all available isolates, i.e. where the input antimicrobial has an
available result (S, I or R). For multiple antimicrobials, these are
isolates with a result for any of the antimicrobials, or with
`only_all_tested = TRUE`, for all of the antimicrobials (see
*Combination Therapy* below). Their use is equal to `dplyr`'s
[`n_distinct()`](https://dplyr.tidyverse.org/reference/n_distinct.html).
For one antimicrobial, their function is equal to
`count_susceptible(...) + count_resistant(...)`.

The function `count_df()` takes any variable from `data` that has an
[`sir`](https://amr-for-r.org/reference/as.sir.md) class (created with
[`as.sir()`](https://amr-for-r.org/reference/as.sir.md)) and counts the
number of S's, I's and R's. It also supports grouped variables. The
function [`sir_df()`](https://amr-for-r.org/reference/proportion.md)
works exactly like `count_df()`, but adds the percentage of S, I and R.

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
  possible in functions that return results per category: `count_df()`,
  [`proportion_df()`](https://amr-for-r.org/reference/proportion.md),
  [`sir_df()`](https://amr-for-r.org/reference/proportion.md),
  [`ggplot_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`geom_sir()`](https://amr-for-r.org/reference/ggplot_sir.md),
  [`labels_sir_count()`](https://amr-for-r.org/reference/ggplot_sir.md),
  and
  [`bug_drug_combinations()`](https://amr-for-r.org/reference/bug_drug_combinations.md).
  Functions that return a single proportion or count, such as
  [`susceptibility()`](https://amr-for-r.org/reference/proportion.md),
  `count_resistant()`,
  [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md),
  [`wisca()`](https://amr-for-r.org/reference/antibiogram.md), and
  [`format()`](https://rdrr.io/r/base/format.html) on the result of
  [`bug_drug_combinations()`](https://amr-for-r.org/reference/bug_drug_combinations.md),
  require 'I' to be considered either susceptible or resistant, and
  return an error.

Use [`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
and [`resistance()`](https://amr-for-r.org/reference/proportion.md) (or
`count_susceptible()` and `count_resistant()`) to follow the guideline,
or e.g.
[`proportion_SI()`](https://amr-for-r.org/reference/proportion.md) and
[`proportion_R()`](https://amr-for-r.org/reference/proportion.md) to be
explicit regardless of the guideline.

The `combine_SI` argument of the functions that return results per
category,
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md), and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md) is
deprecated in favour of `guideline`. If set, it is followed literally
and takes precedence over `guideline`.

## Combination Therapy

Results are ordered from susceptible to resistant: S (and WT), SDD, I, R
(and NWT, NS). With `method = "most_susceptible"` (in
[`sir_combine()`](https://amr-for-r.org/reference/sir_combine.md)), the
most susceptible result of the antimicrobials is returned, and with
`method = "most_resistant"`, the most resistant result. The value `"NI"`
(non-interpretable) is considered a missing result.

A missing result can only be ignored if it cannot change the conclusion:
a susceptible result for `"most_susceptible"`, and a non-susceptible
result for `"most_resistant"`. Under EUCAST, the 'I' category is
considered susceptible (*susceptible, increased exposure*), and under
CLSI, it is not. With `only_all_tested = TRUE`, any missing result leads
to `NA`.

This is an example for two antimicrobials, Drug A and Drug B, under
EUCAST (the order of the antimicrobials does not matter):

    -------------------------------------------------------------------------
                          Result of method =            Result of method =
                         "most_susceptible"              "most_resistant"
                       -----------------------  -----------------------------
                        only_all_tested =        only_all_tested =
     Drug A    Drug B     FALSE      TRUE          FALSE         TRUE
    --------  --------  ---------  ---------    -----------  -------------
        S        S          S          S             S             S
        S        I          S          S             I             I
        S        R          S          S             R             R
        I        R          I          I             R             R
        R        R          R          R             R             R
        S      <NA>         S        <NA>          <NA>          <NA>
        I      <NA>         I        <NA>          <NA>          <NA>
        R      <NA>       <NA>       <NA>            R           <NA>
      <NA>     <NA>       <NA>       <NA>          <NA>          <NA>
    -------------------------------------------------------------------------

Under CLSI, the combination of I with a missing result is `<NA>` for
`"most_susceptible"` and I for `"most_resistant"`.

The `count_*()` and
[`proportion_*()`](https://amr-for-r.org/reference/proportion.md)
functions apply these rules when more than one variable is used. With
`only_all_tested = FALSE` (default), the functions for susceptibility,
such as
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md) and
`count_susceptible()`, count an isolate as susceptible if any of the
antimicrobials is susceptible (`"most_susceptible"`), and the functions
for resistance, such as
[`resistance()`](https://amr-for-r.org/reference/proportion.md) and
`count_resistant()`, count an isolate as resistant if any of the
antimicrobials is resistant (`"most_resistant"`). An isolate is
considered tested if it has such a result, or if it was tested for all
antimicrobials. Therefore, for combinations with
`only_all_tested = FALSE`:

        count_S()    +   count_I()    +   count_R()    >= count_all()
      proportion_S() + proportion_I() + proportion_R() >= 1

With `only_all_tested = TRUE`, all these functions count the result of
the combination therapy (`"most_susceptible"`) of isolates tested for
all antimicrobials, so that:

        count_S()    +   count_I()    +   count_R()    = count_all()
      proportion_S() + proportion_I() + proportion_R() = 1

In [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md), the result
of the combination therapy (`"most_susceptible"`) is used, the 'I'
category is considered susceptible or not following the `guideline`
argument ('SDD' is always considered susceptible), and the result is
then reduced to S or R.

## See also

[`proportion_*`](https://amr-for-r.org/reference/proportion.md) to
calculate microbial resistance and susceptibility.

## Examples

``` r
# example_isolates is a data set available in the AMR package.
# run ?example_isolates for more info.

# base R ------------------------------------------------------------
count_resistant(example_isolates$AMX) # counts "R"
#> ℹ `count_resistant()` assumes the EUCAST guideline and thus considers the 'I'
#>   category susceptible. Use `guideline = "CLSI"` to consider it resistant. Set
#>   the `guideline` argument, or set the `AMR_guideline` option to change the
#>   default (see `?AMR-options`). This message is shown once per session.
#> [1] 804
count_susceptible(example_isolates$AMX) # counts "S" and "I"
#> ℹ `count_susceptible()` assumes the EUCAST guideline and thus considers the 'I'
#>   category susceptible. Use `guideline = "CLSI"` to consider it resistant. Set
#>   the `guideline` argument, or set the `AMR_guideline` option to change the
#>   default (see `?AMR-options`). This message is shown once per session.
#> [1] 546
count_all(example_isolates$AMX) # counts "S", "I" and "R"
#> [1] 1350

# be more specific
count_S(example_isolates$AMX)
#> [1] 543
count_SI(example_isolates$AMX)
#> [1] 546
count_I(example_isolates$AMX)
#> [1] 3
count_IR(example_isolates$AMX)
#> [1] 807
count_R(example_isolates$AMX)
#> [1] 804

# Count all available isolates
count_all(example_isolates$AMX)
#> [1] 1350
n_sir(example_isolates$AMX)
#> [1] 1350

# n_sir() is an alias of count_all().
# Since it counts all available isolates, you can
# calculate back to count e.g. susceptible isolates.
# These results are the same:
count_susceptible(example_isolates$AMX)
#> [1] 546
susceptibility(example_isolates$AMX) * n_sir(example_isolates$AMX)
#> ℹ `susceptibility()` assumes the EUCAST guideline and thus considers the 'I'
#>   category susceptible. Use `guideline = "CLSI"` to consider it resistant. Set
#>   the `guideline` argument, or set the `AMR_guideline` option to change the
#>   default (see `?AMR-options`). This message is shown once per session.
#> [1] 546

# dplyr -------------------------------------------------------------
# \donttest{
if (require("dplyr")) {
  example_isolates %>%
    group_by(ward) %>%
    summarise(
      R = count_R(CIP),
      I = count_I(CIP),
      S = count_S(CIP),
      n1 = count_all(CIP), # the actual total; sum of all three
      n2 = n_sir(CIP), # same - analogous to n_distinct
      total = n()
    ) # NOT the number of tested isolates!

  # Number of available isolates for a whole antibiotic class
  # (i.e., in this data set columns GEN, TOB, AMK, KAN)
  example_isolates %>%
    group_by(ward) %>%
    summarise(across(aminoglycosides(), n_sir))

  # Count co-resistance between amoxicillin/clav acid and gentamicin,
  # so we can see that combination therapy does a lot more than mono therapy.
  # Please mind that `susceptibility()` calculates percentages right away instead.
  example_isolates %>% count_susceptible(AMC) # 1433
  example_isolates %>% count_all(AMC) # 1879

  example_isolates %>% count_susceptible(GEN) # 1399
  example_isolates %>% count_all(GEN) # 1855

  example_isolates %>% count_susceptible(AMC, GEN) # 1764
  example_isolates %>% count_all(AMC, GEN) # 1936

  # Get number of S+I vs. R immediately of selected columns
  example_isolates %>%
    select(AMX, CIP) %>%
    count_df(translate = FALSE)

  # It also supports grouping variables
  example_isolates %>%
    select(ward, AMX, CIP) %>%
    group_by(ward) %>%
    count_df(translate = FALSE)
}
#> ℹ For `aminoglycosides()` using columns GEN (gentamicin), TOB (tobramycin), AMK
#>   (amikacin), and KAN (kanamycin)
#> ℹ `count_df()` assumes the EUCAST guideline and thus considers the 'I' category
#>   susceptible (S+SDD+I vs. R). Use `guideline = "CLSI"` for S+SDD vs. I+R, or
#>   `guideline = NULL` to keep S, I, and R separate. Set the `guideline`
#>   argument, or set the `AMR_guideline` option to change the default (see
#>   `?AMR-options`). This message is shown once per session.
#> # A tibble: 12 × 4
#>    ward       antibiotic interpretation value
#>    <chr>      <chr>      <ord>          <int>
#>  1 Clinical   AMX        SI               357
#>  2 Clinical   AMX        R                487
#>  3 Clinical   CIP        SI               741
#>  4 Clinical   CIP        R                128
#>  5 ICU        AMX        SI               158
#>  6 ICU        AMX        R                270
#>  7 ICU        CIP        SI               362
#>  8 ICU        CIP        R                 85
#>  9 Outpatient AMX        SI                31
#> 10 Outpatient AMX        R                 47
#> 11 Outpatient CIP        SI                78
#> 12 Outpatient CIP        R                 15
# }
```
