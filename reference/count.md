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
count_resistant(..., only_all_tested = FALSE)

count_susceptible(..., only_all_tested = FALSE)

count_S(..., only_all_tested = FALSE)

count_SI(..., only_all_tested = FALSE)

count_I(..., only_all_tested = FALSE)

count_IR(..., only_all_tested = FALSE)

count_R(..., only_all_tested = FALSE)

count_all(..., only_all_tested = FALSE)

n_sir(..., only_all_tested = FALSE)

count_df(data, translate_ab = "name", language = get_AMR_locale(),
  combine_SI = TRUE)
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
  *Combination Therapy* below.

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

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  all values of S, SDD, and I must be merged into one, so the output
  only consists of S+SDD+I vs. R (susceptible vs. resistant) - the
  default is `TRUE`.

## Value

An [integer](https://rdrr.io/r/base/integer.html)

## Details

These functions are meant to count isolates. Use the
[`resistance()`](https://amr-for-r.org/reference/proportion.md)/[`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
functions to calculate microbial resistance/susceptibility.

The function `count_resistant()` is equal to the function `count_R()`.
The function `count_susceptible()` is equal to the function
`count_SI()`.

The function `n_sir()` is an alias of `count_all()`. They can be used to
count all available isolates, i.e. where all input antimicrobials have
an available result (S, I or R). Their use is equal to
[`n_distinct()`](https://dplyr.tidyverse.org/reference/n_distinct.html).
Their function is equal to
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

This AMR package follows insight; use
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
(equal to
[`proportion_SI()`](https://amr-for-r.org/reference/proportion.md)) to
determine antimicrobial susceptibility and `count_susceptible()` (equal
to `count_SI()`) to count susceptible isolates.

## Combination Therapy

When using more than one variable for `...` (= combination therapy), use
`only_all_tested` to only count isolates that are tested for all
antimicrobials/variables that you test them for. See this example for
two antimicrobials, Drug A and Drug B, about how
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
works to calculate the %SI:

    --------------------------------------------------------------------
                        only_all_tested = FALSE  only_all_tested = TRUE
                        -----------------------  -----------------------
     Drug A    Drug B   considered   considered  considered   considered
                        susceptible    tested    susceptible    tested
    --------  --------  -----------  ----------  -----------  ----------
     S or I    S or I        X            X           X            X
       R       S or I        X            X           X            X
      <NA>     S or I        X            X           -            -
     S or I      R           X            X           X            X
       R         R           -            X           -            X
      <NA>       R           -            -           -            -
     S or I     <NA>         X            X           -            -
       R        <NA>         -            -           -            -
      <NA>      <NA>         -            -           -            -
    --------------------------------------------------------------------

Please note that, in combination therapies, for `only_all_tested = TRUE`
applies that:

        count_S()    +   count_I()    +   count_R()    = count_all()
      proportion_S() + proportion_I() + proportion_R() = 1

and that, in combination therapies, for `only_all_tested = FALSE`
applies that:

        count_S()    +   count_I()    +   count_R()    >= count_all()
      proportion_S() + proportion_I() + proportion_R() >= 1

Using `only_all_tested` has no impact when only using one antibiotic as
input.

## See also

[`proportion_*`](https://amr-for-r.org/reference/proportion.md) to
calculate microbial resistance and susceptibility.

## Examples

``` r
# example_isolates is a data set available in the AMR package.
# run ?example_isolates for more info.

# base R ------------------------------------------------------------
count_resistant(example_isolates$AMX) # counts "R"
#> [1] 804
count_susceptible(example_isolates$AMX) # counts "S" and "I"
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
#> ℹ For `aminoglycosides()` using columns 'GEN' (gentamicin), 'TOB'
#>   (tobramycin), 'AMK' (amikacin), and 'KAN' (kanamycin)
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
