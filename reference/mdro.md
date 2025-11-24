# Determine Multidrug-Resistant Organisms (MDRO)

Determine which isolates are multidrug-resistant organisms (MDRO)
according to international, national, or custom guidelines.

## Usage

``` r
mdro(x = NULL, guideline = "CMI 2012", col_mo = NULL, esbl = NA,
  carbapenemase = NA, mecA = NA, mecC = NA, vanA = NA, vanB = NA,
  info = interactive(), pct_required_classes = 0.5, combine_SI = TRUE,
  verbose = FALSE, only_sir_columns = any(is.sir(x)), ...)

brmo(x = NULL, only_sir_columns = any(is.sir(x)), ...)

mrgn(x = NULL, only_sir_columns = any(is.sir(x)), verbose = FALSE, ...)

mdr_tb(x = NULL, only_sir_columns = any(is.sir(x)), verbose = FALSE, ...)

mdr_cmi2012(x = NULL, only_sir_columns = any(is.sir(x)), verbose = FALSE,
  ...)

eucast_exceptional_phenotypes(x = NULL, only_sir_columns = any(is.sir(x)),
  verbose = FALSE, ...)
```

## Arguments

- x:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) with
  antimicrobials columns, like `AMX` or `amox`. Can be left blank for
  automatic determination.

- guideline:

  A specific guideline to follow, see sections *Supported international
  / national guidelines* and *Using Custom Guidelines* below. When left
  empty, the publication by Magiorakos *et al.* (see below) will be
  followed.

- col_mo:

  Column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)) - the default
  is the first column of class
  [`mo`](https://amr-for-r.org/reference/as.mo.md). Values will be
  coerced using [`as.mo()`](https://amr-for-r.org/reference/as.mo.md).

- esbl:

  [logical](https://rdrr.io/r/base/logical.html) values, or a column
  name containing logical values, indicating the presence of an ESBL
  gene (or production of its proteins).

- carbapenemase:

  [logical](https://rdrr.io/r/base/logical.html) values, or a column
  name containing logical values, indicating the presence of a
  carbapenemase gene (or production of its proteins).

- mecA:

  [logical](https://rdrr.io/r/base/logical.html) values, or a column
  name containing logical values, indicating the presence of a *mecA*
  gene (or production of its proteins).

- mecC:

  [logical](https://rdrr.io/r/base/logical.html) values, or a column
  name containing logical values, indicating the presence of a *mecC*
  gene (or production of its proteins).

- vanA:

  [logical](https://rdrr.io/r/base/logical.html) values, or a column
  name containing logical values, indicating the presence of a *vanA*
  gene (or production of its proteins).

- vanB:

  [logical](https://rdrr.io/r/base/logical.html) values, or a column
  name containing logical values, indicating the presence of a *vanB*
  gene (or production of its proteins).

- info:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  progress should be printed to the console - the default is only print
  while in interactive sessions.

- pct_required_classes:

  Minimal required percentage of antimicrobial classes that must be
  available per isolate, rounded down. For example, with the default
  guideline, 17 antimicrobial classes must be available for *S. aureus*.
  Setting this `pct_required_classes` argument to `0.5` (default) means
  that for every *S. aureus* isolate at least 8 different classes must
  be available. Any lower number of available classes will return `NA`
  for that isolate.

- combine_SI:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  all values of S and I must be merged into one, so resistance is only
  considered when isolates are R, not I. As this is the default
  behaviour of the `mdro()` function, it follows the redefinition by
  EUCAST about the interpretation of I (increased exposure) in 2019, see
  section 'Interpretation of S, I and R' below. When using
  `combine_SI = FALSE`, resistance is considered when isolates are R or
  I.

- verbose:

  A [logical](https://rdrr.io/r/base/logical.html) to turn Verbose mode
  on and off (default is off). In Verbose mode, the function returns a
  data set with the MDRO results in logbook form with extensive info
  about which isolates would be MDRO-positive, or why they are not.

- only_sir_columns:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only antimicrobial columns must be included that were transformed to
  class [sir](https://amr-for-r.org/reference/as.sir.md) on beforehand.
  Defaults to `FALSE` if no columns of `x` have a class
  [sir](https://amr-for-r.org/reference/as.sir.md).

- ...:

  Column names of antimicrobials. To automatically detect antimicrobial
  column names, do not provide any named arguments;
  [`guess_ab_col()`](https://amr-for-r.org/reference/guess_ab_col.md)
  will then be used for detection. To manually specify a column, provide
  its name (case-insensitive) as an argument, e.g.
  `AMX = "amoxicillin"`. To skip a specific antimicrobial, set it to
  `NULL`, e.g. `TIC = NULL` to exclude ticarcillin. If a manually
  defined column does not exist in the data, it will be skipped with a
  warning.

## Value

- If `verbose` is set to `TRUE`:  
  A [data.frame](https://rdrr.io/r/base/data.frame.html) containing
  columns `row_number`, `microorganism`, `MDRO`, `reason`,
  `all_nonsusceptible_columns`, `guideline`

- CMI 2012 paper - function `mdr_cmi2012()` or `mdro()`:  
  Ordered [factor](https://rdrr.io/r/base/factor.html) with levels
  `Negative` \< `Multi-drug-resistant (MDR)` \<
  `Extensively drug-resistant (XDR)` \< `Pandrug-resistant (PDR)`

- TB guideline - function `mdr_tb()` or `mdro(..., guideline = "TB")`:  
  Ordered [factor](https://rdrr.io/r/base/factor.html) with levels
  `Negative` \< `Mono-resistant` \< `Poly-resistant` \<
  `Multi-drug-resistant` \< `Extensively drug-resistant`

- German guideline - function `mrgn()` or
  `mdro(..., guideline = "MRGN")`:  
  Ordered [factor](https://rdrr.io/r/base/factor.html) with levels
  `Negative` \< `3MRGN` \< `4MRGN`

- Everything else, except for custom guidelines:  
  Ordered [factor](https://rdrr.io/r/base/factor.html) with levels
  `Negative` \< `Positive, unconfirmed` \< `Positive`. The value
  `"Positive, unconfirmed"` means that, according to the guideline, it
  is not entirely sure if the isolate is multi-drug resistant and this
  should be confirmed with additional (e.g. genotypic) tests

## Details

These functions are context-aware. This means that the `x` argument can
be left blank if used inside a
[data.frame](https://rdrr.io/r/base/data.frame.html) call, see
*Examples*.

For the `pct_required_classes` argument, values above 1 will be divided
by 100. This is to support both fractions (`0.75` or `3/4`) and
percentages (`75`).

**Note:** Every test that involves the Enterobacteriaceae family, will
internally be performed using its newly named *order* Enterobacterales,
since the Enterobacteriaceae family has been taxonomically reclassified
by Adeolu *et al.* in 2016. Before that, Enterobacteriaceae was the only
family under the Enterobacteriales (with an i) order. All species under
the old Enterobacteriaceae family are still under the new
Enterobacterales (without an i) order, but divided into multiple
families. The way tests are performed now by this `mdro()` function
makes sure that results from before 2016 and after 2016 are identical.

### Supported International / National Guidelines

Please suggest to implement guidelines by [letting us
know](https://github.com/msberends/AMR/issues/new?template=2-feature-request.yml&title=Add%20new%20MDRO%20guideline).

Currently supported guidelines are (case-insensitive):

- `guideline = "CMI 2012"` (default)

  Magiorakos AP, Srinivasan A *et al.* "Multidrug-resistant, extensively
  drug-resistant and pandrug-resistant bacteria: an international expert
  proposal for interim standard definitions for acquired resistance."
  Clinical Microbiology and Infection (2012)
  ([doi:10.1111/j.1469-0691.2011.03570.x](https://doi.org/10.1111/j.1469-0691.2011.03570.x)
  )

- `guideline = "EUCAST 3.3"` (or simply `guideline = "EUCAST"`)

  The European international guideline - EUCAST Expert Rules Version 3.3
  "Intrinsic Resistance and Unusual Phenotypes"
  ([link](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2021/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.3_20211018.pdf))

  Also:

  - `guideline = "EUCAST 3.2"`

    The former European international guideline - EUCAST Expert Rules
    Version 3.2 "Intrinsic Resistance and Unusual Phenotypes"
    ([link](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2020/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.2_20200225.pdf))

  - `guideline = "EUCAST 3.1"`

    The former European international guideline - EUCAST Expert Rules
    Version 3.1 "Intrinsic Resistance and Exceptional Phenotypes Tables"
    ([link](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf))

- `guideline = "TB"`

  The international guideline for multi-drug resistant tuberculosis -
  World Health Organization "Companion handbook to the WHO guidelines
  for the programmatic management of drug-resistant tuberculosis"
  ([link](https://www.who.int/publications/i/item/9789241548809))

- `guideline = "MRGN"`

  The German national guideline - Mueller et al. (2015) Antimicrobial
  Resistance and Infection Control 4:7;
  [doi:10.1186/s13756-015-0047-6](https://doi.org/10.1186/s13756-015-0047-6)

- `guideline = "BRMO 2024"` (or simply `guideline = "BRMO"`)

  The Dutch national guideline - Samenwerkingverband Richtlijnen
  Infectiepreventie (SRI) (2024) "Bijzonder Resistente Micro-Organismen
  (BRMO)" ([link](https://www.sri-richtlijnen.nl/brmo))

  Also:

  - `guideline = "BRMO 2017"`

    The former Dutch national guideline - Werkgroep Infectiepreventie
    (WIP), RIVM, last revision as of 2017: "Bijzonder Resistente
    Micro-Organismen (BRMO)"

### Using Custom Guidelines

Using a custom MDRO guideline is of importance if you have custom rules
to determine MDROs in your hospital, e.g., rules that are dependent on
ward, state of contact isolation or other variables in your data.

Custom guidelines can be set with the
[`custom_mdro_guideline()`](https://amr-for-r.org/reference/custom_mdro_guideline.md)
function.

## Interpretation of SIR

In 2019, the European Committee on Antimicrobial Susceptibility Testing
(EUCAST) has decided to change the definitions of susceptibility testing
categories S, I, and R (<https://www.eucast.org/newsiandr>).

This AMR package follows insight; use
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
(equal to
[`proportion_SI()`](https://amr-for-r.org/reference/proportion.md)) to
determine antimicrobial susceptibility and
[`count_susceptible()`](https://amr-for-r.org/reference/count.md) (equal
to [`count_SI()`](https://amr-for-r.org/reference/count.md)) to count
susceptible isolates.

## See also

[`custom_mdro_guideline()`](https://amr-for-r.org/reference/custom_mdro_guideline.md)

## Examples

``` r
out <- mdro(example_isolates)
#> Warning: in `mdro()`: NA introduced for isolates where the available percentage of
#> antimicrobial classes was below 50% (set with `pct_required_classes`)
str(out)
#>  Ord.factor w/ 4 levels "Negative"<"Multi-drug-resistant (MDR)"<..: NA NA 1 1 1 1 NA NA 1 1 ...
table(out)
#> out
#>                         Negative       Multi-drug-resistant (MDR) 
#>                             1617                              128 
#> Extensively drug-resistant (XDR)          Pandrug-resistant (PDR) 
#>                                0                                0 

out <- mdro(example_isolates, guideline = "EUCAST 3.3")
table(out)
#> out
#>              Negative Positive, unconfirmed              Positive 
#>                  1994                     0                     6 

# \donttest{
if (require("dplyr")) {
  # no need to define `x` when used inside dplyr verbs:
  example_isolates %>%
    mutate(MDRO = mdro()) %>%
    count(MDRO)
}
#> Warning: There was 1 warning in `mutate()`.
#> ℹ In argument: `MDRO = mdro()`.
#> Caused by warning:
#> ! in `mdro()`: NA introduced for isolates where the available percentage of
#> antimicrobial classes was below 50% (set with `pct_required_classes`)
#> # A tibble: 3 × 2
#>   MDRO                           n
#>   <ord>                      <int>
#> 1 Negative                    1617
#> 2 Multi-drug-resistant (MDR)   128
#> 3 NA                           255
# }
```
