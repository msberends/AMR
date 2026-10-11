# Combine SIR Results of Multiple Antimicrobials

Combine the SIR results of two or more antimicrobials into one SIR
result per isolate, e.g. to create a variable for the combination
therapy of piperacillin/tazobactam with gentamicin. The result is the
best result of the antimicrobials, as a combination therapy covers the
isolate if any of its antimicrobials is active. See *Combination
Therapy* below for the exact rules.

The same rules are used by
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md) for
combination antibiograms, and by the
[`count_*()`](https://amr-for-r.org/reference/count.md) and
[`proportion_*()`](https://amr-for-r.org/reference/proportion.md)
functions when more than one variable is used, so the result of
`sir_combine()` can be used in any further analysis.

## Usage

``` r
sir_combine(
  ...,
  only_all_tested = FALSE,
  guideline = getOption("AMR_guideline", "EUCAST")
)
```

## Arguments

- ...:

  Two or more vectors (or columns) with antimicrobial interpretations,
  or a data frame with such columns. They will be transformed internally
  with [`as.sir()`](https://amr-for-r.org/reference/as.sir.md) if
  needed. All vectors must be of the same length.

- only_all_tested:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate that
  isolates must be tested for all antimicrobials. If `TRUE`, the result
  is `NA` for isolates with a missing result for any of the
  antimicrobials. If `FALSE` (default), a missing result only leads to
  `NA` if a result for the missing antimicrobial could have changed the
  result of the combination, see *Combination Therapy* below.

- guideline:

  Either `"EUCAST"` (default) or `"CLSI"`, to determine whether the 'I'
  category is considered susceptible. This only affects isolates with a
  missing result for any of the antimicrobials when
  `only_all_tested = FALSE`, see *Combination Therapy* below. The
  default can be changed with the `AMR_guideline` option, see
  [AMR-options](https://amr-for-r.org/reference/AMR-options.md).

## Value

A [`sir`](https://amr-for-r.org/reference/as.sir.md) vector with the
same length as the input.

## Combination Therapy

A combination of antimicrobials is considered as good as its most active
antimicrobial: `"S"` is preferred over `"WT"`, and `"R"` over `"NWT"`
and `"NS"`. The value `"NI"` (non-interpretable) is considered a missing
result.

A missing result can only be ignored if the result of the combination is
susceptible regardless of the missing result, which depends on the
`only_all_tested` argument and on the guideline: under EUCAST, the 'I'
category is considered susceptible (*susceptible, increased exposure*),
and under CLSI, it is not.

This is an example for two antimicrobials, Drug A and Drug B (the order
of the antimicrobials does not matter):

    -------------------------------------------------------------------------
                                Result of the combination
                        -----------------------------------------------------
     Drug A    Drug B   only_all_tested = FALSE       only_all_tested = TRUE
    --------  --------  ----------------------------  ----------------------
        S     S, I, R   S                             S
        I       I, R    I                             I
        R         R     R                             R
        S       <NA>    S                             <NA>
        I       <NA>    I (EUCAST), <NA> (CLSI)       <NA>
        R       <NA>    <NA>                          <NA>
      <NA>      <NA>    <NA>                          <NA>
    -------------------------------------------------------------------------

The [`count_*()`](https://amr-for-r.org/reference/count.md) and
[`proportion_*()`](https://amr-for-r.org/reference/proportion.md)
functions, such as
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md) and
[`resistance()`](https://amr-for-r.org/reference/proportion.md), apply
exactly these rules when more than one variable is used: an isolate is
considered tested if the result of the combination is not `NA`, and is
then counted by the result of the combination. Consequently,
[`count_susceptible()`](https://amr-for-r.org/reference/count.md) `+`
[`count_resistant()`](https://amr-for-r.org/reference/count.md) `=`
[`count_all()`](https://amr-for-r.org/reference/count.md), and
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md) `+`
[`resistance()`](https://amr-for-r.org/reference/proportion.md) `= 1`,
for both values of `only_all_tested`. The guideline is taken from their
`guideline` argument if they have one, and otherwise from the
`AMR_guideline` option.

In [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md), the 'I' and
'SDD' categories are considered susceptible if `combine_SI = TRUE`
(default), and the result of the combination is then reduced to S or R.

## See also

[`count()`](https://amr-for-r.org/reference/count.md) and
[`proportion()`](https://amr-for-r.org/reference/proportion.md) for
counting and calculating susceptibility of (combinations of)
antimicrobials, and
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) for
combination antibiograms.

## Examples

``` r
drug_a <- as.sir(c("S", "R", "R", "I", "I", NA, "R", "SDD"))
drug_b <- as.sir(c("R", "S", "R", "R", NA, NA, NA, "I"))

sir_combine(drug_a, drug_b)
#> Class <sir>
#> [1] S    S    R    I    I    <NA> <NA> SDD 
sir_combine(drug_a, drug_b, guideline = "CLSI")
#> Class <sir>
#> [1] S    S    R    I    <NA> <NA> <NA> SDD 
sir_combine(drug_a, drug_b, only_all_tested = TRUE)
#> Class <sir>
#> [1] S    S    R    I    <NA> <NA> <NA> SDD 

# the count and proportion functions use the same rules
susceptibility(drug_a, drug_b, minimum = 0)
#> Note that `proportion_SI()` will also include dose-dependent susceptibility,
#> "SDD". This note will be shown once for this session.
#> [1] 0.8333333
susceptibility(sir_combine(drug_a, drug_b), minimum = 0)
#> [1] 0.8333333

# \donttest{
if (require("dplyr")) {
  # create a new variable for piperacillin/tazobactam + gentamicin
  example_isolates %>%
    filter(!is.na(TZP)) %>%
    mutate(TZP_GEN = sir_combine(TZP, GEN)) %>%
    select(mo, TZP, GEN, TZP_GEN)
}
#> # A tibble: 1,001 × 4
#>    mo           TZP   GEN   TZP_GEN
#>    <mo>         <sir> <sir> <sir>  
#>  1 B_STRPT_SNGN   S     R     S    
#>  2 B_STRPT_GRPB   S     R     S    
#>  3 B_ESCHR_COLI   R     S     S    
#>  4 B_ESCHR_COLI   R     S     S    
#>  5 B_STRPT_PNMN   S     R     S    
#>  6 B_STRPT_PNMN   S     R     S    
#>  7 B_STRPT_ANGN   S     R     S    
#>  8 B_STRPT_ANGN   S     R     S    
#>  9 B_STRPT_PNMN   S     R     S    
#> 10 B_STRPT_SLVR   S     R     S    
#> # ℹ 991 more rows
if (require("dplyr")) {
  # combine any number of antimicrobials
  example_isolates %>%
    mutate(any_aminoglycoside = sir_combine(GEN, TOB, AMK, KAN)) %>%
    count(any_aminoglycoside)
}
#> # A tibble: 4 × 2
#>   any_aminoglycoside     n
#>   <sir>              <int>
#> 1   S                 1383
#> 2   I                   24
#> 3   R                  427
#> 4   NA                 166
if (require("dplyr")) {
  # use the new variable in further analyses
  # THOUGH CONSIDER wisca() FOR RELIABLE ESTIMATES
  example_isolates %>%
    mutate(TZP_GEN = sir_combine(TZP, GEN)) %>%
    group_by(ward) %>%
    summarise(
      TZP_mono = susceptibility(TZP),
      GEN_mono = susceptibility(GEN),
      TZP_GEN = susceptibility(TZP_GEN),
      # or directly, without combining first:
      TZP_GEN2 = susceptibility(TZP, GEN)
    )
}
#> # A tibble: 3 × 5
#>   ward       TZP_mono GEN_mono TZP_GEN TZP_GEN2
#>   <chr>         <dbl>    <dbl>   <dbl>    <dbl>
#> 1 Clinical      0.898    0.771   0.986    0.986
#> 2 ICU           0.819    0.710   0.962    0.962
#> 3 Outpatient    0.894    0.8     0.991    0.991
# }
```
