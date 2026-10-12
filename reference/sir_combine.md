# Combine SIR Results of Multiple Antimicrobials

Combine the SIR results of two or more antimicrobials into one SIR
result per isolate. This can answer two different questions, set with
the `method` argument:

- `"most_susceptible"` (default): *Is at least one of these
  antimicrobials expected to work?* This is the result of a combination
  therapy, e.g. piperacillin/tazobactam with gentamicin, which covers
  the isolate if any of its antimicrobials is active.

- `"most_resistant"`: *Is resistance demonstrated to any of these
  antimicrobials?* This is the result for a group of antimicrobials,
  e.g. resistance to third-generation cephalosporins if an isolate is
  resistant to ceftazidime or ceftriaxone, as used in surveillance
  reports.

See *Combination Therapy* below for the exact rules. The same rules are
used by
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md) for
combination antibiograms, and by the
[`count_*()`](https://amr-for-r.org/reference/count.md) and
[`proportion_*()`](https://amr-for-r.org/reference/proportion.md)
functions when more than one variable is used.

## Usage

``` r
sir_combine(
  ...,
  method = c("most_susceptible", "most_resistant"),
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

- method:

  Either `"most_susceptible"` (default) to return the most susceptible
  result, as for combination therapy, or `"most_resistant"` to return
  the most resistant result, as for resistance within a group of
  antimicrobials. See *Combination Therapy* below.

- only_all_tested:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate that
  isolates must be tested for all antimicrobials. If `TRUE`, the result
  is `NA` for isolates with a missing result for any of the
  antimicrobials. If `FALSE` (default), a missing result only leads to
  `NA` if a result for the missing antimicrobial could have changed the
  result, see *Combination Therapy* below.

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

Results are ordered from susceptible to resistant: S (and WT), SDD, I, R
(and NWT, NS). With `method = "most_susceptible"` (in `sir_combine()`),
the most susceptible result of the antimicrobials is returned, and with
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

The [`count_*()`](https://amr-for-r.org/reference/count.md) and
[`proportion_*()`](https://amr-for-r.org/reference/proportion.md)
functions apply these rules when more than one variable is used. With
`only_all_tested = FALSE` (default), the functions for susceptibility,
such as
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md) and
[`count_susceptible()`](https://amr-for-r.org/reference/count.md), count
an isolate as susceptible if any of the antimicrobials is susceptible
(`"most_susceptible"`), and the functions for resistance, such as
[`resistance()`](https://amr-for-r.org/reference/proportion.md) and
[`count_resistant()`](https://amr-for-r.org/reference/count.md), count
an isolate as resistant if any of the antimicrobials is resistant
(`"most_resistant"`). An isolate is considered tested if it has such a
result, or if it was tested for all antimicrobials. Therefore, for
combinations with `only_all_tested = FALSE`:

        count_S()    +   count_I()    +   count_R()    >= count_all()
      proportion_S() + proportion_I() + proportion_R() >= 1

With `only_all_tested = TRUE`, all these functions count the result of
the combination therapy (`"most_susceptible"`) of isolates tested for
all antimicrobials, so that:

        count_S()    +   count_I()    +   count_R()    = count_all()
      proportion_S() + proportion_I() + proportion_R() = 1

In [`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md), the result
of the combination therapy (`"most_susceptible"`) is used, the 'I' and
'SDD' categories are considered susceptible if `combine_SI = TRUE`
(default), and the result is then reduced to S or R.

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

# combination therapy: is at least one of them expected to work?
sir_combine(drug_a, drug_b)
#> Class <sir>
#> [1] S    S    R    I    I    <NA> <NA> SDD 
sir_combine(drug_a, drug_b, guideline = "CLSI")
#> Class <sir>
#> [1] S    S    R    I    <NA> <NA> <NA> SDD 
sir_combine(drug_a, drug_b, only_all_tested = TRUE)
#> Class <sir>
#> [1] S    S    R    I    <NA> <NA> <NA> SDD 

# resistance within a group: is resistance demonstrated to any of them?
sir_combine(drug_a, drug_b, method = "most_resistant")
#> Class <sir>
#> [1] R    R    R    R    <NA> <NA> R    I   

# the count and proportion functions use the same rules
susceptibility(drug_a, drug_b, minimum = 0)
#> Note that `proportion_SI()` will also include dose-dependent susceptibility,
#> "SDD". This note will be shown once for this session.
#> [1] 0.8333333
susceptibility(sir_combine(drug_a, drug_b), minimum = 0)
#> [1] 0.8333333
resistance(drug_a, drug_b, minimum = 0)
#> [1] 0.8333333
resistance(sir_combine(drug_a, drug_b, method = "most_resistant"), minimum = 0)
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
  # resistance to third-generation cephalosporins
  example_isolates %>%
    mutate(cephalosporins_3rd = sir_combine(CAZ, CRO, method = "most_resistant")) %>%
    count(cephalosporins_3rd)
}
#> # A tibble: 4 × 2
#>   cephalosporins_3rd     n
#>   <sir>              <int>
#> 1   S                  573
#> 2   I                    1
#> 3   R                 1236
#> 4   NA                 190
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
