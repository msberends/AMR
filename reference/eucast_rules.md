# Apply EUCAST Rules

Apply rules from clinical breakpoints notes and expected resistant
phenotypes as defined by the European Committee on Antimicrobial
Susceptibility Testing (EUCAST, <https://www.eucast.org>), see *Source*.
Use `eucast_dosage()` to get a
[data.frame](https://rdrr.io/r/base/data.frame.html) with advised
dosages of a certain bug-drug combination, which is based on the
[dosage](https://amr-for-r.org/reference/dosage.md) data set.

To improve the interpretation of the antibiogram before EUCAST rules are
applied, some non-EUCAST rules can applied at default, see *Details*.

## Usage

``` r
eucast_rules(x, col_mo = NULL, info = interactive(),
  rules = getOption("AMR_eucastrules", default = c("breakpoints",
  "expected_phenotypes")), verbose = FALSE, version_breakpoints = 15,
  version_expected_phenotypes = 1.2, version_expertrules = 3.3,
  ampc_cephalosporin_resistance = NA, only_sir_columns = any(is.sir(x)),
  custom_rules = NULL, overwrite = FALSE, ...)

eucast_dosage(ab, administration = "iv", version_breakpoints = 15)
```

## Source

- EUCAST Expert Rules. Version 2.0, 2012.  
  Leclercq et al. **EUCAST expert rules in antimicrobial susceptibility
  testing.** *Clin Microbiol Infect.* 2013;19(2):141-60;
  [doi:10.1111/j.1469-0691.2011.03703.x](https://doi.org/10.1111/j.1469-0691.2011.03703.x)

- EUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes
  Tables. Version 3.1, 2016.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf)

- EUCAST Intrinsic Resistance and Unusual Phenotypes. Version 3.2, 2020.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2020/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.2_20200225.pdf)

- EUCAST Intrinsic Resistance and Unusual Phenotypes. Version 3.3, 2021.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2021/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.3_20211018.pdf)

- EUCAST Breakpoint tables for interpretation of MICs and zone
  diameters. Version 9.0, 2019.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_9.0_Breakpoint_Tables.xlsx)

- EUCAST Breakpoint tables for interpretation of MICs and zone
  diameters. Version 10.0, 2020.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_10.0_Breakpoint_Tables.xlsx)

- EUCAST Breakpoint tables for interpretation of MICs and zone
  diameters. Version 11.0, 2021.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_11.0_Breakpoint_Tables.xlsx)

- EUCAST Breakpoint tables for interpretation of MICs and zone
  diameters. Version 12.0, 2022.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_12.0_Breakpoint_Tables.xlsx)

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

- info:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  progress should be printed to the console - the default is only print
  while in interactive sessions.

- rules:

  A [character](https://rdrr.io/r/base/character.html) vector that
  specifies which rules should be applied. Must be one or more of
  `"breakpoints"`, `"expected_phenotypes"`, `"expert"`, `"other"`,
  `"custom"`, `"all"`, and defaults to
  `c("breakpoints", "expected_phenotypes")`. The default value can be
  set to another value using the package option
  [`AMR_eucastrules`](https://amr-for-r.org/reference/AMR-options.md):
  `options(AMR_eucastrules = "all")`. If using `"custom"`, be sure to
  fill in argument `custom_rules` too. Custom rules can be created with
  [`custom_eucast_rules()`](https://amr-for-r.org/reference/custom_eucast_rules.md).

- verbose:

  A [logical](https://rdrr.io/r/base/logical.html) to turn Verbose mode
  on and off (default is off). In Verbose mode, the function does not
  apply rules to the data, but instead returns a data set in logbook
  form with extensive info about which rows and columns would be
  effected and in which way. Using Verbose mode takes a lot more time.

- version_breakpoints:

  The version number to use for the EUCAST Clinical Breakpoints
  guideline. Can be "15.0", "14.0", "13.1", "12.0", "11.0", or "10.0".

- version_expected_phenotypes:

  The version number to use for the EUCAST Expected Phenotypes. Can be
  "1.2".

- version_expertrules:

  The version number to use for the EUCAST Expert Rules and Intrinsic
  Resistance guideline. Can be "3.3", "3.2", or "3.1".

- ampc_cephalosporin_resistance:

  (only applies when `rules` contains `"expert"` or `"all"`) a
  [character](https://rdrr.io/r/base/character.html) value that should
  be applied to cefotaxime, ceftriaxone and ceftazidime for AmpC
  de-repressed cephalosporin-resistant mutants - the default is `NA`.
  Currently only works when `version_expertrules` is `3.2` and higher;
  these versions of '*EUCAST Expert Rules on Enterobacterales*' state
  that results of cefotaxime, ceftriaxone and ceftazidime should be
  reported with a note, or results should be suppressed (emptied) for
  these three drugs. A value of `NA` (the default) for this argument
  will remove results for these three drugs, while e.g. a value of `"R"`
  will make the results for these drugs resistant. Use `NULL` or `FALSE`
  to not alter results for these three drugs of AmpC de-repressed
  cephalosporin-resistant mutants. Using `TRUE` is equal to using
  `"R"`.  
  For *EUCAST Expert Rules* v3.2, this rule applies to: *Citrobacter
  braakii*, *Citrobacter freundii*, *Citrobacter gillenii*, *Citrobacter
  murliniae*, *Citrobacter rodenticum*, *Citrobacter sedlakii*,
  *Citrobacter werkmanii*, *Citrobacter youngae*, *Enterobacter*,
  *Hafnia alvei*, *Klebsiella aerogenes*, *Morganella morganii*,
  *Providencia*, and *Serratia*.

- only_sir_columns:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only antimicrobial columns must be included that were transformed to
  class [sir](https://amr-for-r.org/reference/as.sir.md) on beforehand.
  Defaults to `FALSE` if no columns of `x` have a class
  [sir](https://amr-for-r.org/reference/as.sir.md).

- custom_rules:

  Custom rules to apply, created with
  [`custom_eucast_rules()`](https://amr-for-r.org/reference/custom_eucast_rules.md).

- overwrite:

  A [logical](https://rdrr.io/r/base/logical.html) indicating whether to
  overwrite existing SIR values (default: `FALSE`). When `FALSE`, only
  non-SIR values are modified (i.e., any value that is not already S, I
  or R). To ensure compliance with EUCAST guidelines, **this should
  remain** `FALSE`, as EUCAST notes often state that an organism "should
  be tested for susceptibility to individual agents or be reported
  resistant".

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

- ab:

  Any (vector of) text that can be coerced to a valid antimicrobial drug
  code with [`as.ab()`](https://amr-for-r.org/reference/as.ab.md).

- administration:

  Route of administration, either "", "im", "iv", or "oral".

## Value

The input of `x`, possibly with edited values of antimicrobials. Or, if
`verbose = TRUE`, a [data.frame](https://rdrr.io/r/base/data.frame.html)
with all original and new values of the affected bug-drug combinations.

## Details

**Note:** This function does not translate MIC values to SIR values. Use
[`as.sir()`](https://amr-for-r.org/reference/as.sir.md) for that.  
**Note:** When ampicillin (AMP, J01CA01) is not available but
amoxicillin (AMX, J01CA04) is, the latter will be used for all rules
where there is a dependency on ampicillin. These drugs are
interchangeable when it comes to expression of antimicrobial
resistance.  

The file containing all EUCAST rules is located here:
<https://github.com/msberends/AMR/blob/main/data-raw/eucast_rules.tsv>.
**Note:** Old taxonomic names are replaced with the current taxonomy
where applicable. For example, *Ochrobactrum anthropi* was renamed to
*Brucella anthropi* in 2020; the original EUCAST rules v3.1 and v3.2 did
not yet contain this new taxonomic name. The `AMR` package contains the
full microbial taxonomy updated until June 24th, 2024, see
[microorganisms](https://amr-for-r.org/reference/microorganisms.md).

### Custom Rules

Custom rules can be created using
[`custom_eucast_rules()`](https://amr-for-r.org/reference/custom_eucast_rules.md),
e.g.:

    x <- custom_eucast_rules(AMC == "R" & genus == "Klebsiella" ~ aminopenicillins == "R",
                             AMC == "I" & genus == "Klebsiella" ~ aminopenicillins == "I")

    eucast_rules(example_isolates, rules = "custom", custom_rules = x)

### 'Other' Rules

Before further processing, two non-EUCAST rules about drug combinations
can be applied to improve the efficacy of the EUCAST rules, and the
reliability of your data (analysis). These rules are:

1.  A drug **with** enzyme inhibitor will be set to S if the same drug
    **without** enzyme inhibitor is S

2.  A drug **without** enzyme inhibitor will be set to R if the same
    drug **with** enzyme inhibitor is R

Important examples include amoxicillin and amoxicillin/clavulanic acid,
and trimethoprim and trimethoprim/sulfamethoxazole. Needless to say, for
these rules to work, both drugs must be available in the data set.

Since these rules are not officially approved by EUCAST, they are not
applied at default. To use these rules, include `"other"` to the `rules`
argument, or use `eucast_rules(..., rules = "all")`. You can also set
the package option
[`AMR_eucastrules`](https://amr-for-r.org/reference/AMR-options.md),
i.e. run `options(AMR_eucastrules = "all")`.

## Download Our Reference Data

All reference data sets in the AMR package - including information on
microorganisms, antimicrobials, and clinical breakpoints - are freely
available for download in multiple formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, and Stata.

For maximum compatibility, we also provide machine-readable,
tab-separated plain text files suitable for use in any software,
including laboratory information systems.

Visit [our website for direct download
links](https://amr-for-r.org/articles/datasets.html), or explore the
actual files in [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw/datasets).

## Examples

``` r
# \donttest{
a <- data.frame(
  mo = c(
    "Staphylococcus aureus",
    "Enterococcus faecalis",
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Pseudomonas aeruginosa"
  ),
  VAN = "-", # Vancomycin
  AMX = "-", # Amoxicillin
  COL = "-", # Colistin
  CAZ = "-", # Ceftazidime
  CXM = "-", # Cefuroxime
  PEN = "S", # Benzylpenicillin
  FOX = "S", # Cefoxitin
  stringsAsFactors = FALSE
)

head(a)
#>                       mo VAN AMX COL CAZ CXM PEN FOX
#> 1  Staphylococcus aureus   -   -   -   -   -   S   S
#> 2  Enterococcus faecalis   -   -   -   -   -   S   S
#> 3       Escherichia coli   -   -   -   -   -   S   S
#> 4  Klebsiella pneumoniae   -   -   -   -   -   S   S
#> 5 Pseudomonas aeruginosa   -   -   -   -   -   S   S


# apply EUCAST rules: some results wil be changed
b <- eucast_rules(a, overwrite = TRUE)
#> Warning: in `eucast_rules()`: not all columns with antimicrobial results are of
#> class 'sir'. Transform them on beforehand, with e.g.:
#>   - a %>% as.sir(CXM:AMX)
#>   - a %>% mutate_if(is_sir_eligible, as.sir)
#>   - a %>% mutate(across(where(is_sir_eligible), as.sir))

head(b)
#>                       mo VAN AMX COL CAZ CXM PEN FOX
#> 1  Staphylococcus aureus   -   S   R   R   S   S   S
#> 2  Enterococcus faecalis   -   -   R   R   R   S   R
#> 3       Escherichia coli   R   -   -   -   -   R   S
#> 4  Klebsiella pneumoniae   R   R   -   -   -   R   S
#> 5 Pseudomonas aeruginosa   R   R   -   -   R   R   R


# do not apply EUCAST rules, but rather get a data.frame
# containing all details about the transformations:
c <- eucast_rules(a, overwrite = TRUE, verbose = TRUE)
#> Warning: in `eucast_rules()`: not all columns with antimicrobial results are of
#> class 'sir'. Transform them on beforehand, with e.g.:
#>   - a %>% as.sir(CXM:AMX)
#>   - a %>% mutate_if(is_sir_eligible, as.sir)
#>   - a %>% mutate(across(where(is_sir_eligible), as.sir))
head(c)
#>   row col           mo_fullname old new rule          rule_group
#> 1   1 AMX Staphylococcus aureus   -   S              Breakpoints
#> 2   1 CXM Staphylococcus aureus   -   S              Breakpoints
#> 3   1 CAZ Staphylococcus aureus   -   R      Expected phenotypes
#> 4   1 COL Staphylococcus aureus   -   R      Expected phenotypes
#> 5   2 CAZ Enterococcus faecalis   -   R      Expected phenotypes
#> 6   2 COL Enterococcus faecalis   -   R      Expected phenotypes
#>                                                         rule_name
#> 1                                                  Staphylococcus
#> 2                                                  Staphylococcus
#> 3 Table 4: Expected resistant phenotype in gram-positive bacteria
#> 4 Table 4: Expected resistant phenotype in gram-positive bacteria
#> 5 Table 4: Expected resistant phenotype in gram-positive bacteria
#> 6 Table 4: Expected resistant phenotype in gram-positive bacteria
#>                                         rule_source
#> 1   'EUCAST Clinical Breakpoint Tables' v15.0, 2025
#> 2   'EUCAST Clinical Breakpoint Tables' v15.0, 2025
#> 3 'EUCAST Expected Resistant Phenotypes' v1.2, 2023
#> 4 'EUCAST Expected Resistant Phenotypes' v1.2, 2023
#> 5 'EUCAST Expected Resistant Phenotypes' v1.2, 2023
#> 6 'EUCAST Expected Resistant Phenotypes' v1.2, 2023
# }

# Dosage guidelines:

eucast_dosage(c("tobra", "genta", "cipro"), "iv")
#> ℹ Dosages for antimicrobial drugs, as meant for 'EUCAST Clinical Breakpoint
#>   Tables' v15.0 (2025). This note will be shown once per session.
#> # A tibble: 3 × 5
#>   ab   name          standard_dosage  high_dosage  eucast_version
#>   <ab> <chr>         <chr>            <chr>                 <dbl>
#> 1 TOB  Tobramycin    6-7 mg/kg x 1 iv NA                       15
#> 2 GEN  Gentamicin    6-7 mg/kg x 1 iv NA                       15
#> 3 CIP  Ciprofloxacin 0.4 g x 2 iv     0.4 g x 3 iv             15

eucast_dosage(c("tobra", "genta", "cipro"), "iv", version_breakpoints = 10)
#> # A tibble: 3 × 5
#>   ab   name          standard_dosage high_dosage eucast_version
#>   <ab> <chr>         <chr>           <chr>                <dbl>
#> 1 TOB  Tobramycin    NA              NA                      NA
#> 2 GEN  Gentamicin    NA              NA                      NA
#> 3 CIP  Ciprofloxacin NA              NA                      NA
```
