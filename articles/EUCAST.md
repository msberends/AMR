# Apply EUCAST rules

## Introduction

What are EUCAST rules? The European Committee on Antimicrobial
Susceptibility Testing (EUCAST) states [on their
website](https://www.eucast.org/expert_rules_and_expected_phenotypes):

> *EUCAST expert rules (see below) are a tabulated collection of expert
> knowledge on interpretive rules, expected resistant phenotypes and
> expected susceptible phenotypes which should be applied to
> antimicrobial susceptibility testing in order to reduce testing,
> reduce errors and make appropriate recommendations for reporting
> particular resistances.*

In Europe, a lot of medical microbiological laboratories already apply
these rules ([Brown *et al.*,
2015](https://www.eurosurveillance.org/content/10.2807/1560-7917.ES2015.20.2.21008)).
Our package features their latest insights on expected resistant
phenotypes (v1.2, 2023).

## Examples

These rules can be used to discard improbable bug-drug combinations in
your data. For example, *Klebsiella* produces beta-lactamase that
prevents ampicillin (or amoxicillin) from working against it. In other
words, practically every strain of *Klebsiella* is resistant to
ampicillin.

Sometimes, laboratory data can still contain such strains with
*Klebsiella* being susceptible to ampicillin. This could be because an
antibiogram is available before an identification is available, and the
antibiogram is then not re-interpreted based on the identification. The
[`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md)
function resolves this, by applying the latest ‘EUCAST Expected
Resistant Phenotypes’ guideline:

``` r
oops <- tibble::tibble(
  mo = c(
    "Klebsiella pneumoniae",
    "Escherichia coli"
  ),
  ampicillin = as.sir("S")
)
oops
#> # A tibble: 2 × 2
#>   mo                    ampicillin
#>   <chr>                 <sir>     
#> 1 Klebsiella pneumoniae   S       
#> 2 Escherichia coli        S  

eucast_rules(oops, info = FALSE, overwrite = TRUE)
#> # A tibble: 2 × 2
#>   mo                    ampicillin
#>   <chr>                 <sir>     
#> 1 Klebsiella pneumoniae   R       
#> 2 Escherichia coli        S  
```

A more convenient function is
[`mo_is_intrinsic_resistant()`](https://amr-for-r.org/reference/mo_property.md)
that uses the same guideline, but allows to check for one or more
specific microorganisms or antimicrobials:

``` r
mo_is_intrinsic_resistant(
  c("Klebsiella pneumoniae", "Escherichia coli"),
  "ampicillin"
)
#> [1]  TRUE FALSE

mo_is_intrinsic_resistant(
  "Klebsiella pneumoniae",
  c("ampicillin", "kanamycin")
)
#> [1]  TRUE FALSE
```

EUCAST rules can not only be used for correction, they can also be used
for filling in known resistance and susceptibility based on results of
other antimicrobials drugs. This process is called *interpretive
reading*, and is basically a form of imputation:

``` r
data <- tibble::tibble(
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
  PEN = "S", # Benzylenicillin
  FOX = "S"  # Cefoxitin
)
```

``` r
data
```

| mo                     | VAN | AMX | COL | CAZ | CXM | PEN | FOX |
|:-----------------------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Staphylococcus aureus  | \-  | \-  | \-  | \-  | \-  |  S  |  S  |
| Enterococcus faecalis  | \-  | \-  | \-  | \-  | \-  |  S  |  S  |
| Escherichia coli       | \-  | \-  | \-  | \-  | \-  |  S  |  S  |
| Klebsiella pneumoniae  | \-  | \-  | \-  | \-  | \-  |  S  |  S  |
| Pseudomonas aeruginosa | \-  | \-  | \-  | \-  | \-  |  S  |  S  |

``` r
eucast_rules(data, overwrite = TRUE)
```

| mo                     | VAN | AMX | COL | CAZ | CXM | PEN | FOX |
|:-----------------------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Staphylococcus aureus  | \-  |  S  |  R  |  R  |  S  |  S  |  S  |
| Enterococcus faecalis  | \-  | \-  |  R  |  R  |  R  |  S  |  R  |
| Escherichia coli       |  R  | \-  | \-  | \-  | \-  |  R  |  S  |
| Klebsiella pneumoniae  |  R  |  R  | \-  | \-  | \-  |  R  |  S  |
| Pseudomonas aeruginosa |  R  |  R  | \-  | \-  |  R  |  R  |  R  |
