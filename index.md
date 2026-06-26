# The `AMR` Package for R

- Provides an **all-in-one solution** for antimicrobial resistance (AMR)
  data analysis in a One Health approach
- **Peer-reviewed**, used in over 175 countries, available in 28
  languages
- Generates **antibiograms** - WISCA for empiric coverage estimates, or
  traditional/syndromic for AMR surveillance
- Provides the **full microbiological taxonomy** of ~97 000 distinct
  species and extensive info of ~620 antimicrobial drugs
- Applies **CLSI 2011-2026** and **EUCAST 2011-2026** clinical and
  veterinary breakpoints, and ECOFFs, for MIC and disk zone
  interpretation
- Corrects for duplicate isolates, **calculates** and **predicts** AMR
  per antimicrobial class
- Integrates with **WHONET**, ATC, **EARS-Net**, PubChem, **LOINC**,
  **SNOMED CT**, and **NCBI**
- 100% free of costs and dependencies, highly suitable for places with
  **limited resources**

> Available for Python too! [Click
> here](https://amr-for-r.org/articles/AMR_for_Python.md) to read more.

[amr-for-r.org](https://amr-for-r.org/)

[doi.org/10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)

[![](./endorsement_clsi_eucast.jpg)](https://amr-for-r.org/reference/clinical_breakpoints.html#response-from-clsi-and-eucast)

------------------------------------------------------------------------

## Introduction

The `AMR` package is a peer-reviewed, [free and open-source](#copyright)
R package with [zero
dependencies](https://en.wikipedia.org/wiki/Dependency_hell) to simplify
the analysis and prediction of Antimicrobial Resistance (AMR) and to
work with microbial and antimicrobial data and properties, by using
evidence-based methods. **Our aim is to provide a standard** for clean
and reproducible AMR data analysis, that can therefore empower
epidemiological analyses to continuously enable surveillance and
treatment evaluation in any setting. We are a team of [many different
researchers](https://amr-for-r.org/authors.md) from around the globe to
make this a successful and durable project!

This work was published in the Journal of Statistical Software (Volume
104(3); [DOI
10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)) and
formed the basis of two PhD theses ([DOI
10.33612/diss.177417131](https://doi.org/10.33612/diss.177417131) and
[DOI 10.33612/diss.192486375](https://doi.org/10.33612/diss.192486375)).

After installing this package, R knows [**~97 000 distinct microbial
species**](https://amr-for-r.org/reference/microorganisms.md) (updated
mei 2026) and all [**~620 antimicrobial and antiviral
drugs**](https://amr-for-r.org/reference/antimicrobials.md) by name and
code (including ATC, EARS-Net, ASIARS-Net, PubChem, LOINC and SNOMED
CT), and knows all about valid SIR and MIC values. The integral clinical
breakpoint guidelines from CLSI 2011-2026 and EUCAST 2011-2026 are
included, even with epidemiological cut-off (ECOFF) values. It supports
and can read any data format, including WHONET data. This package works
on Windows, macOS and Linux with all versions of R since R-3.0 (April
2013). **It was designed to work in any setting, including those with
very limited resources**. It was created for both routine data analysis
and academic research at the Faculty of Medical Sciences of the
[University of Groningen](https://www.rug.nl) and the [University
Medical Center Groningen](https://www.umcg.nl).

### Used in over 175 countries, available in 28 languages

[![](./countries.png)](https://amr-for-r.org/countries_large.png)

Since its first public release in early 2018, this R package has been
used in almost all countries in the world. Click the map to enlarge and
to see the country names.

With the help of contributors from all corners of the world, the `AMR`
package is available in ![](lang_en.svg) English, ![](lang_ar.svg)
Arabic, ![](lang_bn.svg) Bengali, ![](lang_zh.svg) Chinese,
![](lang_cs.svg) Czech, ![](lang_da.svg) Danish, ![](lang_nl.svg) Dutch,
![](lang_fi.svg) Finnish, ![](lang_fr.svg) French, ![](lang_de.svg)
German, ![](lang_el.svg) Greek, ![](lang_hi.svg) Hindi, ![](lang_id.svg)
Indonesian, ![](lang_it.svg) Italian, ![](lang_ja.svg) Japanese,
![](lang_ko.svg) Korean, ![](lang_no.svg) Norwegian, ![](lang_pl.svg)
Polish, ![](lang_pt.svg) Portuguese, ![](lang_ro.svg) Romanian,
![](lang_ru.svg) Russian, ![](lang_es.svg) Spanish, ![](lang_sw.svg)
Swahili, ![](lang_sv.svg) Swedish, ![](lang_tr.svg) Turkish,
![](lang_uk.svg) Ukrainian, ![](lang_ur.svg) Urdu, and ![](lang_vi.svg)
Vietnamese. Antimicrobial drug (group) names and colloquial
microorganism names are provided in these languages.

## Practical examples

### Filtering and selecting data

One of the most powerful functions of this package, aside from
calculating and plotting AMR, is selecting and filtering based on
antimicrobial columns. This can be done using the so-called
[antimicrobial
selectors](https://amr-for-r.org/reference/antimicrobial_selectors.html),
which work in base R, `dplyr` and `data.table`.

``` r

# AMR works great with dplyr, but it's not required or neccesary
library(AMR)
library(dplyr, warn.conflicts = FALSE)

example_isolates %>%
  mutate(bacteria = mo_fullname()) %>%
  # filtering functions for microorganisms:
  filter(mo_is_gram_negative(),
         mo_is_intrinsic_resistant(ab = "cefotax")) %>%
  # antimicrobial selectors:
  select(bacteria,
         aminoglycosides(),
         carbapenems())
#> ℹ Using column mo as input for `mo_fullname()`
#> ℹ Using column mo as input for `mo_is_gram_negative()`
#> ℹ Using column mo as input for `mo_is_intrinsic_resistant()`
#> ℹ Determining intrinsic resistance based on 'EUCAST Expected
#>   Resistant Phenotypes' v1.2 (2023). This note will be shown
#>   once per session.
#> ℹ For `aminoglycosides()` using columns GEN (gentamicin), TOB
#>   (tobramycin), AMK (amikacin), and KAN (kanamycin)
#> ℹ For `carbapenems()` using columns IPM (imipenem) and MEM
#>   (meropenem)
#> # A tibble: 35 × 7
#>    bacteria                     GEN   TOB   AMK   KAN   IPM   MEM  
#>    <chr>                        <sir> <sir> <sir> <sir> <sir> <sir>
#>  1 Pseudomonas aeruginosa       I     S     NA    R     S     NA   
#>  2 Pseudomonas aeruginosa       I     S     NA    R     S     NA   
#>  3 Pseudomonas aeruginosa       I     S     NA    R     S     NA   
#>  4 Pseudomonas aeruginosa       S     S     S     R     NA    S    
#>  5 Pseudomonas aeruginosa       S     S     S     R     S     S    
#>  6 Pseudomonas aeruginosa       S     S     S     R     S     S    
#>  7 Stenotrophomonas maltophilia R     R     R     R     R     R    
#>  8 Pseudomonas aeruginosa       S     S     S     R     NA    S    
#>  9 Pseudomonas aeruginosa       S     S     S     R     NA    S    
#> 10 Pseudomonas aeruginosa       S     S     S     R     S     S    
#> # ℹ 25 more rows
```

With only having defined a row filter on Gram-negative bacteria with
intrinsic resistance to cefotaxime
([`mo_is_gram_negative()`](https://amr-for-r.org/reference/mo_property.md)
and
[`mo_is_intrinsic_resistant()`](https://amr-for-r.org/reference/mo_property.md))
and a column selection on two antibiotic groups
([`aminoglycosides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
and
[`carbapenems()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)),
the reference data about [all
microorganisms](https://amr-for-r.org/reference/microorganisms.md) and
[all antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
in the `AMR` package make sure you get what you meant.

### Generating antibiograms

The `AMR` package supports four types of antibiograms, with support for
28 languages. If used inside [R Markdown](https://rmarkdown.rstudio.com)
or [Quarto](https://quarto.org), the table will be printed in the right
output format automatically (such as markdown, LaTeX, HTML, etc.).

**For empirical therapy guidance (i.e., coverage estimates), use WISCA**
(Weighted-Incidence Syndromic Combination Antibiogram). When a clinician
starts empirical treatment, the causative pathogen is unknown. The
relevant question is not *“what percentage of E. coli is susceptible?”*
but *“what is the probability that this regimen will cover whatever
pathogen is causing the infection?”*. WISCA answers that question
directly, weighting susceptibility by pathogen incidence and providing
credible intervals via Bayesian simulation. See
[`vignette("WISCA")`](https://amr-for-r.org/articles/WISCA.md) for the
full explanation.

``` r

wisca(example_isolates,
      antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
      minimum = 10) # Recommended threshold: >=30
#> Warning: invalid microorganism code, NA generated
```

| Piperacillin/tazobactam | Piperacillin/tazobactam + Gentamicin | Piperacillin/tazobactam + Tobramycin |
|:---|:---|:---|
| 70% (64.7-75.2%) | 93.6% (92.2-95.1%) | 89.8% (87-92.5%) |

WISCA supports stratification by any clinical variable, so you can
generate syndrome-specific or ward-specific coverage estimates:

``` r

wisca(example_isolates,
      antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"),
      syndromic_group = "ward",
      minimum = 10) # Recommended threshold: >=30
#> Warning: invalid microorganism code, NA generated
```

| Syndromic Group | Piperacillin/tazobactam | Piperacillin/tazobactam + Gentamicin | Piperacillin/tazobactam + Tobramycin |
|:---|:---|:---|:---|
| Clinical | 74.6% (68.6-80.6%) | 93.7% (92.1-95.1%) | 90.4% (87-93.1%) |
| ICU | 57% (48.6-65.7%) | 86.8% (83.6-89.8%) | 82.9% (78.1-87.3%) |
| Outpatient | 56.9% (45.9-68.2%) | 76.7% (70.6-82.3%) | 68% (57.6-77.2%) |

**For AMR surveillance**, traditional antibiograms remain the right tool
for tracking resistance per species over time:

``` r

antibiogram(example_isolates,
            mo_transform = "gramstain",
            antimicrobials = c("AMC", carbapenems(), "TZP"))
#> ℹ For `carbapenems()` using columns IPM (imipenem) and MEM
#>   (meropenem)
```

| Pathogen | Amoxicillin/clavulanic acid | Imipenem | Meropenem | Piperacillin/tazobactam |
|:---|:---|:---|:---|:---|
| Gram-negative | 76% (73-79%,N=726) | 99% (98-100%,N=631) | 100% (99-100%,N=626) | 88% (85-91%,N=641) |
| Gram-positive | 76% (74-79%,N=1138) | 81% (75-85%,N=257) | 77% (70-82%,N=203) | 86% (82-89%,N=345) |

Combination antibiograms show the additional coverage gained by adding a
second agent, stratified by species:

``` r

antibiogram(example_isolates,
            mo_transform = "gramstain",
            antimicrobials = c("TZP", "TZP+TOB", "TZP+GEN"))
```

| Pathogen | Piperacillin/tazobactam | Piperacillin/tazobactam + Gentamicin | Piperacillin/tazobactam + Tobramycin |
|:---|:---|:---|:---|
| Gram-negative | 88% (85-91%,N=641) | 99% (97-99%,N=691) | 98% (97-99%,N=693) |
| Gram-positive | 86% (82-89%,N=345) | 98% (96-98%,N=1044) | 95% (93-97%,N=550) |

Like many other functions in this package,
[`antibiogram()`](https://amr-for-r.org/reference/antibiogram.md) and
[`wisca()`](https://amr-for-r.org/reference/antibiogram.md) come with
support for 28 languages that are often detected automatically based on
system language:

``` r

antibiogram(example_isolates,
            antimicrobials = c("cipro", "tobra", "genta"), # any arbitrary name or code will work
            mo_transform = "gramstain",
            ab_transform = "name",
            language = "uk") # Ukrainian
```

| Збудник       | Гентаміцин          | Тобраміцин         | Ципрофлоксацин     |
|:--------------|:--------------------|:-------------------|:-------------------|
| Грамнегативні | 96% (95-98%,N=684)  | 96% (94-97%,N=686) | 91% (88-93%,N=684) |
| Грампозитивні | 63% (60-66%,N=1170) | 34% (31-38%,N=665) | 77% (74-80%,N=724) |

### Interpreting and plotting MIC and SIR values

The `AMR` package allows interpretation of MIC and disk diffusion values
based on CLSI and EUCAST. Moreover, the `ggplot2` package is extended
with new scale functions, to allow plotting of log2-distributed MIC
values and SIR values.

``` r

library(ggplot2)
library(AMR)

# generate some random values
some_mic_values <- random_mic(size = 100)
some_groups <- sample(LETTERS[1:5], 20, replace = TRUE)
interpretation <- as.sir(some_mic_values,
                         guideline = "EUCAST 2024",
                         mo = "E. coli", # or any code or name resembling a known species
                         ab = "Cipro")   # or any code or name resembling an antibiotic

# create the plot
ggplot(data.frame(mic = some_mic_values,
                  group = some_groups,
                  sir = interpretation),
       aes(x = group, y = mic, colour = sir)) +
  theme_minimal() +
  geom_boxplot(fill = NA, colour = "grey30") +
  geom_jitter(width = 0.25) +
  
  # NEW scale function: plot MIC values to x, y, colour or fill
  scale_y_mic() +
  
  # NEW scale function: write out S/I/R in any of the 20 supported languages
  #                     and set colourblind-friendly colours
  scale_colour_sir()
```

[![](./plot_readme.png)](https://amr-for-r.org/reference/plotting.md "Plotting Helpers for AMR Data Analysis")

### Calculating resistance per group

For a manual approach, you can use the `resistance` or
[`susceptibility()`](https://amr-for-r.org/reference/proportion.md)
function:

``` r

example_isolates %>%
  # group by ward:
  group_by(ward) %>%
  # calculate AMR using resistance() for gentamicin and tobramycin
  # and get their 95% confidence intervals using sir_confidence_interval():
  summarise(across(c(GEN, TOB),
                   list(total_R = resistance,
                        conf_int = function(x) sir_confidence_interval(x, collapse = "-"))))
#> ℹ `resistance()` assumes the EUCAST guideline and thus
#>   considers the 'I' category susceptible. Set the `guideline`
#>   argument or the `AMR_guideline` option to either "CLSI" or
#>   "EUCAST", see `?AMR-options`.
#> ℹ This message will be shown once per session.
#> # A tibble: 3 × 5
#>   ward       GEN_total_R GEN_conf_int TOB_total_R TOB_conf_int
#>   <chr>            <dbl> <chr>              <dbl> <chr>       
#> 1 Clinical         0.229 0.205-0.254        0.315 0.284-0.347 
#> 2 ICU              0.290 0.253-0.33         0.400 0.353-0.449 
#> 3 Outpatient       0.2   0.131-0.285        0.368 0.254-0.493
```

Or use [antimicrobial
selectors](https://amr-for-r.org/reference/antimicrobial_selectors.html)
to select a series of antibiotic columns:

``` r

library(AMR)
library(dplyr)

out <- example_isolates %>%
  # group by ward:
  group_by(ward) %>%
  # calculate AMR using resistance(), over all aminoglycosides and polymyxins:
  summarise(across(c(aminoglycosides(), polymyxins()),
            resistance))
#> ℹ For `aminoglycosides()` using columns GEN (gentamicin), TOB
#>   (tobramycin), AMK (amikacin), and KAN (kanamycin)
#> ℹ For `polymyxins()` using column COL (colistin)
#> Warning: There was 1 warning in `summarise()`.
#> ℹ In argument: `across(c(aminoglycosides(), polymyxins()),
#>   resistance)`.
#> ℹ In group 3: `ward = "Outpatient"`.
#> Caused by warning:
#> ! Introducing NA: only 23 results available for KAN in group:
#> ward = "Outpatient" (whilst `minimum = 30`).
out
#> # A tibble: 3 × 6
#>   ward         GEN   TOB   AMK   KAN   COL
#>   <chr>      <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 Clinical   0.229 0.315 0.626     1 0.780
#> 2 ICU        0.290 0.400 0.662     1 0.857
#> 3 Outpatient 0.2   0.368 0.605    NA 0.889
```

``` r

# transform the antibiotic columns to names:
out %>% set_ab_names()
#> # A tibble: 3 × 6
#>   ward       gentamicin tobramycin amikacin kanamycin colistin
#>   <chr>           <dbl>      <dbl>    <dbl>     <dbl>    <dbl>
#> 1 Clinical        0.229      0.315    0.626         1    0.780
#> 2 ICU             0.290      0.400    0.662         1    0.857
#> 3 Outpatient      0.2        0.368    0.605        NA    0.889
```

``` r

# transform the antibiotic column to ATC codes:
out %>% set_ab_names(property = "atc")
#> # A tibble: 3 × 6
#>   ward       J01GB03 J01GB01 J01GB06 J01GB04 J01XB01
#>   <chr>        <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1 Clinical     0.229   0.315   0.626       1   0.780
#> 2 ICU          0.290   0.400   0.662       1   0.857
#> 3 Outpatient   0.2     0.368   0.605      NA   0.889
```

## What else can you do with this package?

This package was intended as a comprehensive toolbox for integrated AMR
data analysis. This package can be used for:

- Reference for the taxonomy of microorganisms, since the package
  contains all microbial (sub)species from the List of Prokaryotic names
  with Standing in Nomenclature ([LPSN](https://lpsn.dsmz.de)) and the
  Global Biodiversity Information Facility
  ([GBIF](https://www.gbif.org))
  ([manual](https://amr-for-r.org/reference/mo_property.md))
- Interpreting raw MIC and disk diffusion values, based on any CLSI or
  EUCAST guideline ([manual](https://amr-for-r.org/reference/as.sir.md))
- Retrieving antimicrobial drug names, doses and forms of administration
  from clinical health care records
  ([manual](https://amr-for-r.org/reference/ab_from_text.md))
- Determining first isolates to be used for AMR data analysis
  ([manual](https://amr-for-r.org/reference/first_isolate.md))
- Calculating antimicrobial resistance
  ([tutorial](https://amr-for-r.org/articles/AMR.md))
- Determining multi-drug resistance (MDR) / multi-drug resistant
  organisms (MDRO) ([tutorial](https://amr-for-r.org/reference/mdro.md))
- Calculating (empirical) susceptibility of both mono therapy and
  combination therapies
  ([tutorial](https://amr-for-r.org/articles/AMR.md))
- Apply AMR functions in predictive modelling
  ([tutorial](https://amr-for-r.org/articles/AMR_with_tidymodels.md))
- Getting properties for any microorganism (like Gram stain, species,
  genus or family)
  ([manual](https://amr-for-r.org/reference/mo_property.md))
- Getting properties for any antimicrobial (like name, code of
  EARS-Net/ATC/LOINC/PubChem, defined daily dose or trade name)
  ([manual](https://amr-for-r.org/reference/ab_property.md))
- Plotting antimicrobial resistance
  ([tutorial](https://amr-for-r.org/articles/AMR.md))
- Applying EUCAST expert rules
  ([manual](https://amr-for-r.org/reference/eucast_rules.md))
- Getting SNOMED codes of a microorganism, or getting properties of a
  microorganism based on a SNOMED code
  ([manual](https://amr-for-r.org/reference/mo_property.md))
- Getting LOINC codes of an antibiotic, or getting properties of an
  antibiotic based on a LOINC code
  ([manual](https://amr-for-r.org/reference/ab_property.md))
- Machine reading the EUCAST and CLSI guidelines from 2011-2021 to
  translate MIC values and disk diffusion diameters to SIR
  ([link](https://amr-for-r.org/articles/datasets.md))
- Principal component analysis for AMR
  ([tutorial](https://amr-for-r.org/articles/PCA.md))

## Get this package

### Latest official version

[![CRAN](https://www.r-pkg.org/badges/version-ago/AMR)](https://cran.r-project.org/package=AMR)
[![CRANlogs](https://cranlogs.r-pkg.org/badges/grand-total/AMR)](https://cran.r-project.org/package=AMR)

This package is available [here on the official R network
(CRAN)](https://cran.r-project.org/package=AMR). Install this package in
R from CRAN by using the command:

``` r

install.packages("AMR")
```

It will be downloaded and installed automatically. For RStudio, click on
the menu *Tools* \> *Install Packages…* and then type in “AMR” and press
Install.

**Note:** Not all functions on this website may be available in this
latest release. To use all functions and data sets mentioned on this
website, install the latest beta version.

### Latest beta version

[![check-old](https://github.com/msberends/AMR/actions/workflows/check-old-tinytest.yaml/badge.svg?branch=main)](https://github.com/msberends/AMR/actions/workflows/check-old-tinytest.yaml)
[![check-recent](https://github.com/msberends/AMR/actions/workflows/check-current-testthat.yaml/badge.svg?branch=main)](https://github.com/msberends/AMR/actions/workflows/check-current-testthat.yaml)
[![CodeFactor](https://www.codefactor.io/repository/github/msberends/amr/badge)](https://www.codefactor.io/repository/github/msberends/amr)
[![Codecov](https://codecov.io/gh/msberends/AMR/branch/main/graph/badge.svg)](https://codecov.io/gh/msberends/AMR?branch=main)

Please read our [Developer Guideline
here](https://github.com/msberends/AMR/wiki/Developer-Guideline).

To install the latest and unpublished beta version:

``` r

install.packages("AMR", repos = "beta.amr-for-r.org")

# if this does not work, try to install directly from GitHub using the 'remotes' package:
remotes::install_github("msberends/AMR")
```

## Get started

To find out how to conduct AMR data analysis, please [continue reading
here to get started](https://amr-for-r.org/articles/AMR.md) or click a
link in the [‘How to’ menu](https://amr-for-r.org/articles/).

## Partners

The initial development of this package was part of, related to, or made
possible by the following non-profit organisations and initiatives:

[![](./logo_rug.svg)](https://www.rug.nl "University of Groningen")
[![](./logo_umcg.svg)](https://www.umcg.nl "University Medical Center Groningen")
[![](./logo_certe.svg)](https://www.certe.nl "Certe Medical Diagnostics and Advice Foundation")
[![](./logo_eh1h.png)](https://www.deutschland-nederland.eu "EurHealth-1-Health")
[![](./logo_interreg.png)](https://www.deutschland-nederland.eu "INTERREG")

## Copyright

This R package is free, open-source software and licensed under the [GNU
General Public License v2.0
(GPL-2)](https://amr-for-r.org/LICENSE-text.md). In a nutshell, this
means that this package:

- May be used for commercial purposes

- May be used for private purposes

- May **not** be used for patent purposes

- May be modified, although:

  - Modifications **must** be released under the same license when
    distributing the package
  - Changes made to the code **must** be documented

- May be distributed, although:

  - Source code **must** be made available when the package is
    distributed
  - A copy of the license and copyright notice **must** be included with
    the package.

- Comes with a LIMITATION of liability

- Comes with NO warranty
