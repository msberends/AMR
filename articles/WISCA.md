# Estimating Empirical Coverage with WISCA

> This explainer was largely written by our [AMR for R
> Assistant](https://chat.amr-for-r.org), a ChatGPT manually-trained
> model able to answer any question about the `AMR` package.

## Introduction

Clinical guidelines for empirical antimicrobial therapy require
*probabilistic reasoning*: what is the chance that a regimen will cover
the likely infecting organisms, before culture results are available?

This is the purpose of **WISCA**, or **Weighted-Incidence Syndromic
Combination Antibiogram**.

WISCA is a Bayesian approach that integrates:

- **Pathogen prevalence** (how often each species causes the syndrome),
- **Regimen susceptibility** (how often a regimen works *if* the
  pathogen is known),

to estimate the **overall empirical coverage** of antimicrobial
regimens, with quantified uncertainty.

This vignette explains how WISCA works, why it is useful, and how to
apply it using the `AMR` package.

## Why traditional antibiograms fall short

A standard antibiogram gives you:

    Species → Antibiotic → Susceptibility %

But clinicians don’t know the species *a priori*. They need to choose a
regimen that covers the **likely pathogens**, without knowing which one
is present.

Traditional antibiograms calculate the susceptibility % as just the
number of resistant isolates divided by the total number of tested
isolates. Therefore, traditional antibiograms:

- Fragment information by organism,
- Do not weight by real-world prevalence,
- Do not account for combination therapy or sample size,
- Do not provide uncertainty.

## The idea of WISCA

WISCA asks:

> “What is the **probability** that this regimen **will cover** the
> pathogen, given the syndrome?”

This means combining two things:

- **Incidence** of each pathogen in the syndrome,
- **Susceptibility** of each pathogen to the regimen.

We can write this as:

$$\text{Coverage} = \sum\limits_{i}\left( \text{Incidence}_{i} \times \text{Susceptibility}_{i} \right)$$

For example, suppose:

- *E. coli* causes 60% of cases, and 90% of *E. coli* are susceptible to
  a drug.
- *Klebsiella* causes 40% of cases, and 70% of *Klebsiella* are
  susceptible.

Then:

$$\text{Coverage} = (0.6 \times 0.9) + (0.4 \times 0.7) = 0.82$$

But in real data, incidence and susceptibility are **estimated from
samples**, so they carry uncertainty. WISCA models this
**probabilistically**, using conjugate Bayesian distributions.

## The Bayesian engine behind WISCA

### Pathogen incidence

Let:

- $K$ be the number of pathogens,
- $\alpha = (1,1,\ldots,1)$ be a **Dirichlet** prior (uniform),
- $n = \left( n_{1},\ldots,n_{K} \right)$ be the observed counts per
  species.

Then the posterior incidence is:

$$p \sim \text{Dirichlet}\left( \alpha_{1} + n_{1},\ldots,\alpha_{K} + n_{K} \right)$$

To simulate from this, we use:

$$x_{i} \sim \text{Gamma}\left( \alpha_{i} + n_{i},\ 1 \right),\quad p_{i} = \frac{x_{i}}{\sum\limits_{j = 1}^{K}x_{j}}$$

### Susceptibility

Each pathogen–regimen pair has a prior and data:

- Prior: $\text{Beta}\left( \alpha_{0},\beta_{0} \right)$, with default
  $\alpha_{0} = \beta_{0} = 1$
- Data: $S$ susceptible out of $N$ tested

The $S$ category could also include values SDD (susceptible,
dose-dependent) and I (intermediate \[CLSI\], or susceptible, increased
exposure \[EUCAST\]).

Then the posterior is:

$$\theta \sim \text{Beta}\left( \alpha_{0} + S,\ \beta_{0} + N - S \right)$$

### Final coverage estimate

Putting it together:

1.  Simulate pathogen incidence: $\mathbf{p} \sim \text{Dirichlet}$
2.  Simulate susceptibility:
    $\theta_{i} \sim \text{Beta}\left( 1 + S_{i},\ 1 + R_{i} \right)$
3.  Combine:

$$\text{Coverage} = \sum\limits_{i = 1}^{K}p_{i} \cdot \theta_{i}$$

Repeat this simulation (e.g. 1000×) and summarise:

- **Mean** = expected coverage
- **Quantiles** = credible interval

## Practical use in the `AMR` package

### Prepare data and simulate synthetic syndrome

``` r
library(AMR)
data <- example_isolates

# Structure of our data
data
#> # A tibble: 2,000 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-03 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  3 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  4 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  5 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  6 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  7 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  8 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  9 2002-01-16 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#> 10 2002-01-17 858515     79 F      ICU      B_STPHY_EPDR   R     NA    S     NA 
#> # ℹ 1,990 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …

# Add a fake syndrome column
data$syndrome <- ifelse(data$mo %like% "coli", "UTI", "No UTI")
```

### Basic WISCA antibiogram

``` r
wisca(data,
      antimicrobials = c("AMC", "CIP", "GEN"))
```

| Amoxicillin/clavulanic acid | Ciprofloxacin    | Gentamicin         |
|:----------------------------|:-----------------|:-------------------|
| 73.7% (71.7-75.8%)          | 77% (74.3-79.4%) | 72.8% (70.7-74.8%) |

### Use combination regimens

``` r
wisca(data,
      antimicrobials = c("AMC", "AMC + CIP", "AMC + GEN"))
```

| Amoxicillin/clavulanic acid | Amoxicillin/clavulanic acid + Ciprofloxacin | Amoxicillin/clavulanic acid + Gentamicin |
|:----------------------------|:--------------------------------------------|:-----------------------------------------|
| 73.8% (71.8-75.7%)          | 87.5% (85.9-89%)                            | 89.7% (88.2-91.1%)                       |

### Stratify by syndrome

``` r
wisca(data,
      antimicrobials = c("AMC", "AMC + CIP", "AMC + GEN"),
      syndromic_group = "syndrome")
```

| Syndromic Group | Amoxicillin/clavulanic acid | Amoxicillin/clavulanic acid + Ciprofloxacin | Amoxicillin/clavulanic acid + Gentamicin |
|:----------------|:----------------------------|:--------------------------------------------|:-----------------------------------------|
| No UTI          | 70.1% (67.8-72.3%)          | 85.2% (83.1-87.2%)                          | 87.1% (85.3-88.7%)                       |
| UTI             | 80.9% (77.7-83.8%)          | 88.2% (85.7-90.5%)                          | 90.9% (88.7-93%)                         |

The `AMR` package is available in 28 languages, which can all be used
for the [`wisca()`](https://amr-for-r.org/reference/antibiogram.md)
function too:

``` r
wisca(data,
      antimicrobials = c("AMC", "AMC + CIP", "AMC + GEN"),
      syndromic_group = gsub("UTI", "UCI", data$syndrome),
      language = "Spanish")
```

| Grupo sindrómico | Amoxicillin ácido clavulánico | Amoxicillin ácido clavulánico + Ciprofloxacin | Amoxicillin ácido clavulánico + Gentamicina |
|:-----------------|:------------------------------|:----------------------------------------------|:--------------------------------------------|
| No UCI           | 70% (67.8-72.4%)              | 85.3% (83.3-87.2%)                            | 87% (85.3-88.8%)                            |
| UCI              | 80.9% (77.7-83.9%)            | 88.2% (85.5-90.6%)                            | 90.9% (88.7-93%)                            |

## Sensible defaults, which can be customised

- `simulations = 1000`: number of Monte Carlo draws
- `conf_interval = 0.95`: coverage interval width
- `combine_SI = TRUE`: count “I” and “SDD” as susceptible

## Limitations

- It assumes your data are representative
- No adjustment for patient-level covariates, although these could be
  passed onto the `syndromic_group` argument
- WISCA does not model resistance over time, you might want to use
  `tidymodels` for that, for which we [wrote a basic
  introduction](https://amr-for-r.org/articles/AMR_with_tidymodels.html)

## Summary

WISCA enables:

- Empirical regimen comparison,
- Syndrome-specific coverage estimation,
- Fully probabilistic interpretation.

It is available in the `AMR` package via either:

``` r
wisca(...)

antibiogram(..., wisca = TRUE)
```

## Reference

Bielicki, JA, et al. (2016). *Selecting appropriate empirical antibiotic
regimens for paediatric bloodstream infections: application of a
Bayesian decision model to local and pooled antimicrobial resistance
surveillance data.* **J Antimicrob Chemother**. 71(3):794-802.
<https://doi.org/10.1093/jac/dkv397>
