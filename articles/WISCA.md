# Estimating Empirical Coverage with WISCA

## Why WISCA?

When a clinician starts empirical antimicrobial therapy, the causative
pathogen is unknown. The question they need answered is not *“what
proportion of* E. coli *is susceptible to ciprofloxacin?“* but rather
*“what is the probability that this regimen will adequately cover
whatever pathogen turns out to be causing my patient’s infection?”*

The traditional cumulative antibiogram, as standardised by CLSI M39,
cannot answer that question. It presents susceptibility percentages per
species per antibiotic, but:

- **It fragments information by organism.** The clinician must mentally
  combine susceptibility rates across multiple species, weighting by how
  often each species causes the syndrome, a calculation nobody does at
  the bedside.
- **It ignores pathogen incidence.** A species that causes 2% of
  infections is given the same visual weight as one that causes 60%.
- **It does not evaluate combination regimens.** Much empirical therapy
  consists of two or more agents, but the traditional antibiogram only
  shows monotherapy per organism.
- **It provides no measure of uncertainty.** A reported “90%
  susceptible” based on 50 isolates has a 95% confidence interval of
  roughly 78-97% (Clopper-Pearson), yet the antibiogram presents it as a
  point estimate without context.

**WISCA** (Weighted-Incidence Syndromic Combination Antibiogram)
resolves all four limitations. It estimates the probability that a
regimen will provide adequate empirical coverage for a given infection
syndrome, weighted by local pathogen incidence, with full uncertainty
quantification via Bayesian inference.

The concept was introduced by Hebert *et al.* (2012), who demonstrated
that traditional antibiogram susceptibility rates could be misleading:
ciprofloxacin appeared 84% effective against *E. coli* in the
traditional antibiogram, but WISCA revealed only 62% coverage for UTI
and 37% for abdominal infections, because enterococci (intrinsically
resistant) and other species contribute substantially to these
syndromes. Randhawa *et al.* (2014) showed that WISCA-guided regimen
selection could improve time-to-adequate-coverage on the ICU by over
40%. Bielicki *et al.* (2016) introduced the Bayesian framework now used
in this package, enabling credible intervals and multi-centre pooling.
Cook *et al.* (2022) applied it globally across 52 hospitals in 23
countries.

## The idea

WISCA asks:

> “What is the **probability** that this regimen **will cover** the
> pathogen, given the syndrome?”

This means combining two quantities:

- **Pathogen incidence** in the syndrome (how often each species causes
  it),
- **Susceptibility** of each pathogen to the regimen.

We can write this as:

``` math
\text{Coverage} = \sum_i (\text{Incidence}_i \times \text{Susceptibility}_i)
```

For example, suppose in your hospital:

- *E. coli* causes 60% of UTIs, and 90% of *E. coli* are susceptible to
  a drug.
- *Klebsiella* causes 40% of UTIs, and 70% of *Klebsiella* are
  susceptible.

Then:

``` math
\text{Coverage} = (0.6 \times 0.9) + (0.4 \times 0.7) = 0.82
```

That 82% is a far more clinically meaningful number than the
species-level “90% of *E. coli*” and “70% of *Klebsiella*” reported
separately in a traditional antibiogram, because it directly answers the
question the clinician actually faces.

But in real data, both incidence and susceptibility are **estimated from
finite samples**, so they carry uncertainty. A sample of 50 isolates is
not a census. WISCA models this uncertainty **probabilistically**, using
conjugate Bayesian distributions.

## The Bayesian engine

### Pathogen incidence

Let:

- $`K`$ be the number of pathogens,
- $`\boldsymbol{\alpha} = (1, 1, \ldots, 1)`$ be a $`\text{Dirichlet}`$
  prior (uniform, non-informative),
- $`\boldsymbol{n} = (n_1, \ldots, n_K)`$ be the observed isolate counts
  per species.

Then the posterior incidence is:

``` math
\boldsymbol{p} \sim \text{Dirichlet}(\alpha_1 + n_1, \ldots, \alpha_K + n_K)
```

To simulate from this, we use:

``` math
x_i \sim \text{Gamma}(\alpha_i + n_i,\ 1), \quad p_i = \frac{x_i}{\sum_{j=1}^{K} x_j}
```

The Dirichlet is the conjugate prior for multinomial data. With the
non-informative prior $`\text{Dirichlet}(1, 1, \ldots, 1)`$, the
posterior is dominated by the data once sample sizes are reasonable.
With small samples, the posterior is appropriately more diffuse,
reflecting genuine uncertainty, and the resulting credible intervals
will be wider.

### Susceptibility

Each pathogen-regimen pair has a prior and observed data:

- Default prior: $`\text{Beta}(0.5, 0.5)`$ (Jeffreys prior)
- Intrinsically resistant pairs: $`\text{Beta}(1, 9999)`$, forcing
  near-zero susceptibility regardless of observed data (based on EUCAST
  Expected Resistant Phenotypes)
- Data: $`S`$ susceptible out of $`N`$ tested

The $`S`$ category could also include values SDD (susceptible,
dose-dependent) and I (intermediate \[CLSI\], or susceptible, increased
exposure \[EUCAST\]).

Then the posterior is:

``` math
\theta \sim \text{Beta}(\alpha_0 + S,\ \beta_0 + N - S)
```

### Final coverage estimate

Putting it together:

1.  Simulate pathogen incidence:
    $`\boldsymbol{p} \sim \text{Dirichlet}`$
2.  Simulate susceptibility:
    $`\theta_i \sim \text{Beta}(\alpha_0 + S_i,\ \beta_0 + N_i - S_i)`$
3.  Combine:

``` math
\text{Coverage} = \sum_{i=1}^{K} p_i \cdot \theta_i
```

Repeat this simulation (e.g., 1000 times) and summarise:

- **Mean** = expected coverage
- **Quantiles** = credible interval (95% by default)

Because each simulation draws from the full posterior, the resulting
distribution of coverage estimates naturally captures the joint
uncertainty in both pathogen incidence and susceptibility. The credible
interval tells you how confident you can be in the coverage estimate,
something a traditional antibiogram never provides.

## When to use WISCA vs. traditional antibiograms

| Goal                                  | Recommended approach      |
|---------------------------------------|---------------------------|
| Guide empirical therapy decisions     | **WISCA**                 |
| Compare regimens for a syndrome       | **WISCA**                 |
| Evaluate combination regimens         | **WISCA**                 |
| Antimicrobial stewardship (A-team)    | **WISCA**                 |
| Track resistance trends per species   | Traditional / Combination |
| AMR surveillance reporting            | Traditional / Syndromic   |
| Understand species-level epidemiology | Traditional               |

In short: if the end goal involves a *patient* who does not yet have a
culture result, WISCA is the appropriate tool. If the end goal is
*surveillance* of resistance at the species level, the traditional
antibiogram remains fit for purpose.

## Practical use in the `AMR` package

### Prepare data

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

# Add a synthetic syndrome column for demonstration
data$syndrome <- ifelse(data$mo %like% "coli", "UTI", "Non-UTI")
```

### Basic WISCA

``` r

wisca(data,
  antimicrobials = c("AMC", "CIP", "GEN")
)
```

| Amoxicillin/clavulanic acid | Ciprofloxacin      | Gentamicin         |
|:----------------------------|:-------------------|:-------------------|
| 74.2% (72.1-76.1%)          | 78.4% (75.6-81.1%) | 72.5% (70.4-74.6%) |

### Use combination regimens

Combination regimens are specified with a `+` separator. WISCA evaluates
whether *at least one* agent in the combination covers the pathogen:

``` r

wisca(data,
  antimicrobials = c("AMC", "AMC + CIP", "AMC + GEN")
)
```

| Amoxicillin/clavulanic acid | Amoxicillin/clavulanic acid + Ciprofloxacin | Amoxicillin/clavulanic acid + Gentamicin |
|:---|:---|:---|
| 74.2% (72.2-76.1%) | 88.8% (87.2-90.4%) | 90.8% (89.4-92.2%) |

### Stratify by syndrome

Use `syndromic_group` to produce separate WISCA estimates per clinical
stratum. You can pass a column name or any expression:

``` r

wisca(data,
  antimicrobials = c("AMC", "AMC + CIP", "AMC + GEN"),
  syndromic_group = "syndrome"
)
```

| Syndromic Group | Amoxicillin/clavulanic acid | Amoxicillin/clavulanic acid + Ciprofloxacin | Amoxicillin/clavulanic acid + Gentamicin |
|:---|:---|:---|:---|
| Non-UTI | 70.3% (67.9-72.7%) | 86.8% (84.9-88.7%) | 88.4% (86.4-90.2%) |
| UTI | 80.3% (77-83.3%) | 88.4% (85.7-90.8%) | 91% (88.3-93.3%) |

The `AMR` package is available in 28 languages, which can all be used
for the [`wisca()`](https://amr-for-r.org/reference/antibiogram.md)
function too:

``` r

wisca(data,
  antimicrobials = c("AMC", "AMC + CIP", "AMC + GEN"),
  syndromic_group = gsub("UTI", "UCI", data$syndrome),
  language = "Spanish"
)
```

| Grupo sindrómico | Amoxicilina/ácido clavulánico | Amoxicilina/ácido clavulánico + Ciprofloxacina | Amoxicilina/ácido clavulánico + Gentamicina |
|:---|:---|:---|:---|
| Non-UCI | 70.4% (68-72.8%) | 86.7% (84.6-88.7%) | 88.5% (86.5-90.2%) |
| UCI | 80.3% (77.2-83.5%) | 88.4% (85.5-90.8%) | 91% (88.4-93.1%) |

### Interpreting the output

Each row shows the estimated empirical coverage for a regimen, with a
95% credible interval. When comparing regimens:

- **Overlapping credible intervals** mean there is no statistically
  significant difference in coverage. If a narrower-spectrum regimen
  overlaps with a broader one, the narrower-spectrum option can be
  preferred on stewardship grounds.
- **Non-overlapping credible intervals** indicate a clinically
  meaningful difference in coverage.

## Sensible defaults, which can be customised

- `simulations = 1000`: number of Monte Carlo draws
- `conf_interval = 0.95`: coverage interval width
- `combine_SI = TRUE`: count “I” and “SDD” as susceptible

## Practical considerations

- **First isolates only**: always deduplicate using
  [`first_isolate()`](https://amr-for-r.org/reference/first_isolate.md)
  before running WISCA. Repeat isolates introduce bias.
- **Pathogen selection**: consider filtering with
  [`top_n_microorganisms()`](https://amr-for-r.org/reference/top_n_microorganisms.md).
  Including rare contaminants (e.g. CoNS without clinical context) can
  distort estimates and may artificially lower coverage (Cook *et al.*,
  2022).
- **Sample size**: coverage estimates become reliable with approximately
  100+ isolates. For smaller datasets, consider pooling data from
  multiple sites, but only after verifying that pathogen distributions
  are sufficiently similar (Bielicki *et al.*, 2016).
- **Culture request bias**: WISCA is only as good as the data it is
  based on. If cultures are selectively requested (e.g. only after
  treatment failure), the dataset will be biased towards resistant
  isolates. A robust culture policy is essential for reliable estimates.

## Limitations

- It assumes your data are representative of the patient population you
  are treating
- No direct adjustment for patient-level covariates, although these can
  be passed onto the `syndromic_group` argument for stratification
- WISCA does not model resistance trends over time; for that, you might
  want to use `tidymodels`, for which we [wrote a basic
  introduction](https://amr-for-r.org/articles/AMR_with_tidymodels.html)

## Summary

WISCA enables:

- **Empirical regimen comparison**, answering the clinician’s actual
  question
- **Syndrome-specific coverage estimation**, stratifiable by any
  clinical variable
- **Fully probabilistic interpretation**, with credible intervals that
  honestly communicate uncertainty

It is available in the `AMR` package via either:

``` r

wisca(...)

antibiogram(..., wisca = TRUE)
```

## References

1.  Hebert C, Ridgway J, Vekhter B, Brown EC, Weber SG, Robicsek A.
    Demonstration of the weighted-incidence syndromic combination
    antibiogram: an empiric prescribing decision aid. *Infect Control
    Hosp Epidemiol.* 2012;33(4):381-388.
    <https://doi.org/10.1086/664768>
2.  Randhawa V, Sarwar S, Walker S, Elligsen M, Palmay L, Daneman N.
    Weighted-incidence syndromic combination antibiograms to guide
    empiric treatment of critical care infections: a retrospective
    cohort study. *Crit Care.* 2014;18(3):R112.
    <https://doi.org/10.1186/cc13901>
3.  Bielicki JA, Sharland M, Johnson AP, Henderson KL, Cromwell DA.
    Selecting appropriate empirical antibiotic regimens for paediatric
    bloodstream infections: application of a Bayesian decision model to
    local and pooled antimicrobial resistance surveillance data. *J
    Antimicrob Chemother.* 2016;71(3):794-802.
    <https://doi.org/10.1093/jac/dkv397>
4.  Cook A, Sharland M, Yau Y, Bielicki J. Improving empiric antibiotic
    prescribing in pediatric bloodstream infections: a potential
    application of weighted-incidence syndromic combination antibiograms
    (WISCA). *Expert Rev Anti Infect Ther.* 2022;20(3):445-456.
    <https://doi.org/10.1080/14787210.2021.1967145>
