
<!-- README.md is generated from README.Rmd; please edit that file. -->

# The `AMR` Package for R

Please visit our comprehensive package website <https://amr-for-r.org>
to read more about this package, including many examples and tutorials.

Overview:

- Provides an **all-in-one solution** for antimicrobial resistance (AMR)
  data analysis in a One Health approach
- Peer-reviewed, used in over 175 countries, available in 28 languages
- Generates **antibiograms** - traditional, combined, syndromic, and
  even WISCA
- Provides the **full microbiological taxonomy** of ~79 000 distinct
  species and extensive info of ~620 antimicrobial drugs
- Applies **CLSI 2011-2025** and **EUCAST 2011-2025** clinical and
  veterinary breakpoints, and ECOFFs, for MIC and disk zone
  interpretation
- Corrects for duplicate isolates, **calculates** and **predicts** AMR
  per antimicrobial class
- Integrates with **WHONET**, ATC, **EARS-Net**, PubChem, **LOINC**,
  **SNOMED CT**, and **NCBI**
- 100% free of costs and dependencies, highly suitable for places with
  **limited resources**

------------------------------------------------------------------------

The `AMR` package is a peer-reviewed, free and open-source R package
with zero dependencies to simplify the analysis and prediction of
Antimicrobial Resistance (AMR) and to work with microbial and
antimicrobial data and properties, by using evidence-based methods.
**Our aim is to provide a standard** for clean and reproducible AMR data
analysis, that can therefore empower epidemiological analyses to
continuously enable surveillance and treatment evaluation in any
setting.

The `AMR` package supports and can read any data format, including
WHONET data. This package works on Windows, macOS and Linux with all
versions of R since R-3.0 (April 2013). **It was designed to work in any
setting, including those with very limited resources**. It was created
for both routine data analysis and academic research at the Faculty of
Medical Sciences of the [University of Groningen](https://www.rug.nl)
and the [University Medical Center Groningen](https://www.umcg.nl).

------------------------------------------------------------------------

### How to get this package

To install the latest ‘release’ version from CRAN:

``` r
install.packages("AMR")
```

To install the latest ‘beta’ version:

``` r
install.packages("AMR", repos = "beta.amr-for-r.org")

# if this does not work, try to install directly from GitHub using the 'remotes' package:
remotes::install_github("msberends/AMR")
```

------------------------------------------------------------------------

<small> This AMR package for R is free, open-source software and
licensed under the [GNU General Public License v2.0
(GPL-2)](https://amr-for-r.org/LICENSE-text.html). These requirements
are consequently legally binding: modifications must be released under
the same license when distributing the package, changes made to the code
must be documented, source code must be made available when the package
is distributed, and a copy of the license and copyright notice must be
included with the package. </small>
