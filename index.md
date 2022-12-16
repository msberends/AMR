# The `AMR` Package for R <a href="https://msberends.github.io/AMR/"><img src="./logo.svg" align="right" height="139" /></a>

* Works on Windows, macOS and Linux with **all versions of R** since R-3.0
* Provides the **full microbiological taxonomy** and data on **all antimicrobial drugs**
* Applies all recent **CLSI and EUCAST clinical breakpoints** for MICs and disk zones
* Corrects for duplicate isolates, **calculates and predicts AMR** per antibiotic class
* Integrates with **WHONET**, ATC, **EARS-Net**, PubChem, **LOINC** and **SNOMED CT**
* Completely **dependency-free**, highly suitable for places with **limited resources**

<div style="display: flex;">
<p style="text-align:left; width: 50%;"><small><a href="https://msberends.github.io/AMR/">https://msberends.github.io/AMR</a></small></p>
<p style="text-align:right; width: 50%;"><small><a href="https://doi.org/10.18637/jss.v104.i03">https://doi.org/10.18637/jss.v104.i03</a></small></p>
</div>

----

### Introduction

The `AMR` package is a [free and open-source](#copyright) R package with [zero dependencies](https://en.wikipedia.org/wiki/Dependency_hell) to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial data and properties, by using evidence-based methods. **Our aim is to provide a standard** for clean and reproducible AMR data analysis, that can therefore empower epidemiological analyses to continuously enable surveillance and treatment evaluation in any setting. 

This work was published in the Journal of Statistical Software (Volume 104(3); [DOI 10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)) and formed the basis of two PhD theses ([DOI 10.33612/diss.177417131](https://doi.org/10.33612/diss.177417131) and [DOI 10.33612/diss.192486375](https://doi.org/10.33612/diss.192486375)).

After installing this package, R knows [**~52,000 distinct microbial species**](./reference/microorganisms.html) (updated December 2022) and all [**~600 antibiotic, antimycotic and antiviral drugs**](./reference/antibiotics.html) by name and code (including ATC, EARS-Net, ASIARS-Net, PubChem, LOINC and SNOMED CT), and knows all about valid R/SI and MIC values. The integral breakpoint guidelines from CLSI and EUCAST are included from the last 10 years. It supports and can read any data format, including WHONET data. This package works on Windows, macOS and Linux with all versions of R since R-3.0 (April 2013). **It was designed to work in any setting, including those with very limited resources**. It was created for both routine data analysis and academic research at the Faculty of Medical Sciences of the [University of Groningen](https://www.rug.nl), in collaboration with non-profit organisations [Certe Medical Diagnostics and Advice Foundation](https://www.certe.nl) and [University Medical Center Groningen](https://www.umcg.nl), and is being [actively and durably maintained](./news) by two public healthcare organisations in the Netherlands.

##### Used in 175 countries, translated to 16 languages

<a href="./countries_large.png"><img src="./countries.png" target="_blank" align="right" style="max-width: 300px;" /></a>

Since its first public release in early 2018, this R package has been used in almost all countries in the world. Click the map to enlarge and to see the country names.

The `AMR` package is available in <img src="lang_en.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> English, <img src="lang_zh.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Chinese, <img src="lang_da.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Danish, <img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Dutch, <img src="lang_fr.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> French, <img src="lang_de.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> German, <img src="lang_el.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Greek, <img src="lang_it.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Italian, <img src="lang_ja.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Japanese, <img src="lang_pl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Polish, <img src="lang_pt.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Portuguese, <img src="lang_ru.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Russian, <img src="lang_es.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Spanish, <img src="lang_sv.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Swedish, <img src="lang_tr.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Turkish, and <img src="lang_uk.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Ukrainian. Antimicrobial drug (group) names and colloquial microorganism names are provided in these languages.

### Practical examples

#### Filtering and selecting data

```r
# AMR works great with dplyr, but it's not required or neccesary
library(AMR)
library(dplyr)

example_isolates %>%
  mutate(bacteria = mo_fullname()) %>%
  filter(mo_is_gram_negative(),
         mo_is_intrinsic_resistant(ab = "cefotax")) %>%
  select(bacteria,
         aminoglycosides(),
         carbapenems())
```

With only having defined a row filter on Gram-negative bacteria with intrinsic resistance to cefotaxime (`mo_is_gram_negative()` and `mo_is_intrinsic_resistant()`) and a column selection on two antibiotic groups (`aminoglycosides()` and `carbapenems()`), the reference data about [all microorganisms](./reference/microorganisms.html) and [all antibiotics](./reference/antibiotics.html) in the `AMR` package make sure you get what you meant:

|bacteria                       | GEN | TOB | AMK | KAN | IPM | MEM |
|:------------------------------|:---:|:---:|:---:|:---:|:---:|:---:|
|*Pseudomonas aeruginosa*       |  I  |  S  |     |  R  |  S  |     |
|*Pseudomonas aeruginosa*       |  I  |  S  |     |  R  |  S  |     |
|*Pseudomonas aeruginosa*       |  I  |  S  |     |  R  |  S  |     |
|*Pseudomonas aeruginosa*       |  S  |  S  |  S  |  R  |     |  S  |
|*Pseudomonas aeruginosa*       |  S  |  S  |  S  |  R  |  S  |  S  |
|*Pseudomonas aeruginosa*       |  S  |  S  |  S  |  R  |  S  |  S  |
|*Stenotrophomonas maltophilia* |  R  |  R  |  R  |  R  |  R  |  R  |
|*Pseudomonas aeruginosa*       |  S  |  S  |  S  |  R  |     |  S  |
|*Pseudomonas aeruginosa*       |  S  |  S  |  S  |  R  |     |  S  |
|*Pseudomonas aeruginosa*       |  S  |  S  |  S  |  R  |  S  |  S  |

A base R equivalent would be:

```r
example_isolates$bacteria <- mo_fullname(example_isolates$mo)
example_isolates[which(mo_is_gram_negative() &
                         mo_is_intrinsic_resistant(ab = "cefotax")),
                 c("bacteria", aminoglycosides(), carbapenems())]
```

This base R snippet will work in any version of R since April 2013 (R-3.0).

#### Calculating resistance per group

```r
library(AMR)
library(dplyr)

out <- example_isolates %>%
  # group by ward:
  group_by(ward) %>%
  # calculate AMR using resistance(), over all aminoglycosides
  # and polymyxins:
  summarise(across(c(aminoglycosides(), polymyxins()),
            resistance))
out
```

| ward       |   GEN |   TOB |   AMK |   KAN |   COL |
|:-----------|------:|------:|------:|------:|------:|
| Clinical   | 0.229 | 0.315 | 0.626 |     1 | 0.780 |
| ICU        | 0.290 | 0.400 | 0.662 |     1 | 0.857 |
| Outpatient | 0.200 | 0.368 | 0.605 |    NA | 0.889 |

```r
# transform the antibiotic columns to names:
out %>% set_ab_names()
```

| ward       | gentamicin | tobramycin | amikacin | kanamycin | colistin  |
|:-----------|-----------:|-----------:|----------|----------:|----------:|
| Clinical   | 0.229      | 0.315      | 0.626    |     1     | 0.780     |
| ICU        | 0.290      | 0.400      | 0.662    |     1     | 0.857     |
| Outpatient | 0.200      | 0.368      | 0.605    |    NA     | 0.889     |

```r
# transform the antibiotic column to ATC codes:
out %>% set_ab_names(property = "atc")
```

| ward       |  J01GB03   |   J01GB01  |  J01GB06 |  J01GB04  |  J01XB01  |
|:-----------|-----------:|-----------:|----------|----------:|----------:|
| Clinical   | 0.229      | 0.315      | 0.626    |     1     | 0.780     |
| ICU        | 0.290      | 0.400      | 0.662    |     1     | 0.857     |
| Outpatient | 0.200      | 0.368      | 0.605    |    NA     | 0.889     |

### What else can you do with this package?

This package was intended as a comprehensive toolbox for integrated AMR data analysis. This package can be used for:

  * Reference for the taxonomy of microorganisms, since the package contains all microbial (sub)species from the List of Prokaryotic names with Standing in Nomenclature ([LPSN]((https://lpsn.dsmz.de))) and the Global Biodiversity Information Facility ([GBIF](https://www.gbif.org)) ([manual](./reference/mo_property.html))
  * Interpreting raw MIC and disk diffusion values, based on any CLSI or EUCAST guideline from the last 10 years ([manual](./reference/as.rsi.html))
  * Retrieving antimicrobial drug names, doses and forms of administration from clinical health care records ([manual](./reference/ab_from_text.html))
  * Determining first isolates to be used for AMR data analysis ([manual](./reference/first_isolate.html))
  * Calculating antimicrobial resistance ([tutorial](./articles/AMR.html))
  * Determining multi-drug resistance (MDR) / multi-drug resistant organisms (MDRO) ([tutorial](./articles/MDR.html))
  * Calculating (empirical) susceptibility of both mono therapy and combination therapies ([tutorial](./articles/AMR.html))
  * Predicting future antimicrobial resistance using regression models ([tutorial](./articles/resistance_predict.html))
  * Getting properties for any microorganism (like Gram stain, species, genus or family) ([manual](./reference/mo_property.html))
  * Getting properties for any antibiotic (like name, code of EARS-Net/ATC/LOINC/PubChem, defined daily dose or trade name) ([manual](./reference/ab_property.html))
  * Plotting antimicrobial resistance ([tutorial](./articles/AMR.html))
  * Applying EUCAST expert rules ([manual](./reference/eucast_rules.html))
  * Getting SNOMED codes of a microorganism, or getting properties of a microorganism based on a SNOMED code ([manual](./reference/mo_property.html))
  * Getting LOINC codes of an antibiotic, or getting properties of an antibiotic based on a LOINC code ([manual](./reference/ab_property.html))
  * Machine reading the EUCAST and CLSI guidelines from 2011-2021 to translate MIC values and disk diffusion diameters to R/SI ([link](./articles/datasets.html))
  * Principal component analysis for AMR ([tutorial](./articles/PCA.html))

### Get this package

#### Latest official version

[![CRAN](https://www.r-pkg.org/badges/version-ago/AMR)](https://cran.r-project.org/package=AMR)
[![CRANlogs](https://cranlogs.r-pkg.org/badges/grand-total/AMR)](https://cran.r-project.org/package=AMR)

This package is available [here on the official R network (CRAN)](https://cran.r-project.org/package=AMR). Install this package in R from CRAN by using the command:

```r
install.packages("AMR")
```

It will be downloaded and installed automatically. For RStudio, click on the menu *Tools* > *Install Packages...* and then type in "AMR" and press <kbd>Install</kbd>.

**Note:** Not all functions on this website may be available in this latest release. To use all functions and data sets mentioned on this website, install the latest development version.

#### Latest development version

[![R-code-check](https://github.com/msberends/AMR/workflows/R-code-check/badge.svg?branch=main)](https://codecov.io/gh/msberends/AMR?branch=main)
[![CodeFactor](https://www.codefactor.io/repository/github/msberends/amr/badge)](https://www.codefactor.io/repository/github/msberends/amr)
[![Codecov](https://codecov.io/gh/msberends/AMR/branch/main/graph/badge.svg)](https://codecov.io/gh/msberends/AMR?branch=main)

Please read our [Developer Guideline here](https://github.com/msberends/AMR/wiki/Developer-Guideline).

The latest and unpublished development version can be installed from GitHub in two ways:

1. Manually, using:

   ```r
   install.packages("remotes") # if you haven't already
   remotes::install_github("msberends/AMR")
   ```
   
2. Automatically, using the [rOpenSci R-universe platform](https://ropensci.org/r-universe/), by adding [our R-universe address](https://msberends.r-universe.dev) to your list of repositories ('repos'):

   ```r
   options(repos = c(getOption("repos"),
                     msberends = "https://msberends.r-universe.dev"))
   ```
   
   After this, you can install and update this `AMR` package like any official release (e.g., using `install.packages("AMR")` or in RStudio via *Tools* > *Check for Package Updates...*).

### Get started

To find out how to conduct AMR data analysis, please [continue reading here to get started](./articles/AMR.html) or click a link in the ['How to' menu](https://msberends.github.io/AMR/articles/).

### Partners

The development of this package is part of, related to, or made possible by the following non-profit organisations and initiatives:

<div align="center">
  <a href="https://www.rug.nl" title="University of Groningen"><img src="./logo_rug.svg" style="max-width: 200px;"></a>
  <a href="https://www.umcg.nl" title="University Medical Center Groningen"><img src="./logo_umcg.svg" style="max-width: 200px;"></a>
  <a href="https://www.certe.nl" title="Certe Medical Diagnostics and Advice Foundation"><img src="./logo_certe.svg" style="max-width: 200px;"></a>
  <a href="https://www.deutschland-nederland.eu" title="EurHealth-1-Health"><img src="./logo_eh1h.png" style="max-width: 200px;"></a>
  <a href="https://www.deutschland-nederland.eu" title="INTERREG"><img src="./logo_interreg.png" style="max-width: 200px;"></a>
</div>

### Copyright

This R package is free, open-source software and licensed under the [GNU General Public License v2.0 (GPL-2)](./LICENSE-text.html). In a nutshell, this means that this package:

- May be used for commercial purposes

- May be used for private purposes

- May **not** be used for patent purposes

- May be modified, although:

  - Modifications **must** be released under the same license when distributing the package
  - Changes made to the code **must** be documented

- May be distributed, although:

  - Source code **must** be made available when the package is distributed
  - A copy of the license and copyright notice **must** be included with the package.

- Comes with a LIMITATION of liability

- Comes with NO warranty
