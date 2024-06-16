# The `AMR` Package for R <a href="https://msberends.github.io/AMR/"><img src="./logo.svg" align="right" height="139" /></a>

* Provides an **all-in-one solution** for antimicrobial resistance (AMR) data analysis in a One Health approach
* Used in over 175 countries, available in 20 languages
* Generates **antibiograms** - traditional, combined, syndromic, and even WISCA
* Provides the **full microbiological taxonomy** and extensive info on **all antimicrobial drugs**
* Applies all recent **CLSI** and **EUCAST** clinical and veterinary breakpoints for MICs and disk zones
* Corrects for duplicate isolates, **calculates** and **predicts** AMR per antibiotic class
* Integrates with **WHONET**, ATC, **EARS-Net**, PubChem, **LOINC**, **SNOMED CT**, and **NCBI**
* 100% free of costs and dependencies, highly suitable for places with **limited resources**

<div style="display: flex; font-size: 0.8em;">
<p style="text-align:left; width: 50%;"><small><a href="https://msberends.github.io/AMR/">https://msberends.github.io/AMR</a></small></p>
<p style="text-align:right; width: 50%;"><small><a href="https://doi.org/10.18637/jss.v104.i03" target="_blank">https://doi.org/10.18637/jss.v104.i03</a></small></p>
</div>

<a href="./reference/clinical_breakpoints.html#response-from-clsi-and-eucast"><img src="./endorsement_clsi_eucast.jpg" class="endorse_img" align="right" height="120" /></a>

----

### Introduction

The `AMR` package is a [free and open-source](#copyright) R package with [zero dependencies](https://en.wikipedia.org/wiki/Dependency_hell) to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial data and properties, by using evidence-based methods. **Our aim is to provide a standard** for clean and reproducible AMR data analysis, that can therefore empower epidemiological analyses to continuously enable surveillance and treatment evaluation in any setting. [Many different researchers](./authors.html) from around the globe are continually helping us to make this a successful and durable project!

This work was published in the Journal of Statistical Software (Volume 104(3); [DOI 10.18637/jss.v104.i03](https://doi.org/10.18637/jss.v104.i03)) and formed the basis of two PhD theses ([DOI 10.33612/diss.177417131](https://doi.org/10.33612/diss.177417131) and [DOI 10.33612/diss.192486375](https://doi.org/10.33612/diss.192486375)).

After installing this package, R knows [**~52,000 distinct microbial species**](./reference/microorganisms.html) (updated December 2022) and all [**~600 antibiotic, antimycotic and antiviral drugs**](./reference/antibiotics.html) by name and code (including ATC, EARS-Net, ASIARS-Net, PubChem, LOINC and SNOMED CT), and knows all about valid SIR and MIC values. The integral clinical breakpoint guidelines from CLSI and EUCAST are included, even with epidemiological cut-off (ECOFF) values. It supports and can read any data format, including WHONET data. This package works on Windows, macOS and Linux with all versions of R since R-3.0 (April 2013). **It was designed to work in any setting, including those with very limited resources**. It was created for both routine data analysis and academic research at the Faculty of Medical Sciences of the [University of Groningen](https://www.rug.nl), in collaboration with non-profit organisations [Certe Medical Diagnostics and Advice Foundation](https://www.certe.nl) and [University Medical Center Groningen](https://www.umcg.nl).

##### Used in over 175 countries, available in 20 languages

<a href="./countries_large.png" target="_blank"><img src="./countries.png" align="right" style="max-width: 300px;" /></a>

Since its first public release in early 2018, this R package has been used in almost all countries in the world. Click the map to enlarge and to see the country names.

With the help of contributors from all corners of the world, the `AMR` package is available in <img src="lang_en.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> English, <img src="lang_cs.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Czech, <img src="lang_zh.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Chinese, <img src="lang_da.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Danish, <img src="lang_nl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Dutch, <img src="lang_fi.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Finnish, <img src="lang_fr.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> French, <img src="lang_de.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> German, <img src="lang_el.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Greek, <img src="lang_it.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Italian, <img src="lang_ja.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Japanese, <img src="lang_no.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Norwegian, <img src="lang_pl.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Polish, <img src="lang_pt.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Portuguese, <img src="lang_ro.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Romanian, <img src="lang_ru.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Russian, <img src="lang_es.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Spanish, <img src="lang_sv.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Swedish, <img src="lang_tr.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Turkish, and <img src="lang_uk.svg" style="height: 13px !important; border: 1px solid #cccccc; vertical-align: initial !important;"> Ukrainian. Antimicrobial drug (group) names and colloquial microorganism names are provided in these languages.

### Practical examples

#### Filtering and selecting data

One of the most powerful functions of this package, aside from calculating and plotting AMR, is selecting and filtering based on antibiotic columns. This can be done using the so-called [antibiotic class selectors](https://msberends.github.io/AMR/reference/antibiotic_class_selectors.html) that work in base R, `dplyr` and `data.table`:

```r
# AMR works great with dplyr, but it's not required or neccesary
library(AMR)
library(dplyr)

example_isolates %>%
  mutate(bacteria = mo_fullname()) %>%
  # filtering functions for microorganisms:
  filter(mo_is_gram_negative(),
         mo_is_intrinsic_resistant(ab = "cefotax")) %>%
  # antibiotic selectors:
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
library(AMR)
example_isolates$bacteria <- mo_fullname(example_isolates$mo)
example_isolates[which(mo_is_gram_negative() &
                         mo_is_intrinsic_resistant(ab = "cefotax")),
                 c("bacteria", aminoglycosides(), carbapenems())]
```

This base R code will work in any version of R since April 2013 (R-3.0). Moreover, this code works identically with the `data.table` package, only by starting with:

```r
example_isolates <- data.table::as.data.table(example_isolates)
```

#### Generating antibiograms

The `AMR` package supports generating traditional, combined, syndromic, and even weighted-incidence syndromic combination antibiograms (WISCA).

If used inside R Markdown or Quarto, the table will be printed in the right output format automatically (such as markdown, LaTeX, HTML, etc.).

```r
antibiogram(example_isolates,
            antibiotics = c(aminoglycosides(), carbapenems()))
```

|Pathogen (N min-max)     | AMK| GEN| IPM| KAN| MEM| TOB|
|:------------------------|---:|---:|---:|---:|---:|---:|
|CoNS (43-309)            |   0|  86|  52|   0|  52|  22|
|*E. coli* (0-462)        | 100|  98| 100|    | 100|  97|
|*E. faecalis* (0-39)     |   0|   0| 100|   0|    |   0|
|*K. pneumoniae* (0-58)   |    |  90| 100|    | 100|  90|
|*P. aeruginosa* (17-30)  |    | 100|    |   0|    | 100|
|*P. mirabilis* (0-34)    |    |  94|  94|    |    |  94|
|*S. aureus* (2-233)      |    |  99|    |    |    |  98|
|*S. epidermidis* (8-163) |   0|  79|    |   0|    |  51|
|*S. hominis* (3-80)      |    |  92|    |    |    |  85|
|*S. pneumoniae* (11-117) |   0|   0|    |   0|    |   0|

In combination antibiograms, it is clear that combined antibiotics yield higher empiric coverage:

```r
antibiogram(example_isolates,
            antibiotics = c("TZP", "TZP+TOB", "TZP+GEN"),
            mo_transform = "gramstain")
```

|Pathogen (N min-max)     | TZP| TZP + GEN| TZP + TOB|
|:------------------------|---:|---------:|---------:|
|Gram-negative (641-693)  |  88|        99|        98|
|Gram-positive (345-1044) |  86|        98|        95|

Like many other functions in this package, `antibiogram()` comes with support for 20 languages that are often detected automatically based on system language:

```r
antibiogram(example_isolates,
            antibiotics = c("cipro", "tobra", "genta"), # any arbitrary name or code will work
            mo_transform = "gramstain",
            ab_transform = "name",
            language = "uk") # Ukrainian
```

|Збудник (N min-max)      | Гентаміцин| Тобраміцин| Ципрофлоксацин|
|:------------------------|----------:|----------:|--------------:|
|Грамнегативні (684-686)  |         96|         96|             91|
|Грампозитивні (665-1170) |         63|         34|             77|


#### Calculating resistance per group

For a manual approach, you can use the `resistance` or `susceptibility()` function:

```r
example_isolates %>%
  # group by ward:
  group_by(ward) %>%
  # calculate AMR using resistance() for gentamicin and tobramycin
  # and get their 95% confidence intervals using sir_confidence_interval():
  summarise(across(c(GEN, TOB),
                   list(total_R = resistance,
                        conf_int = function(x) sir_confidence_interval(x, collapse = "-"))))
```

|ward       | GEN_total_R|GEN_conf_int | TOB_total_R|TOB_conf_int |
|:---------:|:----------:|:-----------:|:----------:|:-----------:|
|Clinical   |   0.229    |0.205-0.254  |   0.315    |0.284-0.347  |
|ICU        |   0.290    |0.253-0.330  |   0.400    |0.353-0.449  |
|Outpatient |   0.200    |0.131-0.285  |   0.368    |0.254-0.493  |

Or use [antibiotic class selectors](https://msberends.github.io/AMR/reference/antibiotic_class_selectors.html) to select a series of antibiotic columns:

```r
library(AMR)
library(dplyr)

out <- example_isolates %>%
  # group by ward:
  group_by(ward) %>%
  # calculate AMR using resistance(), over all aminoglycosides and polymyxins:
  summarise(across(c(aminoglycosides(), polymyxins()),
            resistance))
out
```

| ward       |   GEN |   TOB |   AMK |   KAN |   COL |
|:-----------|------:|------:|------:|------:|------:|
| Clinical   | 0.229 | 0.315 | 0.626 |     1 | 0.780 |
| ICU        | 0.290 | 0.400 | 0.662 |     1 | 0.857 |
| Outpatient | 0.200 | 0.368 | 0.605 |       | 0.889 |

```r
# transform the antibiotic columns to names:
out %>% set_ab_names()
```

| ward       | gentamicin | tobramycin | amikacin | kanamycin | colistin  |
|:-----------|-----------:|-----------:|----------|----------:|----------:|
| Clinical   | 0.229      | 0.315      | 0.626    |     1     | 0.780     |
| ICU        | 0.290      | 0.400      | 0.662    |     1     | 0.857     |
| Outpatient | 0.200      | 0.368      | 0.605    |           | 0.889     |

```r
# transform the antibiotic column to ATC codes:
out %>% set_ab_names(property = "atc")
```

| ward       |  J01GB03   |   J01GB01  |  J01GB06 |  J01GB04  |  J01XB01  |
|:-----------|-----------:|-----------:|----------|----------:|----------:|
| Clinical   | 0.229      | 0.315      | 0.626    |     1     | 0.780     |
| ICU        | 0.290      | 0.400      | 0.662    |     1     | 0.857     |
| Outpatient | 0.200      | 0.368      | 0.605    |           | 0.889     |

### What else can you do with this package?

This package was intended as a comprehensive toolbox for integrated AMR data analysis. This package can be used for:

  * Reference for the taxonomy of microorganisms, since the package contains all microbial (sub)species from the List of Prokaryotic names with Standing in Nomenclature ([LPSN]((https://lpsn.dsmz.de))) and the Global Biodiversity Information Facility ([GBIF](https://www.gbif.org)) ([manual](./reference/mo_property.html))
  * Interpreting raw MIC and disk diffusion values, based on any CLSI or EUCAST guideline ([manual](./reference/as.sir.html))
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
  * Machine reading the EUCAST and CLSI guidelines from 2011-2021 to translate MIC values and disk diffusion diameters to SIR ([link](./articles/datasets.html))
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

[![check-old](https://github.com/msberends/AMR/actions/workflows/check-old.yaml/badge.svg?branch=main)](https://github.com/msberends/AMR/actions/workflows/check-old.yaml?query=branch%3Amain)
[![check-recent](https://github.com/msberends/AMR/actions/workflows/check-recent.yaml/badge.svg?branch=main)](https://github.com/msberends/AMR/actions/workflows/check-recent.yaml?query=branch%3Amain)
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
