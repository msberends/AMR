# `AMR` (for R) <img src="./logo.svg" align="right"/>

> Update: The latest [EUCAST guideline for intrinsic resistance](https://www.eucast.org/expert_rules_and_intrinsic_resistance/) (v3.3, October 2021) is now supported, the CLSI 2021 interpretation guideline is now supported, and our taxonomy tables have been updated as well (LPSN, 5 October 2021).

### What is `AMR` (for R)?

`AMR` is a free, open-source and independent [R package](https://www.r-project.org) (see [Copyright](#copyright)) to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial data and properties, by using evidence-based methods. **Our aim is to provide a standard** for clean and reproducible AMR data analysis, that can therefore empower epidemiological analyses to continuously enable surveillance and treatment evaluation in any setting.

After installing this package, R knows [**~71,000 distinct microbial species**](./reference/microorganisms.html) and all [**~570 antibiotic, antimycotic and antiviral drugs**](./reference/antibiotics.html) by name and code (including ATC, WHONET/EARS-Net, PubChem, LOINC and SNOMED CT), and knows all about valid R/SI and MIC values. It supports any data format, including WHONET/EARS-Net data. 

The `AMR` package is available in <img src="lang_da.svg" style="height: 11px !important; vertical-align: initial !important;"> Danish, <img src="lang_nl.svg" style="height: 11px !important; vertical-align: initial !important;"> Dutch, <img src="lang_en.svg" style="height: 11px !important; vertical-align: initial !important;"> English, <img src="lang_fr.svg" style="height: 11px !important; vertical-align: initial !important;"> French, <img src="lang_de.svg" style="height: 11px !important; vertical-align: initial !important;"> German, <img src="lang_it.svg" style="height: 11px !important; vertical-align: initial !important;"> Italian, <img src="lang_pt.svg" style="height: 11px !important; vertical-align: initial !important;"> Portuguese, <img src="lang_ru.svg" style="height: 11px !important; vertical-align: initial !important;"> Russian, <img src="lang_es.svg" style="height: 11px !important; vertical-align: initial !important;"> Spanish and <img src="lang_sv.svg" style="height: 11px !important; vertical-align: initial !important;"> Swedish. Antimicrobial drug (group) names and colloquial microorganism names are provided in these languages.

This package is [fully independent of any other R package](https://en.wikipedia.org/wiki/Dependency_hell) and works on Windows, macOS and Linux with all versions of R since R-3.0 (April 2013). **It was designed to work in any setting, including those with very limited resources**. It was created for both routine data analysis and academic research at the Faculty of Medical Sciences of the [University of Groningen](https://www.rug.nl), in collaboration with non-profit organisations [Certe Medical Diagnostics and Advice Foundation](https://www.certe.nl) and [University Medical Center Groningen](https://www.umcg.nl). This R package formed the basis of two PhD theses ([DOI 10.33612/diss.177417131](https://doi.org/10.33612/diss.177417131) and [DOI 10.33612/diss.192486375](https://doi.org/10.33612/diss.192486375)) but is [actively and durably maintained](./news) by two public healthcare organisations in the Netherlands.

<div class="main-content" style="display: inline-block;">
  <p>
    <a href="./countries_large.png" target="_blank"><img src="./countries.png" class="countries_map"></a>
    <strong>Used in 175 countries</strong><br>
    Since its first public release in early 2018, this R package has been used in almost all countries in the world. Click the map to enlarge and to see the country names.</p>
</div>

##### With `AMR` (for R), there's always a knowledgeable microbiologist by your side!

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

A base R equivalent would be, giving the exact same results:

```r
example_isolates$bacteria <- mo_fullname(example_isolates$mo)
example_isolates[which(mo_is_gram_negative() &
                         mo_is_intrinsic_resistant(ab = "cefotax")),
                 c("bacteria", aminoglycosides(), carbapenems())]
```

#### Partners

The development of this package is part of, related to, or made possible by:

<div align="center">
  <a href="https://www.rug.nl" title="University of Groningen"><img src="./logo_rug.png" class="partner_logo"></a>
  <a href="https://www.umcg.nl" title="University Medical Center Groningen"><img src="./logo_umcg.png" class="partner_logo"></a>
  <a href="https://www.certe.nl" title="Certe Medical Diagnostics and Advice Foundation"><img src="./logo_certe.png" class="partner_logo"></a>
  <a href="https://www.deutschland-nederland.eu" title="EurHealth-1-Health"><img src="./logo_eh1h.png" class="partner_logo"></a>
  <a href="https://www.deutschland-nederland.eu" title="INTERREG"><img src="./logo_interreg.png" class="partner_logo"></a>
</div>

### What can you do with this package?

This package can be used for:

  * Reference for the taxonomy of microorganisms, since the package contains all microbial (sub)species from the [Catalogue of Life](http://www.catalogueoflife.org) and [List of Prokaryotic names with Standing in Nomenclature](https://lpsn.dsmz.de) ([manual](./reference/mo_property.html))
  * Interpreting raw MIC and disk diffusion values, based on the latest CLSI or EUCAST guidelines ([manual](./reference/as.rsi.html))
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

#### Latest released version

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

You can also download the latest build from our repository: <https://github.com/msberends/AMR/raw/main/data-raw/AMR_latest.tar.gz>

### Get started

To find out how to conduct AMR data analysis, please [continue reading here to get started](./articles/AMR.html) or click a link in the ['How to' menu](https://msberends.github.io/AMR/articles/).

### Short introduction

#### Microbial (taxonomic) reference data

This package contains the complete taxonomic tree of almost all ~70,000 microorganisms from the authoritative and comprehensive Catalogue of Life (CoL, [www.catalogueoflife.org](http://www.catalogueoflife.org)), supplemented by data from the List of Prokaryotic names with Standing in Nomenclature (LPSN, [lpsn.dsmz.de](https://lpsn.dsmz.de)). This supplementation is needed until the [CoL+ project](https://github.com/Sp2000/colplus) is finished, which we await. With `catalogue_of_life_version()` can be checked which version of the CoL is included in this package. 

Read more about which data from the Catalogue of Life [in our manual](./reference/catalogue_of_life.html).

#### Antimicrobial reference data

This package contains **all ~550 antibiotic, antimycotic and antiviral drugs** and their Anatomical Therapeutic Chemical (ATC) codes, ATC groups and Defined Daily Dose (DDD, oral and IV) from the World Health Organization Collaborating Centre for Drug Statistics Methodology (WHOCC, https://www.whocc.no) and the [Pharmaceuticals Community Register of the European Commission](https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm).

**NOTE: The WHOCC copyright does not allow use for commercial purposes, unlike any other info from this package. See https://www.whocc.no/copyright_disclaimer/.**

Read more about the data from WHOCC [in our manual](./reference/WHOCC.html).

#### WHONET / EARS-Net

We support WHONET and EARS-Net data. Exported files from WHONET can be imported into R and can be analysed easily using this package. For education purposes, we created an [example data set `WHONET`](./reference/WHONET.html) with the exact same structure as a WHONET export file. Furthermore, this package also contains a [data set antibiotics](./reference/antibiotics.html) with all EARS-Net antibiotic abbreviations, and knows almost all WHONET abbreviations for microorganisms. When using WHONET data as input for analysis, all input parameters will be set automatically.

Read our tutorial about [how to work with WHONET data here](./articles/WHONET.html).

#### Overview of functions

The `AMR` package basically does four important things:

1. It **cleanses existing data** by providing new *classes* for microoganisms, antibiotics and antimicrobial results (both S/I/R and MIC). By installing this package, you teach R everything about microbiology that is needed for analysis. These functions all use intelligent rules to guess results that you would expect:

   * Use `as.mo()` to get a microbial ID. The IDs are human readable for the trained eye - the ID of *Klebsiella pneumoniae* is "B_KLBSL_PNMN" (B stands for Bacteria) and the ID of *S. aureus* is "B_STPHY_AURS". The function takes almost any text as input that looks like the name or code of a microorganism like "E. coli", "esco" or "esccol" and tries to find expected results using intelligent rules combined with the included Catalogue of Life data set. It only takes milliseconds to find results, please see our [benchmarks](./articles/benchmarks.html). Moreover, it can group *Staphylococci* into coagulase negative and positive (CoNS and CoPS, see [source](./reference/as.mo.html#source)) and can categorise *Streptococci* into Lancefield groups (like beta-haemolytic *Streptococcus* Group B, [source](./reference/as.mo.html#source)).
   * Use `as.ab()` to get an antibiotic ID. Like microbial IDs, these IDs are also human readable based on those used by EARS-Net. For example, the ID of amoxicillin is `AMX` and the ID of gentamicin is `GEN`. The `as.ab()` function also uses intelligent rules to find results like accepting misspelling, trade names and abbrevations used in many laboratory systems. For instance, the values "Furabid", "Furadantin", "nitro" all return the ID of Nitrofurantoine. To accomplish this, the package contains a database with most LIS codes, official names, trade names, ATC codes, defined daily doses (DDD) and drug categories of antibiotics.
   * Use `as.rsi()` to get antibiotic interpretations based on raw MIC values (in mg/L) or disk diffusion values (in mm), or transform existing values to valid antimicrobial results. It produces just S, I or R based on your input and warns about invalid values. Even values like "<=0.002; S" (combined MIC/RSI) will result in "S".
   * Use `as.mic()` to cleanse your MIC values. It produces a so-called factor (called *ordinal* in SPSS) with valid MIC values as levels. A value like "<=0.002; S" (combined MIC/RSI) will result in "<=0.002".
   
2. It **enhances existing data** and **adds new data** from data sets included in this package.

   * Use `eucast_rules()` to apply [EUCAST expert rules to isolates](https://www.eucast.org/expert_rules_and_intrinsic_resistance/) (not the translation from MIC to R/SI values, use `as.rsi()` for that).
   * Use `first_isolate()` to identify the first isolates of every patient [using guidelines from the CLSI](https://clsi.org/standards/products/microbiology/documents/m39/) (Clinical and Laboratory Standards Institute).
     * You can also identify first *weighted* isolates of every patient, an adjusted version of the CLSI guideline. This takes into account key antibiotics of every strain and compares them.
   * Use `mdro()` to determine which micro-organisms are multi-drug resistant organisms (MDRO). It supports a variety of international guidelines, such as the MDR-paper by Magiorakos *et al.* (2012, [PMID 21793988](https://www.ncbi.nlm.nih.gov/pubmed/?term=21793988)), the exceptional phenotype definitions of EUCAST and the WHO guideline on multi-drug resistant TB. It also supports the national guidelines of the Netherlands and Germany.
   * The [data set microorganisms](./reference/microorganisms.html) contains the complete taxonomic tree of ~70,000 microorganisms. Furthermore, some colloquial names and all Gram stains are available, which enables resistance analysis of e.g. different antibiotics per Gram stain. The package also contains functions to look up values in this data set like `mo_genus()`, `mo_family()`, `mo_gramstain()` or even `mo_phylum()`. Use `mo_snomed()` to look up any SNOMED CT code associated with a microorganism. As all these function use `as.mo()` internally, they also use the same intelligent rules for determination. For example, `mo_genus("MRSA")` and `mo_genus("S. aureus")` will both return `"Staphylococcus"`. They also come with support for German, Danish, Dutch, Spanish, Italian, French and Portuguese. These functions can be used to add new variables to your data.
   * The [data set antibiotics](./reference/antibiotics.html) contains ~450 antimicrobial drugs with their EARS-Net code, ATC code, PubChem compound ID, LOINC code, official name, common LIS codes and DDDs of both oral and parenteral administration. It also contains all (thousands of) trade names found in PubChem. Use functions like `ab_name()`, `ab_group()`, `ab_atc()`, `ab_loinc()` and `ab_tradenames()` to look up values. The `ab_*` functions use `as.ab()` internally so they support the same intelligent rules to guess the most probable result. For example, `ab_name("Fluclox")`, `ab_name("Floxapen")` and `ab_name("J01CF05")` will all return `"Flucloxacillin"`. These functions can again be used to add new variables to your data.

3. It **analyses the data** with convenient functions that use well-known methods.

   * Calculate the microbial susceptibility or resistance (and even co-resistance) with the `susceptibility()` and `resistance()` functions, or be even more specific with the `proportion_R()`, `proportion_IR()`, `proportion_I()`, `proportion_SI()` and `proportion_S()` functions. Similarly, the *number* of isolates can be determined with the `count_resistant()`, `count_susceptible()` and `count_all()` functions. All these functions can be used with the `dplyr` package (e.g. in conjunction with `summarise()`)
   * Plot AMR results with `geom_rsi()`, a function made for the `ggplot2` package
   * Predict antimicrobial resistance for the nextcoming years using logistic regression models with the `resistance_predict()` function

4. It **teaches the user** how to use all the above actions.

   * Aside from this website with many tutorials, the package itself contains extensive help pages with many examples for all functions.
   * The package also contains example data sets:
     * The [`example_isolates` data set](./reference/example_isolates.html). This data set contains 2,000 microbial isolates with their full antibiograms. It reflects reality and can be used to practice AMR data analysis.
     * The [`WHONET` data set](./reference/WHONET.html). This data set only contains fake data, but with the exact same structure as files exported by WHONET. Read more about WHONET [on its tutorial page](./articles/WHONET.html).

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
