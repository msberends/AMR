# `AMR` (for R) <img src="./logo.png" align="right" height="120px" />

> *18 October 2019*  
> **METHODS PAPER PREPRINTED**  
> A methods paper about this package has been preprinted at bioRxiv (DOI: 10.1101/810622). It was **updated on 18 December 2019** and in parallel sent to a journal. Please click [here for the paper on bioRxiv's publishers page](https://doi.org/10.1101/810622).

### What is `AMR` (for R)?

*(<help title="Too Long, Didn't Read">TLDR</help> - to find out how to conduct AMR analysis, please [continue reading here to get started](./articles/AMR.html).*

`AMR` is a free and open-source [R package](https://www.r-project.org) to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial data and properties, by using evidence-based methods. **Our aim is to provide a standard** for clean and reproducible antimicrobial resistance data analysis, that can therefore empower epidemiological analyses to continuously enable surveillance and treatment evaluation in any setting.

After installing this package, R knows [**~70,000 distinct microbial species**](./reference/microorganisms.html) and all [**~550 antibiotic, antimycotic and antiviral drugs**](./reference/antibiotics.html) by name and code, and knows all about valid RSI and MIC values. It supports any data format, including WHONET/EARS-Net data. 

We created this package for both routine analysis and academic research (as part of our PhD theses) at the Faculty of Medical Sciences of the University of Groningen, the Netherlands, and the Medical Microbiology & Infection Prevention (MMBI) department of the University Medical Center Groningen (UMCG). This R package is [actively maintained](./news) and is free software (see [Copyright](#copyright)).

<div class="main-content"> 
  <p>
    <a href="./countries_large.png" target="_blank"><img src="./countries.png" class="countries_map"></a>
    <strong>Used in almost 80 countries</strong><br>
    Since its first public release in early 2018, this package has been downloaded over 25,000 times from 79 countries <small>(as of December 2019, <a href="https://cran-logs.rstudio.com" target="_blank">CRAN logs</a>)</small>. Click the map to enlarge.</p><br><br>
</div>

#### Partners

The development of this package is part of, related to, or made possible by:

<div align="center">
  <a href="https://www.rug.nl" title="University of Groningen"><img src="./logo_rug.png" class="partner_logo"></a>
  <a href="https://www.umcg.nl" title="University Medical Center Groningen"><img src="./logo_umcg.png" class="partner_logo"></a>
  <a href="https://www.certe.nl" title="Certe Medical Diagnostics and Advice"><img src="./logo_certe.png" class="partner_logo"></a>
  <a href="http://www.eurhealth-1health.eu" title="EurHealth-1-Health"><img src="./logo_eh1h.png" class="partner_logo"></a>
  <a href="https://www.deutschland-nederland.eu" title="INTERREG"><img src="./logo_interreg.png" class="partner_logo"></a>
</div>

### What can you do with this package?

This package can be used for:

  * Reference for the taxonomy of microorganisms, since the package contains all microbial (sub)species from the [Catalogue of Life](http://www.catalogueoflife.org) ([manual](./reference/mo_property.html))
  * Interpreting raw MIC and disk diffusion values, based on the latest CLSI or EUCAST guidelines ([manual](./reference/as.rsi.html))
  * Determining first isolates to be used for AMR analysis ([manual](./reference/first_isolate.html))
  * Calculating antimicrobial resistance ([tutorial](./articles/AMR.html))
  * Determining multi-drug resistance (MDR) / multi-drug resistant organisms (MDRO) ([tutorial](./articles/MDR.html))
  * Calculating (empirical) susceptibility of both mono therapy and combination therapies ([tutorial](./articles/AMR.html))
  * Predicting future antimicrobial resistance using regression models ([tutorial](./articles/resistance_predict.html))
  * Getting properties for any microorganism (like Gram stain, species, genus or family) ([manual](./reference/mo_property.html))
  * Getting properties for any antibiotic (like name, EARS-Net code, ATC code, PubChem code, defined daily dose or trade name) ([manual](./reference/ab_property.html))
  * Plotting antimicrobial resistance ([tutorial](./articles/AMR.html))
  * Applying EUCAST expert rules ([manual](./reference/eucast_rules.html))

This package is ready-to-use for a professional environment by specialists in the following fields:

Medical Microbiology

  * Epidemiologists (both clinical microbiological and research)
  * Research Microbiologists
  * Biomedical Researchers
  * Research Pharmacologists
  * Data Scientists / Data Analysts
  
Veterinary Microbiology

  * Research Veterinarians
  * Veterinary Epidemiologists

Microbial Ecology

  * Soil Microbiologists
  * Extremophile Researchers
  * Astrobiologists

Developers

  * Package developers for R 
  * Software developers
  * Web application / Shiny developers

### Get this package

#### Latest released version

This package is available [here on the official R network (CRAN)](https://cran.r-project.org/package=AMR), which has a peer-reviewed submission process. Install this package in R from CRAN by using the command:

```r
install.packages("AMR")
```

It will be downloaded and installed automatically. For RStudio, click on the menu *Tools* > *Install Packages...* and then type in "AMR" and press <kbd>Install</kbd>.

**Note:** Not all functions on this website may be available in this latest release. To use all functions and data sets mentioned on this website, install the latest development version.

#### Latest development version

The latest and unpublished development version can be installed with (**precaution: may be unstable**):
```r
install.packages("devtools")
devtools::install_gitlab("msberends/AMR")
```

### Get started

To find out how to conduct AMR analysis, please [continue reading here to get started](./articles/AMR.html) or click the links in the 'How to' menu.

### Short introduction

#### Microbial (taxonomic) reference data

This package contains the complete taxonomic tree of almost all ~70,000 microorganisms from the authoritative and comprehensive Catalogue of Life (CoL, [www.catalogueoflife.org](http://www.catalogueoflife.org)). With `catalogue_of_life_version()` can be checked which version of the CoL is included in this package.

Read more about which data from the Catalogue of Life [in our manual](./reference/catalogue_of_life.html).

#### Antimicrobial reference data

This package contains **all ~550 antibiotic, antimycotic and antiviral drugs** and their Anatomical Therapeutic Chemical (ATC) codes, ATC groups and Defined Daily Dose (DDD, oral and IV) from the World Health Organization Collaborating Centre for Drug Statistics Methodology (WHOCC, https://www.whocc.no) and the [Pharmaceuticals Community Register of the European Commission](http://ec.europa.eu/health/documents/community-register/html/atc.htm).

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

   * Use `eucast_rules()` to apply [EUCAST expert rules to isolates](http://www.eucast.org/expert_rules_and_intrinsic_resistance/) (not the translation from MIC to RSI values, use `as.rsi()` for that).
   * Use `first_isolate()` to identify the first isolates of every patient [using guidelines from the CLSI](https://clsi.org/standards/products/microbiology/documents/m39/) (Clinical and Laboratory Standards Institute).
     * You can also identify first *weighted* isolates of every patient, an adjusted version of the CLSI guideline. This takes into account key antibiotics of every strain and compares them.
   * Use `mdro()` to determine which micro-organisms are multi-drug resistant organisms (MDRO). It supports a variety of international guidelines, such as the MDR-paper by Magiorakos *et al.* (2012, [PMID 21793988](https://www.ncbi.nlm.nih.gov/pubmed/?term=21793988)), the exceptional phenotype definitions of EUCAST and the WHO guideline on multi-drug resistant TB. It also supports the national guidelines of the Netherlands and Germany.
   * The [data set microorganisms](./reference/microorganisms.html) contains the complete taxonomic tree of ~70,000 microorganisms. Furthermore, some colloquial names and all Gram stains are available, which enables resistance analysis of e.g. different antibiotics per Gram stain. The package also contains functions to look up values in this data set like `mo_genus()`, `mo_family()`, `mo_gramstain()` or even `mo_phylum()`. As they use `as.mo()` internally, they also use the same intelligent rules for determination. For example, `mo_genus("MRSA")` and `mo_genus("S. aureus")` will both return `"Staphylococcus"`. They also come with support for German, Dutch, Spanish, Italian, French and Portuguese. These functions can be used to add new variables to your data.
   * The [data set antibiotics](./reference/antibiotics.html) contains ~450 antimicrobial drugs with their EARS-Net code, ATC code, PubChem compound ID, official name, common LIS codes and DDDs of both oral and parenteral administration. It also contains all (thousands of) trade names found in PubChem. Use functions like `ab_name()`, `ab_group()`, `ab_atc()` and `ab_tradenames()` to look up values. The `ab_*` functions use `as.ab()` internally so they support the same intelligent rules to guess the most probable result. For example, `ab_name("Fluclox")`, `ab_name("Floxapen")` and `ab_name("J01CF05")` will all return `"Flucloxacillin"`. These functions can again be used to add new variables to your data.

3. It **analyses the data** with convenient functions that use well-known methods.

   * Calculate the microbial susceptibility or resistance (and even co-resistance) with the `susceptibility()` and `resistance()` functions, or be even more specific with the `proportion_R()`, `proportion_IR()`, `proportion_I()`, `proportion_SI()` and `proportion_S()` functions. Similarly, the *number* of isolates can be determined with the `count_resistant()`, `count_susceptible()` and `count_all()` functions. All these functions can be used with the `dplyr` package (e.g. in conjunction with `summarise()`)
   * Plot AMR results with `geom_rsi()`, a function made for the `ggplot2` package
   * Predict antimicrobial resistance for the nextcoming years using logistic regression models with the `resistance_predict()` function

4. It **teaches the user** how to use all the above actions.

   * Aside from this website with many tutorials, the package itself contains extensive help pages with many examples for all functions.
   * The package also contains example data sets:
     * The [`example_isolates` data set](./reference/example_isolates.html). This data set contains 2,000 microbial isolates with their full antibiograms. It reflects reality and can be used to practice AMR analysis.
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
