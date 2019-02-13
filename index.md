# `AMR` (for R) <img src="./logo.png" align="right" height="120px" />

*(<help title="Too Long, Didn't Read">TLDR</help> - to find out how to conduct AMR analysis, please [continue reading here to get started](./articles/AMR.html).*

----

`AMR` is a free and open-source [R package](https://www.r-project.org) to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and to work with microbial and antimicrobial properties by using evidence-based methods. It supports any table format, including WHONET/EARS-Net data.

After installing this package, R knows almost all ~20,000 microorganisms and ~500 antibiotics by name and code, and knows all about valid RSI and MIC values.

We created this package for both academic research and routine analysis at the Faculty of Medical Sciences of the University of Groningen and the Medical Microbiology & Infection Prevention (MMBI) department of the University Medical Center Groningen (UMCG).
This R package is actively maintained and free software; you can freely use and distribute it for both personal and commercial (but **not** patent) purposes under the terms of the GNU General Public License version 2.0 (GPL-2), as published by the Free Software Foundation. Read the full license [here](./LICENSE-text.html).

This package can be used for:

  * Calculating antimicrobial resistance
  * Calculating empirical susceptibility of both mono therapy and combination therapy
  * Predicting future antimicrobial resistance using regression models
  * Getting properties for any microorganism (like Gram stain, species, genus or family)
  * Getting properties for any antibiotic (like name, ATC code, defined daily dose or trade name)
  * Plotting antimicrobial resistance
  * Determining first isolates to be used for AMR analysis
  * Applying EUCAST rules
  * Determining multi-drug resistant organisms (MDRO)
  * Descriptive statistics: frequency tables, kurtosis and skewness

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

## Get this package

#### Latest released version

This package is available [on the official R network (CRAN)](https://cran.r-project.org/package=AMR), which has a peer-reviewed submission process. Install this package in R with:

```r
install.packages("AMR")
```

It will be downloaded and installed automatically. For RStudio, click on the menu *Tools* > *Install Packages...* and then type in "AMR" and press <kbd>Install</kbd>.

#### Latest development version

The latest and unpublished development version can be installed with (precaution: may be unstable):
```r
install.packages("devtools")
devtools::install_gitlab("msberends/AMR")
```

## Get started

To find out how to conduct AMR analysis, please [continue reading here to get started](./articles/AMR.html) or click the links in the 'How to' menu.

## Short introduction

#### WHONET / EARS-Net

<img src="./whonet.png">

We support WHONET and EARS-Net data. Exported files from WHONET can be imported into R and can be analysed easily using this package. For education purposes, we created an [example data set `WHONET`](./reference/WHONET.html) with the exact same structure as a WHONET export file. Furthermore, this package also contains a [data set `antibiotics`](./reference/antibiotics.html) with all EARS-Net antibiotic abbreviations, and knows almost all WHONET abbreviations for microorganisms. When using WHONET data as input for analysis, all input parameters will be set automatically.

Read our tutorial about [how to work with WHONET data here](./articles/WHONET.html).

#### Antimicrobial reference data

<div><img src="reference/figures/logo_who.png" height="75px" class="logo_img"><p class="logo_txt">WHO Collaborating Centre for Drug Statistics Methodology</p></div>

This package contains **all ~500 antimicrobial drugs** and their Anatomical Therapeutic Chemical (ATC) codes, ATC groups and Defined Daily Dose (DDD) from the World Health Organization Collaborating Centre for Drug Statistics Methodology (WHOCC, https://www.whocc.no) and the [Pharmaceuticals Community Register of the European Commission](http://ec.europa.eu/health/documents/community-register/html/atc.htm).

Read more about the data from WHOCC [in our manual](./reference/WHOCC.html).

#### Microbial (taxonomic) reference data

<img src="man/figures/logo_itis.jpg" height="60px">

This package contains the **complete microbial taxonomic data** (with all nine taxonomic ranks - from kingdom to subspecies) from the publicly available Integrated Taxonomic Information System (ITIS, https://www.itis.gov). 

All ~20,000 (sub)species from **the taxonomic kingdoms Bacteria, Fungi and Protozoa are included in this package**, as well as all their ~2,500 previously accepted names known to ITIS. Furthermore, the responsible authors and year of publication are available. This allows users to use authoritative taxonomic information for their data analysis on any microorganism, not only human pathogens. It also helps to quickly determine the Gram stain of bacteria, since ITIS honours the taxonomic branching order of bacterial phyla according to Cavalier-Smith (2002), which defines that all bacteria are classified into either subkingdom Negibacteria or subkingdom Posibacteria. 

Read more about the data from ITIS [in our manual](./reference/ITIS.html).

#### Overview of functions

The `AMR` package basically does four important things:

1. It **cleanses existing data** by providing new *classes* for microoganisms, antibiotics and antimicrobial results (both S/I/R and MIC). By installing this package, you teach R everything about microbiology that is needed for analysis. These functions all use artificial intelligence to guess results that you would expect:

   * Use `as.mo()` to get an ID of a microorganism. The IDs are human readable for the trained eye - the ID of *Klebsiella pneumoniae* is "B_KLBSL_PNE" (B stands for Bacteria) and the ID of *S. aureus* is "B_STPHY_AUR". The function takes almost any text as input that looks like the name or code of a microorganism like "E. coli", "esco" or "esccol" and tries to find expected results using artificial intelligence (AI) on the included ITIS data set, consisting of almost 20,000 microorganisms. It is *very* fast, please see our [benchmarks](./articles/benchmarks.html). Moreover, it can group *Staphylococci* into coagulase negative and positive (CoNS and CoPS, see [source](./reference/as.mo.html#source)) and can categorise *Streptococci* into Lancefield groups (like beta-haemolytic *Streptococcus* Group B, [source](./reference/as.mo.html#source)).
   * Use `as.rsi()` to transform values to valid antimicrobial results. It produces just S, I or R based on your input and warns about invalid values. Even values like "<=0.002; S" (combined MIC/RSI) will result in "S".
   * Use `as.mic()` to cleanse your MIC values. It produces a so-called factor (called *ordinal* in SPSS) with valid MIC values as levels. A value like "<=0.002; S" (combined MIC/RSI) will result in "<=0.002".
   * Use `as.atc()` to get the ATC code of an antibiotic as defined by the WHO. This package contains a database with most LIS codes, official names, DDDs and even trade names of antibiotics. For example, the values "Furabid", "Furadantin", "nitro" all return the ATC code of Nitrofurantoine.
   
2. It **enhances existing data** and **adds new data** from data sets included in this package.

   * Use `eucast_rules()` to apply [EUCAST expert rules to isolates](http://www.eucast.org/expert_rules_and_intrinsic_resistance/).
   * Use `first_isolate()` to identify the first isolates of every patient [using guidelines from the CLSI](https://clsi.org/standards/products/microbiology/documents/m39/) (Clinical and Laboratory Standards Institute).
     * You can also identify first *weighted* isolates of every patient, an adjusted version of the CLSI guideline. This takes into account key antibiotics of every strain and compares them.
   * Use `mdro()` (abbreviation of Multi Drug Resistant Organisms) to check your isolates for exceptional resistance with country-specific guidelines or EUCAST rules. Currently, national guidelines for Germany and the Netherlands are supported.
   * The [data set `microorganisms`](./reference/microorganisms.html) contains the complete taxonomic tree of almost 20,000 microorganisms (bacteria, fungi/yeasts and protozoa). Furthermore, the colloquial name and Gram stain are available, which enables resistance analysis of e.g. different antibiotics per Gram stain. The package also contains functions to look up values in this data set like `mo_genus()`, `mo_family()`, `mo_gramstain()` or even `mo_phylum()`. As they use `as.mo()` internally, they also use artificial intelligence. For example, `mo_genus("MRSA")` and `mo_genus("S. aureus")` will both return `"Staphylococcus"`. They also come with support for German, Dutch, Spanish, Italian, French and Portuguese. These functions can be used to add new variables to your data.
   * The [data set `antibiotics`](./reference/antibiotics.html) contains almost 500 antimicrobial drugs with their ATC code, EARS-Net code, common LIS codes, official name, trivial name and DDD of both oral and parenteral administration. It also contains hundreds of trade names. Use functions like `atc_name()` and `atc_tradenames()` to look up values. The `atc_*` functions use `as.atc()` internally so they support AI to guess your expected result. For example, `atc_name("Fluclox")`, `atc_name("Floxapen")` and `atc_name("J01CF05")` will all return `"Flucloxacillin"`. These functions can again be used to add new variables to your data.

3. It **analyses the data** with convenient functions that use well-known methods.

   * Calculate the resistance (and even co-resistance) of microbial isolates with the `portion_R()`, `portion_IR()`, `portion_I()`, `portion_SI()` and `portion_S()` functions. Similarly, the *number* of isolates can be determined with the `count_R()`, `count_IR()`, `count_I()`, `count_SI()` and `count_S()` functions. All these functions can be used with the `dplyr` package (e.g. in conjunction with `summarise()`)
   * Plot AMR results with `geom_rsi()`, a function made for the `ggplot2` package
   * Predict antimicrobial resistance for the nextcoming years using logistic regression models with the `resistance_predict()` function
   * Conduct descriptive statistics to enhance base R: calculate `kurtosis()`, `skewness()` and create frequency tables with `freq()`

4. It **teaches the user** how to use all the above actions.

   * Aside from this website with many tutorials, the package itself contains extensive help pages with many examples for all functions.
   * The package also contains example data sets:
     * The [`septic_patients` data set](.reference/septic_patients.html). This data set contains:
       * 2,000 blood culture isolates from anonymised septic patients between 2001 and 2017 in the Northern Netherlands
       * Results of 40 antibiotics (each antibiotic in its own column) with a total ~40,000 antimicrobial results
       * Real and genuine data
     * The [`WHONET` data set](.reference/WHONET.html). This data set only contains fake data, but with the exact same structure as files exported by WHONET. Read more about WHONET [on its tutorial page](./articles/WHONET.html).
    

#### Partners

The development of this package is part of, related to, or made possible by:

<a href="https://www.rug.nl"><img src="./logo_rug.png" class="partner_logo"></a>
<a href="https://www.umcg.nl"><img src="./logo_umcg.png" class="partner_logo"></a>
<a href="https://www.certe.nl"><img src="./logo_certe.png" class="partner_logo"></a>
<a href="http://www.eurhealth-1health.eu"><img src="./logo_eh1h.png" class="partner_logo"></a>
<a href="http://www.eurhealth-1health.eu"><img src="./logo_interreg.png" class="partner_logo"></a>
