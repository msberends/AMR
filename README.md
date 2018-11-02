# `AMR` <img src="man/figures/logo_amr.png" align="right" height="120px" />
### An [R package](https://www.r-project.org) to simplify the analysis and prediction of Antimicrobial Resistance (AMR) and work with antibiotic properties by using evidence-based methods.

This R package was created for academic research by PhD students of the Faculty of Medical Sciences of the [University of Groningen](https://www.rug.nl) and the Medical Microbiology & Infection Prevention (MMBI) department of the [University Medical Center Groningen (UMCG)](https://www.umcg.nl).

:arrow_forward: Get it with `install.packages("AMR")` or see below for other possibilities.

:arrow_forward: Read the [changelog here](https://gitlab.com/msberends/AMR/blob/master/NEWS.md).

## Authors
Matthijs S. Berends <a href="https://orcid.org/0000-0001-7620-1800"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,2,a</sup>,
Christian F. Luz <a href="https://orcid.org/0000-0001-5809-5995"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,a</sup>,
Erwin E.A. Hassing<sup>2</sup>,
Corinna Glasner <a href="https://orcid.org/0000-0003-1241-1328"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,b</sup>,
Alex W. Friedrich <a href="https://orcid.org/0000-0003-4881-038X"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,b</sup>,
Bhanu Sinha <a href="https://orcid.org/0000-0003-1634-0010"><img src="https://cran.r-project.org/web/orcid.svg" height="16px"></a> <sup>1,b</sup>
  
<sup>1</sup> Department of Medical Microbiology, University of Groningen, University Medical Center Groningen, Groningen, the Netherlands - [rug.nl](http://www.rug.nl) [umcg.nl](http://www.umcg.nl)<br>
<sup>2</sup> Certe Medical Diagnostics & Advice, Groningen, the Netherlands - [certe.nl](http://www.certe.nl)<br>
<sup>a</sup> R package author and thesis dissertant<br>
<sup>b</sup> Thesis advisor

<a href="https://www.rug.nl"><img src="man/figures/logo_rug.png" height="60px"></a>
<a href="https://www.umcg.nl"><img src="man/figures/logo_umcg.png" height="60px"></a>
<a href="https://www.certe.nl"><img src="man/figures/logo_certe.png" height="60px"></a>
<a href="http://www.eurhealth-1health.eu"><img src="man/figures/logo_eh1h.png" height="60px"></a>
<a href="http://www.eurhealth-1health.eu"><img src="man/figures/logo_interreg.png" height="60px"></a>

## Contents
* [Why this package?](#why-this-package)
  * [ITIS](#itis)
* [How to get it?](#how-to-get-it)
  * [Install from CRAN](#install-from-cran)
  * [Install from Zenodo](#install-from-zenodo)
  * [Install from GitLab](#install-from-gitlab)
* [How to use it?](#how-to-use-it)
  * [New classes](#new-classes)
  * [Overwrite/force resistance based on EUCAST rules](#overwriteforce-resistance-based-on-eucast-rules)
  * [Other (microbial) epidemiological functions](#other-microbial-epidemiological-functions)
  * [Frequency tables](#frequency-tables)
  * [Data sets included in package](#data-sets-included-in-package)
* [Benchmarks](#benchmarks)
* [Copyright](#copyright)

## Why this package?
This R package was intended **to make microbial epidemiology easier**. Most functions contain extensive help pages to get started.

The `AMR` package basically does four important things:

1. It **cleanses existing data**, by transforming it to reproducible and profound *classes*, making the most efficient use of R. These functions all use artificial intelligence to guess results that you would expect:

   * Use `as.mo` to get an ID of a microorganism. The IDs are human readable for the trained eye - the ID of *Klebsiella pneumoniae* is "B_KLBSL_PNE" (B stands for Bacteria) and the ID of *S. aureus* is "B_STPHY_AUR". The function takes almost any text as input that looks like the name or code of a microorganism like "E. coli", "esco" and "esccol". Even `as.mo("MRSA")` will return the ID of *S. aureus*. Moreover, it can group all coagulase negative and positive *Staphylococci*, and can transform *Streptococci* into Lancefield groups. To find bacteria based on your input, it uses Artificial Intelligence to look up values in the included ITIS data, consisting of more than 18,000 microorganisms.
   * Use `as.rsi` to transform values to valid antimicrobial results. It produces just S, I or R based on your input and warns about invalid values. Even values like "<=0.002; S" (combined MIC/RSI) will result in "S".
   * Use `as.mic` to cleanse your MIC values. It produces a so-called factor (called *ordinal* in SPSS) with valid MIC values as levels. A value like "<=0.002; S" (combined MIC/RSI) will result in "<=0.002".
   * Use `as.atc` to get the ATC code of an antibiotic as defined by the WHO. This package contains a database with most LIS codes, official names, DDDs and even trade names of antibiotics. For example, the values "Furabid", "Furadantin", "nitro" all return the ATC code of Nitrofurantoine.
   
2. It **enhances existing data** and **adds new data** from data sets included in this package.

   * Use `EUCAST_rules` to apply [EUCAST expert rules to isolates](http://www.eucast.org/expert_rules_and_intrinsic_resistance/).
   * Use `first_isolate` to identify the first isolates of every patient [using guidelines from the CLSI](https://clsi.org/standards/products/microbiology/documents/m39/) (Clinical and Laboratory Standards Institute).
     * You can also identify first *weighted* isolates of every patient, an adjusted version of the CLSI guideline. This takes into account key antibiotics of every strain and compares them.
   * Use `MDRO` (abbreviation of Multi Drug Resistant Organisms) to check your isolates for exceptional resistance with country-specific guidelines or EUCAST rules. Currently, national guidelines for Germany and the Netherlands are supported.
   * The data set `microorganisms` contains the complete taxonomic tree of more than 18,000 microorganisms (bacteria, fungi/yeasts and protozoa). Furthermore, the colloquial name and Gram stain are available, which enables resistance analysis of e.g. different antibiotics per Gram stain. The package also contains functions to look up values in this data set like `mo_genus`, `mo_family`, `mo_gramstain` or even `mo_phylum`. As they use `as.mo` internally, they also use artificial intelligence. For example, `mo_genus("MRSA")` and `mo_genus("S. aureus")` will both return `"Staphylococcus"`. They also come with support for German, Dutch, French, Italian, Spanish and Portuguese. These functions can be used to add new variables to your data.
   * The data set `antibiotics` contains the ATC code, LIS codes, official name, trivial name and DDD of both oral and parenteral administration. It also contains a total of 298 trade names. Use functions like `ab_name` and `ab_tradenames` to look up values. The `ab_*` functions use `as.atc` internally so they support AI to guess your expected result. For example, `ab_name("Fluclox")`, `ab_name("Floxapen")` and `ab_name("J01CF05")` will all return `"Flucloxacillin"`. These functions can again be used to add new variables to your data.

3. It **analyses the data** with convenient functions that use well-known methods.

   * Calculate the resistance (and even co-resistance) of microbial isolates with the `portion_R`, `portion_IR`, `portion_I`, `portion_SI` and `portion_S` functions. Similarly, the *number* of isolates can be determined with the `count_R`, `count_IR`, `count_I`, `count_SI` and `count_S` functions. All these functions can be used [with the `dplyr` package](https://dplyr.tidyverse.org/#usage) (e.g. in conjunction with [`summarise`](https://dplyr.tidyverse.org/reference/summarise.html))
   * Plot AMR results with `geom_rsi`, a function made for the `ggplot2` package
   * Predict antimicrobial resistance for the nextcoming years using logistic regression models with the `resistance_predict` function
   * Conduct descriptive statistics to enhance base R: calculate kurtosis, skewness and create frequency tables

4. It **teaches the user** how to use all the above actions.

   * The package contains extensive help pages with many examples.
   * It also contains an example data set called `septic_patients`. This data set contains:
     * 2,000 blood culture isolates from anonymised septic patients between 2001 and 2017 in the Northern Netherlands
     * Results of 40 antibiotics (each antibiotic in its own column) with a total of 38,414 antimicrobial results
     * Real and genuine data
     
### ITIS
<img src="man/figures/itis_logo.jpg" height="100px">

This package contains the **complete microbial taxonomic data** (with all  seven taxonomic ranks - from subkingdom to subspecies) from the publicly available Integrated Taxonomic Information System (ITIS, https://www.itis.gov). 

All (sub)species from the **taxonomic kingdoms Bacteria, Fungi and Protozoa are included in this package**, as well as all previously accepted names known to ITIS. Furthermore, the responsible authors and year of publication are available. This allows users to use authoritative taxonomic information for their data analysis on any microorganism, not only human pathogens. It also helps to quickly determine the Gram stain of bacteria, since all bacteria are classified into subkingdom Negibacteria or Posibacteria.

ITIS is a partnership of U.S., Canadian, and Mexican agencies and taxonomic specialists.

**Get a note when a species was renamed**
```r
mo_shortname("Chlamydia psittaci")
# Note: 'Chlamydia psittaci' (Page, 1968) was renamed 'Chlamydophila psittaci' (Everett et al., 1999)
# [1] "C. psittaci"
```

**Get any property from the entire taxonomic tree for all included species**
```r
mo_class("E. coli")
# [1] "Gammaproteobacteria"

mo_family("E. coli")
# [1] "Enterobacteriaceae"

mo_subkingdom("E. coli")
# [1] "Negibacteria"

mo_gramstain("E. coli") # based on subkingdom
# [1] "Gram negative"

mo_ref("E. coli")
# [1] "Castellani and Chalmers, 1919"
```

**Do not get mistaken - the package only includes microorganisms**
```r
mo_phylum("C. elegans")
# [1] "Cyanobacteria"                   # Bacteria?!
mo_fullname("C. elegans")
# [1] "Chroococcus limneticus elegans"  # Because a microorganism was found 
```

## How to get it?
All stable versions of this package [are published on CRAN](http://cran.r-project.org/package=AMR), the official R network with a peer-reviewed submission process.

### Install from CRAN
[![CRAN_Badge](https://www.r-pkg.org/badges/version/AMR)](http://cran.r-project.org/package=AMR) [![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/AMR)](http://cran.r-project.org/package=AMR)

(Note: Downloads measured only by [cran.rstudio.com](https://cran.rstudio.com/package=AMR), this excludes e.g. the official [cran.r-project.org](https://cran.r-project.org/package=AMR))

- <img src="http://www.rstudio.com/favicon.ico" alt="RStudio favicon" height="20px"> Install using [RStudio](http://www.rstudio.com) (recommended):
  - Click on `Tools` and then `Install Packages...`
  - Type in `AMR` and press <kbd>Install</kbd>

- <img src="https://cran.r-project.org/favicon.ico" alt="R favicon" height="20px"> Install in R directly:
  - `install.packages("AMR")`

### Install from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1305355.svg)](https://doi.org/10.5281/zenodo.1305355)

This package was also published on Zenodo (stable releases only): https://doi.org/10.5281/zenodo.1305355

### Install from GitLab

This is the latest **development version**. Although it may contain bugfixes and even new functions compared to the latest released version on CRAN, it is also subject to change and may be unstable or behave unexpectedly. Always consider this a beta version. All below 'badges' should be green:

Development Test | Result | Reference
--- | :---: | ---
All functions checked on Linux | [![pipeline status](https://gitlab.com/msberends/AMR/badges/master/pipeline.svg)](https://gitlab.com/msberends/AMR/commits/master) | GitLab CI [[ref 1]](https://gitlab.com/msberends/AMR/pipelines) 
All functions checked on Windows | [![AppVeyor_Build](https://ci.appveyor.com/api/projects/status/gitlab/msberends/AMR?branch=master&svg=true)](https://ci.appveyor.com/project/msberends/amr-svxon) | Appveyor Systems Inc. [[ref 2]](https://ci.appveyor.com/project/msberends/amr-svxon)
Percentage of syntax lines checked | [![Code_Coverage](https://codecov.io/gl/msberends/AMR/branch/master/graph/badge.svg)](https://codecov.io/gl/msberends/AMR) | Codecov LLC [[ref 3]](https://codecov.io/gl/msberends/AMR)

If so, try it with:
```r
install.packages("devtools") 
devtools::install_gitlab("msberends/AMR")
```

## How to use it?
```r
# Call it with:
library(AMR)

# For a list of functions:
help(package = "AMR")
```

### New classes
This package contains two new S3 classes: `mic` for MIC values (e.g. from Vitek or Phoenix) and `rsi` for antimicrobial drug interpretations (i.e. S, I and R). Both are actually ordered factors under the hood (an MIC of `2` being higher than `<=1` but lower than `>=32`, and for class `rsi` factors are ordered as `S < I < R`). 
Both classes have extensions for existing generic functions like `print`, `summary` and `plot`.

These functions also try to coerce valid values.

#### RSI
The `septic_patients` data set comes with antimicrobial results of more than 40 different drugs. For example, columns `amox` and `cipr` contain results of amoxicillin and ciprofloxacin, respectively.
```r
summary(septic_patients[, c("amox", "cipr")])
#      amox          cipr     
#  Mode  :rsi    Mode  :rsi   
#  <NA>  :1002   <NA>  :596   
#  Sum S :336    Sum S :1108  
#  Sum IR:662    Sum IR:296   
#  -Sum R:659    -Sum R:227   
#  -Sum I:3      -Sum I:69  
```

You can use the `plot` function from base R:
```r
plot(septic_patients$cipr)
```

![example_1_rsi](man/figures/rsi_example1.png)

Or use the `ggplot2` and `dplyr` packages to create more appealing plots:

```r
library(dplyr)
library(ggplot2)

septic_patients %>%
  select(amox, nitr, fosf, trim, cipr) %>%
  ggplot_rsi()
```

![example_2_rsi](man/figures/rsi_example2.png)

Adjust it with any parameter you know from the `ggplot2` package:

```r
septic_patients %>%
  select(amox, nitr, fosf, trim, cipr) %>%
  ggplot_rsi(datalabels = FALSE, 
             width = 0.5, colour = "purple", size = 1, linetype = 2, alpha = 0.5)
```

![example_3_rsi](man/figures/rsi_example3.png)

It also supports grouping variables. Let's say we want to compare resistance of drugs against Urine Tract Infections (UTI) between hospitals A to D (variable `hospital_id`):

```r
septic_patients %>%
  select(hospital_id, amox, nitr, fosf, trim, cipr) %>%
  group_by(hospital_id) %>%
  ggplot_rsi(x = "hospital_id",
             facet = "Antibiotic",
             nrow = 1,
             datalabels = FALSE) +
  labs(title = "AMR of Anti-UTI Drugs Per Hospital",
       x = "Hospital")
```

![example_4_rsi](man/figures/rsi_example4.png)

You could use this to group on anything in your plots: Gram stain, age (group), genus, geographic location, et cetera.

Is there a significant difference between hospital A and D when it comes to Fosfomycin?
```r
check_A_and_D <- septic_patients %>%
  filter(hospital_id %in% c("A", "D")) %>% # filter on only hospitals A and D
  select(hospital_id, fosf) %>%            # select the hospitals and fosfomycin
  group_by(hospital_id) %>%
  count_df(combine_IR = TRUE) %>%          # count all isolates per group (hospital_id)
  tidyr::spread(hospital_id, Value) %>%    # transform output so A and D are columns
  select(A, D) %>%                         # and select these only
  as.matrix()                              # transform to good old matrix for fisher.test

check_A_and_D
#       A  D
# [1,] 24 33
# [2,] 25 77
```

Total sum is lower than 1,000 so we'd prefer a [Fisher's exact test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test), not a [*G*-test](https://en.wikipedia.org/wiki/G-test) (or its formerly used equivalent, the famous [Chi<sup>2</sup> test](https://en.wikipedia.org/wiki/Chi-squared_test)):
```r
fisher.test(check_A_and_D)
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  check_A_and_D
# p-value = 0.03104
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.054283 4.735995
# sample estimates:
# odds ratio 
#   2.228006 
```

Well, there you go!

#### MIC

```r
# Transform values to new class
mic_data <- as.mic(c(">=32", "1.0", "8", "<=0.128", "8", "16", "16"))

summary(mic_data)
#  Mode:mic      
#  <NA>:0        
#  Min.:<=0.128  
#  Max.:>=32 

plot(mic_data)
```
![example_mic](man/figures/mic_example.png)


### Overwrite/force resistance based on EUCAST rules
This is also called *interpretive reading*.
```r
a <- data.frame(mo = c("Staphylococcus aureus",
                       "Enterococcus faecalis",
                       "Escherichia coli",
                       "Klebsiella pneumoniae",
                       "Pseudomonas aeruginosa"),
                vanc = "-",       # Vancomycin
                amox = "-",       # Amoxicillin
                coli = "-",       # Colistin
                cfta = "-",       # Ceftazidime
                cfur = "-",       # Cefuroxime
                peni = "S",       # Benzylpenicillin
                cfox = "S",       # Cefoxitin
                stringsAsFactors = FALSE)
                
a
#                       mo vanc amox coli cfta cfur peni cfox
# 1  Staphylococcus aureus    -    -    -    -    -    S    S
# 2  Enterococcus faecalis    -    -    -    -    -    S    S
# 3       Escherichia coli    -    -    -    -    -    S    S
# 4  Klebsiella pneumoniae    -    -    -    -    -    S    S
# 5 Pseudomonas aeruginosa    -    -    -    -    -    S    S

b <- EUCAST_rules(a) # 18 results are forced as R or S

b
#                       mo vanc amox coli cfta cfur peni cfox
# 1  Staphylococcus aureus    -    S    R    R    S    S    S
# 2  Enterococcus faecalis    -    -    R    R    R    S    R
# 3       Escherichia coli    R    -    -    -    -    R    S
# 4  Klebsiella pneumoniae    R    R    -    -    -    R    S
# 5 Pseudomonas aeruginosa    R    R    -    -    R    R    R
```

Bacteria IDs can be retrieved with the `guess_mo` function. It uses any type of info about a microorganism as input. For example, all these will return value `B_STPHY_AUR`, the ID of *S. aureus*:
```r
guess_mo("stau")
guess_mo("STAU")
guess_mo("staaur")
guess_mo("S. aureus")
guess_mo("S aureus")
guess_mo("Staphylococcus aureus")
guess_mo("MRSA") # Methicillin Resistant S. aureus
guess_mo("MSSA") # Methicillin Susceptible S. aureus
guess_mo("VISA") # Vancomycin Intermediate S. aureus
guess_mo("VRSA") # Vancomycin Resistant S. aureus
```

### Other (microbial) epidemiological functions

```r
# G-test to replace Chi squared test
g.test(...)

# Determine key antibiotic based on bacteria ID
key_antibiotics(...)

# Selection of first isolates of any patient
first_isolate(...)

# Calculate resistance levels of antibiotics, can be used with `summarise` (dplyr)
rsi(...)
# Predict resistance levels of antibiotics
rsi_predict(...)

# Get name of antibiotic by ATC code
abname(...)
abname("J01CR02", from = "atc", to = "umcg") # "AMCL"
```

### Frequency tables
Base R lacks a simple function to create frequency tables. We created such a function that works with almost all data types: `freq` (or `frequency_tbl`). It can be used in two ways:
```r
# Like base R:
freq(mydata$myvariable)

# And like tidyverse:
mydata %>% freq(myvariable)
```

Frequency are of course sorted by count at default:
```r
septic_patients %>% freq(hospital_id)
# Class:     factor (numeric)
# Length:    2000 (of which NA: 0 = 0.00%)
# Unique:    4
# 
#      Item    Count   Percent   Cum. Count   Cum. Percent
# ---  -----  ------  --------  -----------  -------------
# 1    D         762     38.1%          762          38.1%
# 2    B         663     33.1%         1425          71.2%
# 3    A         321     16.1%         1746          87.3%
# 4    C         254     12.7%         2000         100.0%
```

This can be changed with the `sort.count` parameter:
```r
septic_patients %>% freq(hospital_id, sort.count = FALSE)
# Class:     factor (numeric)
# Length:    2000 (of which NA: 0 = 0.00%)
# Unique:    4
# 
#      Item    Count   Percent   Cum. Count   Cum. Percent
# ---  -----  ------  --------  -----------  -------------
# 1    A         321     16.1%          321          16.1%
# 2    B         663     33.1%          984          49.2%
# 3    C         254     12.7%         1238          61.9%
# 4    D         762     38.1%         2000         100.0%
```

For numeric values, some extra descriptive statistics will be calculated:
```r
freq(runif(n = 10, min = 1, max = 5))
# Frequency table  
# Class:     numeric
# Length:    10 (of which NA: 0 = 0.00%)
# Unique:    10
# 
# Mean:      3.1
# Std. dev.: 1.3 (CV: 0.43, MAD: 1.8)
# Five-Num:  1.3 | 1.7 | 3.2 | 4.3 | 5.0 (IQR: 2.6, CQV: 0.43)
# Outliers:  0
# 
#           Item   Count   Percent   Cum. Count   Cum. Percent
# ---  ---------  ------  --------  -----------  -------------
# 1     1.271079       1     10.0%            1          10.0%
# 2     1.333975       1     10.0%            2          20.0%
# 3     1.714946       1     10.0%            3          30.0%
# 4     2.751871       1     10.0%            4          40.0%
# 5     3.090140       1     10.0%            5          50.0%
# 6     3.260850       1     10.0%            6          60.0%
# 7     3.824105       1     10.0%            7          70.0%
# 8     4.278028       1     10.0%            8          80.0%
# 9     4.436265       1     10.0%            9          90.0%
# 10    4.996694       1     10.0%           10         100.0%
# 
# Warning message:
# All observations are unique. 
```
Learn more about this function with:
```r
?freq
```

### Data sets included in package
Data sets to work with antibiotics and bacteria properties.
```r
# Data set with complete taxonomic trees from ITIS, containing of 
# the three kingdoms Bacteria, Fungi and Protozoa
microorganisms     # data.frame: 18,833 x 15
microorganisms.old # data.frame: 2,383 x 4

# Data set with ATC antibiotics codes, official names, trade names 
# and DDDs (oral and parenteral)
antibiotics       #  data.frame: 423 x 18

# Data set with 2000 random blood culture isolates from anonymised
# septic patients between 2001 and 2017 in 5 Dutch hospitals
septic_patients    # data.frame: 2,000 x 49

```

## Benchmarks

One of the most important features of this package is the complete microbial taxonomic database, supplied by ITIS (https://www.itis.gov). We created a function `as.mo` that transforms any user input value to a valid microbial ID by using AI (Artificial Intelligence) and based on the taxonomic tree of ITIS. 

Using the `microbenchmark` package, we can review the calculation performance of this function.

```r
library(microbenchmark)
```

In the next test, we try to 'coerce' different input values for *Staphylococcus aureus*. The actual result is the same every time: it returns its MO code `B_STPHY_AUR` (*B* stands for *Bacteria*, the taxonomic kingdom). 

But the calculation time differs a lot. Here, the AI effect can be reviewed best:

```r
microbenchmark(A = as.mo("stau"),
               B = as.mo("staaur"),
               C = as.mo("S. aureus"),
               D = as.mo("S.  aureus"),
               E = as.mo("STAAUR"),
               F = as.mo("Staphylococcus aureus"),
               G = as.mo("B_STPHY_AUR"),
               times = 10,
               unit = "ms")
# Unit: milliseconds
#  expr       min        lq       mean    median        uq        max neval
#     A 38.864859 38.923316 42.5410391 39.172790 39.394955  70.512389    10
#     B 13.912175 14.002899 14.1044062 14.084962 14.254467  14.281845    10
#     C 11.492663 11.555520 76.6953055 11.652670 11.864149 662.026786    10
#     D 11.616702 11.683261 12.1807189 11.873159 12.142327  14.761724    10
#     E 13.761108 14.012048 14.1360584 14.106509 14.293229  14.547522    10
#     F  6.743735  6.785151  6.8962407  6.871335  7.000961   7.158383    10
#     G  0.119220  0.137030  0.1411503  0.142512  0.145061   0.176909    10
```

In the table above, all measurements are shown in milliseconds (thousands of seconds), tested on a quite regular Linux server from 2007 (Core 2 Duo 2.7 GHz, 2 GB DDR2 RAM). A value of 6.9 milliseconds means it will roughly determine 144 different (unique) input values per second. It case of 39.2 milliseconds, this is only 26 input values per second. The more an input value resembles a full name (like C, D and F), the faster the result will be found. In case of G, the input is already a valid MO code, so it only almost takes no time at all (0.0001 seconds on our server).

To achieve this speed, the `as.mo` function also takes into account the prevalence of human pathogenic microorganisms. The downside is of course that less prevalent microorganisms will be determined far less faster. See this example for the ID of *Burkholderia nodosa* (`B_BRKHL_NOD`):

```r
microbenchmark(A = as.mo("buno"),
               B = as.mo("burnod"),
               C = as.mo("B. nodosa"),
               D = as.mo("B.  nodosa"),
               E = as.mo("BURNOD"),
               F = as.mo("Burkholderia nodosa"),
               G = as.mo("B_BRKHL_NOD"),
               times = 10,
               unit = "ms")
# Unit: milliseconds
#  expr        min         lq        mean      median         uq        max neval
#     A 124.175427 124.474837 125.8610536 125.3750560 126.160945 131.485994    10
#     B 154.249713 155.364729 160.9077032 156.8738940 157.136183 197.315105    10
#     C  66.066571  66.162393  66.5538611  66.4488130  66.698077  67.623404    10
#     D  86.747693  86.918665  90.7831016  87.8149725  89.440982 116.767991    10
#     E 154.863827 155.208563 162.6535954 158.4062465 168.593785 187.378088    10
#     F  32.427028  32.638648  32.9929454  32.7860475  32.992813  34.674241    10
#     G   0.213155   0.216578   0.2369226   0.2338985   0.253734   0.285581    10
```

That takes up to 11 times as much time! A value of 158.4 milliseconds means it can only determine ~6 different input values per second. We can conclude that looking up arbitrary codes of less prevalent microorganisms is the worst way to go, in terms of calculation performance.

To relieve this pitfall and further improve performance, two important calculations take almost no time at all: **repetive results** and **already precalculated results**.

Let's set up 25,000 entries of `"Staphylococcus aureus"` and check its speed:
```r
repetive_results <- rep("Staphylococcus aureus", 25000)
microbenchmark(F = as.mo(repetive_results),
               times = 10,
               unit = "ms")
# Unit: milliseconds
#  expr      min       lq     mean   median       uq      max neval
#     F 12.24381 12.34707 13.84736 12.37689 12.43266 40.36833   100
```

So transforming 25,000 times (!) `"Staphylococcus aureus"` only takes 6 ms (0.006 seconds) more than transforming it once. You only lose time on your unique input values.

What about precalculated results? This package also contains helper functions for specific microbial properties, for example `mo_fullname`. It returns the full microbial name (genus, species and possibly subspecies) and uses `as.mo` internally. If the input is however an already precalculated result, it almost doesn't take any time at all (see 'C' below):

```r
microbenchmark(A = mo_fullname("B_STPHY_AUR"),
               B = mo_fullname("S. aureus"),
               C = mo_fullname("Staphylococcus aureus"),
               times = 10,
               unit = "ms")
# Unit: milliseconds
#  expr       min        lq       mean     median        uq       max neval
#     A 11.364086 11.460537 11.5104799 11.4795330 11.524860 11.818263    10
#     B 11.976454 12.012352 12.1704592 12.0853020 12.210004 12.881737    10
#     C  0.095823  0.102528  0.1167754  0.1153785  0.132629  0.140661    10
```

So going from `mo_fullname("Staphylococcus aureus")` to `"Staphylococcus aureus"` takes 0.0001 seconds - it doesn't even start calculating *if the result would be the same as the expected resulting value*. That goes for all helper functions:

```r
microbenchmark(A = mo_species("aureus"),
               B = mo_genus("Staphylococcus"),
               C = mo_fullname("Staphylococcus aureus"),
               D = mo_family("Staphylococcaceae"),
               E = mo_order("Bacillales"),
               F = mo_class("Bacilli"),
               G = mo_phylum("Firmicutes"),
               H = mo_subkingdom("Posibacteria"),
               times = 10,
               unit = "ms")
# Unit: milliseconds
#  expr      min       lq      mean    median       uq      max neval
#     A 0.096801 0.120966 0.1264836 0.1262045 0.135773 0.158192    10
#     B 0.102807 0.123899 0.1258339 0.1286835 0.132420 0.143245    10
#     C 0.122503 0.128299 0.1374623 0.1292070 0.139683 0.187315    10
#     D 0.087372 0.093239 0.1053774 0.1026330 0.113633 0.128299    10
#     E 0.084020 0.098617 0.1124383 0.1094420 0.113423 0.178515    10
#     F 0.080667 0.085346 0.1068579 0.1128295 0.115030 0.133537    10
#     G 0.087443 0.090026 0.1030171 0.0995250 0.106369 0.152325    10
#     H 0.084648 0.103156 0.1058313 0.1095120 0.112864 0.117265    10
```

Of course, when running `mo_phylum("Firmicutes")` the function has zero knowledge about the actual microorganism, namely *S. aureus*. But since the result would be `"Firmicutes"` too, there is no point in calculating the result. And because this package 'knows' all phyla of all known microorganisms (according to ITIS), it can just return the initial value immediately.

## Copyright

This R package is licensed under the [GNU General Public License (GPL) v2.0](https://gitlab.com/msberends/AMR/blob/master/LICENSE). In a nutshell, this means that this package:

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
