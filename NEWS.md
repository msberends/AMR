# AMR 0.7.1.9067

### Breaking
* Determination of first isolates now **excludes** all 'unknown' microorganisms at default, i.e. microbial code `"UNKNOWN"`. They can be included with the new parameter `include_unknown`:
  ```r
  first_isolate(..., include_unknown = TRUE)
  ```
  For WHONET users, this means that all records/isolates with organism code `"con"` (*contamination*) will be excluded at default, since `as.mo("con") = "UNKNOWN"`. The function always shows a note with the number of 'unknown' microorganisms that were included or excluded.
* For code consistency, classes `ab` and `mo` will now be preserved in any subsetting or assignment. For the sake of data integrity, this means that invalid assignments will now result in `NA`:
  ```r
  # how it works in base R:
  x <- factor("A")
  x[1] <- "B"
  #> Warning message:
  #> invalid factor level, NA generated
  
  # how it now works similarly for classes 'mo' and 'ab':
  x <- as.mo("E. coli")
  x[1] <- "testvalue"
  #> Warning message:
  #> invalid microbial code, NA generated
  ```
  This is important, because a value like `"testvalue"` could never be understood by e.g. `mo_name()`, although the class would suggest a valid microbial code.
* Function `freq()` has moved to a new package, [`clean`](https://github.com/msberends/clean) ([CRAN link](https://cran.r-project.org/package=clean)), since creating frequency tables actually does not fit the scope of this package. The `freq()` function still works, since it is re-exported from the `clean` package (which will be installed automatically upon updating this `AMR` package).

### New
* Function `bug_drug_combinations()` to quickly get a `data.frame` with the antimicrobial resistance of any bug-drug combination in a data set:
  ```r
  x <- bug_drug_combinations(septic_patients)
  x
  #>      ab          mo   S  I   R total
  #> 1   AMC B_ESCHR_COL 332 74  61   467
  #> 2   AMC B_KLBSL_PNE  49  3   6    58
  #> 3   AMC B_PROTS_MIR  28  7   1    36
  #> 4   AMC B_PSDMN_AER   0  0  30    30
  #> 5   AMC B_STPHY_AUR 234  0   1   235
  ```
  You can format this to a printable format, ready for reporting or exporting to e.g. Excel with the base R `format()` function:
  ```r
  format(x, combine_IR = FALSE)
  ```
* Additional way to calculate co-resistance, i.e. when using multiple antimicrobials as input for `portion_*` functions or `count_*` functions. This can be used to determine the empiric susceptibily of a combination therapy. A new parameter `only_all_tested` (**which defaults to `FALSE`**) replaces the old `also_single_tested` and can be used to select one of the two methods to count isolates and calculate portions. The difference can be seen in this example table (which is also on the `portion` and `count` help pages), where the %SI is being determined:

  ```r
  # --------------------------------------------------------------------
  #                     only_all_tested = FALSE  only_all_tested = TRUE
  #                     -----------------------  -----------------------
  #  Drug A    Drug B   include as  include as   include as  include as
  #                     numerator   denominator  numerator   denominator
  # --------  --------  ----------  -----------  ----------  -----------
  #  S or I    S or I       X            X            X            X
  #    R       S or I       X            X            X            X
  #   <NA>     S or I       X            X            -            -
  #  S or I      R          X            X            X            X
  #    R         R          -            X            -            X
  #   <NA>       R          -            -            -            -
  #  S or I     <NA>        X            X            -            -
  #    R        <NA>        -            -            -            -
  #   <NA>      <NA>        -            -            -            -
  # --------------------------------------------------------------------
  ```
  
  Since this is a major change, usage of the old `also_single_tested` will throw an informative error that it has been replaced by `only_all_tested`.
* `tibble` printing support for classes `rsi`, `mic`, `disk`, `ab` `mo`. When using `tibble`s containing antimicrobial columns, values `S` will print in green, values `I` will print in yellow and values `R` will print in red. Microbial IDs (class `mo`) will emphasise on the genus and species, not on the kingdom.
  ```r
  # (run this on your own console, as this page does not support colour printing)
  library(dplyr)
  septic_patients %>%
    select(mo:AMC) %>% 
    as_tibble()
  ```

### Changed
* Renamed data set `septic_patients` to `example_isolates`
* Function `eucast_rules()`:
  * Fixed a bug for *Yersinia pseudotuberculosis*
  * Added more informative errors and warnings
  * Printed info now distinguishes between added and changes values
  * Using Verbose mode (i.e. `eucast_rules(..., verbose = TRUE)`) returns more informative and readable output
  * Using factors as input now adds missing factors levels when the function changes antibiotic results
* Improved the internal auto-guessing function for determining antibiotics in your data set (`AMR:::get_column_abx()`)
* Removed class `atc` - using `as.atc()` is now deprecated in favour of `ab_atc()` and this will return a character, not the `atc` class anymore
* Removed deprecated functions `abname()`, `ab_official()`, `atc_name()`, `atc_official()`, `atc_property()`, `atc_tradenames()`, `atc_trivial_nl()`
* Fix and speed improvement for `mo_shortname()`
* Algorithm improvements for `as.mo()` (by which some additions were made to the `microorganisms` data set:
  * Big improvement for misspelled input
  * These new trivial names known to the field are now understood: meningococcus, gonococcus, pneumococcus
  * Updated to the latest taxonomic data (updated to August 2019, from the International Journal of Systematic and Evolutionary Microbiology
  * Added support for Viridans Group Streptococci (VGS) and Milleri Group Streptococci (MGS)
  * Added support for 5,000 new fungi
  * Added support for unknown yeasts and fungi
* Fix for using `mo_*` functions where the coercion uncertainties and failures would not be available through `mo_uncertainties()` and `mo_failures()` anymore
* Deprecated the `country` parameter of `mdro()` in favour of the already existing `guideline` parameter to support multiple guidelines within one country
* The `name` of `RIF` is now Rifampicin instead of Rifampin
* The `antibiotics` data set is now sorted by name and all cephalosporins now have their generation between brackets
* Speed improvement for `guess_ab_col()` which is now 30 times faster for antibiotic abbreviations
* Improved `filter_ab_class()` to be more reliable and to support 5th generation cephalosporins
* Function `availability()` now uses `portion_R()` instead of `portion_IR()`, to comply with EUCAST insights

#### Other
* Added Prof Dr Casper Albers as doctoral advisor and Dr Bart Meijer, Dr Dennis Souverein and Annick Lenglet as contributors

# AMR 0.7.1

#### New
* Function `rsi_df()` to transform a `data.frame` to a data set containing only the microbial interpretation (S, I, R), the antibiotic, the percentage of S/I/R and the number of available isolates. This is a convenient combination of the existing functions `count_df()` and `portion_df()` to immediately show resistance percentages and number of available isolates:
  ```r
  septic_patients %>%
    select(AMX, CIP) %>%
    rsi_df()
  #      antibiotic  interpretation      value  isolates
  # 1   Amoxicillin              SI  0.4442636       546
  # 2   Amoxicillin               R  0.5557364       683
  # 3 Ciprofloxacin              SI  0.8381831      1181
  # 4 Ciprofloxacin               R  0.1618169       228
  ```
* Support for all scientifically published pathotypes of *E. coli* to date (that we could find). Supported are: 

  * AIEC (Adherent-Invasive *E. coli*) 
  * ATEC (Atypical Entero-pathogenic *E. coli*) 
  * DAEC (Diffusely Adhering *E. coli*) 
  * EAEC (Entero-Aggresive *E. coli*) 
  * EHEC (Entero-Haemorrhagic *E. coli*) 
  * EIEC (Entero-Invasive *E. coli*) 
  * EPEC (Entero-Pathogenic *E. coli*) 
  * ETEC (Entero-Toxigenic *E. coli*) 
  * NMEC (Neonatal Meningitis‐causing *E. coli*) 
  * STEC (Shiga-toxin producing *E. coli*) 
  * UPEC (Uropathogenic *E. coli*)
  
  All these lead to the microbial ID of *E. coli*:
  ```r
  as.mo("UPEC")
  # B_ESCHR_COL
  mo_name("UPEC")
  # "Escherichia coli"
  mo_gramstain("EHEC")
  # "Gram-negative"
  ```
* Function `mo_info()` as an analogy to `ab_info()`. The `mo_info()` prints a list with the full taxonomy, authors, and the URL to the online database of a microorganism
* Function `mo_synonyms()` to get all previously accepted taxonomic names of a microorganism

#### Changed
* Column names of output `count_df()` and `portion_df()` are now lowercase
* Fixed bug in translation of microorganism names
* Fixed bug in determining taxonomic kingdoms
* Algorithm improvements for `as.ab()` and `as.mo()` to understand even more severely misspelled input
* Function `as.ab()` now allows spaces for coercing antibiotics names
* Added `ggplot2` methods for automatically determining the scale type of classes `mo` and `ab`
* Added names of object in the header in frequency tables, even when using pipes
* Prevented `"bacteria"` from getting coerced by `as.ab()` because Bacterial is a brand name of trimethoprim (TMP)
* Fixed a bug where setting an antibiotic would not work for `eucast_rules()` and `mdro()`
* Fixed a EUCAST rule for Staphylococci, where amikacin resistance would not be inferred from tobramycin
* Removed `latest_annual_release` from the `catalogue_of_life_version()` function
* Removed antibiotic code `PVM1` from the `antibiotics` data set as this was a duplicate of `PME`
* Fixed bug where not all old taxonomic names would be printed, when using a vector as input for `as.mo()`
* Manually added *Trichomonas vaginalis* from the kingdom of Protozoa, which is missing from the Catalogue of Life
* Small improvements to `plot()` and `barplot()` for MIC and RSI classes
* Allow Catalogue of Life IDs to be coerced by `as.mo()`

#### Other
* Fixed a note thrown by CRAN tests

# AMR 0.7.0

#### New
* Support for translation of disk diffusion and MIC values to RSI values (i.e. antimicrobial interpretations). Supported guidelines are EUCAST (2011 to 2019) and CLSI (2011 to 2019). Use `as.rsi()` on an MIC value (created with `as.mic()`), a disk diffusion value (created with the new `as.disk()`) or on a complete date set containing columns with MIC or disk diffusion values.
* Function `mo_name()` as alias of `mo_fullname()`
* Added guidelines of the WHO to determine multi-drug resistance (MDR) for TB (`mdr_tb()`) and added a new vignette about MDR. Read this tutorial [here on our website](https://msberends.gitlab.io/AMR/articles/MDR.html).

#### Changed
* Fixed a critical bug in `first_isolate()` where missing species would lead to incorrect FALSEs. This bug was not present in AMR v0.5.0, but was in v0.6.0 and v0.6.1.
* Fixed a bug in `eucast_rules()` where antibiotics from WHONET software would not be recognised
* Completely reworked the `antibiotics` data set:
  * All entries now have 3 different identifiers:
    * Column `ab` contains a human readable EARS-Net code, used by ECDC and WHO/WHONET - this is the primary identifier used in this package
    * Column `atc` contains the ATC code, used by WHO/WHOCC
    * Column `cid` contains the CID code (Compound ID), used by PubChem
  * Based on the Compound ID, almost 5,000 official brand names have been added from many different countries
  * All references to antibiotics in our package now use EARS-Net codes, like `AMX` for amoxicillin
  * Functions `atc_certe`, `ab_umcg` and `atc_trivial_nl` have been removed
  * All `atc_*` functions are superceded by `ab_*` functions
  * All output will be translated by using an included translation file which [can be viewed here](https://gitlab.com/msberends/AMR/blob/master/data-raw/translations.tsv).
    
    Please [create an issue in one of our repositories](https://gitlab.com/msberends/AMR/issues/new?issue[title]=Translation%20suggestion) if you want additions in this file.
* Improvements to plotting AMR results with `ggplot_rsi()`:
  * New parameter `colours` to set the bar colours
  * New parameters `title`, `subtitle`, `caption`, `x.title` and `y.title` to set titles and axis descriptions
* Improved intelligence of looking up antibiotic columns in a data set using `guess_ab_col()`
* Added ~5,000 more old taxonomic names to the `microorganisms.old` data set, which leads to better results finding when using the `as.mo()` function
* This package now honours the new EUCAST insight (2019) that S and I are but classified as susceptible, where I is defined as 'increased exposure' and not 'intermediate' anymore. For functions like `portion_df()` and `count_df()` this means that their new parameter `combine_SI` is TRUE at default. Our plotting function `ggplot_rsi()` also reflects this change since it uses `count_df()` internally.
* The `age()` function gained a new parameter `exact` to determine ages with decimals
* Removed deprecated functions `guess_mo()`, `guess_atc()`, `EUCAST_rules()`, `interpretive_reading()`, `rsi()`
* Frequency tables (`freq()`):
  * speed improvement for microbial IDs
  * fixed factor level names for R Markdown
  * when all values are unique it now shows a message instead of a warning
  * support for boxplots:
    ```r
    septic_patients %>% 
      freq(age) %>% 
      boxplot()
    # grouped boxplots:
    septic_patients %>% 
      group_by(hospital_id) %>% 
      freq(age) %>%
      boxplot()
    ```
* Removed all hardcoded EUCAST rules and replaced them with a new reference file which [can be viewed here](https://gitlab.com/msberends/AMR/blob/master/data-raw/eucast_rules.tsv).
  
  Please [create an issue in one of our repositories](https://gitlab.com/msberends/AMR/issues/new?issue[title]=EUCAST%20edit) if you want changes in this file.
* Added ceftazidim intrinsic resistance to *Streptococci*
* Changed default settings for `age_groups()`, to let groups of fives and tens end with 100+ instead of 120+
* Fix for `freq()` for when all values are `NA`
* Fix for `first_isolate()` for when dates are missing
* Improved speed of `guess_ab_col()`
* Function `as.mo()` now gently interprets any number of whitespace characters (like tabs) as one space
* Function `as.mo()` now returns `UNKNOWN` for `"con"` (WHONET ID of 'contamination') and returns `NA` for `"xxx"`(WHONET ID of 'no growth')
* Small algorithm fix for `as.mo()`
* Removed viruses from data set `microorganisms.codes` and cleaned it up
* Fix for `mo_shortname()` where species would not be determined correctly

#### Other
* Support for R 3.6.0 and later by providing support for [staged install](https://developer.r-project.org/Blog/public/2019/02/14/staged-install/index.html)

# AMR 0.6.1

#### Changed
* Fixed a critical bug when using `eucast_rules()` with `verbose = TRUE`
* Coercion of microbial IDs are now written to the package namespace instead of the user's home folder, to comply with the CRAN policy

# AMR 0.6.0

**New website!**

We've got a new website: [https://msberends.gitlab.io/AMR](https://msberends.gitlab.io/AMR/) (built with the great [`pkgdown`](https://pkgdown.r-lib.org/))

* Contains the complete manual of this package and all of its functions with an explanation of their parameters
* Contains a comprehensive tutorial about how to conduct antimicrobial resistance analysis, import data from WHONET or SPSS and many more.

#### New
* **BREAKING**: removed deprecated functions, parameters and references to 'bactid'. Use `as.mo()` to identify an MO code.
* Catalogue of Life as a new taxonomic source for data about microorganisms, which also contains all ITIS data we used previously. The `microorganisms` data set now contains:
  * All ~55,000 (sub)species from the kingdoms of Archaea, Bacteria and Protozoa
  * All ~3,000 (sub)species from these orders of the kingdom of Fungi: Eurotiales, Onygenales, Pneumocystales, Saccharomycetales and Schizosaccharomycetales (covering at least like all species of *Aspergillus*, *Candida*, *Pneumocystis*, *Saccharomyces* and *Trichophyton*)
  * All ~2,000 (sub)species from ~100 other relevant genera, from the kingdoms of Animalia and Plantae (like *Strongyloides* and *Taenia*)
  * All ~15,000 previously accepted names of included (sub)species that have been taxonomically renamed
  * The responsible author(s) and year of scientific publication
  
    This data is updated annually - check the included version with the new function `catalogue_of_life_version()`.
  * Due to this change, some `mo` codes changed (e.g. *Streptococcus* changed from `B_STRPTC` to `B_STRPT`). A translation table is  used internally to support older microorganism IDs, so users will not notice this difference.
  * New function `mo_rank()` for the taxonomic rank (genus, species, infraspecies, etc.)
  * New function `mo_url()` to get the direct URL of a species from the Catalogue of Life
* Support for data from [WHONET](https://whonet.org/) and [EARS-Net](https://ecdc.europa.eu/en/about-us/partnerships-and-networks/disease-and-laboratory-networks/ears-net) (European Antimicrobial Resistance Surveillance Network):
  * Exported files from WHONET can be read and used in this package. For functions like `first_isolate()` and `eucast_rules()`, all parameters will be filled in automatically.
  * This package now knows all antibiotic abbrevations by EARS-Net (which are also being used by WHONET) - the `antibiotics` data set now contains a column `ears_net`.
  * The function `as.mo()` now knows all WHONET species abbreviations too, because almost 2,000 microbial abbreviations were added to the `microorganisms.codes` data set.
* New filters for antimicrobial classes. Use these functions to filter isolates on results in one of more antibiotics from a specific class:
  ```r
  filter_aminoglycosides()
  filter_carbapenems()
  filter_cephalosporins()
  filter_1st_cephalosporins()
  filter_2nd_cephalosporins()
  filter_3rd_cephalosporins()
  filter_4th_cephalosporins()
  filter_fluoroquinolones()
  filter_glycopeptides()
  filter_macrolides()
  filter_tetracyclines()
  ```
  The `antibiotics` data set will be searched, after which the input data will be checked for column names with a value in any abbreviations, codes or official names found in the `antibiotics` data set.
  For example:
  ```r
  septic_patients %>% filter_glycopeptides(result = "R")
  # Filtering on glycopeptide antibacterials: any of `vanc` or `teic` is R
  septic_patients %>% filter_glycopeptides(result = "R", scope = "all")
  # Filtering on glycopeptide antibacterials: all of `vanc` and `teic` is R
  ```
* All `ab_*` functions are deprecated and replaced by `atc_*` functions:
  ```r
  ab_property -> atc_property()
  ab_name -> atc_name()
  ab_official -> atc_official()
  ab_trivial_nl -> atc_trivial_nl()
  ab_certe -> atc_certe()
  ab_umcg -> atc_umcg()
  ab_tradenames -> atc_tradenames()
  ```
  These functions use `as.atc()` internally. The old `atc_property` has been renamed `atc_online_property()`. This is done for two reasons: firstly, not all ATC codes are of antibiotics (ab) but can also be of antivirals or antifungals. Secondly, the input must have class `atc` or must be coerable to this class. Properties of these classes should start with the same class name, analogous to `as.mo()` and e.g. `mo_genus`.
* New functions `set_mo_source()` and `get_mo_source()` to use your own predefined MO codes as input for `as.mo()` and consequently all `mo_*` functions
* Support for the upcoming [`dplyr`](https://dplyr.tidyverse.org) version 0.8.0
* New function `guess_ab_col()` to find an antibiotic column in a table
* New function `mo_failures()` to review values that could not be coerced to a valid MO code, using `as.mo()`. This latter function will now only show a maximum of 10 uncoerced values and will refer to `mo_failures()`.
* New function `mo_uncertainties()` to review values that could be coerced to a valid MO code using `as.mo()`, but with uncertainty.
* New function `mo_renamed()` to get a list of all returned values from `as.mo()` that have had taxonomic renaming
* New function `age()` to calculate the (patients) age in years
* New function `age_groups()` to split ages into custom or predefined groups (like children or elderly). This allows for easier demographic antimicrobial resistance analysis per age group.
* New function `ggplot_rsi_predict()` as well as the base R `plot()` function can now be used for resistance prediction calculated with `resistance_predict()`:
  ```r
  x <- resistance_predict(septic_patients, col_ab = "amox")
  plot(x)
  ggplot_rsi_predict(x)
  ```
* Functions `filter_first_isolate()` and `filter_first_weighted_isolate()` to shorten and fasten filtering on data sets with antimicrobial results, e.g.:
  ```r
  septic_patients %>% filter_first_isolate(...)
  # or
  filter_first_isolate(septic_patients, ...)
  ```
  is equal to:
  ```r
  septic_patients %>%
    mutate(only_firsts = first_isolate(septic_patients, ...)) %>%
    filter(only_firsts == TRUE) %>%
    select(-only_firsts)
  ```
* New function `availability()` to check the number of available (non-empty) results in a `data.frame`
* New vignettes about how to conduct AMR analysis, predict antimicrobial resistance, use the *G*-test and more. These are also available (and even easier readable) on our website: https://msberends.gitlab.io/AMR.

#### Changed
* Function `eucast_rules()`:
  * Updated EUCAST Clinical breakpoints to [version 9.0 of 1 January 2019](http://www.eucast.org/clinical_breakpoints/), the data set `septic_patients` now reflects these changes
  * Fixed a critical bug where some rules that depend on previous applied rules would not be applied adequately
  * Emphasised in manual that penicillin is meant as benzylpenicillin (ATC [J01CE01](https://www.whocc.no/atc_ddd_index/?code=J01CE01))
  * New info is returned when running this function, stating exactly what has been changed or added. Use `eucast_rules(..., verbose = TRUE)` to get a data set with all changed per bug and drug combination.
* Removed data sets `microorganisms.oldDT`, `microorganisms.prevDT`, `microorganisms.unprevDT` and `microorganismsDT` since they were no longer needed and only contained info already available in the `microorganisms` data set
* Added 65 antibiotics to the `antibiotics` data set, from the [Pharmaceuticals Community Register](http://ec.europa.eu/health/documents/community-register/html/atc.htm) of the European Commission
* Removed columns `atc_group1_nl` and `atc_group2_nl` from the `antibiotics` data set
* Functions `atc_ddd()` and `atc_groups()` have been renamed `atc_online_ddd()` and `atc_online_groups()`. The old functions are deprecated and will be removed in a future version.
* Function `guess_mo()` is now deprecated in favour of `as.mo()` and will be removed in future versions
* Function `guess_atc()` is now deprecated in favour of `as.atc()` and will be removed in future versions
* Improvements for `as.mo()`:
  * Now handles incorrect spelling, like `i` instead of `y` and `f` instead of `ph`:
    ```r
    # mo_fullname() uses as.mo() internally
    
    mo_fullname("Sthafilokockus aaureuz")
    #> [1] "Staphylococcus aureus"
    
    mo_fullname("S. klossi")
    #> [1] "Staphylococcus kloosii"
    ```
  * Uncertainty of the algorithm is now divided into four levels, 0 to 3, where the default `allow_uncertain = TRUE` is equal to uncertainty level 2. Run `?as.mo` for more info about these levels.
    ```r
    # equal:
    as.mo(..., allow_uncertain = TRUE)
    as.mo(..., allow_uncertain = 2)
    
    # also equal:
    as.mo(..., allow_uncertain = FALSE)
    as.mo(..., allow_uncertain = 0)
    ```
    Using `as.mo(..., allow_uncertain = 3)` could lead to very unreliable results.
  * Implemented the latest publication of Becker *et al.* (2019), for categorising coagulase-negative *Staphylococci*
  * All microbial IDs that found are now saved to a local file `~/.Rhistory_mo`. Use the new function `clean_mo_history()` to delete this file, which resets the algorithms.
  * Incoercible results will now be considered 'unknown', MO code `UNKNOWN`. On foreign systems, properties of these will be translated to all languages already previously supported: German, Dutch, French, Italian, Spanish and Portuguese:
    ```r
    mo_genus("qwerty", language = "es")
    # Warning: 
    # one unique value (^= 100.0%) could not be coerced and is considered 'unknown': "qwerty". Use mo_failures() to review it.
    #> [1] "(género desconocido)"
    ```
  * Fix for vector containing only empty values
  * Finds better results when input is in other languages
  * Better handling for subspecies
  * Better handling for *Salmonellae*, especially the 'city like' serovars like *Salmonella London*
  * Understanding of highly virulent *E. coli* strains like EIEC, EPEC and STEC
  * There will be looked for uncertain results at default - these results will be returned with an informative warning
  * Manual (help page) now contains more info about the algorithms
  * Progress bar will be shown when it takes more than 3 seconds to get results
  * Support for formatted console text
  * Console will return the percentage of uncoercable input
* Function `first_isolate()`:
  * Fixed a bug where distances between dates would not be calculated right - in the `septic_patients` data set this yielded a difference of 0.15% more isolates
  * Will now use a column named like "patid" for the patient ID (parameter `col_patientid`), when this parameter was left blank
  * Will now use a column named like "key(...)ab" or "key(...)antibiotics" for the key antibiotics (parameter `col_keyantibiotics()`), when this parameter was left blank
  * Removed parameter `output_logical`, the function will now always return a logical value
  * Renamed parameter `filter_specimen` to `specimen_group`, although using `filter_specimen` will still work
* A note to the manual pages of the `portion` functions, that low counts can influence the outcome and that the `portion` functions may camouflage this, since they only return the portion (albeit being dependent on the `minimum` parameter)
* Merged data sets `microorganisms.certe` and `microorganisms.umcg` into `microorganisms.codes`
* Function `mo_taxonomy()` now contains the kingdom too
* Reduce false positives for `is.rsi.eligible()` using the new `threshold` parameter
* New colours for `scale_rsi_colours()`
* Summaries of class `mo` will now return the top 3 and the unique count, e.g. using `summary(mo)`
* Small text updates to summaries of class `rsi` and `mic`
* Function `as.rsi()`:
  * Now gives a warning when inputting MIC values
  * Now accepts high and low resistance: `"HIGH S"` will return `S`
* Frequency tables (`freq()` function):
  * Support for tidyverse quasiquotation! Now you can create frequency tables of function outcomes:
    ```r
    # Determine genus of microorganisms (mo) in `septic_patients` data set:
    # OLD WAY
    septic_patients %>%
      mutate(genus = mo_genus(mo)) %>%
      freq(genus)
    # NEW WAY
    septic_patients %>% 
      freq(mo_genus(mo))
    
    # Even supports grouping variables:
    septic_patients %>%
      group_by(gender) %>% 
      freq(mo_genus(mo))
    ```
  * Header info is now available as a list, with the `header` function
  * The parameter `header` is now set to `TRUE` at default, even for markdown
  * Added header info for class `mo` to show unique count of families, genera and species
  * Now honours the `decimal.mark` setting, which just like `format` defaults to `getOption("OutDec")`
  * The new `big.mark` parameter will at default be `","` when `decimal.mark = "."` and `"."` otherwise
  * Fix for header text where all observations are `NA`
  * New parameter `droplevels` to exclude empty factor levels when input is a factor
  * Factor levels will be in header when present in input data (maximum of 5)
  * Fix for using `select()` on frequency tables
* Function `scale_y_percent()` now contains the `limits` parameter
* Automatic parameter filling for `mdro()`, `key_antibiotics()` and `eucast_rules()`
* Updated examples for resistance prediction (`resistance_predict()` function)
* Fix for `as.mic()` to support more values ending in (several) zeroes
* if using different lengths of pattern and x in `%like%`, it will now return the call

#### Other
* Updated licence text to emphasise GPL 2.0 and that this is an R package.

# AMR 0.5.0

#### New
* Repository moved to GitLab: https://gitlab.com/msberends/AMR
* Function `count_all` to get all available isolates (that like all `portion_*` and `count_*` functions also supports `summarise` and `group_by`), the old `n_rsi` is now an alias of `count_all`
* Function `get_locale` to determine language for language-dependent output for some `mo_*` functions. This is now the default value for their `language` parameter, by which the system language will be used at default.
* Data sets `microorganismsDT`, `microorganisms.prevDT`, `microorganisms.unprevDT` and `microorganisms.oldDT` to improve the speed of `as.mo`. They are for reference only, since they are primarily for internal use of `as.mo`.
* Function `read.4D` to read from the 4D database of the MMB department of the UMCG
* Functions `mo_authors` and `mo_year` to get specific values about the scientific reference of a taxonomic entry

#### Changed
* Functions `MDRO`, `BRMO`, `MRGN` and `EUCAST_exceptional_phenotypes` were renamed to `mdro`, `brmo`, `mrgn` and `eucast_exceptional_phenotypes`
* `EUCAST_rules` was renamed to `eucast_rules`, the old function still exists as a deprecated function
* Big changes to the `eucast_rules` function:
  * Now also applies rules from the EUCAST 'Breakpoint tables for bacteria', version 8.1, 2018, http://www.eucast.org/clinical_breakpoints/ (see Source of the function)
  * New parameter `rules` to specify which rules should be applied (expert rules, breakpoints, others or all)
  * New parameter `verbose` which can be set to `TRUE` to get very specific messages about which columns and rows were affected
  * Better error handling when rules cannot be applied (i.e. new values could not be inserted)
  * The number of affected values will now only be measured once per row/column combination
  * Data set `septic_patients` now reflects these changes
  * Added parameter `pipe` for piperacillin (J01CA12), also to the `mdro` function
  * Small fixes to EUCAST clinical breakpoint rules
* Added column `kingdom` to the microorganisms data set, and function `mo_kingdom` to look up values
* Tremendous speed improvement for `as.mo` (and subsequently all `mo_*` functions), as empty values wil be ignored *a priori*
* Fewer than 3 characters as input for `as.mo` will return NA
* Function `as.mo` (and all `mo_*` wrappers) now supports genus abbreviations with "species" attached
  ```r
  as.mo("E. species")        # B_ESCHR
  mo_fullname("E. spp.")     # "Escherichia species"
  as.mo("S. spp")            # B_STPHY
  mo_fullname("S. species")  # "Staphylococcus species"
  ```
* Added parameter `combine_IR` (TRUE/FALSE) to functions `portion_df` and `count_df`, to indicate that all values of I and R must be merged into one, so the output only consists of S vs. IR (susceptible vs. non-susceptible)
* Fix for `portion_*(..., as_percent = TRUE)` when minimal number of isolates would not be met
* Added parameter `also_single_tested` for `portion_*` and `count_*` functions to also include cases where not all antibiotics were tested but at least one of the tested antibiotics includes the target antimicribial interpretation, see `?portion`
* Using `portion_*` functions now throws a warning when total available isolate is below parameter `minimum`
* Functions `as.mo`, `as.rsi`, `as.mic`, `as.atc` and `freq` will not set package name as attribute anymore
* Frequency tables - `freq()`:
  * Support for grouping variables, test with:
    ```r
    septic_patients %>% 
      group_by(hospital_id) %>% 
      freq(gender)
    ```
  * Support for (un)selecting columns:
    ```r
    septic_patients %>% 
      freq(hospital_id) %>% 
      select(-count, -cum_count) # only get item, percent, cum_percent
    ```
  * Check for `hms::is.hms`
  * Now prints in markdown at default in non-interactive sessions
  * No longer adds the factor level column and sorts factors on count again
  * Support for class `difftime`
  * New parameter `na`, to choose which character to print for empty values
  * New parameter `header` to turn the header info off (default when `markdown = TRUE`)
  * New parameter `title` to manually setbthe title of the frequency table
* `first_isolate` now tries to find columns to use as input when parameters are left blank
* Improvements for MDRO algorithm (function `mdro`)
* Data set `septic_patients` is now a `data.frame`, not a tibble anymore
* Removed diacritics from all authors (columns `microorganisms$ref` and `microorganisms.old$ref`) to comply with CRAN policy to only allow ASCII characters
* Fix for `mo_property` not working properly
* Fix for `eucast_rules` where some Streptococci would become ceftazidime R in EUCAST rule 4.5
* Support for named vectors of class `mo`, useful for `top_freq()`
* `ggplot_rsi` and `scale_y_percent` have `breaks` parameter
* AI improvements for `as.mo`:
  * `"CRS"` -> *Stenotrophomonas maltophilia*
  * `"CRSM"` -> *Stenotrophomonas maltophilia*
  * `"MSSA"` -> *Staphylococcus aureus*
  * `"MSSE"` -> *Staphylococcus epidermidis*
* Fix for `join` functions
* Speed improvement for `is.rsi.eligible`, now 15-20 times faster
* In `g.test`, when `sum(x)` is below 1000 or any of the expected values is below 5, Fisher's Exact Test will be suggested
* `ab_name` will try to fall back on `as.atc` when no results are found
* Removed the addin to view data sets
* Percentages will now will rounded more logically (e.g. in `freq` function)

#### Other
* New dependency on package `crayon`, to support formatted text in the console
* Dependency `tidyr` is now mandatory (went to `Import` field) since `portion_df` and `count_df` rely on it
* Updated vignettes to comply with README


# AMR 0.4.0

#### New
* The data set `microorganisms` now contains **all microbial taxonomic data from ITIS** (kingdoms Bacteria, Fungi and Protozoa), the Integrated Taxonomy Information System, available via https://itis.gov. The data set now contains more than 18,000 microorganisms with all known bacteria, fungi and protozoa according ITIS with genus, species, subspecies, family, order, class, phylum and subkingdom. The new data set `microorganisms.old` contains all previously known taxonomic names from those kingdoms.
* New functions based on the existing function `mo_property`:
  * Taxonomic names: `mo_phylum`, `mo_class`, `mo_order`, `mo_family`, `mo_genus`, `mo_species`, `mo_subspecies`
  * Semantic names: `mo_fullname`, `mo_shortname`
  * Microbial properties: `mo_type`, `mo_gramstain`
  * Author and year: `mo_ref`
  
  They also come with support for German, Dutch, French, Italian, Spanish and Portuguese:
  ```r
  mo_gramstain("E. coli")
  # [1] "Gram negative"
  mo_gramstain("E. coli", language = "de") # German
  # [1] "Gramnegativ"
  mo_gramstain("E. coli", language = "es") # Spanish
  # [1] "Gram negativo"
  mo_fullname("S. group A", language = "pt") # Portuguese
  # [1] "Streptococcus grupo A"
  ```
  
  Furthermore, former taxonomic names will give a note about the current taxonomic name:
  ```r
  mo_gramstain("Esc blattae")
  # Note: 'Escherichia blattae' (Burgess et al., 1973) was renamed 'Shimwellia blattae' (Priest and Barker, 2010)
  # [1] "Gram negative"
  ```
* Functions `count_R`, `count_IR`, `count_I`, `count_SI` and `count_S` to selectively count resistant or susceptible isolates
  * Extra function `count_df` (which works like `portion_df`) to get all counts of S, I and R of a data set with antibiotic columns, with support for grouped variables
* Function `is.rsi.eligible` to check for columns that have valid antimicrobial results, but do not have the `rsi` class yet. Transform the columns of your raw data with: `data %>% mutate_if(is.rsi.eligible, as.rsi)`
* Functions `as.mo` and `is.mo` as replacements for `as.bactid` and `is.bactid` (since the `microoganisms` data set not only contains bacteria). These last two functions are deprecated and will be removed in a future release. The `as.mo` function determines microbial IDs using intelligent rules:
  ```r
  as.mo("E. coli")
  # [1] B_ESCHR_COL
  as.mo("MRSA")
  # [1] B_STPHY_AUR
  as.mo("S group A")
  # [1] B_STRPTC_GRA
  ```
  And with great speed too - on a quite regular Linux server from 2007 it takes us less than 0.02 seconds to transform 25,000 items:
  ```r
  thousands_of_E_colis <- rep("E. coli", 25000)
  microbenchmark::microbenchmark(as.mo(thousands_of_E_colis), unit = "s")
  # Unit: seconds
  #         min       median         max  neval
  #  0.01817717  0.01843957  0.03878077    100
  ```
* Added parameter `reference_df` for `as.mo`, so users can supply their own microbial IDs, name or codes as a reference table
* Renamed all previous references to `bactid` to `mo`, like:
  * Column names inputs of `EUCAST_rules`, `first_isolate` and `key_antibiotics`
  * Column names of datasets `microorganisms` and `septic_patients`
  * All old syntaxes will still work with this version, but will throw warnings
* Function `labels_rsi_count` to print datalabels on a RSI `ggplot2` model
* Functions `as.atc` and `is.atc` to transform/look up antibiotic ATC codes as defined by the WHO. The existing function `guess_atc` is now an alias of `as.atc`.

* Function `ab_property` and its aliases: `ab_name`, `ab_tradenames`, `ab_certe`, `ab_umcg` and `ab_trivial_nl`
* Introduction to AMR as a vignette
* Removed clipboard functions as it violated the CRAN policy
* Renamed `septic_patients$sex` to `septic_patients$gender`

#### Changed
* Added three antimicrobial agents to the `antibiotics` data set: Terbinafine (D01BA02), Rifaximin (A07AA11) and Isoconazole (D01AC05)
* Added 163 trade names to the `antibiotics` data set, it now contains 298 different trade names in total, e.g.:
  ```r
  ab_official("Bactroban")
  # [1] "Mupirocin"
  ab_name(c("Bactroban", "Amoxil", "Zithromax", "Floxapen"))
  # [1] "Mupirocin" "Amoxicillin" "Azithromycin" "Flucloxacillin"
  ab_atc(c("Bactroban", "Amoxil", "Zithromax", "Floxapen"))
  # [1] "R01AX06" "J01CA04" "J01FA10" "J01CF05"
  ```
* For `first_isolate`, rows will be ignored when there's no species available
* Function `ratio` is now deprecated and will be removed in a future release, as it is not really the scope of this package
* Fix for `as.mic` for values ending in zeroes after a real number
* Small fix where *B. fragilis* would not be found in the `microorganisms.umcg` data set
* Added `prevalence` column to the `microorganisms` data set
* Added parameters `minimum` and `as_percent` to `portion_df`
* Support for quasiquotation in the functions series `count_*` and `portions_*`, and `n_rsi`. This allows to check for more than 2 vectors or columns.
  ```r
  septic_patients %>% select(amox, cipr) %>% count_IR()
  # which is the same as:
  septic_patients %>% count_IR(amox, cipr)
  
  septic_patients %>% portion_S(amcl)
  septic_patients %>% portion_S(amcl, gent)
  septic_patients %>% portion_S(amcl, gent, pita)
  ```
* Edited `ggplot_rsi` and `geom_rsi` so they can cope with `count_df`. The new `fun` parameter has value `portion_df` at default, but can be set to `count_df`.
* Fix for `ggplot_rsi` when the `ggplot2` package was not loaded
* Added datalabels function `labels_rsi_count` to `ggplot_rsi`
* Added possibility to set any parameter to `geom_rsi` (and `ggplot_rsi`) so you can set your own preferences
* Fix for joins, where predefined suffices would not be honoured
* Added parameter `quote` to the `freq` function
* Added generic function `diff` for frequency tables
* Added longest en shortest character length in the frequency table (`freq`) header of class `character`
* Support for types (classes) list and matrix for `freq`
  ```r
  my_matrix = with(septic_patients, matrix(c(age, gender), ncol = 2))
  freq(my_matrix)
  ```
  For lists, subsetting is possible:
  ```r
  my_list = list(age = septic_patients$age, gender = septic_patients$gender)
  my_list %>% freq(age)
  my_list %>% freq(gender)
  ```

#### Other
* More unit tests to ensure better integrity of functions

# AMR 0.3.0

#### New
* **BREAKING**: `rsi_df` was removed in favour of new functions `portion_R`, `portion_IR`, `portion_I`, `portion_SI` and `portion_S` to selectively calculate resistance or susceptibility. These functions are 20 to 30 times faster than the old `rsi` function. The old function still works, but is deprecated.
  * New function `portion_df` to get all portions of S, I and R of a data set with antibiotic columns, with support for grouped variables
* **BREAKING**: the methodology for determining first weighted isolates was changed. The antibiotics that are compared between isolates (call *key antibiotics*) to include more first isolates (afterwards called first *weighted* isolates) are now as follows:
  * Universal: amoxicillin, amoxicillin/clavlanic acid, cefuroxime, piperacillin/tazobactam, ciprofloxacin,  trimethoprim/sulfamethoxazole
  * Gram-positive: vancomycin, teicoplanin, tetracycline, erythromycin, oxacillin, rifampicin
  * Gram-negative: gentamicin, tobramycin, colistin, cefotaxime, ceftazidime, meropenem
* Support for `ggplot2`
  * New functions `geom_rsi`, `facet_rsi`, `scale_y_percent`, `scale_rsi_colours` and `theme_rsi`
  * New wrapper function `ggplot_rsi` to apply all above functions on a data set:
    * `septic_patients %>% select(tobr, gent) %>% ggplot_rsi` will show portions of S, I and R immediately in a pretty plot
    * Support for grouped variables, see `?ggplot_rsi`
* Determining bacterial ID:
  * New functions `as.bactid` and `is.bactid` to transform/ look up microbial ID's.
  * The existing function `guess_bactid` is now an alias of `as.bactid`
  * New Becker classification for *Staphylococcus* to categorise them into Coagulase Negative *Staphylococci* (CoNS) and Coagulase Positve *Staphylococci* (CoPS)
  * New Lancefield classification for *Streptococcus* to categorise them into Lancefield groups
* For convience, new descriptive statistical functions `kurtosis` and `skewness` that are lacking in base R - they are generic functions and have support for vectors, data.frames and matrices
* Function `g.test` to perform the Χ<sup>2</sup> distributed [*G*-test](https://en.wikipedia.org/wiki/G-test), which use is the same as `chisq.test`
* ~~Function `ratio` to transform a vector of values to a preset ratio~~
  * ~~For example: `ratio(c(10, 500, 10), ratio = "1:2:1")` would return `130, 260, 130`~~
* Support for Addins menu in RStudio to quickly insert `%in%` or `%like%` (and give them keyboard shortcuts), or to view the datasets that come with this package
* Function `p.symbol` to transform p values to their related symbols: `0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`
* Functions `clipboard_import` and `clipboard_export` as helper functions to quickly copy and paste from/to software like Excel and SPSS. These functions use the `clipr` package, but are a little altered to also support headless Linux servers (so you can use it in RStudio Server)
* New for frequency tables (function `freq`):
  * A vignette to explain its usage
  * Support for `rsi` (antimicrobial resistance) to use as input
  * Support for `table` to use as input: `freq(table(x, y))`
  * Support for existing functions `hist` and `plot` to use a frequency table as input: `hist(freq(df$age))`
  * Support for `as.vector`, `as.data.frame`, `as_tibble` and `format`
  * Support for quasiquotation: `freq(mydata, mycolumn)` is the same as `mydata %>% freq(mycolumn)`
  * Function `top_freq` function to return the top/below *n* items as vector
  * Header of frequency tables now also show Mean Absolute Deviaton (MAD) and Interquartile Range (IQR)
  * Possibility to globally set the default for the amount of items to print, with `options(max.print.freq = n)` where *n* is your preset value

#### Changed
* Improvements for forecasting with `resistance_predict` and added more examples
* More antibiotics added as parameters for EUCAST rules
* Updated version of the `septic_patients` data set to better reflect the reality
* Pretty printing for tibbles removed as it is not really the scope of this package
* Printing of `mic` and `rsi` classes now returns all values - use `freq` to check distributions
* Improved speed of key antibiotics comparison for determining first isolates
* Column names for the `key_antibiotics` function are now generic: 6 for broadspectrum ABs, 6 for Gram-positive specific and 6 for Gram-negative specific ABs
* Speed improvement for the `abname` function
* `%like%` now supports multiple patterns
* Frequency tables are now actual `data.frame`s with altered console printing to make it look like a frequency table. Because of this, the parameter `toConsole` is not longer needed.
* Fix for `freq` where the class of an item would be lost
* Small translational improvements to the `septic_patients` dataset and the column `bactid` now has the new class `"bactid"`
* Small improvements to the `microorganisms` dataset (especially for *Salmonella*) and the column `bactid` now has the new class `"bactid"`
* Combined MIC/RSI values will now be coerced by the `rsi` and `mic` functions:
  * `as.rsi("<=0.002; S")` will return `S`
  * `as.mic("<=0.002; S")` will return `<=0.002`
* Now possible to coerce MIC values with a space between operator and value, i.e. `as.mic("<= 0.002")` now works
* Classes `rsi` and `mic` do not add the attribute `package.version` anymore
* Added `"groups"` option for `atc_property(..., property)`. It will return a vector of the ATC hierarchy as defined by the [WHO](https://www.whocc.no/atc/structure_and_principles/). The new function `atc_groups` is a convenient wrapper around this.
* Build-in host check for `atc_property` as it requires the host set by `url` to be responsive
* Improved `first_isolate` algorithm to exclude isolates where bacteria ID or genus is unavailable
* Fix for warning *hybrid evaluation forced for row_number* ([`924b62`](https://github.com/tidyverse/dplyr/commit/924b62)) from the `dplyr` package v0.7.5 and above
* Support for empty values and for 1 or 2 columns as input for `guess_bactid` (now called `as.bactid`)
  * So `yourdata %>% select(genus, species) %>% as.bactid()` now also works
* Other small fixes

#### Other
* Added integration tests (check if everything works as expected) for all releases of R 3.1 and higher
  * Linux and macOS: https://travis-ci.org/msberends/AMR
  * Windows: https://ci.appveyor.com/project/msberends/amr
* Added thesis advisors to DESCRIPTION file

# AMR 0.2.0

#### New
* Full support for Windows, Linux and macOS
* Full support for old R versions, only R-3.0.0 (April 2013) or later is needed (needed packages may have other dependencies)
* Function `n_rsi` to count cases where antibiotic test results were available, to be used in conjunction with `dplyr::summarise`, see ?rsi
* Function `guess_bactid` to **determine the ID** of a microorganism based on genus/species or known abbreviations like MRSA
* Function `guess_atc` to **determine the ATC** of an antibiotic based on name, trade name, or known abbreviations
* Function `freq` to create **frequency tables**, with additional info in a header
* Function `MDRO` to **determine Multi Drug Resistant Organisms (MDRO)** with support for country-specific guidelines.
  * [Exceptional resistances defined by EUCAST](http://www.eucast.org/expert_rules_and_intrinsic_resistance) are also supported instead of countries alone
  * Functions `BRMO` and `MRGN` are wrappers for Dutch and German guidelines, respectively
* New algorithm to determine weighted isolates, can now be `"points"` or `"keyantibiotics"`, see `?first_isolate`
* New print format for `tibble`s and `data.table`s

#### Changed
* Fixed `rsi` class for vectors that contain only invalid antimicrobial interpretations
* Renamed dataset `ablist` to `antibiotics`
* Renamed dataset `bactlist` to `microorganisms`
* Added common abbreviations and trade names to the `antibiotics` dataset
* Added more microorganisms to the `microorganisms` dataset
* Added analysis examples on help page of dataset `septic_patients`
* Added support for character vector in `join` functions
* Added warnings when a join results in more rows after than before the join
* Altered `%like%` to make it case insensitive
* For parameters of functions `first_isolate` and `EUCAST_rules` column names are now case-insensitive
* Functions `as.rsi` and `as.mic` now add the package name and version as attributes

#### Other
* Expanded `README.md` with more examples
* Added [ORCID](https://orcid.org) of authors to DESCRIPTION file
* Added unit testing with the `testthat` package
* Added build tests for Linux and macOS using Travis CI (https://travis-ci.org/msberends/AMR)
* Added line coverage checking using CodeCov (https://codecov.io/gh/msberends/AMR/tree/master/R)

# AMR 0.1.1

* `EUCAST_rules` applies for amoxicillin even if ampicillin is missing
* Edited column names to comply with GLIMS, the laboratory information system
* Added more valid MIC values
* Renamed 'Daily Defined Dose' to 'Defined Daily Dose'
* Added barplots for `rsi` and `mic` classes

# AMR 0.1.0

* First submission to CRAN.
