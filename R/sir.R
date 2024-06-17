# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, Sinha BNM, Albers CJ, Glasner C    #
# (2022). AMR: An R Package for Working with Antimicrobial Resistance  #
# Data. Journal of Statistical Software, 104(3), 1-31.                 #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in The Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Translate MIC and Disk Diffusion to SIR, or Clean Existing SIR Data
#'
#' @description Clean up existing SIR values, or interpret minimum inhibitory concentration (MIC) values and disk diffusion diameters according to EUCAST or CLSI. [as.sir()] transforms the input to a new class [`sir`], which is an ordered [factor] containing the levels `S`, `SDD`, `I`, `R`, `NI`.
#' 
#' These breakpoints are currently implemented:
#' - For **clinical microbiology**: EUCAST `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "human")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "human")$guideline)))` and CLSI `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "human")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "human")$guideline)))`;
#' - For **veterinary microbiology**: EUCAST `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "animal")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "animal")$guideline)))` and CLSI `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "animal")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "animal")$guideline)))`;
#' - ECOFFs (Epidemiological cut-off values): EUCAST `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "ECOFF")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "ECOFF")$guideline)))` and CLSI `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "ECOFF")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "ECOFF")$guideline)))`.
#' 
#' All breakpoints used for interpretation are available in our [clinical_breakpoints] data set.
#' @rdname as.sir
#' @param x vector of values (for class [`mic`]: MIC values in mg/L, for class [`disk`]: a disk diffusion radius in millimetres)
#' @param mo a vector (or column name) with [character]s that can be coerced to valid microorganism codes with [as.mo()], can be left empty to determine it automatically
#' @param ab a vector (or column name) with [character]s that can be coerced to a valid antimicrobial drug code with [as.ab()]
#' @param uti (Urinary Tract Infection) a vector (or column name) with [logical]s (`TRUE` or `FALSE`) to specify whether a UTI specific interpretation from the guideline should be chosen. For using [as.sir()] on a [data.frame], this can also be a column containing [logical]s or when left blank, the data set will be searched for a column 'specimen', and rows within this column containing 'urin' (such as 'urine', 'urina') will be regarded isolates from a UTI. See *Examples*.
#' @inheritParams first_isolate
#' @param guideline defaults to EUCAST `r  max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST")$guideline)))` (the latest implemented EUCAST guideline in the [AMR::clinical_breakpoints] data set), but can be set with the [package option][AMR-options] [`AMR_guideline`][AMR-options]. Currently supports EUCAST (`r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`) and CLSI (`r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI")$guideline)))`), see *Details*.
#' @param conserve_capped_values a [logical] to indicate that MIC values starting with `">"` (but not `">="`) must always return "R" , and that MIC values starting with `"<"` (but not `"<="`) must always return "S"
#' @param add_intrinsic_resistance *(only useful when using a EUCAST guideline)* a [logical] to indicate whether intrinsic antibiotic resistance must also be considered for applicable bug-drug combinations, meaning that e.g. ampicillin will always return "R" in *Klebsiella* species. Determination is based on the [intrinsic_resistant] data set, that itself is based on `r format_eucast_version_nr(3.3)`.
#' @param include_screening a [logical] to indicate that clinical breakpoints for screening are allowed - the default is `FALSE`. Can also be set with the [package option][AMR-options] [`AMR_include_screening`][AMR-options].
#' @param include_PKPD a [logical] to indicate that PK/PD clinical breakpoints must be applied as a last resort - the default is `TRUE`. Can also be set with the [package option][AMR-options] [`AMR_include_PKPD`][AMR-options].
#' @param breakpoint_type the type of breakpoints to use, either `r vector_or(clinical_breakpoints$type)`. ECOFF stands for Epidemiological Cut-Off values. The default is `"human"`, which can also be set with the [package option][AMR-options] [`AMR_breakpoint_type`][AMR-options]. If `host` is set to values of veterinary species, this will automatically be set to `"animal"`.
#' @param host a vector (or column name) with [character]s to indicate the host. Only useful for veterinary breakpoints, as it requires `breakpoint_type = "animal"`. The values can be any text resembling the animal species, even in any of the `r length(LANGUAGES_SUPPORTED)` supported languages of this package. For foreign languages, be sure to set the language with [set_AMR_locale()] (though it will be automatically guessed based on the system language).
#' @param verbose a [logical] to indicate that all notes should be printed during interpretation of MIC values or disk diffusion values.
#' @param reference_data a [data.frame] to be used for interpretation, which defaults to the [clinical_breakpoints] data set. Changing this argument allows for using own interpretation guidelines. This argument must contain a data set that is equal in structure to the [clinical_breakpoints] data set (same column names and column types). Please note that the `guideline` argument will be ignored when `reference_data` is manually set.
#' @param threshold maximum fraction of invalid antimicrobial interpretations of `x`, see *Examples*
#' @param ... for using on a [data.frame]: names of columns to apply [as.sir()] on (supports tidy selection such as `column1:column4`). Otherwise: arguments passed on to methods.
#' @details
#' *Note: The clinical breakpoints in this package were validated through, and imported from, [WHONET](https://whonet.org). The public use of this `AMR` package has been endorsed by both CLSI and EUCAST. See [clinical_breakpoints] for more information.*
#' 
#' ### How it Works
#'
#' The [as.sir()] function can work in four ways:
#'
#' 1. For **cleaning raw / untransformed data**. The data will be cleaned to only contain valid values, namely: **S** for susceptible, **I** for intermediate or 'susceptible, increased exposure', **R** for resistant, **NI** for non-interpretable, and **SDD** for susceptible dose-dependent. Each of these can be set using a [regular expression][base::regex]. Furthermore, [as.sir()] will try its best to clean with some intelligence. For example, mixed values with SIR interpretations and MIC values such as `"<0.25; S"` will be coerced to `"S"`. Combined interpretations for multiple test methods (as seen in laboratory records) such as `"S; S"` will be coerced to `"S"`, but a value like `"S; I"` will return `NA` with a warning that the input is invalid.
#'
#' 2. For **interpreting minimum inhibitory concentration (MIC) values** according to EUCAST or CLSI. You must clean your MIC values first using [as.mic()], that also gives your columns the new data class [`mic`]. Also, be sure to have a column with microorganism names or codes. It will be found automatically, but can be set manually using the `mo` argument.
#'    * Using `dplyr`, SIR interpretation can be done very easily with either:
#'      ```r
#'      your_data %>% mutate_if(is.mic, as.sir)
#'      your_data %>% mutate(across(where(is.mic), as.sir))
#'      your_data %>% mutate_if(is.mic, as.sir, ab = "column_with_antibiotics", mo = "column_with_microorganisms")
#'      your_data %>% mutate_if(is.mic, as.sir, ab = c("cipro", "ampicillin", ...), mo = c("E. coli", "K. pneumoniae", ...))
#'      
#'      # for veterinary breakpoints, also set `host`:
#'      your_data %>% mutate_if(is.mic, as.sir, host = "column_with_animal_hosts", guideline = "CLSI")
#'      ```
#'    * Operators like "<=" will be stripped before interpretation. When using `conserve_capped_values = TRUE`, an MIC value of e.g. ">2" will always return "R", even if the breakpoint according to the chosen guideline is ">=4". This is to prevent that capped values from raw laboratory data would not be treated conservatively. The default behaviour (`conserve_capped_values = FALSE`) considers ">2" to be lower than ">=4" and might in this case return "S" or "I".
#' 3. For **interpreting disk diffusion diameters** according to EUCAST or CLSI. You must clean your disk zones first using [as.disk()], that also gives your columns the new data class [`disk`]. Also, be sure to have a column with microorganism names or codes. It will be found automatically, but can be set manually using the `mo` argument.
#'    * Using `dplyr`, SIR interpretation can be done very easily with either:
#'      ```r
#'      your_data %>% mutate_if(is.disk, as.sir)
#'      your_data %>% mutate(across(where(is.disk), as.sir))
#'      your_data %>% mutate_if(is.disk, as.sir, ab = "column_with_antibiotics", mo = "column_with_microorganisms")
#'      your_data %>% mutate_if(is.disk, as.sir, ab = c("cipro", "ampicillin", ...), mo = c("E. coli", "K. pneumoniae", ...))
#'      
#'      # for veterinary breakpoints, also set `host`:
#'      your_data %>% mutate_if(is.disk, as.sir, host = "column_with_animal_hosts", guideline = "CLSI")
#'      ```
#' 4. For **interpreting a complete data set**, with automatic determination of MIC values, disk diffusion diameters, microorganism names or codes, and antimicrobial test results. This is done very simply by running `as.sir(your_data)`.
#'
#' **For points 2, 3 and 4: Use [sir_interpretation_history()]** to retrieve a [data.frame] (or [tibble][tibble::tibble()] if the `tibble` package is installed) with all results of the last [as.sir()] call.
#'
#' ### Supported Guidelines
#'
#' For interpreting MIC values as well as disk diffusion diameters, currently implemented guidelines are for **clinical microbiology**: EUCAST `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "human")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "human")$guideline)))` and CLSI `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "human")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "human")$guideline)))`, and for **veterinary microbiology**: EUCAST `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "animal")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST" & type == "animal")$guideline)))` and CLSI `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "animal")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "animal")$guideline)))`.
#'
#' Thus, the `guideline` argument must be set to e.g., ``r paste0('"', subset(AMR::clinical_breakpoints, guideline %like% "EUCAST")$guideline[1], '"')`` or ``r paste0('"', subset(AMR::clinical_breakpoints, guideline %like% "CLSI")$guideline[1], '"')``. By simply using `"EUCAST"` (the default) or `"CLSI"` as input, the latest included version of that guideline will automatically be selected. You can set your own data set using the `reference_data` argument. The `guideline` argument will then be ignored.
#'
#' You can set the default guideline with the [package option][AMR-options] [`AMR_guideline`][AMR-options] (e.g. in your `.Rprofile` file), such as:
#'
#' ```
#'   options(AMR_guideline = "CLSI")
#'   options(AMR_guideline = "CLSI 2018")
#'   options(AMR_guideline = "EUCAST 2020")
#'   # or to reset:
#'   options(AMR_guideline = NULL)
#' ```
#' 
#' For veterinary guidelines, these might be the best options:
#' 
#' ```
#'   options(AMR_guideline = "CLSI")
#'   options(AMR_breakpoint_type = "animal")
#' ```
#'
#' ### After Interpretation
#'
#' After using [as.sir()], you can use the [eucast_rules()] defined by EUCAST to (1) apply inferred susceptibility and resistance based on results of other antimicrobials and (2) apply intrinsic resistance based on taxonomic properties of a microorganism.
#'
#' ### Machine-Readable Clinical Breakpoints
#'
#' The repository of this package [contains a machine-readable version](https://github.com/msberends/AMR/blob/main/data-raw/clinical_breakpoints.txt) of all guidelines. This is a CSV file consisting of `r format(nrow(AMR::clinical_breakpoints), big.mark = " ")` rows and `r ncol(AMR::clinical_breakpoints)` columns. This file is machine-readable, since it contains one row for every unique combination of the test method (MIC or disk diffusion), the antimicrobial drug and the microorganism. **This allows for easy implementation of these rules in laboratory information systems (LIS)**. Note that it only contains interpretation guidelines for humans - interpretation guidelines from CLSI for animals were removed.
#'
#' ### Other
#'
#' The function [is.sir()] detects if the input contains class `sir`. If the input is a [data.frame], it iterates over all columns and returns a [logical] vector.
#' 
#' The base R function [as.double()] can be used to retrieve quantitative values from a `sir` object: `"S"` = 1, `"I"`/`"SDD"` = 2, `"R"` = 3. All other values are rendered `NA` . **Note:** Do not use `as.integer()`, since that (because of how R works internally) will return the factor level indices, and not these aforementioned quantitative values.
#'
#' The function [is_sir_eligible()] returns `TRUE` when a column contains at most 5% invalid antimicrobial interpretations (not S and/or I and/or R and/or NI and/or SDD), and `FALSE` otherwise. The threshold of 5% can be set with the `threshold` argument. If the input is a [data.frame], it iterates over all columns and returns a [logical] vector.
#' @section Interpretation of SIR:
#' In 2019, the European Committee on Antimicrobial Susceptibility Testing (EUCAST) has decided to change the definitions of susceptibility testing categories S, I, and R as shown below (<https://www.eucast.org/newsiandr>):
#'
#' - **S - Susceptible, standard dosing regimen**\cr
#'   A microorganism is categorised as "Susceptible, standard dosing regimen", when there is a high likelihood of therapeutic success using a standard dosing regimen of the agent.
#' - **I - Susceptible, increased exposure** *\cr
#'   A microorganism is categorised as "Susceptible, Increased exposure*" when there is a high likelihood of therapeutic success because exposure to the agent is increased by adjusting the dosing regimen or by its concentration at the site of infection.
#' - **R = Resistant**\cr
#'   A microorganism is categorised as "Resistant" when there is a high likelihood of therapeutic failure even when there is increased exposure.
#'
#'   * *Exposure* is a function of how the mode of administration, dose, dosing interval, infusion time, as well as distribution and excretion of the antimicrobial agent will influence the infecting organism at the site of infection.
#'
#' This AMR package honours this insight. Use [susceptibility()] (equal to [proportion_SI()]) to determine antimicrobial susceptibility and [count_susceptible()] (equal to [count_SI()]) to count susceptible isolates.
#' @return Ordered [factor] with new class `sir`
#' @aliases sir
#' @export
#' @seealso [as.mic()], [as.disk()], [as.mo()]
#' @source
#' For interpretations of minimum inhibitory concentration (MIC) values and disk diffusion diameters:
#'
#' - **CLSI M39: Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data**, `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI")$guideline)))`, *Clinical and Laboratory Standards Institute* (CLSI). <https://clsi.org/standards/products/microbiology/documents/m39/>.
#' - **CLSI M100: Performance Standard for Antimicrobial Susceptibility Testing**, `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type != "animal")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type != "animal")$guideline)))`, *Clinical and Laboratory Standards Institute* (CLSI). <https://clsi.org/standards/products/microbiology/documents/m100/>.
#' - **CLSI VET01: Performance Standards for Antimicrobial Disk and Dilution Susceptibility Tests for Bacteria Isolated From Animals**, `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "animal")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "CLSI" & type == "animal")$guideline)))`, *Clinical and Laboratory Standards Institute* (CLSI). <https://clsi.org/standards/products/veterinary-medicine/documents/vet01//>.
#' - **EUCAST Breakpoint tables for interpretation of MICs and zone diameters**, `r min(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`-`r max(as.integer(gsub("[^0-9]", "", subset(AMR::clinical_breakpoints, guideline %like% "EUCAST")$guideline)))`, *European Committee on Antimicrobial Susceptibility Testing* (EUCAST). <https://www.eucast.org/clinical_breakpoints>.
#' - **WHONET** as a source for machine-reading the clinical breakpoints ([read more here](https://msberends.github.io/AMR/reference/clinical_breakpoints.html#imported-from-whonet)), 1989-`r max(as.integer(gsub("[^0-9]", "", AMR::clinical_breakpoints$guideline)))`, *WHO Collaborating Centre for Surveillance of Antimicrobial Resistance*. <https://whonet.org/>.
#' 
#' @inheritSection AMR Reference Data Publicly Available
#' @examples
#' example_isolates
#' summary(example_isolates) # see all SIR results at a glance
#'
#' # For INTERPRETING disk diffusion and MIC values -----------------------
#' 
#' # example data sets, with combined MIC values and disk zones
#' df_wide <- data.frame(
#'   microorganism = "Escherichia coli",
#'   amoxicillin = as.mic(8),
#'   cipro = as.mic(0.256),
#'   tobra = as.disk(16),
#'   genta = as.disk(18),
#'   ERY = "R"
#' )
#' df_long <- data.frame(
#'   bacteria = rep("Escherichia coli", 4),
#'   antibiotic = c("amoxicillin", "cipro", "tobra", "genta"),
#'   mics = as.mic(c(0.01, 1, 4, 8)),
#'   disks = as.disk(c(6, 10, 14, 18))
#' )
#'
#' \donttest{
#' ## Using dplyr -------------------------------------------------
#' if (require("dplyr")) {
#'   # approaches that all work without additional arguments:
#'   df_wide %>% mutate_if(is.mic, as.sir)
#'   df_wide %>% mutate_if(function(x) is.mic(x) | is.disk(x), as.sir)
#'   df_wide %>% mutate(across(where(is.mic), as.sir))
#'   df_wide %>% mutate_at(vars(amoxicillin:tobra), as.sir)
#'   df_wide %>% mutate(across(amoxicillin:tobra, as.sir))
#'   
#'   # approaches that all work with additional arguments:
#'   df_long %>%
#'     # given a certain data type, e.g. MIC values
#'     mutate_if(is.mic, as.sir,
#'               mo = "bacteria",
#'               ab = "antibiotic",
#'               guideline = "CLSI")
#'   df_long %>%
#'     mutate(across(where(is.mic),
#'                   function(x) as.sir(x,
#'                                      mo = "bacteria",
#'                                      ab = "antibiotic",
#'                                      guideline = "CLSI")))
#'   df_wide %>%
#'     # given certain columns, e.g. from 'cipro' to 'genta'
#'     mutate_at(vars(cipro:genta), as.sir,
#'               mo = "bacteria",
#'               guideline = "CLSI")
#'   df_wide %>%
#'     mutate(across(cipro:genta,
#'                        function(x) as.sir(x,
#'                                           mo = "bacteria",
#'                                           guideline = "CLSI")))
#'                        
#'   # for veterinary breakpoints, add 'host':
#'   df_long$animal_species <- c("cats", "dogs", "horses", "cattle")
#'   df_long %>%
#'     # given a certain data type, e.g. MIC values
#'     mutate_if(is.mic, as.sir,
#'               mo = "bacteria",
#'               ab = "antibiotic",
#'               host = "animal_species",
#'               guideline = "CLSI")
#'   df_long %>%
#'     mutate(across(where(is.mic),
#'                   function(x) as.sir(x,
#'                                      mo = "bacteria",
#'                                      ab = "antibiotic",
#'                                      host = "animal_species",
#'                                      guideline = "CLSI")))
#'   df_wide %>%
#'     mutate_at(vars(cipro:genta), as.sir,
#'               mo = "bacteria",
#'               ab = "antibiotic",
#'               host = "animal_species",
#'               guideline = "CLSI")
#'   df_wide %>%
#'     mutate(across(cipro:genta,
#'                        function(x) as.sir(x,
#'                                           mo = "bacteria",
#'                                           host = "animal_species",
#'                                           guideline = "CLSI")))
#'   
#'   # to include information about urinary tract infections (UTI)
#'   data.frame(mo = "E. coli",
#'              nitrofuratoin = c("<= 2", 32),
#'              from_the_bladder = c(TRUE, FALSE)) %>%
#'     as.sir(uti = "from_the_bladder")
#'
#'   data.frame(mo = "E. coli",
#'              nitrofuratoin = c("<= 2", 32),
#'              specimen = c("urine", "blood")) %>%
#'     as.sir() # automatically determines urine isolates
#'
#'   df_wide %>%
#'     mutate_at(vars(cipro:genta), as.sir, mo = "E. coli", uti = TRUE)
#' }
#'
#'
#' ## Using base R ------------------------------------------------
#'
#' as.sir(df_wide)
#'
#' # return a 'logbook' about the results:
#' sir_interpretation_history()
#'
#' # for single values
#' as.sir(
#'   x = as.mic(2),
#'   mo = as.mo("S. pneumoniae"),
#'   ab = "AMP",
#'   guideline = "EUCAST"
#' )
#'
#' as.sir(
#'   x = as.disk(18),
#'   mo = "Strep pneu", # `mo` will be coerced with as.mo()
#'   ab = "ampicillin", # and `ab` with as.ab()
#'   guideline = "EUCAST"
#' )
#'
#'
#' # For CLEANING existing SIR values ------------------------------------
#'
#' as.sir(c("S", "SDD", "I", "R", "NI", "A", "B", "C"))
#' as.sir("<= 0.002; S") # will return "S"
#' sir_data <- as.sir(c(rep("S", 474), rep("I", 36), rep("R", 370)))
#' is.sir(sir_data)
#' plot(sir_data) # for percentages
#' barplot(sir_data) # for frequencies
#' 
#' # as common in R, you can use as.integer() to return factor indices:
#' as.integer(as.sir(c("S", "SDD", "I", "R", "NI", NA)))
#' # but for computational use, as.double() will return 1 for S, 2 for I/SDD, and 3 for R:
#' as.double(as.sir(c("S", "SDD", "I", "R", "NI", NA)))
#' 
#' # the dplyr way
#' if (require("dplyr")) {
#'   example_isolates %>%
#'     mutate_at(vars(PEN:RIF), as.sir)
#'   # same:
#'   example_isolates %>%
#'     as.sir(PEN:RIF)
#'
#'   # fastest way to transform all columns with already valid AMR results to class `sir`:
#'   example_isolates %>%
#'     mutate_if(is_sir_eligible, as.sir)
#'
#'   # since dplyr 1.0.0, this can also be:
#'   # example_isolates %>%
#'   #   mutate(across(where(is_sir_eligible), as.sir))
#' }
#' }
as.sir <- function(x, ...) {
  UseMethod("as.sir")
}

as_sir_structure <- function(x) {
  structure(factor(as.character(unlist(unname(x))),
                   levels = c("S", "SDD", "I", "R", "NI"),
                   ordered = TRUE),
            class = c("sir", "ordered", "factor"))
}

#' @rdname as.sir
#' @details `NA_sir_` is a missing value of the new `sir` class, analogous to e.g. base \R's [`NA_character_`][base::NA].
#' @format NULL
#' @export
NA_sir_ <- as_sir_structure(NA_character_)

#' @rdname as.sir
#' @export
is.sir <- function(x) {
  if (inherits(x, "data.frame")) {
    unname(vapply(FUN.VALUE = logical(1), x, is.sir))
  } else {
    isTRUE(inherits(x, "sir"))
  }
}

#' @rdname as.sir
#' @export
is_sir_eligible <- function(x, threshold = 0.05) {
  meet_criteria(threshold, allow_class = "numeric", has_length = 1)

  if (inherits(x, "data.frame")) {
    # iterate this function over all columns
    return(unname(vapply(FUN.VALUE = logical(1), x, is_sir_eligible)))
  }

  stop_if(NCOL(x) > 1, "`x` must be a one-dimensional vector.")
  if (any(c(
    "numeric",
    "integer",
    "mo",
    "ab",
    "Date",
    "POSIXt",
    "raw",
    "hms",
    "mic",
    "disk"
  )
  %in% class(x))) {
    # no transformation needed
    return(FALSE)
  } else if (all(x %in% c("S", "SDD", "I", "R", "NI", NA)) & !all(is.na(x))) {
    return(TRUE)
  } else if (!any(c("S", "SDD", "I", "R", "NI") %in% x, na.rm = TRUE) & !all(is.na(x))) {
    return(FALSE)
  } else {
    x <- x[!is.na(x) & !is.null(x) & !x %in% c("", "-", "NULL")]
    if (length(x) == 0) {
      # no other values than empty
      cur_col <- get_current_column()
      if (!is.null(cur_col)) {
        ab <- suppressWarnings(as.ab(cur_col, fast_mode = TRUE, info = FALSE))
        if (!is.na(ab)) {
          # this is a valid antibiotic drug code
          message_(
            "Column '", font_bold(cur_col), "' is SIR eligible (despite only having empty values), since it seems to be ",
            ab_name(ab, language = NULL, tolower = TRUE), " (", ab, ")"
          )
          return(TRUE)
        }
      }
      # all values empty and no antibiotic col name - return FALSE
      return(FALSE)
    }
    # transform all values and see if it meets the set threshold
    checked <- suppressWarnings(as.sir(x))
    outcome <- sum(is.na(checked)) / length(x)
    outcome <= threshold
  }
}

#' @rdname as.sir
#' @export
#' @param S,I,R,NI,SDD a case-independent [regular expression][base::regex] to translate input to this result. This regular expression will be run *after* all non-letters and whitespaces are removed from the input.
# extra param: warn (logical, to never throw a warning)
as.sir.default <- function(x,
                           S = "^(S|U)+$",
                           I = "^(I)+$",
                           R = "^(R)+$",
                           NI = "^(N|NI|V)+$",
                           SDD = "^(SDD|D|H)+$",
                           ...) {
  if (inherits(x, "sir")) {
    return(as_sir_structure(x))
  }

  x.bak <- x
  x <- as.character(x) # this is needed to prevent the vctrs pkg from throwing an error

  if (inherits(x.bak, c("numeric", "integer")) && all(x %in% c(1:3, NA))) {
    # support haven package for importing e.g., from SPSS - it adds the 'labels' attribute
    lbls <- attributes(x.bak)$labels
    if (!is.null(lbls) && all(c("S", "I", "R") %in% names(lbls)) && all(c(1:3) %in% lbls)) {
      x[x.bak == 1] <- names(lbls[lbls == 1])
      x[x.bak == 2] <- names(lbls[lbls == 2])
      x[x.bak == 3] <- names(lbls[lbls == 3])
    } else {
      x[x.bak == 1] <- "S"
      x[x.bak == 2] <- "I"
      x[x.bak == 3] <- "R"
    }
  } else if (inherits(x.bak, "character") && all(x %in% c("1", "2", "3", "S", "I", "R", NA_character_))) {
    x[x.bak == "1"] <- "S"
    x[x.bak == "2"] <- "I"
    x[x.bak == "3"] <- "R"
  } else if (inherits(x.bak, "character") && all(x %in% c("1", "2", "3", "4", "5", "S", "SDD", "I", "R", "NI", NA_character_))) {
    x[x.bak == "1"] <- "S"
    x[x.bak == "2"] <- "SDD"
    x[x.bak == "3"] <- "I"
    x[x.bak == "4"] <- "R"
    x[x.bak == "5"] <- "NI"
  } else if (!all(is.na(x)) && !identical(levels(x), c("S", "SDD", "I", "R", "NI")) && !all(x %in% c("S", "SDD", "I", "R", "NI", NA))) {
    if (all(x %unlike% "(S|I|R)", na.rm = TRUE)) {
      # check if they are actually MICs or disks
      if (all_valid_mics(x)) {
        warning_("in `as.sir()`: the input seems to contain MIC values. You can transform them with `as.mic()` before running `as.sir()` to interpret them.")
      } else if (all_valid_disks(x)) {
        warning_("in `as.sir()`: the input seems to contain disk diffusion values. You can transform them with `as.disk()` before running `as.sir()` to interpret them.")
      }
    }

    # trim leading and trailing spaces, new lines, etc.
    x <- trimws2(as.character(unlist(x)))
    x[x %in% c(NA, "", "-", "NULL")] <- NA_character_
    x.bak <- x

    na_before <- length(x[is.na(x)])

    # correct for translations
    trans_R <- unlist(TRANSLATIONS[
      which(TRANSLATIONS$pattern == "Resistant"),
      LANGUAGES_SUPPORTED[LANGUAGES_SUPPORTED %in% colnames(TRANSLATIONS)]
    ])
    trans_S <- unlist(TRANSLATIONS[
      which(TRANSLATIONS$pattern == "Susceptible"),
      LANGUAGES_SUPPORTED[LANGUAGES_SUPPORTED %in% colnames(TRANSLATIONS)]
    ])
    trans_I <- unlist(TRANSLATIONS[
      which(TRANSLATIONS$pattern %in% c("Incr. exposure", "Susceptible, incr. exp.", "Intermediate")),
      LANGUAGES_SUPPORTED[LANGUAGES_SUPPORTED %in% colnames(TRANSLATIONS)]
    ])
    x <- gsub(paste0(unique(trans_R[!is.na(trans_R)]), collapse = "|"), "R", x, ignore.case = TRUE)
    x <- gsub(paste0(unique(trans_S[!is.na(trans_S)]), collapse = "|"), "S", x, ignore.case = TRUE)
    x <- gsub(paste0(unique(trans_I[!is.na(trans_I)]), collapse = "|"), "I", x, ignore.case = TRUE)
    # replace all English textual input
    x[x %like% "([^a-z]|^)res(is(tant)?)?"] <- "R"
    x[x %like% "([^a-z]|^)sus(cep(tible)?)?"] <- "S"
    x[x %like% "not|non"] <- "NI"
    x[x %like% "([^a-z]|^)int(er(mediate)?)?|incr.*exp"] <- "I"
    x[x %like% "dose"] <- "SDD"
    x <- gsub("[^A-Z]+", "", x, perl = TRUE)
    # apply regexes set by user
    x[x %like% S] <- "S"
    x[x %like% I] <- "I"
    x[x %like% R] <- "R"
    x[x %like% NI] <- "NI"
    x[x %like% SDD] <- "SDD"
    x[!x %in% c("S", "SDD", "I", "R", "NI")] <- NA_character_
    na_after <- length(x[is.na(x) | x == ""])

    if (!isFALSE(list(...)$warn)) { # so as.sir(..., warn = FALSE) will never throw a warning
      if (na_before != na_after) {
        list_missing <- x.bak[is.na(x) & !is.na(x.bak) & x.bak != ""] %pm>%
          unique() %pm>%
          sort() %pm>%
          vector_and(quotes = TRUE)
        cur_col <- get_current_column()
        warning_("in `as.sir()`: ", na_after - na_before, " result",
          ifelse(na_after - na_before > 1, "s", ""),
          ifelse(is.null(cur_col), "", paste0(" in column '", cur_col, "'")),
          " truncated (",
          round(((na_after - na_before) / length(x)) * 100),
          "%) that were invalid antimicrobial interpretations: ",
          list_missing,
          call = FALSE
        )
      }
    }
  }

  as_sir_structure(x)
}

#' @rdname as.sir
#' @export
as.sir.mic <- function(x,
                       mo = NULL,
                       ab = deparse(substitute(x)),
                       guideline = getOption("AMR_guideline", "EUCAST"),
                       uti = NULL,
                       conserve_capped_values = FALSE,
                       add_intrinsic_resistance = FALSE,
                       reference_data = AMR::clinical_breakpoints,
                       include_screening = getOption("AMR_include_screening", FALSE),
                       include_PKPD = getOption("AMR_include_PKPD", TRUE),
                       breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                       host = NULL,
                       verbose = FALSE,
                       ...) {
  as_sir_method(
    method_short = "mic",
    method_long = "MIC values",
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    uti = uti,
    conserve_capped_values = conserve_capped_values,
    add_intrinsic_resistance = add_intrinsic_resistance,
    reference_data = reference_data,
    include_screening = include_screening,
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    host = host,
    verbose = verbose,
    ...
  )
}

#' @rdname as.sir
#' @export
as.sir.disk <- function(x,
                        mo = NULL,
                        ab = deparse(substitute(x)),
                        guideline = getOption("AMR_guideline", "EUCAST"),
                        uti = NULL,
                        add_intrinsic_resistance = FALSE,
                        reference_data = AMR::clinical_breakpoints,
                        include_screening = getOption("AMR_include_screening", FALSE),
                        include_PKPD = getOption("AMR_include_PKPD", TRUE),
                        breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                        host = NULL,
                        verbose = FALSE,
                        ...) {
  as_sir_method(
    method_short = "disk",
    method_long = "disk diffusion zones",
    x = x,
    mo = mo,
    ab = ab,
    guideline = guideline,
    uti = uti,
    conserve_capped_values = FALSE,
    add_intrinsic_resistance = add_intrinsic_resistance,
    reference_data = reference_data,
    include_screening = include_screening,
    include_PKPD = include_PKPD,
    breakpoint_type = breakpoint_type,
    host = NULL,
    verbose = verbose,
    ...
  )
}

#' @rdname as.sir
#' @export
as.sir.data.frame <- function(x,
                              ...,
                              col_mo = NULL,
                              guideline = getOption("AMR_guideline", "EUCAST"),
                              uti = NULL,
                              conserve_capped_values = FALSE,
                              add_intrinsic_resistance = FALSE,
                              reference_data = AMR::clinical_breakpoints,
                              include_screening = getOption("AMR_include_screening", FALSE),
                              include_PKPD = getOption("AMR_include_PKPD", TRUE),
                              breakpoint_type = getOption("AMR_breakpoint_type", "human"),
                              host = NULL,
                              verbose = FALSE) {
  meet_criteria(x, allow_class = "data.frame") # will also check for dimensions > 0
  meet_criteria(col_mo, allow_class = "character", is_in = colnames(x), allow_NULL = TRUE)
  meet_criteria(guideline, allow_class = "character", has_length = 1)
  meet_criteria(uti, allow_class = c("logical", "character"), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(conserve_capped_values, allow_class = "logical", has_length = 1)
  meet_criteria(add_intrinsic_resistance, allow_class = "logical", has_length = 1)
  meet_criteria(reference_data, allow_class = "data.frame")
  meet_criteria(include_screening, allow_class = "logical", has_length = 1)
  meet_criteria(include_PKPD, allow_class = "logical", has_length = 1)
  meet_criteria(breakpoint_type, allow_class = "character", is_in = reference_data$type, has_length = 1)
  meet_criteria(host, allow_class = c("character", "factor"), allow_NULL = TRUE, allow_NA = TRUE)
  meet_criteria(verbose, allow_class = "logical", has_length = 1)
  x.bak <- x
  for (i in seq_len(ncol(x))) {
    # don't keep factors, overwriting them is hard
    if (is.factor(x[, i, drop = TRUE])) {
      x[, i] <- as.character(x[, i, drop = TRUE])
    }
  }

  # -- MO
  col_mo.bak <- col_mo
  if (is.null(col_mo)) {
    col_mo <- search_type_in_df(x = x, type = "mo", info = FALSE)
  }
  
  # -- host
  if (missing(breakpoint_type) && any(host %in% AMR_env$host_preferred_order, na.rm = TRUE)) {
    message_("Assuming `breakpoint_type = \"animal\"` since `host` contains animal species.")
    breakpoint_type <- "animal"
  } else if (any(!convert_host(host) %in% c("human", "ECOFF"), na.rm = TRUE)) {
    message_("Assuming `breakpoint_type = \"animal\"`.")
    breakpoint_type <- "animal"
  }
  if (breakpoint_type == "animal") {
    if (is.null(host)) {
      host <- search_type_in_df(x = x, type = "host", add_col_prefix = FALSE)
    } else if (length(host) == 1 && as.character(host) %in% colnames(x)) {
      host <- x[[as.character(host)]]
    }
  } else {
    host <- breakpoint_type
  }
  
  # -- UTIs
  col_uti <- uti
  if (is.null(col_uti)) {
    col_uti <- search_type_in_df(x = x, type = "uti", add_col_prefix = FALSE)
  }
  if (!is.null(col_uti)) {
    if (is.logical(col_uti)) {
      # already a logical vector as input
      if (length(col_uti) == 1) {
        uti <- rep(col_uti, NROW(x))
      } else {
        uti <- col_uti
      }
    } else {
      # column found, transform to logical
      stop_if(
        length(col_uti) != 1 | !col_uti %in% colnames(x),
        "argument `uti` must be a [logical] vector, of must be a single column name of `x`"
      )
      uti <- as.logical(x[, col_uti, drop = TRUE])
    }
  } else {
    # col_uti is still NULL - look for specimen column and make logicals of the urines
    col_specimen <- suppressMessages(search_type_in_df(x = x, type = "specimen"))
    if (!is.null(col_specimen)) {
      uti <- x[, col_specimen, drop = TRUE] %like% "urin"
      values <- sort(unique(x[uti, col_specimen, drop = TRUE]))
      if (length(values) > 1) {
        plural <- c("s", "", "")
      } else {
        plural <- c("", "s", "a ")
      }
      message_(
        "Assuming value", plural[1], " ",
        vector_and(values, quotes = TRUE),
        " in column '", font_bold(col_specimen),
        "' reflect", plural[2], " ", plural[3], "urinary tract infection", plural[1],
        ".\n  Use `as.sir(uti = FALSE)` to prevent this."
      )
    } else {
      # no data about UTI's found
      uti <- NULL
    }
  }

  i <- 0
  if (tryCatch(length(list(...)) > 0, error = function(e) TRUE)) {
    sel <- colnames(pm_select(x, ...))
  } else {
    sel <- colnames(x)
  }
  if (!is.null(col_mo)) {
    sel <- sel[sel != col_mo]
  }

  ab_cols <- colnames(x)[vapply(FUN.VALUE = logical(1), x, function(y) {
    i <<- i + 1
    check <- is.mic(y) | is.disk(y)
    ab <- colnames(x)[i]
    if (!is.null(col_mo) && ab == col_mo) {
      return(FALSE)
    }
    if (!is.null(col_uti) && ab == col_uti) {
      return(FALSE)
    }
    if (length(sel) == 0 || (length(sel) > 0 && ab %in% sel)) {
      ab_coerced <- suppressWarnings(as.ab(ab))
      if (is.na(ab_coerced) || (length(sel) > 0 & !ab %in% sel)) {
        # not even a valid AB code
        return(FALSE)
      } else {
        return(TRUE)
      }
    } else {
      return(FALSE)
    }
  })]

  stop_if(
    length(ab_cols) == 0,
    "no columns with MIC values, disk zones or antibiotic column names found in this data set. Use as.mic() or as.disk() to transform antimicrobial columns."
  )
  # set type per column
  types <- character(length(ab_cols))
  types[vapply(FUN.VALUE = logical(1), x.bak[, ab_cols, drop = FALSE], is.disk)] <- "disk"
  types[vapply(FUN.VALUE = logical(1), x.bak[, ab_cols, drop = FALSE], is.mic)] <- "mic"
  types[types == "" & vapply(FUN.VALUE = logical(1), x[, ab_cols, drop = FALSE], all_valid_disks)] <- "disk"
  types[types == "" & vapply(FUN.VALUE = logical(1), x[, ab_cols, drop = FALSE], all_valid_mics)] <- "mic"
  types[types == "" & !vapply(FUN.VALUE = logical(1), x.bak[, ab_cols, drop = FALSE], is.sir)] <- "sir"
  if (any(types %in% c("mic", "disk"), na.rm = TRUE)) {
    # now we need an mo column
    stop_if(is.null(col_mo), "`col_mo` must be set")
    # if not null, we already found it, now find again so a message will show
    if (is.null(col_mo.bak)) {
      col_mo <- search_type_in_df(x = x, type = "mo")
    }
    x_mo <- as.mo(x[, col_mo, drop = TRUE])
  }

  for (i in seq_len(length(ab_cols))) {
    if (types[i] == "mic") {
      x[, ab_cols[i]] <- x %pm>%
        pm_pull(ab_cols[i]) %pm>%
        as.character() %pm>%
        as.mic() %pm>%
        as.sir(
          mo = x_mo,
          mo.bak = x[, col_mo, drop = TRUE],
          ab = ab_cols[i],
          guideline = guideline,
          uti = uti,
          conserve_capped_values = conserve_capped_values,
          add_intrinsic_resistance = add_intrinsic_resistance,
          reference_data = reference_data,
          include_screening = include_screening,
          include_PKPD = include_PKPD,
          breakpoint_type = breakpoint_type,
          host = host,
          verbose = verbose,
          is_data.frame = TRUE
        )
    } else if (types[i] == "disk") {
      x[, ab_cols[i]] <- x %pm>%
        pm_pull(ab_cols[i]) %pm>%
        as.character() %pm>%
        as.disk() %pm>%
        as.sir(
          mo = x_mo,
          mo.bak = x[, col_mo, drop = TRUE],
          ab = ab_cols[i],
          guideline = guideline,
          uti = uti,
          add_intrinsic_resistance = add_intrinsic_resistance,
          reference_data = reference_data,
          include_screening = include_screening,
          include_PKPD = include_PKPD,
          breakpoint_type = breakpoint_type,
          host = host,
          verbose = verbose,
          is_data.frame = TRUE
        )
    } else if (types[i] == "sir") {
      show_message <- FALSE
      ab <- ab_cols[i]
      ab_coerced <- suppressWarnings(as.ab(ab))
      if (!all(x[, ab_cols[i], drop = TRUE] %in% c("S", "SDD", "I", "R", "NI", NA), na.rm = TRUE)) {
        show_message <- TRUE
        # only print message if values are not already clean
        message_("Cleaning values in column '", font_bold(ab), "' (",
          ifelse(ab_coerced != toupper(ab), paste0(ab_coerced, ", "), ""),
          ab_name(ab_coerced, tolower = TRUE), ")... ",
          appendLF = FALSE,
          as_note = FALSE
        )
      } else if (!is.sir(x.bak[, ab_cols[i], drop = TRUE])) {
        show_message <- TRUE
        # only print message if class not already set
        message_("Assigning class 'sir' to already clean column '", font_bold(ab), "' (",
          ifelse(ab_coerced != toupper(ab), paste0(ab_coerced, ", "), ""),
          ab_name(ab_coerced, tolower = TRUE, language = NULL), ")... ",
          appendLF = FALSE,
          as_note = FALSE
        )
      }
      x[, ab_cols[i]] <- as.sir.default(x = as.character(x[, ab_cols[i], drop = TRUE]))
      if (show_message == TRUE) {
        message(font_green_bg(" OK "))
      }
    }
  }

  x
}

get_guideline <- function(guideline, reference_data) {
  if (!identical(reference_data, AMR::clinical_breakpoints)) {
    return(guideline)
  }
  guideline_param <- toupper(guideline)
  if (guideline_param %in% c("CLSI", "EUCAST")) {
    guideline_param <- rev(sort(subset(reference_data, guideline %like% guideline_param)$guideline))[1L]
  }
  if (guideline_param %unlike% " ") {
    # like 'EUCAST2020', should be 'EUCAST 2020'
    guideline_param <- gsub("([a-z]+)([0-9]+)", "\\1 \\2", guideline_param, ignore.case = TRUE)
  }

  stop_ifnot(guideline_param %in% reference_data$guideline,
    "invalid guideline: '", guideline,
    "'.\nValid guidelines are: ", vector_and(reference_data$guideline, quotes = TRUE, reverse = TRUE),
    call = FALSE
  )

  guideline_param
}

convert_host <- function(x, lang = get_AMR_locale()) {
  x <- trimws2(tolower(as.character(x)))
  x_out <- rep(NA_character_, length(x))
  x_out[trimws2(tolower(x)) == "human"] <- "human"
  x_out[trimws2(tolower(x)) == "ecoff"] <- "ecoff"
  # this order is based on: clinical_breakpoints |> filter(type == "animal") |> count(host, sort = TRUE)
  x_out[is.na(x_out) & (x %like% "dog|canine" | x %like% translate_AMR("dog|dogs|canine", lang))] <- "dogs"
  x_out[is.na(x_out) & (x %like% "cattle|bovine" | x %like% translate_AMR("cattle|bovine", lang))] <- "cattle"
  x_out[is.na(x_out) & (x %like% "swine|suida(e)?" | x %like% translate_AMR("swine|swines", lang))] <- "swine"
  x_out[is.na(x_out) & (x %like% "cat|feline" | x %like% translate_AMR("cat|cats|feline", lang))] <- "cats"
  x_out[is.na(x_out) & (x %like% "horse|equine" | x %like% translate_AMR("horse|horses|equine", lang))] <- "horse"
  x_out[is.na(x_out) & (x %like% "aqua|fish" | x %like% translate_AMR("aquatic|fish", lang))] <- "aquatic"
  x_out[is.na(x_out) & (x %like% "bird|chicken|poultry|avia" | x %like% translate_AMR("bird|birds|poultry", lang))] <- "poultry"
  # additional animals, not necessarily currently in breakpoint guidelines:
  x_out[is.na(x_out) & (x %like% "camel|camelid" | x %like% translate_AMR("camel|camels|camelid", lang))] <- "camels"
  x_out[is.na(x_out) & (x %like% "deer|cervine" | x %like% translate_AMR("deer|deers|cervine", lang))] <- "deer"
  x_out[is.na(x_out) & (x %like% "donkey|asinine" | x %like% translate_AMR("donkey|donkeys|asinine", lang))] <- "donkeys"
  x_out[is.na(x_out) & (x %like% "ferret|musteline" | x %like% translate_AMR("ferret|ferrets|musteline", lang))] <- "ferrets"
  x_out[is.na(x_out) & (x %like% "goat|caprine" | x %like% translate_AMR("goat|goats|caprine", lang))] <- "goats"
  x_out[is.na(x_out) & (x %like% "guinea pig|caviine" | x %like% translate_AMR("guinea pig|guinea pigs|caviine", lang))] <- "guinea pigs"
  x_out[is.na(x_out) & (x %like% "hamster|cricetine" | x %like% translate_AMR("hamster|hamsters|cricetine", lang))] <- "hamsters"
  x_out[is.na(x_out) & (x %like% "monkey|simian" | x %like% translate_AMR("monkey|monkeys|simian", lang))] <- "monkeys"
  x_out[is.na(x_out) & (x %like% "mouse|murine" | x %like% translate_AMR("mouse|mice|murine", lang))] <- "mice"
  x_out[is.na(x_out) & (x %like% "pig|porcine" | x %like% translate_AMR("pig|pigs|porcine", lang))] <- "pigs"
  x_out[is.na(x_out) & (x %like% "rabbit|leporine" | x %like% translate_AMR("rabbit|rabbits|leporine", lang))] <- "rabbits"
  x_out[is.na(x_out) & (x %like% "rat|ratine" | x %like% translate_AMR("rat|rats|ratine", lang))] <- "rats"
  x_out[is.na(x_out) & (x %like% "sheep|ovine" | x %like% translate_AMR("sheep|sheeps|ovine", lang))] <- "sheep"
  x_out[is.na(x_out) & (x %like% "snake|serpentine" | x %like% translate_AMR("snake|snakes|serpentine", lang))] <- "snakes"
  x_out[is.na(x_out) & (x %like% "turkey|meleagrine" | x %like% translate_AMR("turkey|turkeys|meleagrine", lang))] <- "turkey"
  if (message_not_thrown_before("convert_host", x) && any(is.na(x_out) & !is.na(x))) {
    warning_("The following host(s) are invalid: ", vector_and(x[is.na(x_out) & !is.na(x)]), call = FALSE, immediate = TRUE)
  }
  x_out[x_out == "ecoff"] <- "ECOFF"
  x_out
}

as_sir_method <- function(method_short,
                          method_long,
                          x,
                          mo,
                          ab,
                          guideline,
                          uti,
                          conserve_capped_values,
                          add_intrinsic_resistance,
                          reference_data,
                          include_screening,
                          include_PKPD,
                          breakpoint_type,
                          host,
                          verbose,
                          ...) {
  meet_criteria(x, allow_NA = TRUE, .call_depth = -2)
  meet_criteria(mo, allow_class = c("mo", "character"), has_length = c(1, length(x)), allow_NULL = TRUE, .call_depth = -2)
  meet_criteria(ab, allow_class = c("ab", "character"), has_length = c(1, length(x)), .call_depth = -2)
  meet_criteria(guideline, allow_class = "character", has_length = 1, .call_depth = -2)
  meet_criteria(uti, allow_class = c("logical", "character"), has_length = c(1, length(x)), allow_NULL = TRUE, allow_NA = TRUE, .call_depth = -2)
  meet_criteria(conserve_capped_values, allow_class = "logical", has_length = 1, .call_depth = -2)
  meet_criteria(add_intrinsic_resistance, allow_class = "logical", has_length = 1, .call_depth = -2)
  meet_criteria(reference_data, allow_class = "data.frame", .call_depth = -2)
  meet_criteria(include_screening, allow_class = "logical", has_length = 1, .call_depth = -2)
  meet_criteria(include_PKPD, allow_class = "logical", has_length = 1, .call_depth = -2)
  check_reference_data(reference_data, .call_depth = -2)
  meet_criteria(breakpoint_type, allow_class = "character", is_in = reference_data$type, has_length = 1, .call_depth = -2)
  meet_criteria(host, allow_class = c("character", "factor"), allow_NULL = TRUE, allow_NA = TRUE, .call_depth = -2)
  meet_criteria(verbose, allow_class = "logical", has_length = 1, .call_depth = -2)
  
  # backward compatibilty
  dots <- list(...)
  dots <- dots[which(!names(dots) %in% c("warn", "mo.bak", "is_data.frame"))]
  if (length(dots) != 0) {
    warning_("These arguments in `as.sir()` are no longer used: ", vector_and(names(dots), quotes = "`"), ".", call = FALSE)
  }
  
  guideline_coerced <- get_guideline(guideline, reference_data)
  
  if (message_not_thrown_before("as.sir", "sir_interpretation_history")) {
    message_("Run `sir_interpretation_history()` afterwards to retrieve a logbook with all the details of the breakpoint interpretations.\n\n", add_fn = font_green)
  }
  
  current_df <- tryCatch(get_current_data(NA, 0), error = function(e) NULL)
  
  # get host
  if (breakpoint_type == "animal") {
    if (is.null(host)) {
      host <- AMR_env$host_preferred_order[1]
      if (message_not_thrown_before("as.sir", "host_missing")) {
        message_("Animal hosts not set in `host`, assuming `host = \"", host, "\"`, since these have the highest breakpoint availability.\n\n")  
      }
    }
  } else {
    if (!is.null(host) && !all(toupper(as.character(host)) %in% c("HUMAN", "ECOFF"))) {
      if (message_not_thrown_before("as.sir", "assumed_breakpoint_animal")) {
        message_("Assuming `breakpoint_type = \"animal\"`, since `host` is set.", ifelse(guideline_coerced %like% "EUCAST", " Do you also need to set `guideline = \"CLSI\"`?", ""), "\n\n")
      }
      breakpoint_type <- "animal"
    } else {
      host <- breakpoint_type
    }
  }
  if (!is.null(host) && !all(toupper(as.character(host)) %in% c("HUMAN", "ECOFF"))) {
    if (!is.null(current_df) && length(host) == 1 && host %in% colnames(current_df) && any(current_df[[host]] %like% "[A-Z]", na.rm = TRUE)) {
      host <- current_df[[host]]
    } else if (length(host) != length(x)) {
      # for dplyr's across()
      cur_column_dplyr <- import_fn("cur_column", "dplyr", error_on_fail = FALSE)
      if (!is.null(cur_column_dplyr) && is.data.frame(current_df)) {
        # try to get current column, which will only be available when in across()
        host <- tryCatch(cur_column_dplyr(),
                         error = function(e) host
        )
      }
    }
  }
  host.bak <- host
  host <- convert_host(host)
  if (breakpoint_type == "animal" && message_not_thrown_before("as.sir", "host_preferred_order")) {
    message_("Please note that in the absence of specific veterinary breakpoints for certain animal hosts, breakpoints for dogs, cattle, swine, cats, horse, aquatic, and poultry, in that order, are used as substitutes.\n\n")
  }
  
  # get ab
  if (!is.null(current_df) && length(ab) == 1 && ab %in% colnames(current_df) && any(current_df[[ab]] %like% "[A-Z]", na.rm = TRUE)) {
    ab <- current_df[[ab]]
  } else if (length(ab) != length(x)) {
    # for dplyr's across()
    cur_column_dplyr <- import_fn("cur_column", "dplyr", error_on_fail = FALSE)
    if (!is.null(cur_column_dplyr) && is.data.frame(current_df)) {
      # try to get current column, which will only be available when in across()
      ab <- tryCatch(cur_column_dplyr(),
                     error = function(e) ab
      )
    }
  }
  
  # get mo
  if (!is.null(current_df) && length(mo) == 1 && mo %in% colnames(current_df)) {
    mo_var_found <- paste0(" based on column '", font_bold(mo), "'")
    mo <- current_df[[mo]]
  } else if (length(mo) != length(x)) {
    mo_var_found <- ""
    if (is.null(mo)) {
      tryCatch(
        {
          df <- get_current_data(arg_name = "mo", call = -3) # will return an error if not found
          mo <- NULL
          try(
            {
              mo <- suppressMessages(search_type_in_df(df, "mo", add_col_prefix = FALSE))
            },
            silent = TRUE
          )
          if (!is.null(df) && !is.null(mo) && is.data.frame(df)) {
            mo_var_found <- paste0(" based on column '", font_bold(mo), "'")
            mo <- df[, mo, drop = TRUE]
          }
        },
        error = function(e) {
          mo <- NULL
        }
      )
    }
  } else {
    mo_var_found <- ""
  }
  if (is.null(mo)) {
    stop_("No information was supplied about the microorganisms (missing argument `mo` and no column of class 'mo' found). See ?as.sir.\n\n",
      "To transform certain columns with e.g. mutate(), use `data %>% mutate(across(..., as.sir, mo = x))`, where x is your column with microorganisms.\n",
      "To transform all ", method_long, " in a data set, use `data %>% as.sir()` or `data %>% mutate_if(is.", method_short, ", as.sir)`.",
      call = FALSE
    )
  }
  
  # get uti
  if (!is.null(current_df) && length(uti) == 1 && uti %in% colnames(current_df)) {
    uti <- current_df[[uti]]
  } else if (length(uti) != length(x)) {
    if (is.null(uti)) {
      tryCatch(
        {
          df <- get_current_data(arg_name = "uti", call = -3) # will return an error if not found
          uti <- NULL
          try(
            {
              uti <- suppressMessages(search_type_in_df(df, "uti", add_col_prefix = FALSE))
            },
            silent = TRUE
          )
          if (!is.null(df) && !is.null(uti) && is.data.frame(df)) {
            uti <- df[, uti, drop = TRUE]
          }
        },
        error = function(e) {
          uti <- NULL
        }
      )
    }
  }

  if (length(ab) == 1 && ab %like% paste0("as.", method_short)) {
    stop_("No unambiguous name was supplied about the antibiotic (argument `ab`). See ?as.sir.", call = FALSE)
  }

  ab.bak <- ab
  ab <- suppressWarnings(as.ab(ab))
  if (!is.null(list(...)$mo.bak)) {
    mo.bak <- list(...)$mo.bak
  } else {
    mo.bak <- mo
  }
  # be sure to take current taxonomy, as the 'clinical_breakpoints' data set only contains current taxonomy
  mo <- suppressWarnings(suppressMessages(as.mo(mo, keep_synonyms = FALSE, info = FALSE)))
  if (all(is.na(ab))) {
    message_("Returning NAs for unknown antibiotic: ", vector_and(ab.bak, sort = FALSE, quotes = TRUE),
      ". Rename this column to a valid name or code, and check the output with `as.ab()`.",
      add_fn = font_red,
      as_note = FALSE
    )
    return(as.sir(rep(NA, length(x))))
  }
  if (length(mo) == 1) {
    mo <- rep(mo, length(x))
  }
  if (length(ab) == 1) {
    ab <- rep(ab, length(x))
    ab.bak <- rep(ab.bak, length(ab))
  }
  if (length(host) == 1) {
    host <- rep(host, length(x))
  }
  if (is.null(uti)) {
    uti <- NA
  }
  if (length(uti) == 1) {
    uti <- rep(uti, length(x))
  }
  uti[is.na(uti)] <- FALSE
  if (isTRUE(add_intrinsic_resistance) && guideline_coerced %unlike% "EUCAST") {
    if (message_not_thrown_before("as.sir", "intrinsic")) {
      warning_("in `as.sir()`: using 'add_intrinsic_resistance' is only useful when using EUCAST guidelines, since the rules for intrinsic resistance are based on EUCAST.")
    }
  }
  
  agent_formatted <- paste0("'", font_bold(ab.bak, collapse = NULL), "'")
  agent_name <- ab_name(ab, tolower = TRUE, language = NULL)
  same_ab <- generalise_antibiotic_name(ab) == generalise_antibiotic_name(agent_name)
  same_ab.bak <- generalise_antibiotic_name(ab.bak) == generalise_antibiotic_name(agent_name)
  agent_formatted[same_ab.bak] <- paste0(agent_formatted[same_ab.bak], " (", ab[same_ab.bak], ")")
  agent_formatted[!same_ab.bak & !same_ab] <- paste0(agent_formatted[!same_ab.bak & !same_ab],
                                                    " (", ifelse(ab.bak[!same_ab.bak & !same_ab] == ab[!same_ab.bak & !same_ab],
                                                                 "",
                                                                 paste0(ab[!same_ab.bak & !same_ab], ", ")),
                                                    agent_name[!same_ab.bak & !same_ab],
                                                    ")")
  # this intro text will also be printed in the progress bar if the `progress` package is installed
  intro_txt <- paste0("Interpreting ", method_long, ": ", ifelse(isTRUE(list(...)$is_data.frame), "column ", ""),
                      ifelse(length(unique(agent_formatted)) == 1, unique(agent_formatted), paste0(vector_and(agent_formatted, quotes = FALSE, sort = FALSE))),
                      mo_var_found,
                      ifelse(identical(reference_data, AMR::clinical_breakpoints),
                             paste0(", ", font_bold(guideline_coerced)),
                             ""),
                      "... ")

  method <- method_short

  metadata_mo <- get_mo_uncertainties()

  rise_warning <- FALSE
  rise_notes <- FALSE
  method_coerced <- toupper(method)
  ab_coerced <- as.ab(ab)
  
  if (identical(reference_data, AMR::clinical_breakpoints)) {
    breakpoints <- reference_data %pm>%
      subset(guideline == guideline_coerced & method == method_coerced & ab %in% ab_coerced)
    if (any(ab_coerced == "AMX") && nrow(breakpoints[breakpoints$ab == "AMX", , drop = FALSE]) == 0) {
      ab_coerced[ab_coerced == "AMX"] <- "AMP"
      breakpoints <- reference_data %pm>%
        subset(guideline == guideline_coerced & method == method_coerced & ab %in% ab_coerced)
    }
  } else {
    breakpoints <- reference_data %pm>%
      subset(method == method_coerced & ab %in% ab_coerced)
  }
  
  # create the unique data frame to be filled to save time
  df <- data.frame(
    values = x,
    mo = mo,
    ab = ab,
    result = NA_sir_,
    uti = uti,
    host = host,
    stringsAsFactors = FALSE
  )
  
  if (method == "mic") {
    # when as.sir.mic is called directly
    df$values <- as.mic(df$values)
  } else if (method == "disk") {
    # when as.sir.disk is called directly
    df$values <- as.disk(df$values)
  }
  df_unique <- unique(df[ , c("mo", "ab", "uti", "host"), drop = FALSE])
  
  # get all breakpoints
  breakpoints <- breakpoints %pm>%
    subset(type == breakpoint_type)
  
  if (isFALSE(include_screening)) {
    # remove screening rules from the breakpoints table
    breakpoints <- breakpoints %pm>%
      subset(site %unlike% "screen" & ref_tbl %unlike% "screen")
  }
  if (isFALSE(include_PKPD)) {
    # remove PKPD rules from the breakpoints table
    breakpoints <- breakpoints %pm>%
      subset(mo != "UNKNOWN" & ref_tbl %unlike% "PK.*PD")
  }

  notes <- character(0)

  if (guideline_coerced %like% "EUCAST") {
    any_is_intrinsic_resistant <- FALSE
    add_intrinsic_resistance_to_AMR_env()
  }
  
  if (nrow(df_unique) < 10 || nrow(breakpoints) == 0) {
    # only print intro under 10 items, otherwise progressbar will print this and then it will be printed double
    message_(intro_txt, appendLF = FALSE, as_note = FALSE)
  }
  p <- progress_ticker(n = nrow(df_unique), n_min = 10, title = font_blue(intro_txt), only_bar_percent = TRUE)
  has_progress_bar <- !is.null(import_fn("progress_bar", "progress", error_on_fail = FALSE)) && nrow(df_unique) >= 10
  on.exit(close(p))
  
  if (nrow(breakpoints) == 0) {
    # apparently no breakpoints found
    message(
      paste0(font_rose_bg(" WARNING "), "\n"),
      font_black(paste0("  ", AMR_env$bullet_icon, " No ", guideline_coerced, " ", method_coerced, " breakpoints available for ",
                        suppressMessages(suppressWarnings(ab_name(unique(ab_coerced), language = NULL, tolower = TRUE))),
                        " (", unique(ab_coerced), ")."), collapse = "\n"))
    
    load_mo_uncertainties(metadata_mo)
    return(rep(NA_sir_, nrow(df)))
  }
  
  vectorise_log_entry <- function(x, len) {
    if (length(x) == 1 && len > 1) {
      rep(x, len)
    } else {
      x
    }
  }
  
  # run the rules (df_unique is a row combination per mo/ab/uti/host)
  for (i in seq_len(nrow(df_unique))) {
    p$tick()
    mo_current <- df_unique[i, "mo", drop = TRUE]
    ab_current <- df_unique[i, "ab", drop = TRUE]
    host_current <- df_unique[i, "host", drop = TRUE]
    if (is.na(host_current)) {
      # fall back to human
      host_current <- "human"
    }
    uti_current <- df_unique[i, "uti", drop = TRUE]
    notes_current <- character(0)
    if (isFALSE(uti_current)) {
      # no preference, so no filter on UTIs
      rows <- which(df$mo == mo_current & df$ab == ab_current & df$host == host_current)
    } else {
      rows <- which(df$mo == mo_current & df$ab == ab_current & df$host == host_current & df$uti == uti_current)
    }
    if (length(rows) == 0) {
      notes_current <- c(notes_current, font_red("Returned an empty result, which is unexpected. Are all of `mo`, `ab`, and `host` set and available?"))
      notes <- c(notes, notes_current)
      rise_warning <- TRUE
      next
    }
    values <- df[rows, "values", drop = TRUE]
    new_sir <- rep(NA_sir_, length(rows))

    # find different mo properties, as fast as possible
    mo_current_genus <- AMR_env$MO_lookup$mo[match(AMR_env$MO_lookup$genus[match(mo_current, AMR_env$MO_lookup$mo)], AMR_env$MO_lookup$fullname)]
    mo_current_family <- AMR_env$MO_lookup$mo[match(AMR_env$MO_lookup$family[match(mo_current, AMR_env$MO_lookup$mo)], AMR_env$MO_lookup$fullname)]
    mo_current_order <- AMR_env$MO_lookup$mo[match(AMR_env$MO_lookup$order[match(mo_current, AMR_env$MO_lookup$mo)], AMR_env$MO_lookup$fullname)]
    mo_current_class <- AMR_env$MO_lookup$mo[match(AMR_env$MO_lookup$class[match(mo_current, AMR_env$MO_lookup$mo)], AMR_env$MO_lookup$fullname)]
    mo_current_rank <- AMR_env$MO_lookup$rank[match(mo_current, AMR_env$MO_lookup$mo)]
    mo_current_name <- AMR_env$MO_lookup$fullname[match(mo_current, AMR_env$MO_lookup$mo)]
    if (mo_current %in% AMR::microorganisms.groups$mo) {
      # get the species group (might be more than 1 entry)
      mo_current_species_group <- AMR::microorganisms.groups$mo_group[which(AMR::microorganisms.groups$mo == mo_current)]
    } else {
      mo_current_species_group <- NULL
    }
    mo_current_other <- structure("UNKNOWN", class = c("mo", "character"))
    # formatted for notes
    mo_formatted <- mo_current_name
    if (!mo_current_rank %in% c("kingdom", "phylum", "class", "order")) {
      mo_formatted <- font_italic(mo_formatted, collapse = NULL)
    }
    ab_formatted <- paste0(
      suppressMessages(suppressWarnings(ab_name(ab_current, language = NULL, tolower = TRUE))),
      " (", ab_current, ")"
    )
   
    # gather all available breakpoints for current MO
    breakpoints_current <- breakpoints %pm>%
      subset(ab == ab_current) %pm>%
      subset(mo %in% c(
        mo_current, mo_current_genus, mo_current_family,
        mo_current_order, mo_current_class,
        mo_current_species_group,
        mo_current_other
      ))
    
    if (NROW(breakpoints_current) == 0) {
      AMR_env$sir_interpretation_history <- rbind_AMR(
        AMR_env$sir_interpretation_history,
        # recycling 1 to 2 rows does not always seem to work, which is why vectorise_log_entry() was added
        data.frame(
          datetime = vectorise_log_entry(Sys.time(), length(rows)),
          index = rows,
          ab_given = vectorise_log_entry(ab.bak[match(ab_current, df$ab)][1], length(rows)),
          mo_given = vectorise_log_entry(mo.bak[match(mo_current, df$mo)][1], length(rows)),
          host_given = vectorise_log_entry(host.bak[match(host_current, df$host)][1], length(rows)),
          ab = vectorise_log_entry(ab_current, length(rows)),
          mo = vectorise_log_entry(mo_current, length(rows)),
          host = vectorise_log_entry(host_current, length(rows)),
          method = vectorise_log_entry(method_coerced, length(rows)),
          input = vectorise_log_entry(as.double(values), length(rows)),
          outcome = vectorise_log_entry(NA_sir_, length(rows)),
          notes = vectorise_log_entry("NO BREAKPOINT AVAILABLE", length(rows)),
          guideline = vectorise_log_entry(guideline_coerced, length(rows)),
          ref_table = vectorise_log_entry(NA_character_, length(rows)),
          uti = vectorise_log_entry(uti_current, length(rows)),
          breakpoint_S_R = vectorise_log_entry(NA_character_, length(rows)),
          stringsAsFactors = FALSE
        )
      )
      notes <- c(notes, notes_current)
      next
    }
    
    # set the host index according to most available breakpoints (see R/zzz.R where this is set in the pkg environment)
    breakpoints_current$host_index <- match(breakpoints_current$host, c("human", "ECOFF", AMR_env$host_preferred_order))
    
    # sort on host and taxonomic rank
    # (this will e.g. prefer 'species' breakpoints over 'order' breakpoints)
    if (all(uti_current == FALSE, na.rm = TRUE)) {
      breakpoints_current <- breakpoints_current %pm>%
        #  `uti` is a column in the data set
        # this will put UTI = FALSE first, then UTI = NA, then UTI = TRUE
        pm_mutate(uti_index = ifelse(is.na(uti) & uti == FALSE, 1, 
                                     ifelse(is.na(uti), 2,
                                            3))) %pm>%
        # be as specific as possible (i.e. prefer species over genus):
        pm_arrange(host_index, rank_index, uti_index)
    } else if (all(uti_current == TRUE, na.rm = TRUE)) {
      breakpoints_current <- breakpoints_current %pm>%
        subset(uti == TRUE) %pm>%
        # be as specific as possible (i.e. prefer species over genus):
        pm_arrange(host_index, rank_index)
    }
    
    # veterinary host check
    host_current <- unique(df_unique[i, "host", drop = TRUE])[1]
    breakpoints_current$host_match <- breakpoints_current$host == host_current
    if (breakpoint_type == "animal") {
      if (any(breakpoints_current$host_match == TRUE, na.rm = TRUE)) {
        breakpoints_current <- breakpoints_current %pm>%
          subset(host_match == TRUE)
      } else {
        # no breakpoint found for this host, so sort on mostly available guidelines
        notes_current <- c(notes_current, paste0("Using ", font_bold(breakpoints_current$host[1]), " breakpoints since ", font_bold(host_current), " for ", ab_formatted, " in ", mo_formatted, " are not available."))
      }
    }
 
    # throw messages for different body sites
    site <- breakpoints_current[1L, "site", drop = FALSE] # this is the one we'll take
    if (is.na(site)) {
      site <- paste0("an unspecified body site")
    } else {
      site <- paste0("body site '", site, "'")
    }
    if (nrow(breakpoints_current) == 1 && all(breakpoints_current$uti == TRUE) && any(uti_current %in% c(FALSE, NA)) && message_not_thrown_before("as.sir", "uti", ab_current)) {
      # only UTI breakpoints available
      notes_current <- c(notes_current, paste0("Breakpoints for ", font_bold(ab_formatted), " in ", mo_formatted, " are only available for (uncomplicated) urinary tract infections (UTI); assuming `uti = TRUE`."))
    } else if (nrow(breakpoints_current) > 1 && length(unique(breakpoints_current$site)) > 1 && any(is.na(uti_current)) && all(c(TRUE, FALSE) %in% breakpoints_current$uti, na.rm = TRUE) && message_not_thrown_before("as.sir", "siteUTI", mo_current, ab_current)) {
      # both UTI and Non-UTI breakpoints available
      notes_current <- c(notes_current, paste0("Breakpoints for UTI ", font_bold("and"), " non-UTI available for ", ab_formatted, " in ", mo_formatted, " - assuming ", site, ". Use argument `uti` to set which isolates are from urine. See `?as.sir`."))
      breakpoints_current <- breakpoints_current %pm>%
        pm_filter(uti == FALSE)
    } else if (nrow(breakpoints_current) > 1 && length(unique(breakpoints_current$site)) > 1 && all(breakpoints_current$uti == FALSE, na.rm = TRUE) && message_not_thrown_before("as.sir", "siteOther", mo_current, ab_current)) {
      # breakpoints for multiple body sites available
      notes_current <- c(notes_current, paste0("Multiple breakpoints available for ", font_bold(ab_formatted), " in ", mo_formatted, " - assuming ", site, "."))
    }
    
    # first check if mo is intrinsic resistant
    if (isTRUE(add_intrinsic_resistance) && guideline_coerced %like% "EUCAST" && paste(mo_current, ab_current) %in% AMR_env$intrinsic_resistant) {
      notes_current <- c(notes_current, paste0("Intrinsic resistance applied for ", ab_formatted, " in ", mo_formatted, ""))
      new_sir <- rep(as.sir("R"), length(rows))
    } else if (nrow(breakpoints_current) == 0) {
      # no rules available
      new_sir <- rep(NA_sir_, length(rows))
    } else {
      # then run the rules
      breakpoints_current <- breakpoints_current[1L, , drop = FALSE]

      if (any(breakpoints_current$mo == "UNKNOWN", na.rm = TRUE) | any(breakpoints_current$ref_tbl %like% "PK.*PD", na.rm = TRUE)) {
        notes_current <- c(notes_current, "Some PK/PD breakpoints were applied - use `include_PKPD = FALSE` to prevent this")
      }
      if (any(breakpoints_current$site %like% "screen", na.rm = TRUE) | any(breakpoints_current$ref_tbl %like% "screen", na.rm = TRUE)) {
        notes_current <- c(notes_current, "Some screening breakpoints were applied - use `include_screening = FALSE` to prevent this")
      }

      if (method == "mic") {
        new_sir <- case_when_AMR(
          is.na(values) ~ NA_sir_,
          values <= breakpoints_current$breakpoint_S ~ as.sir("S"),
          guideline_coerced %like% "EUCAST" & values > breakpoints_current$breakpoint_R ~ as.sir("R"),
          guideline_coerced %like% "CLSI" & values >= breakpoints_current$breakpoint_R ~ as.sir("R"),
          # return "I" when breakpoints are in the middle
          !is.na(breakpoints_current$breakpoint_S) & !is.na(breakpoints_current$breakpoint_R) ~ as.sir("I"),
          # and NA otherwise
          TRUE ~ NA_sir_
        )
      } else if (method == "disk") {
        new_sir <- case_when_AMR(
          is.na(values) ~ NA_sir_,
          as.double(values) >= as.double(breakpoints_current$breakpoint_S) ~ as.sir("S"),
          guideline_coerced %like% "EUCAST" & as.double(values) < as.double(breakpoints_current$breakpoint_R) ~ as.sir("R"),
          guideline_coerced %like% "CLSI" & as.double(values) <= as.double(breakpoints_current$breakpoint_R) ~ as.sir("R"),
          # return "I" when breakpoints are in the middle
          !is.na(breakpoints_current$breakpoint_S) & !is.na(breakpoints_current$breakpoint_R) ~ as.sir("I"),
          # and NA otherwise
          TRUE ~ NA_sir_
        )
      }
      
      # write to verbose output
      AMR_env$sir_interpretation_history <- rbind_AMR(
        AMR_env$sir_interpretation_history,
        # recycling 1 to 2 rows does not always seem to work, which is why vectorise_log_entry() was added
        data.frame(
          datetime = vectorise_log_entry(Sys.time(), length(rows)),
          index = rows,
          ab_given = vectorise_log_entry(ab.bak[match(ab_current, df$ab)][1], length(rows)),
          mo_given = vectorise_log_entry(mo.bak[match(mo_current, df$mo)][1], length(rows)),
          host_given = vectorise_log_entry(host.bak[match(host_current, df$host)][1], length(rows)),
          ab = vectorise_log_entry(breakpoints_current[, "ab", drop = TRUE], length(rows)),
          mo = vectorise_log_entry(breakpoints_current[, "mo", drop = TRUE], length(rows)),
          host = vectorise_log_entry(breakpoints_current[, "host", drop = TRUE], length(rows)),
          method = vectorise_log_entry(method_coerced, length(rows)),
          input = vectorise_log_entry(as.double(values), length(rows)),
          outcome = vectorise_log_entry(as.sir(new_sir), length(rows)),
          notes = vectorise_log_entry(paste0(font_stripstyle(notes_current), collapse = "\n"), length(rows)),
          guideline = vectorise_log_entry(guideline_coerced, length(rows)),
          ref_table = vectorise_log_entry(breakpoints_current[, "ref_tbl", drop = TRUE], length(rows)),
          uti = vectorise_log_entry(breakpoints_current[, "uti", drop = TRUE], length(rows)),
          breakpoint_S_R = vectorise_log_entry(paste0(breakpoints_current[, "breakpoint_S", drop = TRUE], "-", breakpoints_current[, "breakpoint_R", drop = TRUE]), length(rows)),
          stringsAsFactors = FALSE
        )
      )
    }

    notes <- c(notes, notes_current)
    df[rows, "result"] <- new_sir
  }
  
  close(p)
  # printing messages
  if (has_progress_bar == TRUE) {
    # the progress bar has overwritten the intro text, so:
    message_(intro_txt, appendLF = FALSE, as_note = FALSE)
  }
  if (length(notes) > 0) {
    if (isTRUE(rise_warning)) {
      message(font_rose_bg(" WARNING "))
    } else {
      message(font_yellow_bg(" NOTE "))
    }
    notes <- unique(notes)
    if (isTRUE(verbose) || length(notes) == 1 || NROW(AMR_env$sir_interpretation_history) == 0) {
      for (i in seq_len(length(notes))) {
        message(word_wrap("  ", AMR_env$bullet_icon, " ", notes[i], add_fn = font_black))
      }
    } else {
      message(word_wrap("  ", AMR_env$bullet_icon, " There were multiple notes. Print or View `sir_interpretation_history()` to examine them, or use `as.sir(..., verbose = TRUE)` next time to directly print them here.", add_fn = font_black))
    }
  } else {
    message(font_green_bg(" OK "))
  }
  
  load_mo_uncertainties(metadata_mo)
  
  df$result
}

#' @rdname as.sir
#' @param clean a [logical] to indicate whether previously stored results should be forgotten after returning the 'logbook' with results
#' @export
sir_interpretation_history <- function(clean = FALSE) {
  meet_criteria(clean, allow_class = "logical", has_length = 1)

  out <- AMR_env$sir_interpretation_history
  out$outcome <- as.sir(out$outcome)
  if (NROW(out) > 0) {
    # sort descending on time
    out <- out[order(format(out$datetime, "%Y%m%d%H%M"), out$index, decreasing = TRUE), , drop = FALSE]
  }
  
  if (isTRUE(clean)) {
    AMR_env$sir_interpretation_history <- AMR_env$sir_interpretation_history[0, , drop = FALSE]
  }
  
  if (pkg_is_available("tibble")) {
    out <- import_fn("as_tibble", "tibble")(out)
  }
  structure(out, class = c("sir_log", class(out)))
}

#' @method print sir_log
#' @export
#' @noRd
print.sir_log <- function(x, ...) {
  if (NROW(x) == 0) {
    message_("No results to print. Run `as.sir()` on MIC values or disk diffusion zones first to print a 'logbook' data set here.")
    return(invisible(NULL))
  }
  class(x) <- class(x)[class(x) != "sir_log"]
  print(x, ...)
}

# will be exported using s3_register() in R/zzz.R
pillar_shaft.sir <- function(x, ...) {
  out <- trimws(format(x))
  if (has_colour()) {
    # colours will anyway not work when has_colour() == FALSE,
    # but then the indentation should also not be applied
    out[is.na(x)] <- font_grey("  NA")
    out[x == "NI"] <- font_grey_bg("  NI ")
    out[x == "S"] <- font_green_bg("  S  ")
    out[x == "I"] <- font_orange_bg("  I  ")
    out[x == "SDD"] <- font_orange_bg(" SDD ")
    if (is_dark()) {
      out[x == "R"] <- font_red_bg("  R  ")
    } else {
      out[x == "R"] <- font_rose_bg("  R  ")
    }
  }
  create_pillar_column(out, align = "left", width = 5)
}

# will be exported using s3_register() in R/zzz.R
type_sum.sir <- function(x, ...) {
  "sir"
}

# will be exported using s3_register() in R/zzz.R
freq.sir <- function(x, ...) {
  x_name <- deparse(substitute(x))
  x_name <- gsub(".*[$]", "", x_name)
  if (x_name %in% c("x", ".")) {
    # try again going through system calls
    x_name <- stats::na.omit(vapply(
      FUN.VALUE = character(1),
      sys.calls(),
      function(call) {
        call_txt <- as.character(call)
        ifelse(call_txt[1] %like% "freq$", call_txt[length(call_txt)], character(0))
      }
    ))[1L]
  }
  ab <- suppressMessages(suppressWarnings(as.ab(x_name)))
  digits <- list(...)$digits
  if (is.null(digits)) {
    digits <- 2
  }
  if (!is.na(ab)) {
    cleaner::freq.default(
      x = x, ...,
      .add_header = list(
        Drug = paste0(ab_name(ab, language = NULL), " (", ab, ", ", paste(ab_atc(ab), collapse = "/"), ")"),
        `Drug group` = ab_group(ab, language = NULL),
        `%SI` = trimws(percentage(susceptibility(x, minimum = 0, as_percent = FALSE),
          digits = digits
        ))
      )
    )
  } else {
    cleaner::freq.default(
      x = x, ...,
      .add_header = list(
        `%SI` = trimws(percentage(susceptibility(x, minimum = 0, as_percent = FALSE),
          digits = digits
        ))
      )
    )
  }
}


# will be exported using s3_register() in R/zzz.R
get_skimmers.sir <- function(column) {
  # get the variable name 'skim_variable'
  name_call <- function(.data) {
    calls <- sys.calls()
    frms <- sys.frames()
    calls_txt <- vapply(calls, function(x) paste(deparse(x), collapse = ""), FUN.VALUE = character(1))
    if (any(calls_txt %like% "skim_variable", na.rm = TRUE)) {
      ind <- which(calls_txt %like% "skim_variable")[1L]
      vars <- tryCatch(eval(parse(text = ".data$skim_variable$sir"), envir = frms[[ind]]),
        error = function(e) NULL
      )
      tryCatch(ab_name(as.character(calls[[length(calls)]][[2]]), language = NULL),
        error = function(e) NA_character_
      )
    } else {
      NA_character_
    }
  }

  skimr::sfl(
    skim_type = "sir",
    ab_name = name_call,
    count_R = count_R,
    count_S = count_susceptible,
    count_I = count_I,
    prop_R = ~ proportion_R(., minimum = 0),
    prop_S = ~ susceptibility(., minimum = 0),
    prop_I = ~ proportion_I(., minimum = 0)
  )
}

#' @method print sir
#' @export
#' @noRd
print.sir <- function(x, ...) {
  cat("Class 'sir'\n")
  print(as.character(x), quote = FALSE)
}


#' @method as.double sir
#' @export
as.double.sir <- function(x, ...) {
  dbls <- rep(NA_real_, length(x))
  dbls[x == "S"] <- 1
  dbls[x %in% c("SDD", "I")] <- 2
  dbls[x == "R"] <- 3
  dbls
}

#' @method droplevels sir
#' @export
#' @noRd
droplevels.sir <- function(x, exclude = if (any(is.na(levels(x)))) NULL else NA, ...) {
  x <- droplevels.factor(x, exclude = exclude, ...)
  class(x) <- c("sir", "ordered", "factor")
  x
}

#' @method summary sir
#' @export
#' @noRd
summary.sir <- function(object, ...) {
  x <- object
  n <- sum(!is.na(x))
  S <- sum(x == "S", na.rm = TRUE)
  SDD <- sum(x == "SDD", na.rm = TRUE)
  I <- sum(x == "I", na.rm = TRUE)
  R <- sum(x == "R", na.rm = TRUE)
  NI <- sum(x == "NI", na.rm = TRUE)
  pad <- function(x) {
    if (is.na(x)) {
      return("??")
    }
    if (x == "0%") {
      x <- " 0.0%"
    }
    if (nchar(x) < 5) {
      x <- paste0(rep(" ", 5 - nchar(x)), x)
    }
    x
  }
  value <- c(
    "Class" = "sir",
    "%S" = paste0(pad(percentage(S / n, digits = 1)), " (n=", S, ")"),
    "%SDD" = paste0(pad(percentage(SDD / n, digits = 1)), " (n=", SDD, ")"),
    "%I" = paste0(pad(percentage(I / n, digits = 1)), " (n=", I, ")"),
    "%R" = paste0(pad(percentage(R / n, digits = 1)), " (n=", R, ")"),
    "%NI" = paste0(pad(percentage(NI / n, digits = 1)), " (n=", NI, ")")
  )
  class(value) <- c("summaryDefault", "table")
  value
}

#' @method [<- sir
#' @export
#' @noRd
"[<-.sir" <- function(i, j, ..., value) {
  value <- as.sir(value)
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @method [[<- sir
#' @export
#' @noRd
"[[<-.sir" <- function(i, j, ..., value) {
  value <- as.sir(value)
  y <- NextMethod()
  attributes(y) <- attributes(i)
  y
}
#' @method c sir
#' @export
#' @noRd
c.sir <- function(...) {
  as.sir(unlist(lapply(list(...), as.character)))
}

#' @method unique sir
#' @export
#' @noRd
unique.sir <- function(x, incomparables = FALSE, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

#' @method rep sir
#' @export
#' @noRd
rep.sir <- function(x, ...) {
  y <- NextMethod()
  attributes(y) <- attributes(x)
  y
}

check_reference_data <- function(reference_data, .call_depth) {
  if (!identical(reference_data, AMR::clinical_breakpoints)) {
    class_sir <- vapply(FUN.VALUE = character(1), AMR::clinical_breakpoints, function(x) paste0("<", class(x), ">", collapse = " and "))
    class_ref <- vapply(FUN.VALUE = character(1), reference_data, function(x) paste0("<", class(x), ">", collapse = " and "))
    if (!all(names(class_sir) == names(class_ref))) {
      stop_("`reference_data` must have the same column names as the 'clinical_breakpoints' data set.", call = .call_depth)
    }
    if (!all(class_sir == class_ref)) {
      stop_("`reference_data` must be the same structure as the 'clinical_breakpoints' data set. Column '", names(class_ref[class_sir != class_ref][1]), "' is of class ", class_ref[class_sir != class_ref][1], ", but should be of class ", class_sir[class_sir != class_ref][1], ".", call = .call_depth)
    }
  }
}
