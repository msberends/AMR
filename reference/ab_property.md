# Get Properties of an Antibiotic

Use these functions to return a specific property of an antibiotic from
the [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
data set. All input values will be evaluated internally with
[`as.ab()`](https://amr-for-r.org/reference/as.ab.md).

## Usage

``` r
ab_name(x, language = get_AMR_locale(), tolower = FALSE, ...)

ab_cid(x, ...)

ab_synonyms(x, ...)

ab_tradenames(x, ...)

ab_group(x, language = get_AMR_locale(), all_groups = FALSE, ...)

ab_atc(x, only_first = FALSE, ...)

ab_atc_group1(x, language = get_AMR_locale(), ...)

ab_atc_group2(x, language = get_AMR_locale(), ...)

ab_loinc(x, ...)

ab_ddd(x, administration = "oral", ...)

ab_ddd_units(x, administration = "oral", ...)

ab_info(x, language = get_AMR_locale(), ...)

ab_url(x, open = FALSE, ...)

ab_property(x, property = "name", language = get_AMR_locale(), ...)

set_ab_names(data, ..., property = "name", language = get_AMR_locale(),
  snake_case = NULL)
```

## Arguments

- x:

  Any (vector of) text that can be coerced to a valid antibiotic drug
  code with [`as.ab()`](https://amr-for-r.org/reference/as.ab.md).

- language:

  Language of the returned text - the default is the current system
  language (see
  [`get_AMR_locale()`](https://amr-for-r.org/reference/translate.md))
  and can also be set with the package option
  [`AMR_locale`](https://amr-for-r.org/reference/AMR-options.md). Use
  `language = NULL` or `language = ""` to prevent translation.

- tolower:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the first [character](https://rdrr.io/r/base/character.html) of every
  output should be transformed to a lower case
  [character](https://rdrr.io/r/base/character.html). This will lead to
  e.g. "polymyxin B" and not "polymyxin b".

- ...:

  In case of `set_ab_names()` and `data` is a
  [data.frame](https://rdrr.io/r/base/data.frame.html): columns to
  select (supports tidy selection such as `column1:column4`), otherwise
  other arguments passed on to
  [`as.ab()`](https://amr-for-r.org/reference/as.ab.md).

- only_first:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only the first ATC code must be returned, with giving preference to
  J0-codes (i.e., the antimicrobial drug group).

- administration:

  Way of administration, either `"oral"` or `"iv"`.

- open:

  Browse the URL using
  [`utils::browseURL()`](https://rdrr.io/r/utils/browseURL.html).

- property:

  One of the column names of one of the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set: `vector_or(colnames(antimicrobials), sort = FALSE)`.

- data:

  A [data.frame](https://rdrr.io/r/base/data.frame.html) of which the
  columns need to be renamed, or a
  [character](https://rdrr.io/r/base/character.html) vector of column
  names.

- snake_case:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the names should be in so-called [snake
  case](https://en.wikipedia.org/wiki/Snake_case): in lower case and all
  spaces/slashes replaced with an underscore (`_`).

## Value

- An [integer](https://rdrr.io/r/base/integer.html) in case of
  `ab_cid()`

- A named [list](https://rdrr.io/r/base/list.html) in case of
  `ab_info()` and multiple `ab_atc()`/`ab_synonyms()`/`ab_tradenames()`

- A [double](https://rdrr.io/r/base/double.html) in case of `ab_ddd()`

- A [data.frame](https://rdrr.io/r/base/data.frame.html) in case of
  `set_ab_names()`

- A [character](https://rdrr.io/r/base/character.html) in all other
  cases

## Details

All output [will be
translated](https://amr-for-r.org/reference/translate.md) where
possible.

The function `ab_url()` will return the direct URL to the official WHO
website. A warning will be returned if the required ATC code is not
available.

The function `set_ab_names()` is a special column renaming function for
[data.frame](https://rdrr.io/r/base/data.frame.html)s. It renames
columns names that resemble antimicrobial drugs. It always makes sure
that the new column names are unique. If `property = "atc"` is set,
preference is given to ATC codes from the J-group.

## Source

World Health Organization (WHO) Collaborating Centre for Drug Statistics
Methodology: <https://atcddd.fhi.no/atc_ddd_index/>

European Commission Public Health PHARMACEUTICALS - COMMUNITY REGISTER:
<https://ec.europa.eu/health/documents/community-register/html/reg_hum_atc.htm>

## Download Our Reference Data

All reference data sets in the AMR package - including information on
microorganisms, antimicrobials, and clinical breakpoints - are freely
available for download in multiple formats: R, MS Excel, Apache Feather,
Apache Parquet, SPSS, and Stata.

For maximum compatibility, we also provide machine-readable,
tab-separated plain text files suitable for use in any software,
including laboratory information systems.

Visit [our website for direct download
links](https://amr-for-r.org/articles/datasets.html), or explore the
actual files in [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw/datasets).

## See also

[antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)

## Examples

``` r
# all properties:
ab_name("AMX")
#> [1] "Amoxicillin"
ab_atc("AMX")
#> [1] "J01CA04"  "QG51AA03" "QJ01CA04"
ab_cid("AMX")
#> [1] 33613
ab_synonyms("AMX")
#>  [1] "acuotricina"     "alfamox"         "alfida"          "amitron"        
#>  [5] "amoclen"         "amodex"          "amoksicillin"    "amolin"         
#>  [9] "amopen"          "amopenixin"      "amophar"         "amoran"         
#> [13] "amoxi"           "amoxicaps"       "amoxicilina"     "amoxicilline"   
#> [17] "amoxicillinum"   "amoxidal"        "amoxiden"        "amoxil"         
#> [21] "amoxillat"       "amoxina"         "amoxine"         "amoxipen"       
#> [25] "amoxivet"        "amoxycillin"     "amoxycillinsalt" "amoxyke"        
#> [29] "anemolin"        "aspenil"         "atoksilin"       "bristamox"      
#> [33] "cemoxin"         "ciblor"          "clamoxyl"        "damoxy"         
#> [37] "danoxillin"      "delacillin"      "demoksil"        "dispermox"      
#> [41] "efpenix"         "eupen"           "flemoxin"        "flemoxine"      
#> [45] "galenamox"       "gramidil"        "hiconcil"        "himinomax"      
#> [49] "histocillin"     "ibiamox"         "imacillin"       "izoltil"        
#> [53] "kentrocyllin"    "lamoxy"          "largopen"        "larotid"        
#> [57] "matasedrin"      "metifarma"       "moksilin"        "moxacin"        
#> [61] "moxal"           "moxaline"        "moxatag"         "neotetranase"   
#> [65] "novabritine"     "ospamox"         "pacetocin"       "pamocil"        
#> [69] "paradroxil"      "pasetocin"       "penamox"         "piramox"        
#> [73] "promoxil"        "quimiopen"       "remoxil"         "riotapen"       
#> [77] "robamox"         "sawacillin"      "siganopen"       "simplamox"      
#> [81] "sintopen"        "sumox"           "topramoxin"      "trifamox"       
#> [85] "trimox"          "unicillin"       "utimox"          "velamox"        
#> [89] "vetramox"        "wymox"           "zamocillin"      "zamocilline"    
#> [93] "zimox"          
ab_tradenames("AMX")
#>  [1] "acuotricina"     "alfamox"         "alfida"          "amitron"        
#>  [5] "amoclen"         "amodex"          "amoksicillin"    "amolin"         
#>  [9] "amopen"          "amopenixin"      "amophar"         "amoran"         
#> [13] "amoxi"           "amoxicaps"       "amoxicilina"     "amoxicilline"   
#> [17] "amoxicillinum"   "amoxidal"        "amoxiden"        "amoxil"         
#> [21] "amoxillat"       "amoxina"         "amoxine"         "amoxipen"       
#> [25] "amoxivet"        "amoxycillin"     "amoxycillinsalt" "amoxyke"        
#> [29] "anemolin"        "aspenil"         "atoksilin"       "bristamox"      
#> [33] "cemoxin"         "ciblor"          "clamoxyl"        "damoxy"         
#> [37] "danoxillin"      "delacillin"      "demoksil"        "dispermox"      
#> [41] "efpenix"         "eupen"           "flemoxin"        "flemoxine"      
#> [45] "galenamox"       "gramidil"        "hiconcil"        "himinomax"      
#> [49] "histocillin"     "ibiamox"         "imacillin"       "izoltil"        
#> [53] "kentrocyllin"    "lamoxy"          "largopen"        "larotid"        
#> [57] "matasedrin"      "metifarma"       "moksilin"        "moxacin"        
#> [61] "moxal"           "moxaline"        "moxatag"         "neotetranase"   
#> [65] "novabritine"     "ospamox"         "pacetocin"       "pamocil"        
#> [69] "paradroxil"      "pasetocin"       "penamox"         "piramox"        
#> [73] "promoxil"        "quimiopen"       "remoxil"         "riotapen"       
#> [77] "robamox"         "sawacillin"      "siganopen"       "simplamox"      
#> [81] "sintopen"        "sumox"           "topramoxin"      "trifamox"       
#> [85] "trimox"          "unicillin"       "utimox"          "velamox"        
#> [89] "vetramox"        "wymox"           "zamocillin"      "zamocilline"    
#> [93] "zimox"          
ab_group("AMX")
#> [1] "Aminopenicillins"
ab_atc_group1("AMX")
#> [1] "Beta-lactam antibacterials, penicillins"
ab_atc_group2("AMX")
#> [1] "Penicillins with extended spectrum"
ab_url("AMX")
#>                                                             Amoxicillin 
#> "https://atcddd.fhi.no/atc_ddd_index//?code=J01CA04&showdescription=no" 

# smart lowercase transformation
ab_name(x = c("AMC", "PLB"))
#> [1] "Amoxicillin/clavulanic acid" "Polymyxin B"                
ab_name(x = c("AMC", "PLB"), tolower = TRUE)
#> [1] "amoxicillin/clavulanic acid" "polymyxin B"                

# defined daily doses (DDD)
ab_ddd("AMX", "oral")
#> [1] 1.5
ab_ddd_units("AMX", "oral")
#> [1] "g"
ab_ddd("AMX", "iv")
#> [1] 3
ab_ddd_units("AMX", "iv")
#> [1] "g"

ab_info("AMX") # all properties as a list
#> $ab
#> [1] "AMX"
#> 
#> $cid
#> [1] 33613
#> 
#> $name
#> [1] "Amoxicillin"
#> 
#> $group
#> [1] "Aminopenicillins"
#> 
#> $atc
#> [1] "J01CA04"  "QG51AA03" "QJ01CA04"
#> 
#> $atc_group1
#> [1] "Beta-lactam antibacterials, penicillins"
#> 
#> $atc_group2
#> [1] "Penicillins with extended spectrum"
#> 
#> $tradenames
#>  [1] "acuotricina"     "alfamox"         "alfida"          "amitron"        
#>  [5] "amoclen"         "amodex"          "amoksicillin"    "amolin"         
#>  [9] "amopen"          "amopenixin"      "amophar"         "amoran"         
#> [13] "amoxi"           "amoxicaps"       "amoxicilina"     "amoxicilline"   
#> [17] "amoxicillinum"   "amoxidal"        "amoxiden"        "amoxil"         
#> [21] "amoxillat"       "amoxina"         "amoxine"         "amoxipen"       
#> [25] "amoxivet"        "amoxycillin"     "amoxycillinsalt" "amoxyke"        
#> [29] "anemolin"        "aspenil"         "atoksilin"       "bristamox"      
#> [33] "cemoxin"         "ciblor"          "clamoxyl"        "damoxy"         
#> [37] "danoxillin"      "delacillin"      "demoksil"        "dispermox"      
#> [41] "efpenix"         "eupen"           "flemoxin"        "flemoxine"      
#> [45] "galenamox"       "gramidil"        "hiconcil"        "himinomax"      
#> [49] "histocillin"     "ibiamox"         "imacillin"       "izoltil"        
#> [53] "kentrocyllin"    "lamoxy"          "largopen"        "larotid"        
#> [57] "matasedrin"      "metifarma"       "moksilin"        "moxacin"        
#> [61] "moxal"           "moxaline"        "moxatag"         "neotetranase"   
#> [65] "novabritine"     "ospamox"         "pacetocin"       "pamocil"        
#> [69] "paradroxil"      "pasetocin"       "penamox"         "piramox"        
#> [73] "promoxil"        "quimiopen"       "remoxil"         "riotapen"       
#> [77] "robamox"         "sawacillin"      "siganopen"       "simplamox"      
#> [81] "sintopen"        "sumox"           "topramoxin"      "trifamox"       
#> [85] "trimox"          "unicillin"       "utimox"          "velamox"        
#> [89] "vetramox"        "wymox"           "zamocillin"      "zamocilline"    
#> [93] "zimox"          
#> 
#> $loinc
#>  [1] "101498-4" "15-8"     "16-6"     "16365-9"  "17-4"     "18-2"    
#>  [7] "18861-5"  "18862-3"  "19-0"     "20-8"     "21-6"     "22-4"    
#> [13] "25274-2"  "25310-4"  "3344-9"   "55614-2"  "55615-9"  "55616-7" 
#> [19] "6976-5"   "6977-3"   "80133-2" 
#> 
#> $ddd
#> $ddd$oral
#> $ddd$oral$amount
#> [1] 1.5
#> 
#> $ddd$oral$units
#> [1] "g"
#> 
#> 
#> $ddd$iv
#> $ddd$iv$amount
#> [1] 3
#> 
#> $ddd$iv$units
#> [1] "g"
#> 
#> 
#> 

# all ab_* functions use as.ab() internally, so you can go from 'any' to 'any':
ab_atc("AMP")
#> [1] "J01CA01"  "QJ01CA01" "QJ51CA01" "QS01AA19" "S01AA19" 
ab_group("J01CA01")
#> [1] "Aminopenicillins"
ab_loinc("ampicillin")
#>  [1] "101477-8" "101478-6" "18864-9"  "18865-6"  "20374-5"  "21066-6" 
#>  [7] "23618-2"  "27-3"     "28-1"     "29-9"     "30-7"     "31-5"    
#> [13] "32-3"     "33-1"     "3355-5"   "33562-0"  "33919-2"  "34-9"    
#> [19] "43883-8"  "43884-6"  "6979-9"   "6980-7"   "87604-5" 
ab_name("21066-6")
#> [1] "Ampicillin"
ab_name(6249)
#> [1] "Ampicillin"
ab_name("J01CA01")
#> [1] "Ampicillin"

# spelling from different languages and dyslexia are no problem
ab_atc("ceftriaxon")
#> [1] "J01DD04"  "QJ01DD04"
ab_atc("cephtriaxone")
#> [1] "J01DD04"  "QJ01DD04"
ab_atc("cephthriaxone")
#> [1] "J01DD04"  "QJ01DD04"
ab_atc("seephthriaaksone")
#> [1] "J01DD04"  "QJ01DD04"

# use set_ab_names() for renaming columns
colnames(example_isolates)
#>  [1] "date"    "patient" "age"     "gender"  "ward"    "mo"      "PEN"    
#>  [8] "OXA"     "FLC"     "AMX"     "AMC"     "AMP"     "TZP"     "CZO"    
#> [15] "FEP"     "CXM"     "FOX"     "CTX"     "CAZ"     "CRO"     "GEN"    
#> [22] "TOB"     "AMK"     "KAN"     "TMP"     "SXT"     "NIT"     "FOS"    
#> [29] "LNZ"     "CIP"     "MFX"     "VAN"     "TEC"     "TCY"     "TGC"    
#> [36] "DOX"     "ERY"     "CLI"     "AZM"     "IPM"     "MEM"     "MTR"    
#> [43] "CHL"     "COL"     "MUP"     "RIF"    
colnames(set_ab_names(example_isolates))
#>  [1] "date"                          "patient"                      
#>  [3] "age"                           "gender"                       
#>  [5] "ward"                          "mo"                           
#>  [7] "benzylpenicillin"              "oxacillin"                    
#>  [9] "flucloxacillin"                "amoxicillin"                  
#> [11] "amoxicillin_clavulanic_acid"   "ampicillin"                   
#> [13] "piperacillin_tazobactam"       "cefazolin"                    
#> [15] "cefepime"                      "cefuroxime"                   
#> [17] "cefoxitin"                     "cefotaxime"                   
#> [19] "ceftazidime"                   "ceftriaxone"                  
#> [21] "gentamicin"                    "tobramycin"                   
#> [23] "amikacin"                      "kanamycin"                    
#> [25] "trimethoprim"                  "trimethoprim_sulfamethoxazole"
#> [27] "nitrofurantoin"                "fosfomycin"                   
#> [29] "linezolid"                     "ciprofloxacin"                
#> [31] "moxifloxacin"                  "vancomycin"                   
#> [33] "teicoplanin"                   "tetracycline"                 
#> [35] "tigecycline"                   "doxycycline"                  
#> [37] "erythromycin"                  "clindamycin"                  
#> [39] "azithromycin"                  "imipenem"                     
#> [41] "meropenem"                     "metronidazole"                
#> [43] "chloramphenicol"               "colistin"                     
#> [45] "mupirocin"                     "rifampicin"                   
colnames(set_ab_names(example_isolates, NIT:VAN))
#>  [1] "date"           "patient"        "age"            "gender"        
#>  [5] "ward"           "mo"             "PEN"            "OXA"           
#>  [9] "FLC"            "AMX"            "AMC"            "AMP"           
#> [13] "TZP"            "CZO"            "FEP"            "CXM"           
#> [17] "FOX"            "CTX"            "CAZ"            "CRO"           
#> [21] "GEN"            "TOB"            "AMK"            "KAN"           
#> [25] "TMP"            "SXT"            "nitrofurantoin" "fosfomycin"    
#> [29] "linezolid"      "ciprofloxacin"  "moxifloxacin"   "vancomycin"    
#> [33] "TEC"            "TCY"            "TGC"            "DOX"           
#> [37] "ERY"            "CLI"            "AZM"            "IPM"           
#> [41] "MEM"            "MTR"            "CHL"            "COL"           
#> [45] "MUP"            "RIF"           
# \donttest{
if (require("dplyr")) {
  example_isolates %>%
    set_ab_names()

  # this does the same:
  example_isolates %>%
    rename_with(set_ab_names)

  # set_ab_names() works with any AB property:
  example_isolates %>%
    set_ab_names(property = "atc")

  example_isolates %>%
    set_ab_names(where(is.sir)) %>%
    colnames()

  example_isolates %>%
    set_ab_names(NIT:VAN) %>%
    colnames()
}
#>  [1] "date"           "patient"        "age"            "gender"        
#>  [5] "ward"           "mo"             "PEN"            "OXA"           
#>  [9] "FLC"            "AMX"            "AMC"            "AMP"           
#> [13] "TZP"            "CZO"            "FEP"            "CXM"           
#> [17] "FOX"            "CTX"            "CAZ"            "CRO"           
#> [21] "GEN"            "TOB"            "AMK"            "KAN"           
#> [25] "TMP"            "SXT"            "nitrofurantoin" "fosfomycin"    
#> [29] "linezolid"      "ciprofloxacin"  "moxifloxacin"   "vancomycin"    
#> [33] "TEC"            "TCY"            "TGC"            "DOX"           
#> [37] "ERY"            "CLI"            "AZM"            "IPM"           
#> [41] "MEM"            "MTR"            "CHL"            "COL"           
#> [45] "MUP"            "RIF"           
# }
```
