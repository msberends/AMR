# Determine Multidrug-Resistant Organisms (MDRO)

Determine which isolates are multidrug-resistant organisms (MDRO)
according to international, national and custom guidelines.

## Usage

``` r
mdro(
  x = NULL,
  guideline = "CMI2012",
  col_mo = NULL,
  info = interactive(),
  pct_required_classes = 0.5,
  combine_SI = TRUE,
  verbose = FALSE,
  only_sir_columns = FALSE,
  ...
)

custom_mdro_guideline(..., as_factor = TRUE)

brmo(x = NULL, only_sir_columns = FALSE, ...)

mrgn(x = NULL, only_sir_columns = FALSE, ...)

mdr_tb(x = NULL, only_sir_columns = FALSE, ...)

mdr_cmi2012(x = NULL, only_sir_columns = FALSE, ...)

eucast_exceptional_phenotypes(x = NULL, only_sir_columns = FALSE, ...)
```

## Source

See the supported guidelines above for the
[list](https://rdrr.io/r/base/list.html) of publications used for this
function.

## Arguments

- x:

  a [data.frame](https://rdrr.io/r/base/data.frame.html) with
  antibiotics columns, like `AMX` or `amox`. Can be left blank for
  automatic determination.

- guideline:

  a specific guideline to follow, see sections *Supported international
  / national guidelines* and *Using Custom Guidelines* below. When left
  empty, the publication by Magiorakos *et al.* (see below) will be
  followed.

- col_mo:

  column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md)) -
  the default is the first column of class
  [`mo`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/2.0.0/reference/as.mo.md).

- info:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  progress should be printed to the console - the default is only print
  while in interactive sessions

- pct_required_classes:

  minimal required percentage of antimicrobial classes that must be
  available per isolate, rounded down. For example, with the default
  guideline, 17 antimicrobial classes must be available for *S. aureus*.
  Setting this `pct_required_classes` argument to `0.5` (default) means
  that for every *S. aureus* isolate at least 8 different classes must
  be available. Any lower number of available classes will return `NA`
  for that isolate.

- combine_SI:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  all values of S and I must be merged into one, so resistance is only
  considered when isolates are R, not I. As this is the default
  behaviour of the `mdro()` function, it follows the redefinition by
  EUCAST about the interpretation of I (increased exposure) in 2019, see
  section 'Interpretation of S, I and R' below. When using
  `combine_SI = FALSE`, resistance is considered when isolates are R or
  I.

- verbose:

  a [logical](https://rdrr.io/r/base/logical.html) to turn Verbose mode
  on and off (default is off). In Verbose mode, the function does not
  return the MDRO results, but instead returns a data set in logbook
  form with extensive info about which isolates would be MDRO-positive,
  or why they are not.

- only_sir_columns:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only antibiotic columns must be detected that were transformed to
  class `sir` (see
  [`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md))
  on beforehand (default is `FALSE`)

- ...:

  in case of `custom_mdro_guideline()`: a set of rules, see section
  *Using Custom Guidelines* below. Otherwise: column name of an
  antibiotic, see section *Antibiotics* below.

- as_factor:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the returned value should be an ordered
  [factor](https://rdrr.io/pkg/data.table/man/fctr.html) (`TRUE`,
  default), or otherwise a
  [character](https://rdrr.io/r/base/character.html) vector

## Value

- CMI 2012 paper - function `mdr_cmi2012()` or `mdro()`:  
  Ordered [factor](https://rdrr.io/pkg/data.table/man/fctr.html) with
  levels `Negative` \< `Multi-drug-resistant (MDR)` \<
  `Extensively drug-resistant (XDR)` \< `Pandrug-resistant (PDR)`

- TB guideline - function `mdr_tb()` or `mdro(..., guideline = "TB")`:  
  Ordered [factor](https://rdrr.io/pkg/data.table/man/fctr.html) with
  levels `Negative` \< `Mono-resistant` \< `Poly-resistant` \<
  `Multi-drug-resistant` \< `Extensively drug-resistant`

- German guideline - function `mrgn()` or
  `mdro(..., guideline = "MRGN")`:  
  Ordered [factor](https://rdrr.io/pkg/data.table/man/fctr.html) with
  levels `Negative` \< `3MRGN` \< `4MRGN`

- Everything else, except for custom guidelines:  
  Ordered [factor](https://rdrr.io/pkg/data.table/man/fctr.html) with
  levels `Negative` \< `Positive, unconfirmed` \< `Positive`. The value
  `"Positive, unconfirmed"` means that, according to the guideline, it
  is not entirely sure if the isolate is multi-drug resistant and this
  should be confirmed with additional (e.g. molecular) tests

## Details

These functions are context-aware. This means that the `x` argument can
be left blank if used inside a
[data.frame](https://rdrr.io/r/base/data.frame.html) call, see
*Examples*.

For the `pct_required_classes` argument, values above 1 will be divided
by 100. This is to support both fractions (`0.75` or `3/4`) and
percentages (`75`).

**Note:** Every test that involves the Enterobacteriaceae family, will
internally be performed using its newly named *order* Enterobacterales,
since the Enterobacteriaceae family has been taxonomically reclassified
by Adeolu *et al.* in 2016. Before that, Enterobacteriaceae was the only
family under the Enterobacteriales (with an i) order. All species under
the old Enterobacteriaceae family are still under the new
Enterobacterales (without an i) order, but divided into multiple
families. The way tests are performed now by this `mdro()` function
makes sure that results from before 2016 and after 2016 are identical.

## Supported International / National Guidelines

Currently supported guidelines are (case-insensitive):

- `guideline = "CMI2012"` (default)

  Magiorakos AP, Srinivasan A *et al.* "Multidrug-resistant, extensively
  drug-resistant and pandrug-resistant bacteria: an international expert
  proposal for interim standard definitions for acquired resistance."
  Clinical Microbiology and Infection (2012)
  ([link](https://www.clinicalmicrobiologyandinfection.com/article/S1198-743X(14)61632-3/fulltext))

- `guideline = "EUCAST3.3"` (or simply `guideline = "EUCAST"`)

  The European international guideline - EUCAST Expert Rules Version 3.3
  "Intrinsic Resistance and Unusual Phenotypes"
  ([link](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2021/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.3_20211018.pdf))

- `guideline = "EUCAST3.2"`

  The European international guideline - EUCAST Expert Rules Version 3.2
  "Intrinsic Resistance and Unusual Phenotypes"
  ([link](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2020/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.2_20200225.pdf))

- `guideline = "EUCAST3.1"`

  The European international guideline - EUCAST Expert Rules Version 3.1
  "Intrinsic Resistance and Exceptional Phenotypes Tables"
  ([link](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf))

- `guideline = "TB"`

  The international guideline for multi-drug resistant tuberculosis -
  World Health Organization "Companion handbook to the WHO guidelines
  for the programmatic management of drug-resistant tuberculosis"
  ([link](https://www.who.int/publications/i/item/9789241548809))

- `guideline = "MRGN"`

  The German national guideline - Mueller et al. (2015) Antimicrobial
  Resistance and Infection Control 4:7;
  [doi:10.1186/s13756-015-0047-6](https://doi.org/10.1186/s13756-015-0047-6)

- `guideline = "BRMO"`

  The Dutch national guideline - Rijksinstituut voor Volksgezondheid en
  Milieu "WIP-richtlijn BRMO (Bijzonder Resistente Micro-Organismen)
  (ZKH)"
  ([link](https://www.rivm.nl/wip-richtlijn-brmo-bijzonder-resistente-micro-organismen-zkh))

Please suggest your own (country-specific) guidelines by letting us
know: <https://github.com/msberends/AMR/issues/new>.

## Using Custom Guidelines

Custom guidelines can be set with the `custom_mdro_guideline()`
function. This is of great importance if you have custom rules to
determine MDROs in your hospital, e.g., rules that are dependent on
ward, state of contact isolation or other variables in your data.

If you are familiar with the
[`case_when()`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html)
function of the `dplyr` package, you will recognise the input method to
set your own rules. Rules must be set using what R considers to be the
'formula notation'. The rule is written *before* the tilde (`~`) and the
consequence of the rule is written *after* the tilde:

    custom <- custom_mdro_guideline(CIP == "R" & age > 60 ~ "Elderly Type A",
                                    ERY == "R" & age > 60 ~ "Elderly Type B")

If a row/an isolate matches the first rule, the value after the first
`~` (in this case *'Elderly Type A'*) will be set as MDRO value.
Otherwise, the second rule will be tried and so on. The number of rules
is unlimited.

You can print the rules set in the console for an overview. Colours will
help reading it if your console supports colours.

    custom
    #> A set of custom MDRO rules:
    #>   1. CIP is "R" and age is higher than 60 -> Elderly Type A
    #>   2. ERY is "R" and age is higher than 60 -> Elderly Type B
    #>   3. Otherwise -> Negative
    #>
    #> Unmatched rows will return NA.

The outcome of the function can be used for the `guideline` argument in
the `mdro()` function:

    x <- mdro(example_isolates,
              guideline = custom)
    table(x)
    #>       Negative Elderly Type A Elderly Type B
    #>           1070            198            732

Rules can also be combined with other custom rules by using
[`c()`](https://rdrr.io/r/base/c.html):

    x <- mdro(example_isolates,
              guideline = c(custom,
                            custom_mdro_guideline(ERY == "R" & age > 50 ~ "Elderly Type C")))
    table(x)
    #>       Negative Elderly Type A Elderly Type B Elderly Type C
    #>            961            198            732            109

The rules set (the `custom` object in this case) could be exported to a
shared file location using
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) if you collaborate
with multiple users. The custom rules set could then be imported using
[`readRDS()`](https://rdrr.io/r/base/readRDS.html).

## Antibiotics

To define antibiotics column names, leave as it is to determine it
automatically with
[`guess_ab_col()`](https://amr-for-r.org/release/2.0.0/reference/guess_ab_col.md)
or input a text (case-insensitive), or use `NULL` to skip a column (e.g.
`TIC = NULL` to skip ticarcillin). Manually defined but non-existing
columns will be skipped with a warning.

The following antibiotics are eligible for the functions
[`eucast_rules()`](https://amr-for-r.org/release/2.0.0/reference/eucast_rules.md)
and `mdro()`. These are shown below in the format 'name
(`antimicrobial ID`, [ATC
code](https://www.whocc.no/atc/structure_and_principles/))', sorted
alphabetically:

Amikacin (`AMK`,
[S01AE08](https://www.whocc.no/atc_ddd_index/?code=J01GB06&showdescription=no)),
amoxicillin (`AMX`,
[J01MA02](https://www.whocc.no/atc_ddd_index/?code=J01CA04&showdescription=no)),
amoxicillin/clavulanic acid (`AMC`,
[J01MA23](https://www.whocc.no/atc_ddd_index/?code=J01CR02&showdescription=no)),
ampicillin (`AMP`,
[J01MA04](https://www.whocc.no/atc_ddd_index/?code=J01CA01&showdescription=no)),
ampicillin/sulbactam (`SAM`,
[J01MA08](https://www.whocc.no/atc_ddd_index/?code=J01CR01&showdescription=no)),
arbekacin (`ARB`,
[J01MA19](https://www.whocc.no/atc_ddd_index/?code=J01GB12&showdescription=no)),
aspoxicillin (`APX`,
[J01MA16](https://www.whocc.no/atc_ddd_index/?code=J01CA19&showdescription=no)),
azidocillin (`AZD`,
[J01MA15](https://www.whocc.no/atc_ddd_index/?code=J01CE04&showdescription=no)),
azithromycin (`AZM`,
[J01MA11](https://www.whocc.no/atc_ddd_index/?code=J01FA10&showdescription=no)),
azlocillin (`AZL`,
[J01MA25](https://www.whocc.no/atc_ddd_index/?code=J01CA09&showdescription=no)),
aztreonam (`ATM`,
[J01MA12](https://www.whocc.no/atc_ddd_index/?code=J01DF01&showdescription=no)),
bacampicillin (`BAM`,
[J01MA24](https://www.whocc.no/atc_ddd_index/?code=J01CA06&showdescription=no)),
bekanamycin (`BEK`,
[J01MA07](https://www.whocc.no/atc_ddd_index/?code=J01GB13&showdescription=no)),
benzathine benzylpenicillin (`BNB`,
[J01MA14](https://www.whocc.no/atc_ddd_index/?code=J01CE08&showdescription=no)),
benzathine phenoxymethylpenicillin (`BNP`,
[D10AF05](https://www.whocc.no/atc_ddd_index/?code=J01CE10&showdescription=no)),
benzylpenicillin (`PEN`,
[J01MA06](https://www.whocc.no/atc_ddd_index/?code=J01CE01&showdescription=no)),
besifloxacin (`BES`,
[J01MA01](https://www.whocc.no/atc_ddd_index/?code=S01AE08&showdescription=no)),
biapenem (`BIA`,
[J01MA18](https://www.whocc.no/atc_ddd_index/?code=J01DH05&showdescription=no)),
carbenicillin (`CRB`,
[J01MA03](https://www.whocc.no/atc_ddd_index/?code=J01CA03&showdescription=no)),
carindacillin (`CRN`,
[J01MA17](https://www.whocc.no/atc_ddd_index/?code=J01CA05&showdescription=no)),
cefacetrile (`CAC`,
[J01MA10](https://www.whocc.no/atc_ddd_index/?code=J01DB10&showdescription=no)),
cefaclor (`CEC`,
[J01MA21](https://www.whocc.no/atc_ddd_index/?code=J01DC04&showdescription=no)),
cefadroxil (`CFR`,
[J01MA09](https://www.whocc.no/atc_ddd_index/?code=J01DB05&showdescription=no)),
cefalexin (`LEX`,
[J01MA05](https://www.whocc.no/atc_ddd_index/?code=J01DB01&showdescription=no)),
cefaloridine (`RID`,
[P01AA05](https://www.whocc.no/atc_ddd_index/?code=J01DB02&showdescription=no)),
cefalotin (`CEP`,
[J01MA22](https://www.whocc.no/atc_ddd_index/?code=J01DB03&showdescription=no)),
cefamandole (`MAN`,
[J01MA13](https://www.whocc.no/atc_ddd_index/?code=J01DC03&showdescription=no)),
cefapirin (`HAP`,
[J01CA01](https://www.whocc.no/atc_ddd_index/?code=J01DB08&showdescription=no)),
cefatrizine (`CTZ`,
[J01CA04](https://www.whocc.no/atc_ddd_index/?code=J01DB07&showdescription=no)),
cefazedone (`CZD`,
[J01CA12](https://www.whocc.no/atc_ddd_index/?code=J01DB06&showdescription=no)),
cefazolin (`CZO`,
[J01CR05](https://www.whocc.no/atc_ddd_index/?code=J01DB04&showdescription=no)),
cefcapene (`CCP`,
[J01CA13](https://www.whocc.no/atc_ddd_index/?code=J01DD17&showdescription=no)),
cefdinir (`CDR`,
[J01AA02](https://www.whocc.no/atc_ddd_index/?code=J01DD15&showdescription=no)),
cefditoren (`DIT`,
[J01FA10](https://www.whocc.no/atc_ddd_index/?code=J01DD16&showdescription=no)),
cefepime (`FEP`,
[J01FA09](https://www.whocc.no/atc_ddd_index/?code=J01DE01&showdescription=no)),
cefetamet (`CAT`,
[J01CR02](https://www.whocc.no/atc_ddd_index/?code=J01DD10&showdescription=no)),
cefixime (`CFM`,
[J01AA08](https://www.whocc.no/atc_ddd_index/?code=J01DD08&showdescription=no)),
cefmenoxime (`CMX`,
[J01FA06](https://www.whocc.no/atc_ddd_index/?code=J01DD05&showdescription=no)),
cefmetazole (`CMZ`,
[J01CF04](https://www.whocc.no/atc_ddd_index/?code=J01DC09&showdescription=no)),
cefodizime (`DIZ`,
[J01CF05](https://www.whocc.no/atc_ddd_index/?code=J01DD09&showdescription=no)),
cefonicid (`CID`,
[J01CR01](https://www.whocc.no/atc_ddd_index/?code=J01DC06&showdescription=no)),
cefoperazone (`CFP`,
[J01CA19](https://www.whocc.no/atc_ddd_index/?code=J01DD12&showdescription=no)),
cefoperazone/sulbactam (`CSL`,
[J01CE04](https://www.whocc.no/atc_ddd_index/?code=J01DD62&showdescription=no)),
ceforanide (`CND`,
[J01CA09](https://www.whocc.no/atc_ddd_index/?code=J01DC11&showdescription=no)),
cefotaxime (`CTX`,
[J01DF01](https://www.whocc.no/atc_ddd_index/?code=J01DD01&showdescription=no)),
cefotaxime/clavulanic acid (`CTC`,
[J01CA06](https://www.whocc.no/atc_ddd_index/?code=J01DD51&showdescription=no)),
cefotetan (`CTT`,
[J01CE08](https://www.whocc.no/atc_ddd_index/?code=J01DC05&showdescription=no)),
cefotiam (`CTF`,
[J01CE10](https://www.whocc.no/atc_ddd_index/?code=J01DC07&showdescription=no)),
cefoxitin (`FOX`,
[J01CE01](https://www.whocc.no/atc_ddd_index/?code=J01DC01&showdescription=no)),
cefozopran (`ZOP`,
[J01CA03](https://www.whocc.no/atc_ddd_index/?code=J01DE03&showdescription=no)),
cefpiramide (`CPM`,
[J01CA05](https://www.whocc.no/atc_ddd_index/?code=J01DD11&showdescription=no)),
cefpirome (`CPO`,
[J01CE07](https://www.whocc.no/atc_ddd_index/?code=J01DE02&showdescription=no)),
cefpodoxime (`CPD`,
[J01CF02](https://www.whocc.no/atc_ddd_index/?code=J01DD13&showdescription=no)),
cefprozil (`CPR`,
[J01CF01](https://www.whocc.no/atc_ddd_index/?code=J01DC10&showdescription=no)),
cefroxadine (`CRD`,
[J01CA07](https://www.whocc.no/atc_ddd_index/?code=J01DB11&showdescription=no)),
cefsulodin (`CFS`,
[J01CA18](https://www.whocc.no/atc_ddd_index/?code=J01DD03&showdescription=no)),
ceftaroline (`CPT`,
[J01CA11](https://www.whocc.no/atc_ddd_index/?code=J01DI02&showdescription=no)),
ceftazidime (`CAZ`,
[J01CA14](https://www.whocc.no/atc_ddd_index/?code=J01DD02&showdescription=no)),
ceftazidime/clavulanic acid (`CCV`,
[J01CF03](https://www.whocc.no/atc_ddd_index/?code=J01DD52&showdescription=no)),
cefteram (`CEM`,
[J01CA10](https://www.whocc.no/atc_ddd_index/?code=J01DD18&showdescription=no)),
ceftezole (`CTL`,
[J01CF06](https://www.whocc.no/atc_ddd_index/?code=J01DB12&showdescription=no)),
ceftibuten (`CTB`,
[J01CE06](https://www.whocc.no/atc_ddd_index/?code=J01DD14&showdescription=no)),
ceftizoxime (`CZX`,
[J01CE05](https://www.whocc.no/atc_ddd_index/?code=J01DD07&showdescription=no)),
ceftobiprole medocaril (`CFM1`,
[J01CE02](https://www.whocc.no/atc_ddd_index/?code=J01DI01&showdescription=no)),
ceftolozane/tazobactam (`CZT`,
[J01CA02](https://www.whocc.no/atc_ddd_index/?code=J01DI54&showdescription=no)),
ceftriaxone (`CRO`,
[J01CA08](https://www.whocc.no/atc_ddd_index/?code=J01DD04&showdescription=no)),
ceftriaxone/beta-lactamase inhibitor (`CEB`,
[J01CE09](https://www.whocc.no/atc_ddd_index/?code=J01DD63&showdescription=no)),
cefuroxime (`CXM`,
[J01CE03](https://www.whocc.no/atc_ddd_index/?code=J01DC02&showdescription=no)),
cephradine (`CED`,
[J01CG01](https://www.whocc.no/atc_ddd_index/?code=J01DB09&showdescription=no)),
chloramphenicol (`CHL`,
[J01CA16](https://www.whocc.no/atc_ddd_index/?code=J01BA01&showdescription=no)),
ciprofloxacin (`CIP`,
[J01CR04](https://www.whocc.no/atc_ddd_index/?code=J01MA02&showdescription=no)),
clarithromycin (`CLR`,
[J01CA15](https://www.whocc.no/atc_ddd_index/?code=J01FA09&showdescription=no)),
clindamycin (`CLI`,
[J01CG02](https://www.whocc.no/atc_ddd_index/?code=J01FF01&showdescription=no)),
clometocillin (`CLM`,
[J01CA17](https://www.whocc.no/atc_ddd_index/?code=J01CE07&showdescription=no)),
cloxacillin (`CLO`,
[J01CR03](https://www.whocc.no/atc_ddd_index/?code=J01CF02&showdescription=no)),
colistin (`COL`,
[J01DB10](https://www.whocc.no/atc_ddd_index/?code=J01XB01&showdescription=no)),
cycloserine (`CYC`,
[J01DC04](https://www.whocc.no/atc_ddd_index/?code=J04AB01&showdescription=no)),
dalbavancin (`DAL`,
[J01DB05](https://www.whocc.no/atc_ddd_index/?code=J01XA04&showdescription=no)),
daptomycin (`DAP`,
[J01DB01](https://www.whocc.no/atc_ddd_index/?code=J01XX09&showdescription=no)),
delafloxacin (`DFX`,
[J01DB02](https://www.whocc.no/atc_ddd_index/?code=J01MA23&showdescription=no)),
dibekacin (`DKB`,
[J01DB03](https://www.whocc.no/atc_ddd_index/?code=J01GB09&showdescription=no)),
dicloxacillin (`DIC`,
[J01DC03](https://www.whocc.no/atc_ddd_index/?code=J01CF01&showdescription=no)),
dirithromycin (`DIR`,
[J01DB08](https://www.whocc.no/atc_ddd_index/?code=J01FA13&showdescription=no)),
doripenem (`DOR`,
[J01DB07](https://www.whocc.no/atc_ddd_index/?code=J01DH04&showdescription=no)),
doxycycline (`DOX`,
[J01DB06](https://www.whocc.no/atc_ddd_index/?code=J01AA02&showdescription=no)),
enoxacin (`ENX`,
[J01DB04](https://www.whocc.no/atc_ddd_index/?code=J01MA04&showdescription=no)),
epicillin (`EPC`,
[J01DD17](https://www.whocc.no/atc_ddd_index/?code=J01CA07&showdescription=no)),
ertapenem (`ETP`,
[J01DD15](https://www.whocc.no/atc_ddd_index/?code=J01DH03&showdescription=no)),
erythromycin (`ERY`,
[J01DD16](https://www.whocc.no/atc_ddd_index/?code=J01FA01&showdescription=no)),
fleroxacin (`FLE`,
[J01DE01](https://www.whocc.no/atc_ddd_index/?code=J01MA08&showdescription=no)),
flucloxacillin (`FLC`,
[J01DD10](https://www.whocc.no/atc_ddd_index/?code=J01CF05&showdescription=no)),
flurithromycin (`FLR1`,
[J01DD08](https://www.whocc.no/atc_ddd_index/?code=J01FA14&showdescription=no)),
fosfomycin (`FOS`,
[J01DD05](https://www.whocc.no/atc_ddd_index/?code=J01XX01&showdescription=no)),
framycetin (`FRM`,
[J01DC09](https://www.whocc.no/atc_ddd_index/?code=D09AA01&showdescription=no)),
fusidic acid (`FUS`,
[J01DD09](https://www.whocc.no/atc_ddd_index/?code=J01XC01&showdescription=no)),
garenoxacin (`GRN`,
[J01DC06](https://www.whocc.no/atc_ddd_index/?code=J01MA19&showdescription=no)),
gatifloxacin (`GAT`,
[J01DD12](https://www.whocc.no/atc_ddd_index/?code=J01MA16&showdescription=no)),
gemifloxacin (`GEM`,
[J01DD62](https://www.whocc.no/atc_ddd_index/?code=J01MA15&showdescription=no)),
gentamicin (`GEN`,
[J01DC11](https://www.whocc.no/atc_ddd_index/?code=J01GB03&showdescription=no)),
grepafloxacin (`GRX`,
[J01DD01](https://www.whocc.no/atc_ddd_index/?code=J01MA11&showdescription=no)),
hetacillin (`HET`,
[J01DD51](https://www.whocc.no/atc_ddd_index/?code=J01CA18&showdescription=no)),
imipenem (`IPM`,
[J01DC05](https://www.whocc.no/atc_ddd_index/?code=J01DH51&showdescription=no)),
imipenem/relebactam (`IMR`,
[J01DC07](https://www.whocc.no/atc_ddd_index/?code=J01DH56&showdescription=no)),
isepamicin (`ISE`,
[J01DC01](https://www.whocc.no/atc_ddd_index/?code=J01GB11&showdescription=no)),
josamycin (`JOS`,
[J01DE03](https://www.whocc.no/atc_ddd_index/?code=J01FA07&showdescription=no)),
kanamycin (`KAN`,
[J01DD11](https://www.whocc.no/atc_ddd_index/?code=J01GB04&showdescription=no)),
lascufloxacin (`LSC`,
[J01DE02](https://www.whocc.no/atc_ddd_index/?code=J01MA25&showdescription=no)),
latamoxef (`LTM`,
[J01DD13](https://www.whocc.no/atc_ddd_index/?code=J01DD06&showdescription=no)),
levofloxacin (`LVX`,
[J01DC10](https://www.whocc.no/atc_ddd_index/?code=J01MA12&showdescription=no)),
levonadifloxacin (`LND`,
[J01DB11](https://www.whocc.no/atc_ddd_index/?code=J01MA24&showdescription=no)),
lincomycin (`LIN`,
[J01DD03](https://www.whocc.no/atc_ddd_index/?code=J01FF02&showdescription=no)),
linezolid (`LNZ`,
[J01DI02](https://www.whocc.no/atc_ddd_index/?code=J01XX08&showdescription=no)),
lomefloxacin (`LOM`,
[J01DD02](https://www.whocc.no/atc_ddd_index/?code=J01MA07&showdescription=no)),
loracarbef (`LOR`,
[J01DD52](https://www.whocc.no/atc_ddd_index/?code=J01DC08&showdescription=no)),
mecillinam (`MEC`,
[J01DD18](https://www.whocc.no/atc_ddd_index/?code=J01CA11&showdescription=no)),
meropenem (`MEM`,
[J01DB12](https://www.whocc.no/atc_ddd_index/?code=J01DH02&showdescription=no)),
meropenem/vaborbactam (`MEV`,
[J01DD14](https://www.whocc.no/atc_ddd_index/?code=J01DH52&showdescription=no)),
metampicillin (`MTM`,
[J01DD07](https://www.whocc.no/atc_ddd_index/?code=J01CA14&showdescription=no)),
meticillin (`MET`,
[J01DI01](https://www.whocc.no/atc_ddd_index/?code=J01CF03&showdescription=no)),
mezlocillin (`MEZ`,
[J01DI54](https://www.whocc.no/atc_ddd_index/?code=J01CA10&showdescription=no)),
micronomicin (`MCR`,
[J01DD04](https://www.whocc.no/atc_ddd_index/?code=S01AA22&showdescription=no)),
midecamycin (`MID`,
[J01DD63](https://www.whocc.no/atc_ddd_index/?code=J01FA03&showdescription=no)),
minocycline (`MNO`,
[J01DC02](https://www.whocc.no/atc_ddd_index/?code=J01AA08&showdescription=no)),
miocamycin (`MCM`,
[J01DB09](https://www.whocc.no/atc_ddd_index/?code=J01FA11&showdescription=no)),
moxifloxacin (`MFX`,
[J01DD06](https://www.whocc.no/atc_ddd_index/?code=J01MA14&showdescription=no)),
nadifloxacin (`NAD`,
[J01DC08](https://www.whocc.no/atc_ddd_index/?code=D10AF05&showdescription=no)),
nafcillin (`NAF`,
[J01DH05](https://www.whocc.no/atc_ddd_index/?code=J01CF06&showdescription=no)),
nalidixic acid (`NAL`,
[J01DH04](https://www.whocc.no/atc_ddd_index/?code=J01MB02&showdescription=no)),
neomycin (`NEO`,
[J01DH03](https://www.whocc.no/atc_ddd_index/?code=J01GB05&showdescription=no)),
netilmicin (`NET`,
[J01DH51](https://www.whocc.no/atc_ddd_index/?code=J01GB07&showdescription=no)),
nitrofurantoin (`NIT`,
[J01DH56](https://www.whocc.no/atc_ddd_index/?code=J01XE01&showdescription=no)),
norfloxacin (`NOR`,
[J01DH02](https://www.whocc.no/atc_ddd_index/?code=J01MA06&showdescription=no)),
ofloxacin (`OFX`,
[J01DH52](https://www.whocc.no/atc_ddd_index/?code=J01MA01&showdescription=no)),
oleandomycin (`OLE`,
[J01DH55](https://www.whocc.no/atc_ddd_index/?code=J01FA05&showdescription=no)),
oritavancin (`ORI`,
[J01DH06](https://www.whocc.no/atc_ddd_index/?code=J01XA05&showdescription=no)),
oxacillin (`OXA`,
[J01XA02](https://www.whocc.no/atc_ddd_index/?code=J01CF04&showdescription=no)),
panipenem (`PAN`,
[J01XA01](https://www.whocc.no/atc_ddd_index/?code=J01DH55&showdescription=no)),
pazufloxacin (`PAZ`,
[J01XC01](https://www.whocc.no/atc_ddd_index/?code=J01MA18&showdescription=no)),
pefloxacin (`PEF`,
[J01FA13](https://www.whocc.no/atc_ddd_index/?code=J01MA03&showdescription=no)),
penamecillin (`PNM`,
[J01FA01](https://www.whocc.no/atc_ddd_index/?code=J01CE06&showdescription=no)),
pheneticillin (`PHE`,
[J01FA14](https://www.whocc.no/atc_ddd_index/?code=J01CE05&showdescription=no)),
phenoxymethylpenicillin (`PHN`,
[J01FA07](https://www.whocc.no/atc_ddd_index/?code=J01CE02&showdescription=no)),
piperacillin (`PIP`,
[J01FA03](https://www.whocc.no/atc_ddd_index/?code=J01CA12&showdescription=no)),
piperacillin/tazobactam (`TZP`,
[J01FA11](https://www.whocc.no/atc_ddd_index/?code=J01CR05&showdescription=no)),
pivampicillin (`PVM`,
[J01FA05](https://www.whocc.no/atc_ddd_index/?code=J01CA02&showdescription=no)),
pivmecillinam (`PME`,
[J01FA12](https://www.whocc.no/atc_ddd_index/?code=J01CA08&showdescription=no)),
plazomicin (`PLZ`,
[J01FA16](https://www.whocc.no/atc_ddd_index/?code=J01GB14&showdescription=no)),
polymyxin B (`PLB`,
[J01FA02](https://www.whocc.no/atc_ddd_index/?code=J01XB02&showdescription=no)),
pristinamycin (`PRI`,
[J01FA15](https://www.whocc.no/atc_ddd_index/?code=J01FG01&showdescription=no)),
procaine benzylpenicillin (`PRB`,
[J01FA08](https://www.whocc.no/atc_ddd_index/?code=J01CE09&showdescription=no)),
propicillin (`PRP`,
[J01FF02](https://www.whocc.no/atc_ddd_index/?code=J01CE03&showdescription=no)),
prulifloxacin (`PRU`,
[J01FG01](https://www.whocc.no/atc_ddd_index/?code=J01MA17&showdescription=no)),
quinupristin/dalfopristin (`QDA`,
[J01FG02](https://www.whocc.no/atc_ddd_index/?code=J01FG02&showdescription=no)),
ribostamycin (`RST`,
[J04AB02](https://www.whocc.no/atc_ddd_index/?code=J01GB10&showdescription=no)),
rifampicin (`RIF`,
[J01XX09](https://www.whocc.no/atc_ddd_index/?code=J04AB02&showdescription=no)),
rokitamycin (`ROK`,
[J01XX08](https://www.whocc.no/atc_ddd_index/?code=J01FA12&showdescription=no)),
roxithromycin (`RXT`,
[J01AA07](https://www.whocc.no/atc_ddd_index/?code=J01FA06&showdescription=no)),
rufloxacin (`RFL`,
[J01XB01](https://www.whocc.no/atc_ddd_index/?code=J01MA10&showdescription=no)),
sisomicin (`SIS`,
[J01XB02](https://www.whocc.no/atc_ddd_index/?code=J01GB08&showdescription=no)),
sitafloxacin (`SIT`,
[J01XE01](https://www.whocc.no/atc_ddd_index/?code=J01MA21&showdescription=no)),
solithromycin (`SOL`,
[J01AA12](https://www.whocc.no/atc_ddd_index/?code=J01FA16&showdescription=no)),
sparfloxacin (`SPX`,
[J01EA01](https://www.whocc.no/atc_ddd_index/?code=J01MA09&showdescription=no)),
spiramycin (`SPI`,
[J01XX01](https://www.whocc.no/atc_ddd_index/?code=J01FA02&showdescription=no)),
streptoduocin (`STR`,
[J01BA01](https://www.whocc.no/atc_ddd_index/?code=J01GA02&showdescription=no)),
streptomycin (`STR1`,
[J01GB06](https://www.whocc.no/atc_ddd_index/?code=J01GA01&showdescription=no)),
sulbactam (`SUL`,
[J01GB12](https://www.whocc.no/atc_ddd_index/?code=J01CG01&showdescription=no)),
sulbenicillin (`SBC`,
[J01GB13](https://www.whocc.no/atc_ddd_index/?code=J01CA16&showdescription=no)),
sulfadiazine (`SDI`,
[J01GB09](https://www.whocc.no/atc_ddd_index/?code=J01EC02&showdescription=no)),
sulfadiazine/trimethoprim (`SLT1`,
[D09AA01](https://www.whocc.no/atc_ddd_index/?code=J01EE02&showdescription=no)),
sulfadimethoxine (`SUD`,
[J01GB03](https://www.whocc.no/atc_ddd_index/?code=J01ED01&showdescription=no)),
sulfadimidine (`SDM`,
[J01GB11](https://www.whocc.no/atc_ddd_index/?code=J01EB03&showdescription=no)),
sulfadimidine/trimethoprim (`SLT2`,
[J01GB04](https://www.whocc.no/atc_ddd_index/?code=J01EE05&showdescription=no)),
sulfafurazole (`SLF`,
[S01AA22](https://www.whocc.no/atc_ddd_index/?code=J01EB05&showdescription=no)),
sulfaisodimidine (`SLF1`,
[J01GB05](https://www.whocc.no/atc_ddd_index/?code=J01EB01&showdescription=no)),
sulfalene (`SLF2`,
[J01GB07](https://www.whocc.no/atc_ddd_index/?code=J01ED02&showdescription=no)),
sulfamazone (`SZO`,
[J01GB14](https://www.whocc.no/atc_ddd_index/?code=J01ED09&showdescription=no)),
sulfamerazine (`SLF3`,
[J01GB10](https://www.whocc.no/atc_ddd_index/?code=J01ED07&showdescription=no)),
sulfamerazine/trimethoprim (`SLT3`,
[J01GB08](https://www.whocc.no/atc_ddd_index/?code=J01EE07&showdescription=no)),
sulfamethizole (`SLF4`,
[J01GA02](https://www.whocc.no/atc_ddd_index/?code=J01EB02&showdescription=no)),
sulfamethoxazole (`SMX`,
[J01GA01](https://www.whocc.no/atc_ddd_index/?code=J01EC01&showdescription=no)),
sulfamethoxypyridazine (`SLF5`,
[J01GB01](https://www.whocc.no/atc_ddd_index/?code=J01ED05&showdescription=no)),
sulfametomidine (`SLF6`,
[J01EE01](https://www.whocc.no/atc_ddd_index/?code=J01ED03&showdescription=no)),
sulfametoxydiazine (`SLF7`,
[J01MB02](https://www.whocc.no/atc_ddd_index/?code=J01ED04&showdescription=no)),
sulfametrole/trimethoprim (`SLT4`,
[J01FF01](https://www.whocc.no/atc_ddd_index/?code=J01EE03&showdescription=no)),
sulfamoxole (`SLF8`,
[J01XA04](https://www.whocc.no/atc_ddd_index/?code=J01EC03&showdescription=no)),
sulfamoxole/trimethoprim (`SLT5`,
[J01XA05](https://www.whocc.no/atc_ddd_index/?code=J01EE04&showdescription=no)),
sulfanilamide (`SLF9`,
[J01XA03](https://www.whocc.no/atc_ddd_index/?code=J01EB06&showdescription=no)),
sulfaperin (`SLF10`,
[J04AB01](https://www.whocc.no/atc_ddd_index/?code=J01ED06&showdescription=no)),
sulfaphenazole (`SLF11`,
[J01XX11](https://www.whocc.no/atc_ddd_index/?code=J01ED08&showdescription=no)),
sulfapyridine (`SLF12`,
[J01EC02](https://www.whocc.no/atc_ddd_index/?code=J01EB04&showdescription=no)),
sulfathiazole (`SUT`,
[J01ED01](https://www.whocc.no/atc_ddd_index/?code=J01EB07&showdescription=no)),
sulfathiourea (`SLF13`,
[J01EB03](https://www.whocc.no/atc_ddd_index/?code=J01EB08&showdescription=no)),
sultamicillin (`SLT6`,
[J01EB05](https://www.whocc.no/atc_ddd_index/?code=J01CR04&showdescription=no)),
talampicillin (`TAL`,
[J01EB01](https://www.whocc.no/atc_ddd_index/?code=J01CA15&showdescription=no)),
tazobactam (`TAZ`,
[J01ED02](https://www.whocc.no/atc_ddd_index/?code=J01CG02&showdescription=no)),
tebipenem (`TBP`,
[J01ED09](https://www.whocc.no/atc_ddd_index/?code=J01DH06&showdescription=no)),
tedizolid (`TZD`,
[J01ED07](https://www.whocc.no/atc_ddd_index/?code=J01XX11&showdescription=no)),
teicoplanin (`TEC`,
[J01EB02](https://www.whocc.no/atc_ddd_index/?code=J01XA02&showdescription=no)),
telavancin (`TLV`,
[J01EC01](https://www.whocc.no/atc_ddd_index/?code=J01XA03&showdescription=no)),
telithromycin (`TLT`,
[J01ED05](https://www.whocc.no/atc_ddd_index/?code=J01FA15&showdescription=no)),
temafloxacin (`TMX`,
[J01ED03](https://www.whocc.no/atc_ddd_index/?code=J01MA05&showdescription=no)),
temocillin (`TEM`,
[J01ED04](https://www.whocc.no/atc_ddd_index/?code=J01CA17&showdescription=no)),
tetracycline (`TCY`,
[J01EC03](https://www.whocc.no/atc_ddd_index/?code=J01AA07&showdescription=no)),
ticarcillin (`TIC`,
[J01EB06](https://www.whocc.no/atc_ddd_index/?code=J01CA13&showdescription=no)),
ticarcillin/clavulanic acid (`TCC`,
[J01ED06](https://www.whocc.no/atc_ddd_index/?code=J01CR03&showdescription=no)),
tigecycline (`TGC`,
[J01ED08](https://www.whocc.no/atc_ddd_index/?code=J01AA12&showdescription=no)),
tilbroquinol (`TBQ`,
[J01EB04](https://www.whocc.no/atc_ddd_index/?code=P01AA05&showdescription=no)),
tobramycin (`TOB`,
[J01EB07](https://www.whocc.no/atc_ddd_index/?code=J01GB01&showdescription=no)),
tosufloxacin (`TFX`,
[J01EB08](https://www.whocc.no/atc_ddd_index/?code=J01MA22&showdescription=no)),
trimethoprim (`TMP`,
[J01EE02](https://www.whocc.no/atc_ddd_index/?code=J01EA01&showdescription=no)),
trimethoprim/sulfamethoxazole (`SXT`,
[J01EE05](https://www.whocc.no/atc_ddd_index/?code=J01EE01&showdescription=no)),
troleandomycin (`TRL`,
[J01EE07](https://www.whocc.no/atc_ddd_index/?code=J01FA08&showdescription=no)),
trovafloxacin (`TVA`,
[J01EE03](https://www.whocc.no/atc_ddd_index/?code=J01MA13&showdescription=no)),
vancomycin (`VAN`,
[J01EE04](https://www.whocc.no/atc_ddd_index/?code=J01XA01&showdescription=no))

## Interpretation of SIR

In 2019, the European Committee on Antimicrobial Susceptibility Testing
(EUCAST) has decided to change the definitions of susceptibility testing
categories S, I, and R as shown below
(<https://www.eucast.org/newsiandr/>):

- **S - Susceptible, standard dosing regimen**  
  A microorganism is categorised as "Susceptible, standard dosing
  regimen", when there is a high likelihood of therapeutic success using
  a standard dosing regimen of the agent.

- **I - Susceptible, increased exposure** *  
  A microorganism is categorised as "Susceptible, Increased exposure*"
  when there is a high likelihood of therapeutic success because
  exposure to the agent is increased by adjusting the dosing regimen or
  by its concentration at the site of infection.

- **R = Resistant**  
  A microorganism is categorised as "Resistant" when there is a high
  likelihood of therapeutic failure even when there is increased
  exposure.

  - *Exposure* is a function of how the mode of administration, dose,
    dosing interval, infusion time, as well as distribution and
    excretion of the antimicrobial agent will influence the infecting
    organism at the site of infection.

This AMR package honours this insight. Use
[`susceptibility()`](https://amr-for-r.org/release/2.0.0/reference/proportion.md)
(equal to
[`proportion_SI()`](https://amr-for-r.org/release/2.0.0/reference/proportion.md))
to determine antimicrobial susceptibility and
[`count_susceptible()`](https://amr-for-r.org/release/2.0.0/reference/count.md)
(equal to
[`count_SI()`](https://amr-for-r.org/release/2.0.0/reference/count.md))
to count susceptible isolates.

## Examples

``` r
out <- mdro(example_isolates, guideline = "EUCAST")
str(out)
table(out)

out <- mdro(example_isolates,
  guideline = custom_mdro_guideline(
    AMX == "R" ~ "Custom MDRO 1",
    VAN == "R" ~ "Custom MDRO 2"
  )
)
table(out)

# \donttest{
if (require("dplyr")) {
  example_isolates %>%
    mdro() %>%
    table()

  # no need to define `x` when used inside dplyr verbs:
  example_isolates %>%
    mutate(MDRO = mdro()) %>%
    pull(MDRO) %>%
    table()
}
# }
```
