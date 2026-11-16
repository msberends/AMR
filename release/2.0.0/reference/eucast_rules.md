# Apply EUCAST Rules

Apply rules for clinical breakpoints and intrinsic resistance as defined
by the European Committee on Antimicrobial Susceptibility Testing
(EUCAST, <https://www.eucast.org>), see *Source*. Use `eucast_dosage()`
to get a [data.frame](https://rdrr.io/r/base/data.frame.html) with
advised dosages of a certain bug-drug combination, which is based on the
[dosage](https://amr-for-r.org/release/2.0.0/reference/dosage.md) data
set.

To improve the interpretation of the antibiogram before EUCAST rules are
applied, some non-EUCAST rules can applied at default, see *Details*.

## Usage

``` r
eucast_rules(
  x,
  col_mo = NULL,
  info = interactive(),
  rules = getOption("AMR_eucastrules", default = c("breakpoints", "expert")),
  verbose = FALSE,
  version_breakpoints = 12,
  version_expertrules = 3.3,
  ampc_cephalosporin_resistance = NA,
  only_sir_columns = FALSE,
  custom_rules = NULL,
  ...
)

eucast_dosage(ab, administration = "iv", version_breakpoints = 12)
```

## Source

- EUCAST Expert Rules. Version 2.0, 2012.  
  Leclercq et al. **EUCAST expert rules in antimicrobial susceptibility
  testing.** *Clin Microbiol Infect.* 2013;19(2):141-60;
  [doi:10.1111/j.1469-0691.2011.03703.x](https://doi.org/10.1111/j.1469-0691.2011.03703.x)

- EUCAST Expert Rules, Intrinsic Resistance and Exceptional Phenotypes
  Tables. Version 3.1, 2016.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/Expert_rules_intrinsic_exceptional_V3.1.pdf)

- EUCAST Intrinsic Resistance and Unusual Phenotypes. Version 3.2, 2020.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2020/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.2_20200225.pdf)

- EUCAST Intrinsic Resistance and Unusual Phenotypes. Version 3.3, 2021.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Expert_Rules/2021/Intrinsic_Resistance_and_Unusual_Phenotypes_Tables_v3.3_20211018.pdf)

- EUCAST Breakpoint tables for interpretation of MICs and zone
  diameters. Version 9.0, 2019.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_9.0_Breakpoint_Tables.xlsx)

- EUCAST Breakpoint tables for interpretation of MICs and zone
  diameters. Version 10.0, 2020.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_10.0_Breakpoint_Tables.xlsx)

- EUCAST Breakpoint tables for interpretation of MICs and zone
  diameters. Version 11.0, 2021.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_11.0_Breakpoint_Tables.xlsx)

- EUCAST Breakpoint tables for interpretation of MICs and zone
  diameters. Version 12.0, 2022.
  [(link)](https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_12.0_Breakpoint_Tables.xlsx)

## Arguments

- x:

  a data set with antibiotic columns, such as `amox`, `AMX` and `AMC`

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

- rules:

  a [character](https://rdrr.io/r/base/character.html) vector that
  specifies which rules should be applied. Must be one or more of
  `"breakpoints"`, `"expert"`, `"other"`, `"custom"`, `"all"`, and
  defaults to `c("breakpoints", "expert")`. The default value can be set
  to another value using the [package
  option](https://amr-for-r.org/release/2.0.0/reference/AMR-options.md)
  [`AMR_eucastrules`](https://amr-for-r.org/release/2.0.0/reference/AMR-options.md):
  `options(AMR_eucastrules = "all")`. If using `"custom"`, be sure to
  fill in argument `custom_rules` too. Custom rules can be created with
  [`custom_eucast_rules()`](https://amr-for-r.org/release/2.0.0/reference/custom_eucast_rules.md).

- verbose:

  a [logical](https://rdrr.io/r/base/logical.html) to turn Verbose mode
  on and off (default is off). In Verbose mode, the function does not
  apply rules to the data, but instead returns a data set in logbook
  form with extensive info about which rows and columns would be
  effected and in which way. Using Verbose mode takes a lot more time.

- version_breakpoints:

  the version number to use for the EUCAST Clinical Breakpoints
  guideline. Can be either "12.0", "11.0", or "10.0".

- version_expertrules:

  the version number to use for the EUCAST Expert Rules and Intrinsic
  Resistance guideline. Can be either "3.3", "3.2", or "3.1".

- ampc_cephalosporin_resistance:

  a [character](https://rdrr.io/r/base/character.html) value that should
  be applied to cefotaxime, ceftriaxone and ceftazidime for AmpC
  de-repressed cephalosporin-resistant mutants - the default is `NA`.
  Currently only works when `version_expertrules` is `3.2` and higher;
  these version of '*EUCAST Expert Rules on Enterobacterales*' state
  that results of cefotaxime, ceftriaxone and ceftazidime should be
  reported with a note, or results should be suppressed (emptied) for
  these three drugs. A value of `NA` (the default) for this argument
  will remove results for these three drugs, while e.g. a value of `"R"`
  will make the results for these drugs resistant. Use `NULL` or `FALSE`
  to not alter results for these three drugs of AmpC de-repressed
  cephalosporin-resistant mutants. Using `TRUE` is equal to using
  `"R"`.  
  For *EUCAST Expert Rules* v3.2, this rule applies to: *Citrobacter
  braakii*, *Citrobacter freundii*, *Citrobacter gillenii*, *Citrobacter
  murliniae*, *Citrobacter rodenticum*, *Citrobacter sedlakii*,
  *Citrobacter werkmanii*, *Citrobacter youngae*, *Enterobacter*,
  *Hafnia alvei*, *Klebsiella aerogenes*, *Morganella morganii*,
  *Providencia*, and *Serratia*.

- only_sir_columns:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only antibiotic columns must be detected that were transformed to
  class `sir` (see
  [`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md))
  on beforehand (default is `FALSE`)

- custom_rules:

  custom rules to apply, created with
  [`custom_eucast_rules()`](https://amr-for-r.org/release/2.0.0/reference/custom_eucast_rules.md)

- ...:

  column name of an antibiotic, see section *Antibiotics* below

- ab:

  any (vector of) text that can be coerced to a valid antibiotic drug
  code with
  [`as.ab()`](https://amr-for-r.org/release/2.0.0/reference/as.ab.md)

- administration:

  route of administration, either "im", "iv", or "oral"

## Value

The input of `x`, possibly with edited values of antibiotics. Or, if
`verbose = TRUE`, a [data.frame](https://rdrr.io/r/base/data.frame.html)
with all original and new values of the affected bug-drug combinations.

## Details

**Note:** This function does not translate MIC values to SIR values. Use
[`as.sir()`](https://amr-for-r.org/release/2.0.0/reference/as.sir.md)
for that.  
**Note:** When ampicillin (AMP, J01CA01) is not available but
amoxicillin (AMX, J01CA04) is, the latter will be used for all rules
where there is a dependency on ampicillin. These drugs are
interchangeable when it comes to expression of antimicrobial
resistance.  

The file containing all EUCAST rules is located here:
<https://github.com/msberends/AMR/blob/main/data-raw/eucast_rules.tsv>.
**Note:** Old taxonomic names are replaced with the current taxonomy
where applicable. For example, *Ochrobactrum anthropi* was renamed to
*Brucella anthropi* in 2020; the original EUCAST rules v3.1 and v3.2 did
not yet contain this new taxonomic name. The `AMR` package contains the
full microbial taxonomy updated until 11 December, 2022, see
[microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md).

### Custom Rules

Custom rules can be created using
[`custom_eucast_rules()`](https://amr-for-r.org/release/2.0.0/reference/custom_eucast_rules.md),
e.g.:

    x <- custom_eucast_rules(AMC == "R" & genus == "Klebsiella" ~ aminopenicillins == "R",
                             AMC == "I" & genus == "Klebsiella" ~ aminopenicillins == "I")

    eucast_rules(example_isolates, rules = "custom", custom_rules = x)

### 'Other' Rules

Before further processing, two non-EUCAST rules about drug combinations
can be applied to improve the efficacy of the EUCAST rules, and the
reliability of your data (analysis). These rules are:

1.  A drug **with** enzyme inhibitor will be set to S if the same drug
    **without** enzyme inhibitor is S

2.  A drug **without** enzyme inhibitor will be set to R if the same
    drug **with** enzyme inhibitor is R

Important examples include amoxicillin and amoxicillin/clavulanic acid,
and trimethoprim and trimethoprim/sulfamethoxazole. Needless to say, for
these rules to work, both drugs must be available in the data set.

Since these rules are not officially approved by EUCAST, they are not
applied at default. To use these rules, include `"other"` to the `rules`
argument, or use `eucast_rules(..., rules = "all")`. You can also set
the [package
option](https://amr-for-r.org/release/2.0.0/reference/AMR-options.md)
[`AMR_eucastrules`](https://amr-for-r.org/release/2.0.0/reference/AMR-options.md),
i.e. run `options(AMR_eucastrules = "all")`.

## Antibiotics

To define antibiotics column names, leave as it is to determine it
automatically with
[`guess_ab_col()`](https://amr-for-r.org/release/2.0.0/reference/guess_ab_col.md)
or input a text (case-insensitive), or use `NULL` to skip a column (e.g.
`TIC = NULL` to skip ticarcillin). Manually defined but non-existing
columns will be skipped with a warning.

The following antibiotics are eligible for the functions
`eucast_rules()` and
[`mdro()`](https://amr-for-r.org/release/2.0.0/reference/mdro.md). These
are shown below in the format 'name (`antimicrobial ID`, [ATC
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

## Reference Data Publicly Available

All data sets in this `AMR` package (about microorganisms, antibiotics,
SIR interpretation, EUCAST rules, etc.) are publicly and freely
available for download in the following formats: R, MS Excel, Apache
Feather, Apache Parquet, SPSS, SAS, and Stata. We also provide
tab-separated plain text files that are machine-readable and suitable
for input in any software program, such as laboratory information
systems. Please visit [our website for the download
links](https://msberends.github.io/AMR/articles/datasets.html). The
actual files are of course available on [our GitHub
repository](https://github.com/msberends/AMR/tree/main/data-raw).

## Examples

``` r
# \donttest{
a <- data.frame(
  mo = c(
    "Staphylococcus aureus",
    "Enterococcus faecalis",
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Pseudomonas aeruginosa"
  ),
  VAN = "-", # Vancomycin
  AMX = "-", # Amoxicillin
  COL = "-", # Colistin
  CAZ = "-", # Ceftazidime
  CXM = "-", # Cefuroxime
  PEN = "S", # Benzylpenicillin
  FOX = "S", # Cefoxitin
  stringsAsFactors = FALSE
)

head(a)


# apply EUCAST rules: some results wil be changed
b <- eucast_rules(a)

head(b)


# do not apply EUCAST rules, but rather get a data.frame
# containing all details about the transformations:
c <- eucast_rules(a, verbose = TRUE)
head(c)
# }

# Dosage guidelines:

eucast_dosage(c("tobra", "genta", "cipro"), "iv")

eucast_dosage(c("tobra", "genta", "cipro"), "iv", version_breakpoints = 10)
```
