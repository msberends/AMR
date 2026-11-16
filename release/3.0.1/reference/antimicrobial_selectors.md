# Antimicrobial Selectors

These functions allow for filtering rows and selecting columns based on
antimicrobial test results that are of a specific antimicrobial class or
group, without the need to define the columns or antimicrobial
abbreviations. They can be used in base R, tidyverse, tidymodels, and
`data.table`.

Simply puy, if you have a column name that resembles an antimicrobial
drug, it will be picked up by any of these functions that matches its
pharmaceutical class, code or name: column names "cefazolin", "kefzol",
"CZO" and "J01DB04" would all be included in the following selection:

    library(dplyr)
    my_data_with_all_these_columns %>%
      select(cephalosporins())

## Usage

``` r
aminoglycosides(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

aminopenicillins(only_sir_columns = FALSE, return_all = TRUE, ...)

antifungals(only_sir_columns = FALSE, return_all = TRUE, ...)

antimycobacterials(only_sir_columns = FALSE, return_all = TRUE, ...)

betalactams(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

betalactams_with_inhibitor(only_sir_columns = FALSE, return_all = TRUE,
  ...)

carbapenems(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

cephalosporins(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

cephalosporins_1st(only_sir_columns = FALSE, return_all = TRUE, ...)

cephalosporins_2nd(only_sir_columns = FALSE, return_all = TRUE, ...)

cephalosporins_3rd(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

cephalosporins_4th(only_sir_columns = FALSE, return_all = TRUE, ...)

cephalosporins_5th(only_sir_columns = FALSE, return_all = TRUE, ...)

fluoroquinolones(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

glycopeptides(only_sir_columns = FALSE, return_all = TRUE, ...)

isoxazolylpenicillins(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

lincosamides(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

lipoglycopeptides(only_sir_columns = FALSE, return_all = TRUE, ...)

macrolides(only_sir_columns = FALSE, return_all = TRUE, ...)

monobactams(only_sir_columns = FALSE, return_all = TRUE, ...)

nitrofurans(only_sir_columns = FALSE, return_all = TRUE, ...)

oxazolidinones(only_sir_columns = FALSE, return_all = TRUE, ...)

penicillins(only_sir_columns = FALSE, return_all = TRUE, ...)

phenicols(only_sir_columns = FALSE, return_all = TRUE, ...)

polymyxins(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

quinolones(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

rifamycins(only_sir_columns = FALSE, return_all = TRUE, ...)

streptogramins(only_sir_columns = FALSE, return_all = TRUE, ...)

sulfonamides(only_sir_columns = FALSE, return_all = TRUE, ...)

tetracyclines(only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

trimethoprims(only_sir_columns = FALSE, return_all = TRUE, ...)

ureidopenicillins(only_sir_columns = FALSE, return_all = TRUE, ...)

amr_class(amr_class, only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

amr_selector(filter, only_sir_columns = FALSE, only_treatable = TRUE,
  return_all = TRUE, ...)

administrable_per_os(only_sir_columns = FALSE, return_all = TRUE, ...)

administrable_iv(only_sir_columns = FALSE, return_all = TRUE, ...)

not_intrinsic_resistant(only_sir_columns = FALSE, col_mo = NULL,
  version_expected_phenotypes = 1.2, ...)
```

## Arguments

- only_sir_columns:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only antimicrobial columns must be included that were transformed to
  class [sir](https://amr-for-r.org/release/3.0.1/reference/as.sir.md)
  on beforehand. Defaults to `FALSE`.

- only_treatable:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  antimicrobial drugs should be excluded that are only for laboratory
  tests (default is `TRUE`), such as gentamicin-high (`GEH`) and
  imipenem/EDTA (`IPE`).

- return_all:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  all matched columns must be returned (default is `TRUE`). With
  `FALSE`, only the first of each unique antimicrobial will be returned,
  e.g. if both columns `"genta"` and `"gentamicin"` exist in the data,
  only the first hit for gentamicin will be returned.

- ...:

  Ignored, only in place to allow future extensions.

- amr_class:

  An antimicrobial class or a part of it, such as `"carba"` and
  `"carbapenems"`. The columns `group`, `atc_group1` and `atc_group2` of
  the
  [antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
  data set will be searched (case-insensitive) for this value.

- filter:

  An [expression](https://rdrr.io/r/base/expression.html) to be
  evaluated in the
  [antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
  data set, such as `name %like% "trim"`.

- col_mo:

  Column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md)) -
  the default is the first column of class
  [`mo`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/3.0.1/reference/as.mo.md).

- version_expected_phenotypes:

  The version number to use for the EUCAST Expected Phenotypes. Can be
  "1.2".

## Value

When used inside selecting or filtering, this returns a
[character](https://rdrr.io/r/base/character.html) vector of column
names, with additional class `"amr_selector"`. When used individually,
this returns an ['ab'
vector](https://amr-for-r.org/release/3.0.1/reference/as.ab.md) with all
possible antimicrobials that the function would be able to select or
filter.

## Details

These functions can be used in data set calls for selecting columns and
filtering rows. They work with base R, the Tidyverse, and `data.table`.
They are heavily inspired by the [Tidyverse selection
helpers](https://tidyselect.r-lib.org/reference/language.html) such as
[`everything()`](https://tidyselect.r-lib.org/reference/everything.html),
but are not limited to `dplyr` verbs. Nonetheless, they are very
convenient to use with `dplyr` functions such as
[`select()`](https://dplyr.tidyverse.org/reference/select.html),
[`filter()`](https://dplyr.tidyverse.org/reference/filter.html) and
[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html),
see *Examples*.

All selectors can also be used in `tidymodels` packages such as `recipe`
and `parsnip`. See for more info [our
tutorial](https://amr-for-r.org/articles/AMR_with_tidymodels.html) on
using antimicrobial selectors for predictive modelling.

All columns in the data in which these functions are called will be
searched for known antimicrobial names, abbreviations, brand names, and
codes (ATC, EARS-Net, WHO, etc.) according to the
[antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
data set. This means that a selector such as `aminoglycosides()` will
pick up column names like 'gen', 'genta', 'J01GB03', 'tobra',
'Tobracin', etc.

The `amr_class()` function can be used to filter/select on a manually
defined antimicrobial class. It searches for results in the
[antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
data set within the columns `group`, `atc_group1` and `atc_group2`.

The `administrable_per_os()` and `administrable_iv()` functions also
rely on the
[antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
data set - antimicrobials will be matched where a DDD (defined daily
dose) for resp. oral and IV treatment is available in the
[antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
data set.

The `amr_selector()` function can be used to internally filter the
[antimicrobials](https://amr-for-r.org/release/3.0.1/reference/antimicrobials.md)
data set on any results, see *Examples*. It allows for filtering on a
(part of) a certain name, and/or a group name or even a minimum of DDDs
for oral treatment. This function yields the highest flexibility, but is
also the least user-friendly, since it requires a hard-coded filter to
set.

The `not_intrinsic_resistant()` function can be used to only select
antimicrobials that pose no intrinsic resistance for the microorganisms
in the data set. For example, if a data set contains only microorganism
codes or names of *E. coli* and *K. pneumoniae* and contains a column
"vancomycin", this column will be removed (or rather, unselected) using
this function. It currently applies ['EUCAST Expected Resistant
Phenotypes'
v1.2](https://www.eucast.org/expert_rules_and_expected_phenotypes)
(2023) to determine intrinsic resistance, using the
[`eucast_rules()`](https://amr-for-r.org/release/3.0.1/reference/eucast_rules.md)
function internally. Because of this determination, this function is
quite slow in terms of performance.

## Full list of supported (antimicrobial) classes

- `aminoglycosides()` can select:  
  amikacin (AMK), amikacin/fosfomycin (AKF), apramycin (APR), arbekacin
  (ARB), astromicin (AST), bekanamycin (BEK), dibekacin (DKB),
  framycetin (FRM), gentamicin (GEN), gentamicin-high (GEH), habekacin
  (HAB), hygromycin (HYG), isepamicin (ISE), kanamycin (KAN),
  kanamycin-high (KAH), kanamycin/cephalexin (KAC), micronomicin (MCR),
  neomycin (NEO), netilmicin (NET), pentisomicin (PIM), plazomicin
  (PLZ), propikacin (PKA), ribostamycin (RST), sisomicin (SIS),
  streptoduocin (STR), streptomycin (STR1), streptomycin-high (STH),
  tobramycin (TOB), and tobramycin-high (TOH)

- `aminopenicillins()` can select:  
  amoxicillin (AMX) and ampicillin (AMP)

- `antifungals()` can select:  
  amorolfine (AMO), amphotericin B (AMB), amphotericin B-high (AMH),
  anidulafungin (ANI), butoconazole (BUT), caspofungin (CAS), ciclopirox
  (CIX), clotrimazole (CTR), econazole (ECO), fluconazole (FLU),
  flucytosine (FCT), fosfluconazole (FFL), griseofulvin (GRI),
  hachimycin (HCH), ibrexafungerp (IBX), isavuconazole (ISV),
  isoconazole (ISO), itraconazole (ITR), ketoconazole (KET), manogepix
  (MGX), micafungin (MIF), miconazole (MCZ), nystatin (NYS),
  oteseconazole (OTE), pimaricin (PMR), posaconazole (POS), rezafungin
  (RZF), ribociclib (RBC), sulconazole (SUC), terbinafine (TRB),
  terconazole (TRC), and voriconazole (VOR)

- `antimycobacterials()` can select:  
  4-aminosalicylic acid (AMA), calcium aminosalicylate (CLA),
  capreomycin (CAP), clofazimine (CLF), delamanid (DLM), enviomycin
  (ENV), ethambutol (ETH), ethambutol/isoniazid (ETI), ethionamide
  (ETI1), isoniazid (INH),
  isoniazid/sulfamethoxazole/trimethoprim/pyridoxine (IST), morinamide
  (MRN), p-aminosalicylic acid (PAS), pretomanid (PMD), protionamide
  (PTH), pyrazinamide (PZA), rifabutin (RIB), rifampicin (RIF),
  rifampicin/ethambutol/isoniazid (REI), rifampicin/isoniazid (RFI),
  rifampicin/pyrazinamide/ethambutol/isoniazid (RPEI),
  rifampicin/pyrazinamide/isoniazid (RPI), rifamycin (RFM), rifapentine
  (RFP), sodium aminosalicylate (SDA), streptomycin/isoniazid (STI),
  terizidone (TRZ), thioacetazone (TAT), thioacetazone/isoniazid (THI1),
  tiocarlide (TCR), and viomycin (VIO)

- `betalactams()` can select:  
  amoxicillin (AMX), amoxicillin/clavulanic acid (AMC),
  amoxicillin/sulbactam (AXS), ampicillin (AMP), ampicillin/sulbactam
  (SAM), apalcillin (APL), aspoxicillin (APX), azidocillin (AZD),
  azlocillin (AZL), aztreonam (ATM), aztreonam/avibactam (AZA),
  aztreonam/nacubactam (ANC), bacampicillin (BAM), benzathine
  benzylpenicillin (BNB), benzathine phenoxymethylpenicillin (BNP),
  benzylpenicillin (PEN), benzylpenicillin screening test (PEN-S),
  biapenem (BIA), carbenicillin (CRB), carindacillin (CRN), carumonam
  (CAR), cefacetrile (CAC), cefaclor (CEC), cefadroxil (CFR), cefalexin
  (LEX), cefaloridine (RID), cefalotin (CEP), cefamandole (MAN),
  cefapirin (HAP), cefatrizine (CTZ), cefazedone (CZD), cefazolin (CZO),
  cefcapene (CCP), cefcapene pivoxil (CCX), cefdinir (CDR), cefditoren
  (DIT), cefditoren pivoxil (DIX), cefepime (FEP), cefepime/amikacin
  (CFA), cefepime/clavulanic acid (CPC), cefepime/enmetazobactam (FPE),
  cefepime/nacubactam (FNC), cefepime/tazobactam (FPT),
  cefepime/zidebactam (FPZ), cefetamet (CAT), cefetamet pivoxil (CPI),
  cefetecol (CCL), cefetrizole (CZL), cefiderocol (FDC), cefixime (CFM),
  cefmenoxime (CMX), cefmetazole (CMZ), cefodizime (DIZ), cefonicid
  (CID), cefoperazone (CFP), cefoperazone/sulbactam (CSL), ceforanide
  (CND), cefoselis (CSE), cefotaxime (CTX), cefotaxime screening test
  (CTX-S), cefotaxime/clavulanic acid (CTC), cefotaxime/sulbactam (CTS),
  cefotetan (CTT), cefotiam (CTF), cefotiam hexetil (CHE), cefovecin
  (FOV), cefoxitin (FOX), cefoxitin screening test (FOX-S), cefozopran
  (ZOP), cefpimizole (CFZ), cefpiramide (CPM), cefpirome (CPO),
  cefpodoxime (CPD), cefpodoxime proxetil (CPX), cefpodoxime/clavulanic
  acid (CDC), cefprozil (CPR), cefquinome (CEQ), cefroxadine (CRD),
  cefsulodin (CFS), cefsumide (CSU), ceftaroline (CPT),
  ceftaroline/avibactam (CPA), ceftazidime (CAZ), ceftazidime/avibactam
  (CZA), ceftazidime/clavulanic acid (CCV), cefteram (CEM), cefteram
  pivoxil (CPL), ceftezole (CTL), ceftibuten (CTB), ceftiofur (TIO),
  ceftizoxime (CZX), ceftizoxime alapivoxil (CZP), ceftobiprole (BPR),
  ceftobiprole medocaril (CFM1), ceftolozane/tazobactam (CZT),
  ceftriaxone (CRO), ceftriaxone/beta-lactamase inhibitor (CEB),
  cefuroxime (CXM), cefuroxime axetil (CXA), cephradine (CED),
  ciclacillin (CIC), clometocillin (CLM), cloxacillin (CLO),
  dicloxacillin (DIC), doripenem (DOR), epicillin (EPC), ertapenem
  (ETP), flucloxacillin (FLC), hetacillin (HET), imipenem (IPM),
  imipenem/EDTA (IPE), imipenem/relebactam (IMR), latamoxef (LTM),
  lenampicillin (LEN), loracarbef (LOR), mecillinam (MEC), meropenem
  (MEM), meropenem/nacubactam (MNC), meropenem/vaborbactam (MEV),
  metampicillin (MTM), meticillin (MET), mezlocillin (MEZ),
  mezlocillin/sulbactam (MSU), nafcillin (NAF), oxacillin (OXA),
  oxacillin screening test (OXA-S), panipenem (PAN), penamecillin (PNM),
  penicillin/novobiocin (PNO), penicillin/sulbactam (PSU), pheneticillin
  (PHE), phenoxymethylpenicillin (PHN), piperacillin (PIP),
  piperacillin/sulbactam (PIS), piperacillin/tazobactam (TZP),
  piridicillin (PRC), pivampicillin (PVM), pivmecillinam (PME), procaine
  benzylpenicillin (PRB), propicillin (PRP), razupenem (RZM), ritipenem
  (RIT), ritipenem acoxil (RIA), sarmoxicillin (SRX), sulbenicillin
  (SBC), sultamicillin (SLT6), talampicillin (TAL), tebipenem (TBP),
  temocillin (TEM), ticarcillin (TIC), ticarcillin/clavulanic acid
  (TCC), and tigemonam (TMN)

- `betalactams_with_inhibitor()` can select:  
  amoxicillin/clavulanic acid (AMC), amoxicillin/sulbactam (AXS),
  ampicillin/sulbactam (SAM), aztreonam/avibactam (AZA),
  aztreonam/nacubactam (ANC), cefepime/amikacin (CFA),
  cefepime/clavulanic acid (CPC), cefepime/enmetazobactam (FPE),
  cefepime/nacubactam (FNC), cefepime/tazobactam (FPT),
  cefepime/zidebactam (FPZ), cefoperazone/sulbactam (CSL),
  cefotaxime/clavulanic acid (CTC), cefotaxime/sulbactam (CTS),
  cefpodoxime/clavulanic acid (CDC), ceftaroline/avibactam (CPA),
  ceftazidime/avibactam (CZA), ceftazidime/clavulanic acid (CCV),
  ceftolozane/tazobactam (CZT), ceftriaxone/beta-lactamase inhibitor
  (CEB), imipenem/relebactam (IMR), meropenem/nacubactam (MNC),
  meropenem/vaborbactam (MEV), mezlocillin/sulbactam (MSU),
  penicillin/novobiocin (PNO), penicillin/sulbactam (PSU),
  piperacillin/sulbactam (PIS), piperacillin/tazobactam (TZP), and
  ticarcillin/clavulanic acid (TCC)

- `carbapenems()` can select:  
  biapenem (BIA), doripenem (DOR), ertapenem (ETP), imipenem (IPM),
  imipenem/EDTA (IPE), imipenem/relebactam (IMR), meropenem (MEM),
  meropenem/nacubactam (MNC), meropenem/vaborbactam (MEV), panipenem
  (PAN), razupenem (RZM), ritipenem (RIT), ritipenem acoxil (RIA), and
  tebipenem (TBP)

- `cephalosporins()` can select:  
  cefacetrile (CAC), cefaclor (CEC), cefadroxil (CFR), cefalexin (LEX),
  cefaloridine (RID), cefalotin (CEP), cefamandole (MAN), cefapirin
  (HAP), cefatrizine (CTZ), cefazedone (CZD), cefazolin (CZO), cefcapene
  (CCP), cefcapene pivoxil (CCX), cefdinir (CDR), cefditoren (DIT),
  cefditoren pivoxil (DIX), cefepime (FEP), cefepime/amikacin (CFA),
  cefepime/clavulanic acid (CPC), cefepime/enmetazobactam (FPE),
  cefepime/nacubactam (FNC), cefepime/tazobactam (FPT),
  cefepime/zidebactam (FPZ), cefetamet (CAT), cefetamet pivoxil (CPI),
  cefetecol (CCL), cefetrizole (CZL), cefiderocol (FDC), cefixime (CFM),
  cefmenoxime (CMX), cefmetazole (CMZ), cefodizime (DIZ), cefonicid
  (CID), cefoperazone (CFP), cefoperazone/sulbactam (CSL), ceforanide
  (CND), cefoselis (CSE), cefotaxime (CTX), cefotaxime screening test
  (CTX-S), cefotaxime/clavulanic acid (CTC), cefotaxime/sulbactam (CTS),
  cefotetan (CTT), cefotiam (CTF), cefotiam hexetil (CHE), cefovecin
  (FOV), cefoxitin (FOX), cefoxitin screening test (FOX-S), cefozopran
  (ZOP), cefpimizole (CFZ), cefpiramide (CPM), cefpirome (CPO),
  cefpodoxime (CPD), cefpodoxime proxetil (CPX), cefpodoxime/clavulanic
  acid (CDC), cefprozil (CPR), cefquinome (CEQ), cefroxadine (CRD),
  cefsulodin (CFS), cefsumide (CSU), ceftaroline (CPT),
  ceftaroline/avibactam (CPA), ceftazidime (CAZ), ceftazidime/avibactam
  (CZA), ceftazidime/clavulanic acid (CCV), cefteram (CEM), cefteram
  pivoxil (CPL), ceftezole (CTL), ceftibuten (CTB), ceftiofur (TIO),
  ceftizoxime (CZX), ceftizoxime alapivoxil (CZP), ceftobiprole (BPR),
  ceftobiprole medocaril (CFM1), ceftolozane/tazobactam (CZT),
  ceftriaxone (CRO), ceftriaxone/beta-lactamase inhibitor (CEB),
  cefuroxime (CXM), cefuroxime axetil (CXA), cephradine (CED), latamoxef
  (LTM), and loracarbef (LOR)

- `cephalosporins_1st()` can select:  
  cefacetrile (CAC), cefadroxil (CFR), cefalexin (LEX), cefaloridine
  (RID), cefalotin (CEP), cefapirin (HAP), cefatrizine (CTZ), cefazedone
  (CZD), cefazolin (CZO), cefroxadine (CRD), ceftezole (CTL), and
  cephradine (CED)

- `cephalosporins_2nd()` can select:  
  cefaclor (CEC), cefamandole (MAN), cefmetazole (CMZ), cefonicid (CID),
  ceforanide (CND), cefotetan (CTT), cefotiam (CTF), cefoxitin (FOX),
  cefoxitin screening test (FOX-S), cefprozil (CPR), cefuroxime (CXM),
  cefuroxime axetil (CXA), and loracarbef (LOR)

- `cephalosporins_3rd()` can select:  
  cefcapene (CCP), cefcapene pivoxil (CCX), cefdinir (CDR), cefditoren
  (DIT), cefditoren pivoxil (DIX), cefetamet (CAT), cefetamet pivoxil
  (CPI), cefixime (CFM), cefmenoxime (CMX), cefodizime (DIZ),
  cefoperazone (CFP), cefoperazone/sulbactam (CSL), cefotaxime (CTX),
  cefotaxime screening test (CTX-S), cefotaxime/clavulanic acid (CTC),
  cefotaxime/sulbactam (CTS), cefotiam hexetil (CHE), cefovecin (FOV),
  cefpimizole (CFZ), cefpiramide (CPM), cefpodoxime (CPD), cefpodoxime
  proxetil (CPX), cefpodoxime/clavulanic acid (CDC), cefsulodin (CFS),
  ceftazidime (CAZ), ceftazidime/avibactam (CZA), ceftazidime/clavulanic
  acid (CCV), cefteram (CEM), cefteram pivoxil (CPL), ceftibuten (CTB),
  ceftiofur (TIO), ceftizoxime (CZX), ceftizoxime alapivoxil (CZP),
  ceftriaxone (CRO), ceftriaxone/beta-lactamase inhibitor (CEB), and
  latamoxef (LTM)

- `cephalosporins_4th()` can select:  
  cefepime (FEP), cefepime/amikacin (CFA), cefepime/clavulanic acid
  (CPC), cefepime/enmetazobactam (FPE), cefepime/nacubactam (FNC),
  cefepime/tazobactam (FPT), cefepime/zidebactam (FPZ), cefetecol (CCL),
  cefoselis (CSE), cefozopran (ZOP), cefpirome (CPO), and cefquinome
  (CEQ)

- `cephalosporins_5th()` can select:  
  ceftaroline (CPT), ceftaroline/avibactam (CPA), ceftobiprole (BPR),
  ceftobiprole medocaril (CFM1), and ceftolozane/tazobactam (CZT)

- `fluoroquinolones()` can select:  
  besifloxacin (BES), ciprofloxacin (CIP), ciprofloxacin/metronidazole
  (CIM), ciprofloxacin/ornidazole (CIO), ciprofloxacin/tinidazole (CIT),
  clinafloxacin (CLX), danofloxacin (DAN), delafloxacin (DFX),
  difloxacin (DIF), enoxacin (ENX), enrofloxacin (ENR), finafloxacin
  (FIN), fleroxacin (FLE), garenoxacin (GRN), gatifloxacin (GAT),
  gemifloxacin (GEM), grepafloxacin (GRX), lascufloxacin (LSC),
  levofloxacin (LVX), levofloxacin/ornidazole (LEO), levonadifloxacin
  (LND), lomefloxacin (LOM), marbofloxacin (MAR), metioxate (MXT),
  miloxacin (MIL), moxifloxacin (MFX), nadifloxacin (NAD), nemonoxacin
  (NEM), nifuroquine (NIF), nitroxoline (NTR), norfloxacin (NOR),
  norfloxacin screening test (NOR-S), norfloxacin/metronidazole (NME),
  norfloxacin/tinidazole (NTI), ofloxacin (OFX), ofloxacin/ornidazole
  (OOR), orbifloxacin (ORB), pazufloxacin (PAZ), pefloxacin (PEF),
  pefloxacin screening test (PEF-S), pradofloxacin (PRA), premafloxacin
  (PRX), prulifloxacin (PRU), rufloxacin (RFL), sarafloxacin (SAR),
  sitafloxacin (SIT), sparfloxacin (SPX), temafloxacin (TMX),
  tilbroquinol (TBQ), tioxacin (TXC), tosufloxacin (TFX), and
  trovafloxacin (TVA)

- `glycopeptides()` can select:  
  avoparcin (AVO), bleomycin (BLM), dalbavancin (DAL), norvancomycin
  (NVA), oritavancin (ORI), ramoplanin (RAM), teicoplanin (TEC),
  teicoplanin-macromethod (TCM), telavancin (TLV), vancomycin (VAN), and
  vancomycin-macromethod (VAM)

- `isoxazolylpenicillins()` can select:  
  cloxacillin (CLO), dicloxacillin (DIC), flucloxacillin (FLC),
  meticillin (MET), oxacillin (OXA), and oxacillin screening test
  (OXA-S)

- `lincosamides()` can select:  
  clindamycin (CLI), lincomycin (LIN), and pirlimycin (PRL)

- `lipoglycopeptides()` can select:  
  dalbavancin (DAL), oritavancin (ORI), and telavancin (TLV)

- `macrolides()` can select:  
  acetylmidecamycin (ACM), acetylspiramycin (ASP), azithromycin (AZM),
  clarithromycin (CLR), dirithromycin (DIR), erythromycin (ERY),
  flurithromycin (FLR1), gamithromycin (GAM), josamycin (JOS),
  kitasamycin (KIT), meleumycin (MEL), midecamycin (MID), miocamycin
  (MCM), nafithromycin (ZWK), oleandomycin (OLE), rokitamycin (ROK),
  roxithromycin (RXT), solithromycin (SOL), spiramycin (SPI),
  telithromycin (TLT), tildipirosin (TIP), tilmicosin (TIL),
  troleandomycin (TRL), tulathromycin (TUL), tylosin (TYL), and
  tylvalosin (TYL1)

- `monobactams()` can select:  
  aztreonam (ATM), aztreonam/avibactam (AZA), aztreonam/nacubactam
  (ANC), carumonam (CAR), and tigemonam (TMN)

- `nitrofurans()` can select:  
  furazidin (FUR), furazolidone (FRZ), nifurtoinol (NFR), nitrofurantoin
  (NIT), and nitrofurazone (NIZ)

- `oxazolidinones()` can select:  
  cadazolid (CDZ), cycloserine (CYC), linezolid (LNZ), tedizolid (TZD),
  and thiacetazone (THA)

- `penicillins()` can select:  
  amoxicillin (AMX), amoxicillin/clavulanic acid (AMC),
  amoxicillin/sulbactam (AXS), ampicillin (AMP), ampicillin/sulbactam
  (SAM), apalcillin (APL), aspoxicillin (APX), azidocillin (AZD),
  azlocillin (AZL), bacampicillin (BAM), benzathine benzylpenicillin
  (BNB), benzathine phenoxymethylpenicillin (BNP), benzylpenicillin
  (PEN), benzylpenicillin screening test (PEN-S), carbenicillin (CRB),
  carindacillin (CRN), ciclacillin (CIC), clometocillin (CLM),
  cloxacillin (CLO), dicloxacillin (DIC), epicillin (EPC),
  flucloxacillin (FLC), hetacillin (HET), lenampicillin (LEN),
  mecillinam (MEC), metampicillin (MTM), meticillin (MET), mezlocillin
  (MEZ), mezlocillin/sulbactam (MSU), nafcillin (NAF), oxacillin (OXA),
  oxacillin screening test (OXA-S), penamecillin (PNM),
  penicillin/novobiocin (PNO), penicillin/sulbactam (PSU), pheneticillin
  (PHE), phenoxymethylpenicillin (PHN), piperacillin (PIP),
  piperacillin/sulbactam (PIS), piperacillin/tazobactam (TZP),
  piridicillin (PRC), pivampicillin (PVM), pivmecillinam (PME), procaine
  benzylpenicillin (PRB), propicillin (PRP), sarmoxicillin (SRX),
  sulbenicillin (SBC), sultamicillin (SLT6), talampicillin (TAL),
  temocillin (TEM), ticarcillin (TIC), and ticarcillin/clavulanic acid
  (TCC)

- `phenicols()` can select:  
  chloramphenicol (CHL), florfenicol (FLR), and thiamphenicol (THI)

- `polymyxins()` can select:  
  colistin (COL), polymyxin B (PLB), and polymyxin B/polysorbate 80
  (POP)

- `quinolones()` can select:  
  besifloxacin (BES), cinoxacin (CIN), ciprofloxacin (CIP),
  ciprofloxacin/metronidazole (CIM), ciprofloxacin/ornidazole (CIO),
  ciprofloxacin/tinidazole (CIT), clinafloxacin (CLX), danofloxacin
  (DAN), delafloxacin (DFX), difloxacin (DIF), enoxacin (ENX),
  enrofloxacin (ENR), finafloxacin (FIN), fleroxacin (FLE), flumequine
  (FLM), garenoxacin (GRN), gatifloxacin (GAT), gemifloxacin (GEM),
  grepafloxacin (GRX), lascufloxacin (LSC), levofloxacin (LVX),
  levofloxacin/ornidazole (LEO), levonadifloxacin (LND), lomefloxacin
  (LOM), marbofloxacin (MAR), metioxate (MXT), miloxacin (MIL),
  moxifloxacin (MFX), nadifloxacin (NAD), nalidixic acid (NAL),
  nalidixic acid screening test (NAL-S), nemonoxacin (NEM), nifuroquine
  (NIF), nitroxoline (NTR), norfloxacin (NOR), norfloxacin screening
  test (NOR-S), norfloxacin/metronidazole (NME), norfloxacin/tinidazole
  (NTI), ofloxacin (OFX), ofloxacin/ornidazole (OOR), orbifloxacin
  (ORB), oxolinic acid (OXO), pazufloxacin (PAZ), pefloxacin (PEF),
  pefloxacin screening test (PEF-S), pipemidic acid (PPA), piromidic
  acid (PIR), pradofloxacin (PRA), premafloxacin (PRX), prulifloxacin
  (PRU), rosoxacin (ROS), rufloxacin (RFL), sarafloxacin (SAR),
  sitafloxacin (SIT), sparfloxacin (SPX), temafloxacin (TMX),
  tilbroquinol (TBQ), tioxacin (TXC), tosufloxacin (TFX), and
  trovafloxacin (TVA)

- `rifamycins()` can select:  
  rifabutin (RIB), rifampicin (RIF), rifampicin/ethambutol/isoniazid
  (REI), rifampicin/isoniazid (RFI),
  rifampicin/pyrazinamide/ethambutol/isoniazid (RPEI),
  rifampicin/pyrazinamide/isoniazid (RPI), rifamycin (RFM), and
  rifapentine (RFP)

- `streptogramins()` can select:  
  pristinamycin (PRI) and quinupristin/dalfopristin (QDA)

- `sulfonamides()` can select:  
  brodimoprim (BDP), sulfadiazine (SDI), sulfadiazine/tetroxoprim (SLT),
  sulfadimethoxine (SUD), sulfadimidine (SDM), sulfafurazole (SLF),
  sulfaisodimidine (SLF1), sulfalene (SLF2), sulfamazone (SZO),
  sulfamerazine (SLF3), sulfamethizole (SLF4), sulfamethoxazole (SMX),
  sulfamethoxypyridazine (SLF5), sulfametomidine (SLF6),
  sulfametoxydiazine (SLF7), sulfamoxole (SLF8), sulfanilamide (SLF9),
  sulfaperin (SLF10), sulfaphenazole (SLF11), sulfapyridine (SLF12),
  sulfathiazole (SUT), and sulfathiourea (SLF13)

- `tetracyclines()` can select:  
  cetocycline (CTO), chlortetracycline (CTE), clomocycline (CLM1),
  demeclocycline (DEM), doxycycline (DOX), eravacycline (ERV),
  lymecycline (LYM), metacycline (MTC), minocycline (MNO), omadacycline
  (OMC), oxytetracycline (OXY), penimepicycline (PNM1), rolitetracycline
  (RLT), sarecycline (SRC), tetracycline (TCY), tetracycline screening
  test (TCY-S), and tigecycline (TGC)

- `trimethoprims()` can select:  
  brodimoprim (BDP), sulfadiazine (SDI), sulfadiazine/tetroxoprim (SLT),
  sulfadiazine/trimethoprim (SLT1), sulfadimethoxine (SUD),
  sulfadimidine (SDM), sulfadimidine/trimethoprim (SLT2), sulfafurazole
  (SLF), sulfaisodimidine (SLF1), sulfalene (SLF2), sulfamazone (SZO),
  sulfamerazine (SLF3), sulfamerazine/trimethoprim (SLT3),
  sulfamethizole (SLF4), sulfamethoxazole (SMX), sulfamethoxypyridazine
  (SLF5), sulfametomidine (SLF6), sulfametoxydiazine (SLF7),
  sulfametrole/trimethoprim (SLT4), sulfamoxole (SLF8),
  sulfamoxole/trimethoprim (SLT5), sulfanilamide (SLF9), sulfaperin
  (SLF10), sulfaphenazole (SLF11), sulfapyridine (SLF12), sulfathiazole
  (SUT), sulfathiourea (SLF13), trimethoprim (TMP), and
  trimethoprim/sulfamethoxazole (SXT)

- `ureidopenicillins()` can select:  
  azlocillin (AZL), mezlocillin (MEZ), piperacillin (PIP), and
  piperacillin/tazobactam (TZP)

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

## Examples

``` r
# `example_isolates` is a data set available in the AMR package.
# See ?example_isolates.
example_isolates


# you can use the selectors separately to retrieve all possible antimicrobials:
carbapenems()


# Though they are primarily intended to use for selections and filters.
# Examples sections below are split into 'dplyr', 'base R', and 'data.table':

# \donttest{
if (FALSE) { # \dontrun{
# dplyr -------------------------------------------------------------------

library(dplyr, warn.conflicts = FALSE)

example_isolates %>% select(carbapenems())

# select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB'
example_isolates %>% select(mo, aminoglycosides())

# you can combine selectors like you are used with tidyverse
# e.g., for betalactams, but not the ones with an enzyme inhibitor:
example_isolates %>% select(betalactams(), -betalactams_with_inhibitor())

# select only antimicrobials with DDDs for oral treatment
example_isolates %>% select(administrable_per_os())

# get AMR for all aminoglycosides e.g., per ward:
example_isolates %>%
  group_by(ward) %>%
  summarise(across(aminoglycosides(),
                   resistance))

# You can combine selectors with '&' to be more specific:
example_isolates %>%
  select(penicillins() & administrable_per_os())

# get AMR for only drugs that matter - no intrinsic resistance:
example_isolates %>%
  filter(mo_genus() %in% c("Escherichia", "Klebsiella")) %>%
  group_by(ward) %>%
  summarise_at(not_intrinsic_resistant(),
               resistance)

# get susceptibility for antimicrobials whose name contains "trim":
example_isolates %>%
  filter(first_isolate()) %>%
  group_by(ward) %>%
  summarise(across(amr_selector(name %like% "trim"), susceptibility))

# this will select columns 'IPM' (imipenem) and 'MEM' (meropenem):
example_isolates %>%
  select(carbapenems())

# this will select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB':
example_isolates %>%
  select(mo, aminoglycosides())

# any() and all() work in dplyr's filter() too:
example_isolates %>%
  filter(
    any(aminoglycosides() == "R"),
    all(cephalosporins_2nd() == "R")
  )

# also works with c():
example_isolates %>%
  filter(any(c(carbapenems(), aminoglycosides()) == "R"))

# not setting any/all will automatically apply all():
example_isolates %>%
  filter(aminoglycosides() == "R")

# this will select columns 'mo' and all antimycobacterial drugs ('RIF'):
example_isolates %>%
  select(mo, amr_class("mycobact"))

# get bug/drug combinations for only glycopeptides in Gram-positives:
example_isolates %>%
  filter(mo_is_gram_positive()) %>%
  select(mo, glycopeptides()) %>%
  bug_drug_combinations() %>%
  format()

data.frame(
  some_column = "some_value",
  J01CA01 = "S"
) %>% # ATC code of ampicillin
  select(penicillins()) # only the 'J01CA01' column will be selected

# with recent versions of dplyr, this is all equal:
x <- example_isolates[carbapenems() == "R", ]
y <- example_isolates %>% filter(carbapenems() == "R")
z <- example_isolates %>% filter(if_all(carbapenems(), ~ .x == "R"))
identical(x, y) && identical(y, z)

} # }
# base R ------------------------------------------------------------------

# select columns 'IPM' (imipenem) and 'MEM' (meropenem)
example_isolates[, carbapenems()]

# select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB'
example_isolates[, c("mo", aminoglycosides())]

# select only antimicrobials with DDDs for oral treatment
example_isolates[, administrable_per_os()]

# filter using any() or all()
example_isolates[any(carbapenems() == "R"), ]
subset(example_isolates, any(carbapenems() == "R"))

# filter on any or all results in the carbapenem columns (i.e., IPM, MEM):
example_isolates[any(carbapenems()), ]
example_isolates[all(carbapenems()), ]

# filter with multiple antimicrobial selectors using c()
example_isolates[all(c(carbapenems(), aminoglycosides()) == "R"), ]

# filter + select in one go: get penicillins in carbapenem-resistant strains
example_isolates[any(carbapenems() == "R"), penicillins()]

# You can combine selectors with '&' to be more specific. For example,
# penicillins() would select benzylpenicillin ('peni G') and
# administrable_per_os() would select erythromycin. Yet, when combined these
# drugs are both omitted since benzylpenicillin is not administrable per os
# and erythromycin is not a penicillin:
example_isolates[, penicillins() & administrable_per_os()]

# amr_selector() applies a filter in the `antimicrobials` data set and is thus
# very flexible. For instance, to select antimicrobials with an oral DDD
# of at least 1 gram:
example_isolates[, amr_selector(oral_ddd > 1 & oral_units == "g")]


# data.table --------------------------------------------------------------

# data.table is supported as well, just use it in the same way as with
# base R, but add `with = FALSE` if using a single AB selector.

if (require("data.table")) {
  dt <- as.data.table(example_isolates)

  # this does not work, it returns column *names*
  dt[, carbapenems()]
}
if (require("data.table")) {
  # so `with = FALSE` is required
  dt[, carbapenems(), with = FALSE]
}

# for multiple selections or AB selectors, `with = FALSE` is not needed:
if (require("data.table")) {
  dt[, c("mo", aminoglycosides())]
}
if (require("data.table")) {
  dt[, c(carbapenems(), aminoglycosides())]
}

# row filters are also supported:
if (require("data.table")) {
  dt[any(carbapenems() == "S"), ]
}
if (require("data.table")) {
  dt[any(carbapenems() == "S"), penicillins(), with = FALSE]
}
# }
```
