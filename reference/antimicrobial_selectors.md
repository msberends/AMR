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
  class [sir](https://amr-for-r.org/reference/as.sir.md) on beforehand.
  Defaults to `FALSE`.

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
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set will be searched (case-insensitive) for this value.

- filter:

  An [expression](https://rdrr.io/r/base/expression.html) to be
  evaluated in the
  [antimicrobials](https://amr-for-r.org/reference/antimicrobials.md)
  data set, such as `name %like% "trim"`.

- col_mo:

  Column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/reference/as.mo.md)) - the default
  is the first column of class
  [`mo`](https://amr-for-r.org/reference/as.mo.md). Values will be
  coerced using [`as.mo()`](https://amr-for-r.org/reference/as.mo.md).

- version_expected_phenotypes:

  The version number to use for the EUCAST Expected Phenotypes. Can be
  "1.2".

## Value

When used inside selecting or filtering, this returns a
[character](https://rdrr.io/r/base/character.html) vector of column
names, with additional class `"amr_selector"`. When used individually,
this returns an ['ab' vector](https://amr-for-r.org/reference/as.ab.md)
with all possible antimicrobials that the function would be able to
select or filter.

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
[antimicrobials](https://amr-for-r.org/reference/antimicrobials.md) data
set. This means that a selector such as `aminoglycosides()` will pick up
column names like 'gen', 'genta', 'J01GB03', 'tobra', 'Tobracin', etc.

The `amr_class()` function can be used to filter/select on a manually
defined antimicrobial class. It searches for results in the
[antimicrobials](https://amr-for-r.org/reference/antimicrobials.md) data
set within the columns `group`, `atc_group1` and `atc_group2`.

The `administrable_per_os()` and `administrable_iv()` functions also
rely on the
[antimicrobials](https://amr-for-r.org/reference/antimicrobials.md) data
set - antimicrobials will be matched where a DDD (defined daily dose)
for resp. oral and IV treatment is available in the
[antimicrobials](https://amr-for-r.org/reference/antimicrobials.md) data
set.

The `amr_selector()` function can be used to internally filter the
[antimicrobials](https://amr-for-r.org/reference/antimicrobials.md) data
set on any results, see *Examples*. It allows for filtering on a (part
of) a certain name, and/or a group name or even a minimum of DDDs for
oral treatment. This function yields the highest flexibility, but is
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
[`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md)
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
  cefepime/nacubactam (FNC), cefepime/taniborbactam (FTA),
  cefepime/tazobactam (FPT), cefepime/zidebactam (FPZ), cefetamet (CAT),
  cefetamet pivoxil (CPI), cefetecol (CCL), cefetrizole (CZL),
  cefiderocol (FDC), cefixime (CFM), cefmenoxime (CMX), cefmetazole
  (CMZ), cefodizime (DIZ), cefonicid (CID), cefoperazone (CFP),
  cefoperazone/sulbactam (CSL), ceforanide (CND), cefoselis (CSE),
  cefotaxime (CTX), cefotaxime screening test (CTX-S),
  cefotaxime/clavulanic acid (CTC), cefotaxime/sulbactam (CTS),
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
  (SBC), sultamicillin (SLT6), talampicillin (TAL), taniborbactam (TAN),
  tebipenem (TBP), temocillin (TEM), ticarcillin (TIC),
  ticarcillin/clavulanic acid (TCC), and tigemonam (TMN)

- `betalactams_with_inhibitor()` can select:  
  amoxicillin/clavulanic acid (AMC), amoxicillin/sulbactam (AXS),
  ampicillin/sulbactam (SAM), aztreonam/avibactam (AZA),
  aztreonam/nacubactam (ANC), cefepime/amikacin (CFA),
  cefepime/clavulanic acid (CPC), cefepime/enmetazobactam (FPE),
  cefepime/nacubactam (FNC), cefepime/taniborbactam (FTA),
  cefepime/tazobactam (FPT), cefepime/zidebactam (FPZ),
  cefoperazone/sulbactam (CSL), cefotaxime/clavulanic acid (CTC),
  cefotaxime/sulbactam (CTS), cefpodoxime/clavulanic acid (CDC),
  ceftaroline/avibactam (CPA), ceftazidime/avibactam (CZA),
  ceftazidime/clavulanic acid (CCV), ceftolozane/tazobactam (CZT),
  ceftriaxone/beta-lactamase inhibitor (CEB), imipenem/relebactam (IMR),
  meropenem/nacubactam (MNC), meropenem/vaborbactam (MEV),
  mezlocillin/sulbactam (MSU), penicillin/novobiocin (PNO),
  penicillin/sulbactam (PSU), piperacillin/sulbactam (PIS),
  piperacillin/tazobactam (TZP), and ticarcillin/clavulanic acid (TCC)

- `carbapenems()` can select:  
  biapenem (BIA), doripenem (DOR), ertapenem (ETP), imipenem (IPM),
  imipenem/EDTA (IPE), imipenem/relebactam (IMR), meropenem (MEM),
  meropenem/nacubactam (MNC), meropenem/vaborbactam (MEV), panipenem
  (PAN), razupenem (RZM), ritipenem (RIT), ritipenem acoxil (RIA),
  taniborbactam (TAN), and tebipenem (TBP)

- `cephalosporins()` can select:  
  cefacetrile (CAC), cefaclor (CEC), cefadroxil (CFR), cefalexin (LEX),
  cefaloridine (RID), cefalotin (CEP), cefamandole (MAN), cefapirin
  (HAP), cefatrizine (CTZ), cefazedone (CZD), cefazolin (CZO), cefcapene
  (CCP), cefcapene pivoxil (CCX), cefdinir (CDR), cefditoren (DIT),
  cefditoren pivoxil (DIX), cefepime (FEP), cefepime/amikacin (CFA),
  cefepime/clavulanic acid (CPC), cefepime/enmetazobactam (FPE),
  cefepime/nacubactam (FNC), cefepime/taniborbactam (FTA),
  cefepime/tazobactam (FPT), cefepime/zidebactam (FPZ), cefetamet (CAT),
  cefetamet pivoxil (CPI), cefetecol (CCL), cefetrizole (CZL),
  cefiderocol (FDC), cefixime (CFM), cefmenoxime (CMX), cefmetazole
  (CMZ), cefodizime (DIZ), cefonicid (CID), cefoperazone (CFP),
  cefoperazone/sulbactam (CSL), ceforanide (CND), cefoselis (CSE),
  cefotaxime (CTX), cefotaxime screening test (CTX-S),
  cefotaxime/clavulanic acid (CTC), cefotaxime/sulbactam (CTS),
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
  cefepime/taniborbactam (FTA), cefepime/tazobactam (FPT),
  cefepime/zidebactam (FPZ), cefetecol (CCL), cefoselis (CSE),
  cefozopran (ZOP), cefpirome (CPO), and cefquinome (CEQ)

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
#> # A tibble: 2,000 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-02 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-03 A77334     65 F      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  3 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  4 2002-01-07 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  5 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  6 2002-01-13 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#>  7 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  8 2002-01-14 462729     78 M      Clinical B_STPHY_AURS   R     NA    S     R  
#>  9 2002-01-16 067927     45 F      ICU      B_STPHY_EPDR   R     NA    R     NA 
#> 10 2002-01-17 858515     79 F      ICU      B_STPHY_EPDR   R     NA    S     NA 
#> # ℹ 1,990 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …


# you can use the selectors separately to retrieve all possible antimicrobials:
carbapenems()
#> ℹ in `carbapenems()`: Imipenem/EDTA (`IPE`) and meropenem/nacubactam
#>   (`MNC`) are not included since `only_treatable = TRUE`.
#> ℹ This 'ab' vector was retrieved using `carbapenems()`, which should
#>   normally be used inside a `dplyr` verb or `data.frame` call, e.g.:
#>   • your_data %>% select(carbapenems())
#>   • your_data %>% select(column_a, column_b, carbapenems())
#>   • your_data %>% filter(any(carbapenems() == "R"))
#>   • your_data[, carbapenems()]
#>   • your_data[, c("column_a", "column_b", carbapenems())]
#> Class 'ab'
#>  [1] BIA DOR ETP IMR IPM MEM MEV PAN RIA RIT RZM TAN TBP


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
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> # A tibble: 2,000 × 2
#>    IPM   MEM  
#>    <sir> <sir>
#>  1   NA    NA 
#>  2   NA    NA 
#>  3   NA    NA 
#>  4   NA    NA 
#>  5   NA    NA 
#>  6   NA    NA 
#>  7   NA    NA 
#>  8   NA    NA 
#>  9   NA    NA 
#> 10   NA    NA 
#> # ℹ 1,990 more rows

# select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB'
example_isolates[, c("mo", aminoglycosides())]
#> ℹ For `aminoglycosides()` using columns 'GEN' (gentamicin), 'TOB'
#>   (tobramycin), 'AMK' (amikacin), and 'KAN' (kanamycin)
#> # A tibble: 2,000 × 5
#>    mo           GEN   TOB   AMK   KAN  
#>    <mo>         <sir> <sir> <sir> <sir>
#>  1 B_ESCHR_COLI   NA    NA    NA    NA 
#>  2 B_ESCHR_COLI   NA    NA    NA    NA 
#>  3 B_STPHY_EPDR   NA    NA    NA    NA 
#>  4 B_STPHY_EPDR   NA    NA    NA    NA 
#>  5 B_STPHY_EPDR   NA    NA    NA    NA 
#>  6 B_STPHY_EPDR   NA    NA    NA    NA 
#>  7 B_STPHY_AURS   NA    S     NA    NA 
#>  8 B_STPHY_AURS   NA    S     NA    NA 
#>  9 B_STPHY_EPDR   NA    NA    NA    NA 
#> 10 B_STPHY_EPDR   NA    NA    NA    NA 
#> # ℹ 1,990 more rows

# select only antimicrobials with DDDs for oral treatment
example_isolates[, administrable_per_os()]
#> ℹ For `administrable_per_os()` using columns 'OXA' (oxacillin), 'FLC'
#>   (flucloxacillin), 'AMX' (amoxicillin), 'AMC' (amoxicillin/clavulanic acid),
#>   'AMP' (ampicillin), 'CXM' (cefuroxime), 'KAN' (kanamycin), 'TMP'
#>   (trimethoprim), 'NIT' (nitrofurantoin), 'FOS' (fosfomycin), 'LNZ'
#>   (linezolid), 'CIP' (ciprofloxacin), 'MFX' (moxifloxacin), 'VAN'
#>   (vancomycin), 'TCY' (tetracycline), 'DOX' (doxycycline), 'ERY'
#>   (erythromycin), 'CLI' (clindamycin), 'AZM' (azithromycin), 'MTR'
#>   (metronidazole), 'CHL' (chloramphenicol), 'COL' (colistin), and 'RIF'
#>   (rifampicin)
#> # A tibble: 2,000 × 23
#>    OXA   FLC   AMX   AMC   AMP   CXM   KAN   TMP   NIT   FOS   LNZ   CIP   MFX  
#>    <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir>
#>  1   NA    NA    NA    I     NA    I     NA    R     NA    NA    R     NA    NA 
#>  2   NA    NA    NA    I     NA    I     NA    R     NA    NA    R     NA    NA 
#>  3   NA    R     NA    NA    NA    R     NA    S     NA    NA    NA    NA    NA 
#>  4   NA    R     NA    NA    NA    R     NA    S     NA    NA    NA    NA    NA 
#>  5   NA    R     NA    NA    NA    R     NA    R     NA    NA    NA    NA    NA 
#>  6   NA    R     NA    NA    NA    R     NA    R     NA    NA    NA    NA    NA 
#>  7   NA    S     R     S     R     S     NA    R     NA    NA    NA    NA    NA 
#>  8   NA    S     R     S     R     S     NA    R     NA    NA    NA    NA    NA 
#>  9   NA    R     NA    NA    NA    R     NA    S     NA    NA    NA    S     NA 
#> 10   NA    S     NA    NA    NA    S     NA    S     NA    NA    NA    S     NA 
#> # ℹ 1,990 more rows
#> # ℹ 10 more variables: VAN <sir>, TCY <sir>, DOX <sir>, ERY <sir>, CLI <sir>,
#> #   AZM <sir>, MTR <sir>, CHL <sir>, COL <sir>, RIF <sir>

# filter using any() or all()
example_isolates[any(carbapenems() == "R"), ]
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> # A tibble: 55 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2004-06-09 529296     69 M      ICU      B_ENTRC_FACM   NA    NA    NA    NA 
#>  2 2004-06-09 529296     69 M      ICU      B_ENTRC_FACM   NA    NA    NA    NA 
#>  3 2004-11-03 D65308     80 F      ICU      B_STNTR_MLTP   R     NA    NA    R  
#>  4 2005-04-21 452212     82 F      ICU      B_ENTRC        NA    NA    NA    NA 
#>  5 2005-04-22 452212     82 F      ICU      B_ENTRC        NA    NA    NA    NA 
#>  6 2005-04-22 452212     82 F      ICU      B_ENTRC_FACM   NA    NA    NA    NA 
#>  7 2007-02-21 8BBC46     61 F      Clinical B_ENTRC_FACM   NA    NA    NA    NA 
#>  8 2007-12-15 401043     72 M      Clinical B_ENTRC_FACM   NA    NA    NA    NA 
#>  9 2008-01-22 1710B8     82 M      Clinical B_PROTS_MRBL   R     NA    NA    NA 
#> 10 2008-01-22 1710B8     82 M      Clinical B_PROTS_MRBL   R     NA    NA    NA 
#> # ℹ 45 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …
subset(example_isolates, any(carbapenems() == "R"))
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> # A tibble: 55 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2004-06-09 529296     69 M      ICU      B_ENTRC_FACM   NA    NA    NA    NA 
#>  2 2004-06-09 529296     69 M      ICU      B_ENTRC_FACM   NA    NA    NA    NA 
#>  3 2004-11-03 D65308     80 F      ICU      B_STNTR_MLTP   R     NA    NA    R  
#>  4 2005-04-21 452212     82 F      ICU      B_ENTRC        NA    NA    NA    NA 
#>  5 2005-04-22 452212     82 F      ICU      B_ENTRC        NA    NA    NA    NA 
#>  6 2005-04-22 452212     82 F      ICU      B_ENTRC_FACM   NA    NA    NA    NA 
#>  7 2007-02-21 8BBC46     61 F      Clinical B_ENTRC_FACM   NA    NA    NA    NA 
#>  8 2007-12-15 401043     72 M      Clinical B_ENTRC_FACM   NA    NA    NA    NA 
#>  9 2008-01-22 1710B8     82 M      Clinical B_PROTS_MRBL   R     NA    NA    NA 
#> 10 2008-01-22 1710B8     82 M      Clinical B_PROTS_MRBL   R     NA    NA    NA 
#> # ℹ 45 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …

# filter on any or all results in the carbapenem columns (i.e., IPM, MEM):
example_isolates[any(carbapenems()), ]
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> ℹ Filtering any of columns 'IPM' and 'MEM' to contain value "S", "I" or "R"
#> # A tibble: 962 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2002-01-19 738003     71 M      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  2 2002-01-19 738003     71 M      Clinical B_ESCHR_COLI   R     NA    NA    NA 
#>  3 2002-01-22 F35553     50 M      ICU      B_PROTS_MRBL   R     NA    NA    NA 
#>  4 2002-01-22 F35553     50 M      ICU      B_PROTS_MRBL   R     NA    NA    NA 
#>  5 2002-02-05 067927     45 F      ICU      B_SERRT_MRCS   R     NA    NA    R  
#>  6 2002-02-05 067927     45 F      ICU      B_SERRT_MRCS   R     NA    NA    R  
#>  7 2002-02-05 067927     45 F      ICU      B_SERRT_MRCS   R     NA    NA    R  
#>  8 2002-02-27 066895     85 F      Clinical B_KLBSL_PNMN   R     NA    NA    R  
#>  9 2002-02-27 066895     85 F      Clinical B_KLBSL_PNMN   R     NA    NA    R  
#> 10 2002-03-08 4FC193     69 M      Clinical B_ESCHR_COLI   R     NA    NA    R  
#> # ℹ 952 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …
example_isolates[all(carbapenems()), ]
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> ℹ Filtering all of columns 'IPM' and 'MEM' to contain value "S", "I" or "R"
#> # A tibble: 756 × 46
#>    date       patient   age gender ward    mo            PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>   <mo>          <sir> <sir> <sir> <sir>
#>  1 2002-04-14 F30196     73 M      Outpat… B_STRPT_GRPB    S     NA    S     S  
#>  2 2003-04-08 114570     74 M      ICU     B_STRPT_PYGN    S     NA    S     S  
#>  3 2003-04-08 114570     74 M      ICU     B_STRPT_GRPA    S     NA    S     S  
#>  4 2003-04-08 114570     74 M      ICU     B_STRPT_GRPA    S     NA    S     S  
#>  5 2003-08-14 F71508      0 F      Clinic… B_STRPT_GRPB    S     NA    S     S  
#>  6 2003-10-16 650870     63 F      ICU     B_ESCHR_COLI    R     NA    NA    R  
#>  7 2003-10-20 F35553     52 M      ICU     B_ENTRBC_CLOC   R     NA    NA    R  
#>  8 2003-10-20 F35553     52 M      ICU     B_ENTRBC_CLOC   R     NA    NA    R  
#>  9 2003-11-04 2FC253     87 F      ICU     B_ESCHR_COLI    R     NA    NA    NA 
#> 10 2003-11-04 2FC253     87 F      ICU     B_ESCHR_COLI    R     NA    NA    NA 
#> # ℹ 746 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …

# filter with multiple antimicrobial selectors using c()
example_isolates[all(c(carbapenems(), aminoglycosides()) == "R"), ]
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> ℹ For `aminoglycosides()` using columns 'GEN' (gentamicin), 'TOB'
#>   (tobramycin), 'AMK' (amikacin), and 'KAN' (kanamycin)
#> # A tibble: 26 × 46
#>    date       patient   age gender ward     mo           PEN   OXA   FLC   AMX  
#>    <date>     <chr>   <dbl> <chr>  <chr>    <mo>         <sir> <sir> <sir> <sir>
#>  1 2004-11-03 D65308     80 F      ICU      B_STNTR_MLTP   R     NA    NA    R  
#>  2 2005-04-22 452212     82 F      ICU      B_ENTRC_FACM   NA    NA    NA    NA 
#>  3 2007-02-21 8BBC46     61 F      Clinical B_ENTRC_FACM   NA    NA    NA    NA 
#>  4 2007-12-15 401043     72 M      Clinical B_ENTRC_FACM   NA    NA    NA    NA 
#>  5 2008-12-06 501361     43 F      Clinical B_STNTR_MLTP   R     NA    NA    R  
#>  6 2011-05-09 207325     82 F      ICU      B_ENTRC_FACM   NA    NA    NA    NA 
#>  7 2012-03-12 582258     80 M      ICU      B_STPHY_CONS   R     R     R     R  
#>  8 2012-05-19 C25552     89 F      Outpati… B_STPHY_CONS   R     R     R     R  
#>  9 2012-07-17 F05015     83 M      ICU      B_STPHY_CONS   R     R     R     R  
#> 10 2012-07-20 404299     66 F      Clinical B_STPHY_CONS   R     R     R     R  
#> # ℹ 16 more rows
#> # ℹ 36 more variables: AMC <sir>, AMP <sir>, TZP <sir>, CZO <sir>, FEP <sir>,
#> #   CXM <sir>, FOX <sir>, CTX <sir>, CAZ <sir>, CRO <sir>, GEN <sir>,
#> #   TOB <sir>, AMK <sir>, KAN <sir>, TMP <sir>, SXT <sir>, NIT <sir>,
#> #   FOS <sir>, LNZ <sir>, CIP <sir>, MFX <sir>, VAN <sir>, TEC <sir>,
#> #   TCY <sir>, TGC <sir>, DOX <sir>, ERY <sir>, CLI <sir>, AZM <sir>,
#> #   IPM <sir>, MEM <sir>, MTR <sir>, CHL <sir>, COL <sir>, MUP <sir>, …

# filter + select in one go: get penicillins in carbapenem-resistant strains
example_isolates[any(carbapenems() == "R"), penicillins()]
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> ℹ For `penicillins()` using columns 'PEN' (benzylpenicillin), 'OXA'
#>   (oxacillin), 'FLC' (flucloxacillin), 'AMX' (amoxicillin), 'AMC'
#>   (amoxicillin/clavulanic acid), 'AMP' (ampicillin), and 'TZP'
#>   (piperacillin/tazobactam)
#> # A tibble: 55 × 7
#>    PEN   OXA   FLC   AMX   AMC   AMP   TZP  
#>    <sir> <sir> <sir> <sir> <sir> <sir> <sir>
#>  1   NA    NA    NA    NA    NA    NA    NA 
#>  2   NA    NA    NA    NA    NA    NA    NA 
#>  3   R     NA    NA    R     R     R     R  
#>  4   NA    NA    NA    NA    NA    NA    R  
#>  5   NA    NA    NA    NA    NA    NA    R  
#>  6   NA    NA    NA    NA    NA    NA    R  
#>  7   NA    NA    NA    NA    NA    NA    R  
#>  8   NA    NA    NA    NA    NA    NA    R  
#>  9   R     NA    NA    NA    S     NA    S  
#> 10   R     NA    NA    NA    S     NA    S  
#> # ℹ 45 more rows

# You can combine selectors with '&' to be more specific. For example,
# penicillins() would select benzylpenicillin ('peni G') and
# administrable_per_os() would select erythromycin. Yet, when combined these
# drugs are both omitted since benzylpenicillin is not administrable per os
# and erythromycin is not a penicillin:
example_isolates[, penicillins() & administrable_per_os()]
#> ℹ For `penicillins()` using columns 'PEN' (benzylpenicillin), 'OXA'
#>   (oxacillin), 'FLC' (flucloxacillin), 'AMX' (amoxicillin), 'AMC'
#>   (amoxicillin/clavulanic acid), 'AMP' (ampicillin), and 'TZP'
#>   (piperacillin/tazobactam)
#> ℹ For `administrable_per_os()` using columns 'OXA' (oxacillin), 'FLC'
#>   (flucloxacillin), 'AMX' (amoxicillin), 'AMC' (amoxicillin/clavulanic acid),
#>   'AMP' (ampicillin), 'CXM' (cefuroxime), 'KAN' (kanamycin), 'TMP'
#>   (trimethoprim), 'NIT' (nitrofurantoin), 'FOS' (fosfomycin), 'LNZ'
#>   (linezolid), 'CIP' (ciprofloxacin), 'MFX' (moxifloxacin), 'VAN'
#>   (vancomycin), 'TCY' (tetracycline), 'DOX' (doxycycline), 'ERY'
#>   (erythromycin), 'CLI' (clindamycin), 'AZM' (azithromycin), 'MTR'
#>   (metronidazole), 'CHL' (chloramphenicol), 'COL' (colistin), and 'RIF'
#>   (rifampicin)
#> # A tibble: 2,000 × 5
#>    OXA   FLC   AMX   AMC   AMP  
#>    <sir> <sir> <sir> <sir> <sir>
#>  1   NA    NA    NA    I     NA 
#>  2   NA    NA    NA    I     NA 
#>  3   NA    R     NA    NA    NA 
#>  4   NA    R     NA    NA    NA 
#>  5   NA    R     NA    NA    NA 
#>  6   NA    R     NA    NA    NA 
#>  7   NA    S     R     S     R  
#>  8   NA    S     R     S     R  
#>  9   NA    R     NA    NA    NA 
#> 10   NA    S     NA    NA    NA 
#> # ℹ 1,990 more rows

# amr_selector() applies a filter in the `antimicrobials` data set and is thus
# very flexible. For instance, to select antimicrobials with an oral DDD
# of at least 1 gram:
example_isolates[, amr_selector(oral_ddd > 1 & oral_units == "g")]
#> ℹ For `amr_selector(oral_ddd > 1 & oral_units == "g")` using columns 'OXA'
#>   (oxacillin), 'FLC' (flucloxacillin), 'AMX' (amoxicillin), 'AMC'
#>   (amoxicillin/clavulanic acid), 'AMP' (ampicillin), 'KAN' (kanamycin), 'FOS'
#>   (fosfomycin), 'LNZ' (linezolid), 'VAN' (vancomycin), 'ERY' (erythromycin),
#>   'CLI' (clindamycin), 'MTR' (metronidazole), and 'CHL' (chloramphenicol)
#> # A tibble: 2,000 × 13
#>    OXA   FLC   AMX   AMC   AMP   KAN   FOS   LNZ   VAN   ERY   CLI   MTR   CHL  
#>    <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir>
#>  1   NA    NA    NA    I     NA    NA    NA    R     R     R     R     NA    NA 
#>  2   NA    NA    NA    I     NA    NA    NA    R     R     R     R     NA    NA 
#>  3   NA    R     NA    NA    NA    NA    NA    NA    S     R     NA    NA    NA 
#>  4   NA    R     NA    NA    NA    NA    NA    NA    S     R     NA    NA    NA 
#>  5   NA    R     NA    NA    NA    NA    NA    NA    S     R     NA    NA    NA 
#>  6   NA    R     NA    NA    NA    NA    NA    NA    S     R     R     NA    NA 
#>  7   NA    S     R     S     R     NA    NA    NA    S     S     NA    NA    NA 
#>  8   NA    S     R     S     R     NA    NA    NA    S     S     NA    NA    NA 
#>  9   NA    R     NA    NA    NA    NA    NA    NA    S     R     NA    NA    NA 
#> 10   NA    S     NA    NA    NA    NA    NA    NA    S     S     NA    NA    NA 
#> # ℹ 1,990 more rows


# data.table --------------------------------------------------------------

# data.table is supported as well, just use it in the same way as with
# base R, but add `with = FALSE` if using a single AB selector.

if (require("data.table")) {
  dt <- as.data.table(example_isolates)

  # this does not work, it returns column *names*
  dt[, carbapenems()]
}
#> Loading required package: data.table
#> 
#> Attaching package: ‘data.table’
#> The following object is masked from ‘package:purrr’:
#> 
#>     transpose
#> The following objects are masked from ‘package:dplyr’:
#> 
#>     between, first, last
#> The following objects are masked from ‘package:AMR’:
#> 
#>     %like%, like
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> Warning: It should never be needed to print an antimicrobial selector class. Are you
#> using data.table? Then add the argument `with = FALSE`, see our examples at
#> `?amr_selector`.
#> Class 'amr_selector'
#> [1] IPM MEM
if (require("data.table")) {
  # so `with = FALSE` is required
  dt[, carbapenems(), with = FALSE]
}
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#>         IPM   MEM
#>       <sir> <sir>
#>    1:  <NA>  <NA>
#>    2:  <NA>  <NA>
#>    3:  <NA>  <NA>
#>    4:  <NA>  <NA>
#>    5:  <NA>  <NA>
#>   ---            
#> 1996:  <NA>  <NA>
#> 1997:     S     S
#> 1998:     S     S
#> 1999:     S     S
#> 2000:     S     S

# for multiple selections or AB selectors, `with = FALSE` is not needed:
if (require("data.table")) {
  dt[, c("mo", aminoglycosides())]
}
#> ℹ For `aminoglycosides()` using columns 'GEN' (gentamicin), 'TOB'
#>   (tobramycin), 'AMK' (amikacin), and 'KAN' (kanamycin)
#>                 mo   GEN   TOB   AMK   KAN
#>               <mo> <sir> <sir> <sir> <sir>
#>    1: B_ESCHR_COLI  <NA>  <NA>  <NA>  <NA>
#>    2: B_ESCHR_COLI  <NA>  <NA>  <NA>  <NA>
#>    3: B_STPHY_EPDR  <NA>  <NA>  <NA>  <NA>
#>    4: B_STPHY_EPDR  <NA>  <NA>  <NA>  <NA>
#>    5: B_STPHY_EPDR  <NA>  <NA>  <NA>  <NA>
#>   ---                                     
#> 1996: B_STRPT_PNMN     R     R     R     R
#> 1997: B_ESCHR_COLI     S     S     S  <NA>
#> 1998: B_STPHY_CONS     S  <NA>  <NA>  <NA>
#> 1999: B_ESCHR_COLI     S     S  <NA>  <NA>
#> 2000: B_KLBSL_PNMN     S     S  <NA>  <NA>
if (require("data.table")) {
  dt[, c(carbapenems(), aminoglycosides())]
}
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> ℹ For `aminoglycosides()` using columns 'GEN' (gentamicin), 'TOB'
#>   (tobramycin), 'AMK' (amikacin), and 'KAN' (kanamycin)
#>         IPM   MEM   GEN   TOB   AMK   KAN
#>       <sir> <sir> <sir> <sir> <sir> <sir>
#>    1:  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>
#>    2:  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>
#>    3:  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>
#>    4:  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>
#>    5:  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>
#>   ---                                    
#> 1996:  <NA>  <NA>     R     R     R     R
#> 1997:     S     S     S     S     S  <NA>
#> 1998:     S     S     S  <NA>  <NA>  <NA>
#> 1999:     S     S     S     S  <NA>  <NA>
#> 2000:     S     S     S     S  <NA>  <NA>

# row filters are also supported:
if (require("data.table")) {
  dt[any(carbapenems() == "S"), ]
}
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#>            date patient   age gender       ward           mo   PEN   OXA   FLC
#>          <Date>  <char> <num> <char>     <char>         <mo> <sir> <sir> <sir>
#>   1: 2002-01-19  738003    71      M   Clinical B_ESCHR_COLI     R  <NA>  <NA>
#>   2: 2002-01-19  738003    71      M   Clinical B_ESCHR_COLI     R  <NA>  <NA>
#>   3: 2002-01-22  F35553    50      M        ICU B_PROTS_MRBL     R  <NA>  <NA>
#>   4: 2002-01-22  F35553    50      M        ICU B_PROTS_MRBL     R  <NA>  <NA>
#>   5: 2002-02-05  067927    45      F        ICU B_SERRT_MRCS     R  <NA>  <NA>
#>  ---                                                                          
#> 905: 2005-04-12  D71461    70      M        ICU B_ESCHR_COLI     R  <NA>  <NA>
#> 906: 2009-11-12  650870    69      F Outpatient B_ESCHR_COLI     R  <NA>  <NA>
#> 907: 2012-06-14  8CBCF2    41      F   Clinical B_STPHY_CONS     R     S     S
#> 908: 2012-10-11  175532    78      M   Clinical B_ESCHR_COLI     R  <NA>  <NA>
#> 909: 2013-11-23  A97263    77      M   Clinical B_KLBSL_PNMN     R  <NA>  <NA>
#>        AMX   AMC   AMP   TZP   CZO   FEP   CXM   FOX   CTX   CAZ   CRO   GEN
#>      <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir>
#>   1:  <NA>     I  <NA>  <NA>  <NA>  <NA>     S  <NA>     S  <NA>     S  <NA>
#>   2:  <NA>     I  <NA>  <NA>  <NA>  <NA>     S  <NA>     S  <NA>     S  <NA>
#>   3:  <NA>     I  <NA>  <NA>  <NA>  <NA>     S  <NA>     S     S     S  <NA>
#>   4:  <NA>     I  <NA>  <NA>  <NA>  <NA>     S  <NA>     S     S     S  <NA>
#>   5:     R     R     R  <NA>     R  <NA>     R     R  <NA>  <NA>  <NA>  <NA>
#>  ---                                                                        
#> 905:     S     S     S     S  <NA>     S     S     S     S     S     S     S
#> 906:     S     S     S     S     S     S     S     S     S     S     S     S
#> 907:  <NA>     S  <NA>  <NA>     S     S     S     S     S     R     S     S
#> 908:     R     S     R     S  <NA>     S     R     R     S     S     S     S
#> 909:     R     S     R     S  <NA>     S     S     S     S     S     S     S
#>        TOB   AMK   KAN   TMP   SXT   NIT   FOS   LNZ   CIP   MFX   VAN   TEC
#>      <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir>
#>   1:     S  <NA>  <NA>     S     S  <NA>  <NA>     R  <NA>  <NA>     R     R
#>   2:     S  <NA>  <NA>     S     S  <NA>  <NA>     R  <NA>  <NA>     R     R
#>   3:  <NA>  <NA>  <NA>     S     S     R  <NA>     R     S  <NA>     R     R
#>   4:  <NA>  <NA>  <NA>     S     S     R  <NA>     R     S  <NA>     R     R
#>   5:  <NA>  <NA>  <NA>     S     S     R  <NA>     R     S  <NA>     R     R
#>  ---                                                                        
#> 905:     S     S  <NA>  <NA>     S     S  <NA>     R     S  <NA>     R     R
#> 906:     S     S  <NA>     S     S     S  <NA>     R     S  <NA>     R     R
#> 907:  <NA>  <NA>  <NA>     S     S  <NA>  <NA>  <NA>     S  <NA>     S  <NA>
#> 908:     S  <NA>  <NA>     R     R     R  <NA>     R     R     R     R     R
#> 909:     S  <NA>  <NA>     S     S     S  <NA>     R     S  <NA>     R     R
#>        TCY   TGC   DOX   ERY   CLI   AZM   IPM   MEM   MTR   CHL   COL   MUP
#>      <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir> <sir>
#>   1:  <NA>  <NA>  <NA>     R     R     R     S  <NA>  <NA>  <NA>  <NA>  <NA>
#>   2:  <NA>  <NA>  <NA>     R     R     R     S  <NA>  <NA>  <NA>  <NA>  <NA>
#>   3:     R     R     R     R     R     R     S  <NA>  <NA>  <NA>     R  <NA>
#>   4:     R     R     R     R     R     R     S  <NA>  <NA>  <NA>     R  <NA>
#>   5:     R     R     R     R     R     R     S  <NA>  <NA>  <NA>     R  <NA>
#>  ---                                                                        
#> 905:  <NA>  <NA>  <NA>     R     R     R     S     S  <NA>  <NA>  <NA>  <NA>
#> 906:  <NA>  <NA>  <NA>     R     R     R     S     S  <NA>  <NA>  <NA>  <NA>
#> 907:  <NA>  <NA>  <NA>     S     S     S     S     S  <NA>  <NA>     R  <NA>
#> 908:  <NA>  <NA>  <NA>     R     R     R     S     S  <NA>  <NA>     S  <NA>
#> 909:  <NA>  <NA>  <NA>     R     R     R     S     S  <NA>  <NA>     S  <NA>
#>        RIF
#>      <sir>
#>   1:     R
#>   2:     R
#>   3:     R
#>   4:     R
#>   5:     R
#>  ---      
#> 905:     R
#> 906:     R
#> 907:  <NA>
#> 908:     R
#> 909:     R
if (require("data.table")) {
  dt[any(carbapenems() == "S"), penicillins(), with = FALSE]
}
#> ℹ For `carbapenems()` using columns 'IPM' (imipenem) and 'MEM' (meropenem)
#> ℹ For `penicillins()` using columns 'PEN' (benzylpenicillin), 'OXA'
#>   (oxacillin), 'FLC' (flucloxacillin), 'AMX' (amoxicillin), 'AMC'
#>   (amoxicillin/clavulanic acid), 'AMP' (ampicillin), and 'TZP'
#>   (piperacillin/tazobactam)
#>        PEN   OXA   FLC   AMX   AMC   AMP   TZP
#>      <sir> <sir> <sir> <sir> <sir> <sir> <sir>
#>   1:     R  <NA>  <NA>  <NA>     I  <NA>  <NA>
#>   2:     R  <NA>  <NA>  <NA>     I  <NA>  <NA>
#>   3:     R  <NA>  <NA>  <NA>     I  <NA>  <NA>
#>   4:     R  <NA>  <NA>  <NA>     I  <NA>  <NA>
#>   5:     R  <NA>  <NA>     R     R     R  <NA>
#>  ---                                          
#> 905:     R  <NA>  <NA>     S     S     S     S
#> 906:     R  <NA>  <NA>     S     S     S     S
#> 907:     R     S     S  <NA>     S  <NA>  <NA>
#> 908:     R  <NA>  <NA>     R     S     R     S
#> 909:     R  <NA>  <NA>     R     S     R     S
# }
```
