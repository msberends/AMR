# Antibiotic Selectors

These functions allow for filtering rows and selecting columns based on
antibiotic test results that are of a specific antibiotic class or group
(according to the
[antibiotics](https://amr-for-r.org/release/2.1.0/reference/antibiotics.md)
data set), without the need to define the columns or antibiotic
abbreviations.

In short, if you have a column name that resembles an antimicrobial
drug, it will be picked up by any of these functions that matches its
pharmaceutical class: "cefazolin", "kefzol", "CZO" and "J01DB04" will
all be picked up by `cephalosporins()`.

## Usage

``` r
ab_class(ab_class, only_sir_columns = FALSE, only_treatable = TRUE, ...)

ab_selector(filter, only_sir_columns = FALSE, only_treatable = TRUE, ...)

aminoglycosides(only_sir_columns = FALSE, only_treatable = TRUE, ...)

aminopenicillins(only_sir_columns = FALSE, ...)

antifungals(only_sir_columns = FALSE, ...)

antimycobacterials(only_sir_columns = FALSE, ...)

betalactams(only_sir_columns = FALSE, only_treatable = TRUE, ...)

carbapenems(only_sir_columns = FALSE, only_treatable = TRUE, ...)

cephalosporins(only_sir_columns = FALSE, ...)

cephalosporins_1st(only_sir_columns = FALSE, ...)

cephalosporins_2nd(only_sir_columns = FALSE, ...)

cephalosporins_3rd(only_sir_columns = FALSE, ...)

cephalosporins_4th(only_sir_columns = FALSE, ...)

cephalosporins_5th(only_sir_columns = FALSE, ...)

fluoroquinolones(only_sir_columns = FALSE, ...)

glycopeptides(only_sir_columns = FALSE, ...)

lincosamides(only_sir_columns = FALSE, ...)

lipoglycopeptides(only_sir_columns = FALSE, ...)

macrolides(only_sir_columns = FALSE, ...)

oxazolidinones(only_sir_columns = FALSE, ...)

penicillins(only_sir_columns = FALSE, ...)

polymyxins(only_sir_columns = FALSE, only_treatable = TRUE, ...)

streptogramins(only_sir_columns = FALSE, ...)

quinolones(only_sir_columns = FALSE, ...)

tetracyclines(only_sir_columns = FALSE, ...)

trimethoprims(only_sir_columns = FALSE, ...)

ureidopenicillins(only_sir_columns = FALSE, ...)

administrable_per_os(only_sir_columns = FALSE, ...)

administrable_iv(only_sir_columns = FALSE, ...)

not_intrinsic_resistant(
  only_sir_columns = FALSE,
  col_mo = NULL,
  version_expertrules = 3.3,
  ...
)
```

## Arguments

- ab_class:

  an antimicrobial class or a part of it, such as `"carba"` and
  `"carbapenems"`. The columns `group`, `atc_group1` and `atc_group2` of
  the
  [antibiotics](https://amr-for-r.org/release/2.1.0/reference/antibiotics.md)
  data set will be searched (case-insensitive) for this value.

- only_sir_columns:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  only columns of class `sir` must be selected (default is `FALSE`), see
  [`as.sir()`](https://amr-for-r.org/release/2.1.0/reference/as.sir.md)

- only_treatable:

  a [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  antimicrobial drugs should be excluded that are only for laboratory
  tests (default is `TRUE`), such as gentamicin-high (`GEH`) and
  imipenem/EDTA (`IPE`)

- ...:

  ignored, only in place to allow future extensions

- filter:

  an [expression](https://rdrr.io/r/base/expression.html) to be
  evaluated in the
  [antibiotics](https://amr-for-r.org/release/2.1.0/reference/antibiotics.md)
  data set, such as `name %like% "trim"`

- col_mo:

  column name of the names or codes of the microorganisms (see
  [`as.mo()`](https://amr-for-r.org/release/2.1.0/reference/as.mo.md)) -
  the default is the first column of class
  [`mo`](https://amr-for-r.org/release/2.1.0/reference/as.mo.md). Values
  will be coerced using
  [`as.mo()`](https://amr-for-r.org/release/2.1.0/reference/as.mo.md).

- version_expertrules:

  the version number to use for the EUCAST Expert Rules and Intrinsic
  Resistance guideline. Can be "3.3", "3.2", or "3.1".

## Value

(internally) a [character](https://rdrr.io/r/base/character.html) vector
of column names, with additional class `"ab_selector"`

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

All columns in the data in which these functions are called will be
searched for known antibiotic names, abbreviations, brand names, and
codes (ATC, EARS-Net, WHO, etc.) according to the
[antibiotics](https://amr-for-r.org/release/2.1.0/reference/antibiotics.md)
data set. This means that a selector such as `aminoglycosides()` will
pick up column names like 'gen', 'genta', 'J01GB03', 'tobra',
'Tobracin', etc.

The `ab_class()` function can be used to filter/select on a manually
defined antibiotic class. It searches for results in the
[antibiotics](https://amr-for-r.org/release/2.1.0/reference/antibiotics.md)
data set within the columns `group`, `atc_group1` and `atc_group2`.

The `ab_selector()` function can be used to internally filter the
[antibiotics](https://amr-for-r.org/release/2.1.0/reference/antibiotics.md)
data set on any results, see *Examples*. It allows for filtering on a
(part of) a certain name, and/or a group name or even a minimum of DDDs
for oral treatment. This function yields the highest flexibility, but is
also the least user-friendly, since it requires a hard-coded filter to
set.

The `administrable_per_os()` and `administrable_iv()` functions also
rely on the
[antibiotics](https://amr-for-r.org/release/2.1.0/reference/antibiotics.md)
data set - antibiotic columns will be matched where a DDD (defined daily
dose) for resp. oral and IV treatment is available in the
[antibiotics](https://amr-for-r.org/release/2.1.0/reference/antibiotics.md)
data set.

The `not_intrinsic_resistant()` function can be used to only select
antibiotic columns that pose no intrinsic resistance for the
microorganisms in the data set. For example, if a data set contains only
microorganism codes or names of *E. coli* and *K. pneumoniae* and
contains a column "vancomycin", this column will be removed (or rather,
unselected) using this function. It currently applies ['EUCAST Expert
Rules' and 'EUCAST Intrinsic Resistance and Unusual Phenotypes'
v3.3](https://www.eucast.org/expert_rules_and_expected_phenotypes)
(2021) to determine intrinsic resistance, using the
[`eucast_rules()`](https://amr-for-r.org/release/2.1.0/reference/eucast_rules.md)
function internally. Because of this determination, this function is
quite slow in terms of performance.

## Full list of supported (antibiotic) classes

- `aminoglycosides()` can select:  
  amikacin (AMK), amikacin/fosfomycin (AKF), amphotericin B-high (AMH),
  apramycin (APR), arbekacin (ARB), astromicin (AST), bekanamycin (BEK),
  dibekacin (DKB), framycetin (FRM), gentamicin (GEN), gentamicin-high
  (GEH), habekacin (HAB), hygromycin (HYG), isepamicin (ISE), kanamycin
  (KAN), kanamycin-high (KAH), kanamycin/cephalexin (KAC), micronomicin
  (MCR), neomycin (NEO), netilmicin (NET), pentisomicin (PIM),
  plazomicin (PLZ), propikacin (PKA), ribostamycin (RST), sisomicin
  (SIS), streptoduocin (STR), streptomycin (STR1), streptomycin-high
  (STH), tobramycin (TOB), and tobramycin-high (TOH)

- `aminopenicillins()` can select:  
  amoxicillin (AMX) and ampicillin (AMP)

- `antifungals()` can select:  
  amphotericin B (AMB), anidulafungin (ANI), butoconazole (BUT),
  caspofungin (CAS), ciclopirox (CIX), clotrimazole (CTR), econazole
  (ECO), fluconazole (FLU), flucytosine (FCT), fosfluconazole (FFL),
  griseofulvin (GRI), hachimycin (HCH), ibrexafungerp (IBX),
  isavuconazole (ISV), isoconazole (ISO), itraconazole (ITR),
  ketoconazole (KET), manogepix (MGX), micafungin (MIF), miconazole
  (MCZ), nystatin (NYS), oteseconazole (OTE), pimaricin (PMR),
  posaconazole (POS), rezafungin (RZF), ribociclib (RBC), sulconazole
  (SUC), terbinafine (TRB), terconazole (TRC), and voriconazole (VOR)

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
  (RFP), simvastatin/fenofibrate (SMF), sodium aminosalicylate (SDA),
  streptomycin/isoniazid (STI), terizidone (TRZ), thioacetazone (TAT),
  thioacetazone/isoniazid (THI1), tiocarlide (TCR), and viomycin (VIO)

- `betalactams()` can select:  
  amoxicillin (AMX), amoxicillin/clavulanic acid (AMC),
  amoxicillin/sulbactam (AXS), ampicillin (AMP), ampicillin/sulbactam
  (SAM), apalcillin (APL), aspoxicillin (APX), avibactam (AVB),
  azidocillin (AZD), azlocillin (AZL), aztreonam (ATM),
  aztreonam/avibactam (AZA), aztreonam/nacubactam (ANC), bacampicillin
  (BAM), benzathine benzylpenicillin (BNB), benzathine
  phenoxymethylpenicillin (BNP), benzylpenicillin (PEN), biapenem (BIA),
  carbenicillin (CRB), carindacillin (CRN), cefacetrile (CAC), cefaclor
  (CEC), cefadroxil (CFR), cefalexin (LEX), cefaloridine (RID),
  cefalotin (CEP), cefamandole (MAN), cefapirin (HAP), cefatrizine
  (CTZ), cefazedone (CZD), cefazolin (CZO), cefcapene (CCP), cefcapene
  pivoxil (CCX), cefdinir (CDR), cefditoren (DIT), cefditoren pivoxil
  (DIX), cefepime (FEP), cefepime/clavulanic acid (CPC),
  cefepime/nacubactam (FNC), cefepime/tazobactam (FPT), cefetamet (CAT),
  cefetamet pivoxil (CPI), cefetecol (CCL), cefetrizole (CZL), cefixime
  (CFM), cefmenoxime (CMX), cefmetazole (CMZ), cefodizime (DIZ),
  cefonicid (CID), cefoperazone (CFP), cefoperazone/sulbactam (CSL),
  ceforanide (CND), cefoselis (CSE), cefotaxime (CTX),
  cefotaxime/clavulanic acid (CTC), cefotaxime/sulbactam (CTS),
  cefotetan (CTT), cefotiam (CTF), cefotiam hexetil (CHE), cefovecin
  (FOV), cefoxitin (FOX), cefoxitin screening (FOX1), cefozopran (ZOP),
  cefpimizole (CFZ), cefpiramide (CPM), cefpirome (CPO), cefpodoxime
  (CPD), cefpodoxime proxetil (CPX), cefpodoxime/clavulanic acid (CDC),
  cefprozil (CPR), cefquinome (CEQ), cefroxadine (CRD), cefsulodin
  (CFS), cefsumide (CSU), ceftaroline (CPT), ceftaroline/avibactam
  (CPA), ceftazidime (CAZ), ceftazidime/avibactam (CZA),
  ceftazidime/clavulanic acid (CCV), cefteram (CEM), cefteram pivoxil
  (CPL), ceftezole (CTL), ceftibuten (CTB), ceftiofur (TIO), ceftizoxime
  (CZX), ceftizoxime alapivoxil (CZP), ceftobiprole (BPR), ceftobiprole
  medocaril (CFM1), ceftolozane/tazobactam (CZT), ceftriaxone (CRO),
  ceftriaxone/beta-lactamase inhibitor (CEB), cefuroxime (CXM),
  cefuroxime axetil (CXA), cephradine (CED), ciclacillin (CIC),
  clometocillin (CLM), cloxacillin (CLO), dicloxacillin (DIC), doripenem
  (DOR), epicillin (EPC), ertapenem (ETP), flucloxacillin (FLC),
  hetacillin (HET), imipenem (IPM), imipenem/EDTA (IPE),
  imipenem/relebactam (IMR), latamoxef (LTM), lenampicillin (LEN),
  loracarbef (LOR), mecillinam (MEC), meropenem (MEM),
  meropenem/nacubactam (MNC), meropenem/vaborbactam (MEV), metampicillin
  (MTM), meticillin (MET), mezlocillin (MEZ), mezlocillin/sulbactam
  (MSU), nacubactam (NAC), nafcillin (NAF), oxacillin (OXA), panipenem
  (PAN), penamecillin (PNM), penicillin/novobiocin (PNO),
  penicillin/sulbactam (PSU), pheneticillin (PHE),
  phenoxymethylpenicillin (PHN), piperacillin (PIP),
  piperacillin/sulbactam (PIS), piperacillin/tazobactam (TZP),
  piridicillin (PRC), pivampicillin (PVM), pivmecillinam (PME), procaine
  benzylpenicillin (PRB), propicillin (PRP), razupenem (RZM), ritipenem
  (RIT), ritipenem acoxil (RIA), sarmoxicillin (SRX), sulbactam (SUL),
  sulbenicillin (SBC), sultamicillin (SLT6), talampicillin (TAL),
  tazobactam (TAZ), tebipenem (TBP), temocillin (TEM), ticarcillin
  (TIC), and ticarcillin/clavulanic acid (TCC)

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
  cefditoren pivoxil (DIX), cefepime (FEP), cefepime/clavulanic acid
  (CPC), cefepime/tazobactam (FPT), cefetamet (CAT), cefetamet pivoxil
  (CPI), cefetecol (CCL), cefetrizole (CZL), cefixime (CFM), cefmenoxime
  (CMX), cefmetazole (CMZ), cefodizime (DIZ), cefonicid (CID),
  cefoperazone (CFP), cefoperazone/sulbactam (CSL), ceforanide (CND),
  cefoselis (CSE), cefotaxime (CTX), cefotaxime/clavulanic acid (CTC),
  cefotaxime/sulbactam (CTS), cefotetan (CTT), cefotiam (CTF), cefotiam
  hexetil (CHE), cefovecin (FOV), cefoxitin (FOX), cefoxitin screening
  (FOX1), cefozopran (ZOP), cefpimizole (CFZ), cefpiramide (CPM),
  cefpirome (CPO), cefpodoxime (CPD), cefpodoxime proxetil (CPX),
  cefpodoxime/clavulanic acid (CDC), cefprozil (CPR), cefquinome (CEQ),
  cefroxadine (CRD), cefsulodin (CFS), cefsumide (CSU), ceftaroline
  (CPT), ceftaroline/avibactam (CPA), ceftazidime (CAZ),
  ceftazidime/avibactam (CZA), ceftazidime/clavulanic acid (CCV),
  cefteram (CEM), cefteram pivoxil (CPL), ceftezole (CTL), ceftibuten
  (CTB), ceftiofur (TIO), ceftizoxime (CZX), ceftizoxime alapivoxil
  (CZP), ceftobiprole (BPR), ceftobiprole medocaril (CFM1),
  ceftolozane/tazobactam (CZT), ceftriaxone (CRO),
  ceftriaxone/beta-lactamase inhibitor (CEB), cefuroxime (CXM),
  cefuroxime axetil (CXA), cephradine (CED), latamoxef (LTM), and
  loracarbef (LOR)

- `cephalosporins_1st()` can select:  
  cefacetrile (CAC), cefadroxil (CFR), cefalexin (LEX), cefaloridine
  (RID), cefalotin (CEP), cefapirin (HAP), cefatrizine (CTZ), cefazedone
  (CZD), cefazolin (CZO), cefroxadine (CRD), ceftezole (CTL), and
  cephradine (CED)

- `cephalosporins_2nd()` can select:  
  cefaclor (CEC), cefamandole (MAN), cefmetazole (CMZ), cefonicid (CID),
  ceforanide (CND), cefotetan (CTT), cefotiam (CTF), cefoxitin (FOX),
  cefoxitin screening (FOX1), cefprozil (CPR), cefuroxime (CXM),
  cefuroxime axetil (CXA), and loracarbef (LOR)

- `cephalosporins_3rd()` can select:  
  cefcapene (CCP), cefcapene pivoxil (CCX), cefdinir (CDR), cefditoren
  (DIT), cefditoren pivoxil (DIX), cefetamet (CAT), cefetamet pivoxil
  (CPI), cefixime (CFM), cefmenoxime (CMX), cefodizime (DIZ),
  cefoperazone (CFP), cefoperazone/sulbactam (CSL), cefotaxime (CTX),
  cefotaxime/clavulanic acid (CTC), cefotaxime/sulbactam (CTS), cefotiam
  hexetil (CHE), cefovecin (FOV), cefpimizole (CFZ), cefpiramide (CPM),
  cefpodoxime (CPD), cefpodoxime proxetil (CPX), cefpodoxime/clavulanic
  acid (CDC), cefsulodin (CFS), ceftazidime (CAZ), ceftazidime/avibactam
  (CZA), ceftazidime/clavulanic acid (CCV), cefteram (CEM), cefteram
  pivoxil (CPL), ceftibuten (CTB), ceftiofur (TIO), ceftizoxime (CZX),
  ceftizoxime alapivoxil (CZP), ceftriaxone (CRO),
  ceftriaxone/beta-lactamase inhibitor (CEB), and latamoxef (LTM)

- `cephalosporins_4th()` can select:  
  cefepime (FEP), cefepime/clavulanic acid (CPC), cefepime/tazobactam
  (FPT), cefetecol (CCL), cefoselis (CSE), cefozopran (ZOP), cefpirome
  (CPO), and cefquinome (CEQ)

- `cephalosporins_5th()` can select:  
  ceftaroline (CPT), ceftaroline/avibactam (CPA), ceftobiprole (BPR),
  ceftobiprole medocaril (CFM1), and ceftolozane/tazobactam (CZT)

- `fluoroquinolones()` can select:  
  besifloxacin (BES), ciprofloxacin (CIP), clinafloxacin (CLX),
  danofloxacin (DAN), delafloxacin (DFX), difloxacin (DIF), enoxacin
  (ENX), enrofloxacin (ENR), finafloxacin (FIN), fleroxacin (FLE),
  garenoxacin (GRN), gatifloxacin (GAT), gemifloxacin (GEM),
  grepafloxacin (GRX), lascufloxacin (LSC), levofloxacin (LVX),
  levonadifloxacin (LND), lomefloxacin (LOM), marbofloxacin (MAR),
  metioxate (MXT), miloxacin (MIL), moxifloxacin (MFX), nadifloxacin
  (NAD), nifuroquine (NIF), norfloxacin (NOR), ofloxacin (OFX),
  orbifloxacin (ORB), pazufloxacin (PAZ), pefloxacin (PEF),
  pradofloxacin (PRA), premafloxacin (PRX), prulifloxacin (PRU),
  rufloxacin (RFL), sarafloxacin (SAR), sitafloxacin (SIT), sparfloxacin
  (SPX), temafloxacin (TMX), tilbroquinol (TBQ), tioxacin (TXC),
  tosufloxacin (TFX), and trovafloxacin (TVA)

- `glycopeptides()` can select:  
  avoparcin (AVO), dalbavancin (DAL), norvancomycin (NVA), oritavancin
  (ORI), ramoplanin (RAM), teicoplanin (TEC), teicoplanin-macromethod
  (TCM), telavancin (TLV), vancomycin (VAN), and vancomycin-macromethod
  (VAM)

- `lincosamides()` can select:  
  acetylmidecamycin (ACM), acetylspiramycin (ASP), clindamycin (CLI),
  gamithromycin (GAM), kitasamycin (KIT), lincomycin (LIN), meleumycin
  (MEL), nafithromycin (ZWK), pirlimycin (PRL), primycin (PRM),
  solithromycin (SOL), tildipirosin (TIP), tilmicosin (TIL),
  tulathromycin (TUL), tylosin (TYL), and tylvalosin (TYL1)

- `lipoglycopeptides()` can select:  
  dalbavancin (DAL), oritavancin (ORI), and telavancin (TLV)

- `macrolides()` can select:  
  acetylmidecamycin (ACM), acetylspiramycin (ASP), azithromycin (AZM),
  clarithromycin (CLR), dirithromycin (DIR), erythromycin (ERY),
  flurithromycin (FLR1), gamithromycin (GAM), josamycin (JOS),
  kitasamycin (KIT), meleumycin (MEL), midecamycin (MID), miocamycin
  (MCM), nafithromycin (ZWK), oleandomycin (OLE), pirlimycin (PRL),
  primycin (PRM), rokitamycin (ROK), roxithromycin (RXT), solithromycin
  (SOL), spiramycin (SPI), telithromycin (TLT), tildipirosin (TIP),
  tilmicosin (TIL), troleandomycin (TRL), tulathromycin (TUL), tylosin
  (TYL), and tylvalosin (TYL1)

- `oxazolidinones()` can select:  
  cadazolid (CDZ), cycloserine (CYC), linezolid (LNZ), tedizolid (TZD),
  and thiacetazone (THA)

- `penicillins()` can select:  
  amoxicillin (AMX), amoxicillin/clavulanic acid (AMC),
  amoxicillin/sulbactam (AXS), ampicillin (AMP), ampicillin/sulbactam
  (SAM), apalcillin (APL), aspoxicillin (APX), avibactam (AVB),
  azidocillin (AZD), azlocillin (AZL), aztreonam (ATM),
  aztreonam/avibactam (AZA), aztreonam/nacubactam (ANC), bacampicillin
  (BAM), benzathine benzylpenicillin (BNB), benzathine
  phenoxymethylpenicillin (BNP), benzylpenicillin (PEN), carbenicillin
  (CRB), carindacillin (CRN), cefepime/nacubactam (FNC), ciclacillin
  (CIC), clometocillin (CLM), cloxacillin (CLO), dicloxacillin (DIC),
  epicillin (EPC), flucloxacillin (FLC), hetacillin (HET), lenampicillin
  (LEN), mecillinam (MEC), metampicillin (MTM), meticillin (MET),
  mezlocillin (MEZ), mezlocillin/sulbactam (MSU), nacubactam (NAC),
  nafcillin (NAF), oxacillin (OXA), penamecillin (PNM),
  penicillin/novobiocin (PNO), penicillin/sulbactam (PSU), pheneticillin
  (PHE), phenoxymethylpenicillin (PHN), piperacillin (PIP),
  piperacillin/sulbactam (PIS), piperacillin/tazobactam (TZP),
  piridicillin (PRC), pivampicillin (PVM), pivmecillinam (PME), procaine
  benzylpenicillin (PRB), propicillin (PRP), sarmoxicillin (SRX),
  sulbactam (SUL), sulbenicillin (SBC), sultamicillin (SLT6),
  talampicillin (TAL), tazobactam (TAZ), temocillin (TEM), ticarcillin
  (TIC), and ticarcillin/clavulanic acid (TCC)

- `polymyxins()` can select:  
  colistin (COL), polymyxin B (PLB), and polymyxin B/polysorbate 80
  (POP)

- `quinolones()` can select:  
  besifloxacin (BES), cinoxacin (CIN), ciprofloxacin (CIP),
  clinafloxacin (CLX), danofloxacin (DAN), delafloxacin (DFX),
  difloxacin (DIF), enoxacin (ENX), enrofloxacin (ENR), finafloxacin
  (FIN), fleroxacin (FLE), flumequine (FLM), garenoxacin (GRN),
  gatifloxacin (GAT), gemifloxacin (GEM), grepafloxacin (GRX),
  lascufloxacin (LSC), levofloxacin (LVX), levonadifloxacin (LND),
  lomefloxacin (LOM), marbofloxacin (MAR), metioxate (MXT), miloxacin
  (MIL), moxifloxacin (MFX), nadifloxacin (NAD), nalidixic acid (NAL),
  nemonoxacin (NEM), nifuroquine (NIF), nitroxoline (NTR), norfloxacin
  (NOR), ofloxacin (OFX), orbifloxacin (ORB), oxolinic acid (OXO),
  pazufloxacin (PAZ), pefloxacin (PEF), pipemidic acid (PPA), piromidic
  acid (PIR), pradofloxacin (PRA), premafloxacin (PRX), prulifloxacin
  (PRU), rosoxacin (ROS), rufloxacin (RFL), sarafloxacin (SAR),
  sitafloxacin (SIT), sparfloxacin (SPX), temafloxacin (TMX),
  tilbroquinol (TBQ), tioxacin (TXC), tosufloxacin (TFX), and
  trovafloxacin (TVA)

- `streptogramins()` can select:  
  pristinamycin (PRI) and quinupristin/dalfopristin (QDA)

- `tetracyclines()` can select:  
  cetocycline (CTO), chlortetracycline (CTE), clomocycline (CLM1),
  demeclocycline (DEM), doxycycline (DOX), eravacycline (ERV),
  lymecycline (LYM), metacycline (MTC), minocycline (MNO), omadacycline
  (OMC), oxytetracycline (OXY), penimepicycline (PNM1), rolitetracycline
  (RLT), sarecycline (SRC), tetracycline (TCY), and tigecycline (TGC)

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
# `example_isolates` is a data set available in the AMR package.
# See ?example_isolates.
example_isolates


# Examples sections below are split into 'base R', 'dplyr', and 'data.table':


# base R ------------------------------------------------------------------

# select columns 'IPM' (imipenem) and 'MEM' (meropenem)
example_isolates[, carbapenems()]

# select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB'
example_isolates[, c("mo", aminoglycosides())]

# select only antibiotic columns with DDDs for oral treatment
example_isolates[, administrable_per_os()]

# filter using any() or all()
example_isolates[any(carbapenems() == "R"), ]
subset(example_isolates, any(carbapenems() == "R"))

# filter on any or all results in the carbapenem columns (i.e., IPM, MEM):
example_isolates[any(carbapenems()), ]
example_isolates[all(carbapenems()), ]

# filter with multiple antibiotic selectors using c()
example_isolates[all(c(carbapenems(), aminoglycosides()) == "R"), ]

# filter + select in one go: get penicillins in carbapenem-resistant strains
example_isolates[any(carbapenems() == "R"), penicillins()]

# You can combine selectors with '&' to be more specific. For example,
# penicillins() would select benzylpenicillin ('peni G') and
# administrable_per_os() would select erythromycin. Yet, when combined these
# drugs are both omitted since benzylpenicillin is not administrable per os
# and erythromycin is not a penicillin:
example_isolates[, penicillins() & administrable_per_os()]

# ab_selector() applies a filter in the `antibiotics` data set and is thus
# very flexible. For instance, to select antibiotic columns with an oral DDD
# of at least 1 gram:
example_isolates[, ab_selector(oral_ddd > 1 & oral_units == "g")]

# \donttest{
# dplyr -------------------------------------------------------------------

if (require("dplyr")) {
  tibble(kefzol = random_sir(5)) %>%
    select(cephalosporins())
}

if (require("dplyr")) {
  # get AMR for all aminoglycosides e.g., per ward:
  example_isolates %>%
    group_by(ward) %>%
    summarise(across(aminoglycosides(), resistance))
}
if (require("dplyr")) {
  # You can combine selectors with '&' to be more specific:
  example_isolates %>%
    select(penicillins() & administrable_per_os())
}
if (require("dplyr")) {
  # get AMR for only drugs that matter - no intrinsic resistance:
  example_isolates %>%
    filter(mo_genus() %in% c("Escherichia", "Klebsiella")) %>%
    group_by(ward) %>%
    summarise(across(not_intrinsic_resistant(), resistance))
}
if (require("dplyr")) {
  # get susceptibility for antibiotics whose name contains "trim":
  example_isolates %>%
    filter(first_isolate()) %>%
    group_by(ward) %>%
    summarise(across(ab_selector(name %like% "trim"), susceptibility))
}
if (require("dplyr")) {
  # this will select columns 'IPM' (imipenem) and 'MEM' (meropenem):
  example_isolates %>%
    select(carbapenems())
}
if (require("dplyr")) {
  # this will select columns 'mo', 'AMK', 'GEN', 'KAN' and 'TOB':
  example_isolates %>%
    select(mo, aminoglycosides())
}
if (require("dplyr")) {
  # any() and all() work in dplyr's filter() too:
  example_isolates %>%
    filter(
      any(aminoglycosides() == "R"),
      all(cephalosporins_2nd() == "R")
    )
}
if (require("dplyr")) {
  # also works with c():
  example_isolates %>%
    filter(any(c(carbapenems(), aminoglycosides()) == "R"))
}
if (require("dplyr")) {
  # not setting any/all will automatically apply all():
  example_isolates %>%
    filter(aminoglycosides() == "R")
}
if (require("dplyr")) {
  # this will select columns 'mo' and all antimycobacterial drugs ('RIF'):
  example_isolates %>%
    select(mo, ab_class("mycobact"))
}
if (require("dplyr")) {
  # get bug/drug combinations for only glycopeptides in Gram-positives:
  example_isolates %>%
    filter(mo_is_gram_positive()) %>%
    select(mo, glycopeptides()) %>%
    bug_drug_combinations() %>%
    format()
}
if (require("dplyr")) {
  data.frame(
    some_column = "some_value",
    J01CA01 = "S"
  ) %>% # ATC code of ampicillin
    select(penicillins()) # only the 'J01CA01' column will be selected
}
if (require("dplyr")) {
  # with recent versions of dplyr, this is all equal:
  x <- example_isolates[carbapenems() == "R", ]
  y <- example_isolates %>% filter(carbapenems() == "R")
  z <- example_isolates %>% filter(if_all(carbapenems(), ~ .x == "R"))
  identical(x, y) && identical(y, z)
}


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
