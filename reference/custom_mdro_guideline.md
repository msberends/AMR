# Define Custom MDRO Guideline

Define custom a MDRO guideline for your organisation or specific
analysis and use the output of this function in
[`mdro()`](https://amr-for-r.org/reference/mdro.md).

## Usage

``` r
custom_mdro_guideline(..., as_factor = TRUE)

# S3 method for class 'custom_mdro_guideline'
c(x, ..., as_factor = NULL)
```

## Arguments

- ...:

  Guideline rules in [formula](https://rdrr.io/r/base/tilde.html)
  notation, see below for instructions, and in *Examples*.

- as_factor:

  A [logical](https://rdrr.io/r/base/logical.html) to indicate whether
  the returned value should be an ordered
  [factor](https://rdatatable.gitlab.io/data.table/reference/fctr.html)
  (`TRUE`, default), or otherwise a
  [character](https://rdrr.io/r/base/character.html) vector. For
  combining rules sets (using [`c()`](https://rdrr.io/r/base/c.html))
  this value will be inherited from the first set at default.

- x:

  Existing custom MDRO rules

## Value

A [list](https://rdrr.io/r/base/list.html) containing the custom rules

## Details

Using a custom MDRO guideline is of importance if you have custom rules
to determine MDROs in your hospital, e.g., rules that are dependent on
ward, state of contact isolation or other variables in your data.

### Basics

If you are familiar with the
[`case_when()`](https://dplyr.tidyverse.org/reference/case_when.html)
function of the `dplyr` package, you will recognise the input method to
set your own rules. Rules must be set using what R considers to be the
'formula notation'. The rule itself is written *before* the tilde (`~`)
and the consequence of the rule is written *after* the tilde:

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
    #>   1. If CIP is R and age is higher than 60 then: Elderly Type A
    #>   2. If ERY is R and age is higher than 60 then: Elderly Type B
    #>   3. Otherwise: Negative

    #> Unmatched rows will return NA.
    #> Results will be of class 'factor', with ordered levels: Negative < Elderly Type A < Elderly Type B

The outcome of the function can be used for the `guideline` argument in
the [`mdro()`](https://amr-for-r.org/reference/mdro.md) function:

    x <- mdro(example_isolates, guideline = custom)
    #> Determining MDROs based on custom rules, resulting in factor levels: Negative < Elderly Type A < Elderly Type B.
    #> - Custom MDRO rule 1: CIP == "R" & age > 60 (198 rows matched)
    #> - Custom MDRO rule 2: ERY == "R" & age > 60 (732 rows matched)
    #> => Found 930 custom defined MDROs out of 2000 isolates (46.5%)

    table(x)
    #> x
    #>       Negative  Elderly Type A  Elderly Type B
    #>           1070             198             732

Rules can also be combined with other custom rules by using
[`c()`](https://rdrr.io/r/base/c.html):

    x <- mdro(example_isolates,
              guideline = c(custom,
                            custom_mdro_guideline(ERY == "R" & age > 50 ~ "Elderly Type C")))
    #> Determining MDROs based on custom rules, resulting in factor levels: Negative < Elderly Type A < Elderly Type B < Elderly Type C.
    #> - Custom MDRO rule 1: CIP == "R" & age > 60 (198 rows matched)
    #> - Custom MDRO rule 2: ERY == "R" & age > 60 (732 rows matched)
    #> - Custom MDRO rule 3: ERY == "R" & age > 50 (109 rows matched)
    #> => Found 1039 custom defined MDROs out of 2000 isolates (52.0%)

    table(x)
    #> x
    #>       Negative  Elderly Type A  Elderly Type B  Elderly Type C
    #>            961             198             732             109

### Sharing rules among multiple users

The rules set (the `custom` object in this case) could be exported to a
shared file location using
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) if you collaborate
with multiple users. The custom rules set could then be imported using
[`readRDS()`](https://rdrr.io/r/base/readRDS.html).

### Usage of multiple antimicrobials and antimicrobial group names

You can define antimicrobial groups instead of single antimicrobials for
the rule itself, which is the part *before* the tilde (~). Use
[`any()`](https://rdrr.io/r/base/any.html) or
[`all()`](https://rdrr.io/r/base/all.html) to specify the scope of the
antimicrobial group:

    custom_mdro_guideline(
      AMX == "R"                       ~ "My MDRO #1",
      any(cephalosporins_2nd() == "R") ~ "My MDRO #2",
      all(glycopeptides() == "R")      ~ "My MDRO #3"
    )

All 38 antimicrobial selectors are supported for use in the rules:

- [`aminoglycosides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  amikacin, amikacin/fosfomycin, apramycin, arbekacin, astromicin,
  bekanamycin, dibekacin, framycetin, gentamicin, gentamicin-high,
  habekacin, hygromycin, isepamicin, kanamycin, kanamycin-high,
  kanamycin/cephalexin, micronomicin, neomycin, netilmicin,
  pentisomicin, plazomicin, propikacin, ribostamycin, sisomicin,
  streptoduocin, streptomycin, streptomycin-high, tobramycin, and
  tobramycin-high

- [`aminopenicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  amoxicillin, amoxicillin/clavulanic acid, and ampicillin

- [`antifungals()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  amorolfine, amphotericin B, amphotericin B-high, anidulafungin,
  butoconazole, caspofungin, ciclopirox, clotrimazole, econazole,
  fluconazole, flucytosine, fosfluconazole, griseofulvin, hachimycin,
  ibrexafungerp, isavuconazole, isoconazole, itraconazole, ketoconazole,
  manogepix, micafungin, miconazole, nystatin, oteseconazole, pimaricin,
  posaconazole, rezafungin, ribociclib, sulconazole, terbinafine,
  terconazole, and voriconazole

- [`antimycobacterials()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  4-aminosalicylic acid, calcium aminosalicylate, capreomycin,
  clofazimine, delamanid, enviomycin, ethambutol, ethambutol/isoniazid,
  ethionamide, isoniazid,
  isoniazid/sulfamethoxazole/trimethoprim/pyridoxine, morinamide,
  p-aminosalicylic acid, pretomanid, protionamide, pyrazinamide,
  rifabutin, rifampicin, rifampicin/ethambutol/isoniazid,
  rifampicin/isoniazid, rifampicin/pyrazinamide/ethambutol/isoniazid,
  rifampicin/pyrazinamide/isoniazid, rifamycin, rifapentine, sodium
  aminosalicylate, streptomycin/isoniazid, terizidone, thioacetazone,
  thioacetazone/isoniazid, tiocarlide, and viomycin

- [`betalactams()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  amoxicillin, amoxicillin/clavulanic acid, amoxicillin/sulbactam,
  ampicillin, ampicillin/sulbactam, apalcillin, aspoxicillin,
  azidocillin, azlocillin, aztreonam, aztreonam/avibactam,
  aztreonam/nacubactam, bacampicillin, benzathine benzylpenicillin,
  benzathine phenoxymethylpenicillin, benzylpenicillin, benzylpenicillin
  screening test, biapenem, carbenicillin, carindacillin, carumonam,
  cefacetrile, cefaclor, cefadroxil, cefalexin, cefaloridine, cefalotin,
  cefamandole, cefapirin, cefatrizine, cefazedone, cefazolin, cefcapene,
  cefcapene pivoxil, cefdinir, cefditoren, cefditoren pivoxil, cefepime,
  cefepime/amikacin, cefepime/clavulanic acid, cefepime/enmetazobactam,
  cefepime/nacubactam, cefepime/taniborbactam, cefepime/tazobactam,
  cefepime/zidebactam, cefetamet, cefetamet pivoxil, cefetecol,
  cefetrizole, cefiderocol, cefixime, cefmenoxime, cefmetazole,
  cefodizime, cefonicid, cefoperazone, cefoperazone/sulbactam,
  ceforanide, cefoselis, cefotaxime, cefotaxime screening test,
  cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefotetan, cefotiam,
  cefotiam hexetil, cefovecin, cefoxitin, cefoxitin screening test,
  cefozopran, cefpimizole, cefpiramide, cefpirome, cefpodoxime,
  cefpodoxime proxetil, cefpodoxime/clavulanic acid, cefprozil,
  cefquinome, cefroxadine, cefsulodin, cefsumide, ceftaroline,
  ceftaroline/avibactam, ceftazidime, ceftazidime/avibactam,
  ceftazidime/clavulanic acid, cefteram, cefteram pivoxil, ceftezole,
  ceftibuten, ceftiofur, ceftizoxime, ceftizoxime alapivoxil,
  ceftobiprole, ceftobiprole medocaril, ceftolozane/tazobactam,
  ceftriaxone, ceftriaxone/beta-lactamase inhibitor, cefuroxime,
  cefuroxime axetil, cephradine, ciclacillin, clometocillin,
  cloxacillin, dicloxacillin, doripenem, epicillin, ertapenem,
  flucloxacillin, hetacillin, imipenem, imipenem/EDTA,
  imipenem/relebactam, latamoxef, lenampicillin, loracarbef, mecillinam,
  meropenem, meropenem/nacubactam, meropenem/vaborbactam, metampicillin,
  meticillin, mezlocillin, mezlocillin/sulbactam, nafcillin, oxacillin,
  oxacillin screening test, panipenem, penamecillin,
  penicillin/novobiocin, penicillin/sulbactam, pheneticillin,
  phenoxymethylpenicillin, piperacillin, piperacillin/sulbactam,
  piperacillin/tazobactam, piridicillin, pivampicillin, pivmecillinam,
  procaine benzylpenicillin, propicillin, razupenem, ritipenem,
  ritipenem acoxil, sarmoxicillin, sulbenicillin, sultamicillin,
  talampicillin, taniborbactam, tebipenem, temocillin, ticarcillin,
  ticarcillin/clavulanic acid, and tigemonam

- [`betalactams_with_inhibitor()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  amoxicillin/clavulanic acid, amoxicillin/sulbactam,
  ampicillin/sulbactam, aztreonam/avibactam, aztreonam/nacubactam,
  cefepime/amikacin, cefepime/clavulanic acid, cefepime/enmetazobactam,
  cefepime/nacubactam, cefepime/taniborbactam, cefepime/tazobactam,
  cefepime/zidebactam, cefoperazone/sulbactam, cefotaxime/clavulanic
  acid, cefotaxime/sulbactam, cefpodoxime/clavulanic acid,
  ceftaroline/avibactam, ceftazidime/avibactam, ceftazidime/clavulanic
  acid, ceftolozane/tazobactam, ceftriaxone/beta-lactamase inhibitor,
  imipenem/relebactam, meropenem/nacubactam, meropenem/vaborbactam,
  mezlocillin/sulbactam, penicillin/novobiocin, penicillin/sulbactam,
  piperacillin/sulbactam, piperacillin/tazobactam, and
  ticarcillin/clavulanic acid

- [`carbapenems()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  biapenem, doripenem, ertapenem, imipenem, imipenem/EDTA,
  imipenem/relebactam, meropenem, meropenem/nacubactam,
  meropenem/vaborbactam, panipenem, razupenem, ritipenem, ritipenem
  acoxil, taniborbactam, and tebipenem

- [`cephalosporins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  cefacetrile, cefaclor, cefadroxil, cefalexin, cefaloridine, cefalotin,
  cefamandole, cefapirin, cefatrizine, cefazedone, cefazolin, cefcapene,
  cefcapene pivoxil, cefdinir, cefditoren, cefditoren pivoxil, cefepime,
  cefepime/amikacin, cefepime/clavulanic acid, cefepime/enmetazobactam,
  cefepime/nacubactam, cefepime/taniborbactam, cefepime/tazobactam,
  cefepime/zidebactam, cefetamet, cefetamet pivoxil, cefetecol,
  cefetrizole, cefiderocol, cefixime, cefmenoxime, cefmetazole,
  cefodizime, cefonicid, cefoperazone, cefoperazone/sulbactam,
  ceforanide, cefoselis, cefotaxime, cefotaxime screening test,
  cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefotetan, cefotiam,
  cefotiam hexetil, cefovecin, cefoxitin, cefoxitin screening test,
  cefozopran, cefpimizole, cefpiramide, cefpirome, cefpodoxime,
  cefpodoxime proxetil, cefpodoxime/clavulanic acid, cefprozil,
  cefquinome, cefroxadine, cefsulodin, cefsumide, ceftaroline,
  ceftaroline/avibactam, ceftazidime, ceftazidime/avibactam,
  ceftazidime/clavulanic acid, cefteram, cefteram pivoxil, ceftezole,
  ceftibuten, ceftiofur, ceftizoxime, ceftizoxime alapivoxil,
  ceftobiprole, ceftobiprole medocaril, ceftolozane/tazobactam,
  ceftriaxone, ceftriaxone/beta-lactamase inhibitor, cefuroxime,
  cefuroxime axetil, cephradine, latamoxef, and loracarbef

- [`cephalosporins_1st()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  cefacetrile, cefadroxil, cefalexin, cefaloridine, cefalotin,
  cefapirin, cefatrizine, cefazedone, cefazolin, cefroxadine, ceftezole,
  and cephradine

- [`cephalosporins_2nd()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  cefaclor, cefamandole, cefmetazole, cefonicid, ceforanide, cefotetan,
  cefotiam, cefoxitin, cefoxitin screening test, cefprozil, cefuroxime,
  cefuroxime axetil, and loracarbef

- [`cephalosporins_3rd()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren
  pivoxil, cefetamet, cefetamet pivoxil, cefixime, cefmenoxime,
  cefodizime, cefoperazone, cefoperazone/sulbactam, cefotaxime,
  cefotaxime screening test, cefotaxime/clavulanic acid,
  cefotaxime/sulbactam, cefotiam hexetil, cefovecin, cefpimizole,
  cefpiramide, cefpodoxime, cefpodoxime proxetil, cefpodoxime/clavulanic
  acid, cefsulodin, ceftazidime, ceftazidime/avibactam,
  ceftazidime/clavulanic acid, cefteram, cefteram pivoxil, ceftibuten,
  ceftiofur, ceftizoxime, ceftizoxime alapivoxil, ceftriaxone,
  ceftriaxone/beta-lactamase inhibitor, and latamoxef

- [`cephalosporins_4th()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  cefepime, cefepime/amikacin, cefepime/clavulanic acid,
  cefepime/enmetazobactam, cefepime/nacubactam, cefepime/taniborbactam,
  cefepime/tazobactam, cefepime/zidebactam, cefetecol, cefoselis,
  cefozopran, cefpirome, and cefquinome

- [`cephalosporins_5th()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  ceftaroline, ceftaroline/avibactam, ceftobiprole, ceftobiprole
  medocaril, and ceftolozane/tazobactam

- [`fluoroquinolones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  besifloxacin, ciprofloxacin, ciprofloxacin/metronidazole,
  ciprofloxacin/ornidazole, ciprofloxacin/tinidazole, clinafloxacin,
  danofloxacin, delafloxacin, difloxacin, enoxacin, enrofloxacin,
  finafloxacin, fleroxacin, garenoxacin, gatifloxacin, gemifloxacin,
  grepafloxacin, lascufloxacin, levofloxacin, levofloxacin/ornidazole,
  levonadifloxacin, lomefloxacin, marbofloxacin, metioxate, miloxacin,
  moxifloxacin, nadifloxacin, nemonoxacin, nifuroquine, nitroxoline,
  norfloxacin, norfloxacin screening test, norfloxacin/metronidazole,
  norfloxacin/tinidazole, ofloxacin, ofloxacin/ornidazole, orbifloxacin,
  pazufloxacin, pefloxacin, pefloxacin screening test, pradofloxacin,
  premafloxacin, prulifloxacin, rufloxacin, sarafloxacin, sitafloxacin,
  sparfloxacin, temafloxacin, tilbroquinol, tioxacin, tosufloxacin, and
  trovafloxacin

- [`glycopeptides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  avoparcin, bleomycin, dalbavancin, norvancomycin, oritavancin,
  ramoplanin, teicoplanin, teicoplanin-macromethod, telavancin,
  vancomycin, and vancomycin-macromethod

- [`isoxazolylpenicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  cloxacillin, dicloxacillin, flucloxacillin, meticillin, oxacillin, and
  oxacillin screening test

- [`lincosamides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  clindamycin, clindamycin inducible screening test, lincomycin, and
  pirlimycin

- [`lipoglycopeptides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  dalbavancin, oritavancin, and telavancin

- [`macrolides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  acetylmidecamycin, acetylspiramycin, azithromycin, clarithromycin,
  clindamycin inducible screening test, dirithromycin, erythromycin,
  flurithromycin, gamithromycin, josamycin, kitasamycin, meleumycin,
  midecamycin, miocamycin, nafithromycin, oleandomycin, pirlimycin,
  primycin, rokitamycin, roxithromycin, solithromycin, spiramycin,
  telithromycin, tildipirosin, tilmicosin, troleandomycin,
  tulathromycin, tylosin, and tylvalosin

- [`monobactams()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  aztreonam, aztreonam/avibactam, aztreonam/nacubactam, carumonam, and
  tigemonam

- [`nitrofurans()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  furazidin, furazolidone, nifurtoinol, nitrofurantoin, and
  nitrofurazone

- [`oxazolidinones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  cadazolid, cycloserine, linezolid, tedizolid, and thiacetazone

- [`penicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  amoxicillin, amoxicillin/clavulanic acid, amoxicillin/sulbactam,
  ampicillin, ampicillin/sulbactam, apalcillin, aspoxicillin,
  azidocillin, azlocillin, bacampicillin, benzathine benzylpenicillin,
  benzathine phenoxymethylpenicillin, benzylpenicillin, benzylpenicillin
  screening test, carbenicillin, carindacillin, ciclacillin,
  clometocillin, cloxacillin, dicloxacillin, epicillin, flucloxacillin,
  hetacillin, lenampicillin, mecillinam, metampicillin, meticillin,
  mezlocillin, mezlocillin/sulbactam, nafcillin, oxacillin, oxacillin
  screening test, penamecillin, penicillin/novobiocin,
  penicillin/sulbactam, pheneticillin, phenoxymethylpenicillin,
  piperacillin, piperacillin/sulbactam, piperacillin/tazobactam,
  piridicillin, pivampicillin, pivmecillinam, procaine benzylpenicillin,
  propicillin, sarmoxicillin, sulbenicillin, sultamicillin,
  talampicillin, temocillin, ticarcillin, and ticarcillin/clavulanic
  acid

- [`phenicols()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  chloramphenicol, florfenicol, and thiamphenicol

- [`phosphonics()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  amikacin/fosfomycin and fosfomycin

- [`polymyxins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  colistin, polymyxin B, and polymyxin B/polysorbate 80

- [`quinolones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  besifloxacin, cinoxacin, ciprofloxacin, ciprofloxacin/metronidazole,
  ciprofloxacin/ornidazole, ciprofloxacin/tinidazole, clinafloxacin,
  danofloxacin, delafloxacin, difloxacin, enoxacin, enrofloxacin,
  finafloxacin, fleroxacin, flumequine, garenoxacin, gatifloxacin,
  gemifloxacin, grepafloxacin, lascufloxacin, levofloxacin,
  levofloxacin/ornidazole, levonadifloxacin, lomefloxacin,
  marbofloxacin, metioxate, miloxacin, moxifloxacin, nadifloxacin,
  nalidixic acid, nalidixic acid screening test, nemonoxacin,
  nifuroquine, nitroxoline, norfloxacin, norfloxacin screening test,
  norfloxacin/metronidazole, norfloxacin/tinidazole, ofloxacin,
  ofloxacin/ornidazole, orbifloxacin, oxolinic acid, ozenoxacin,
  pazufloxacin, pefloxacin, pefloxacin screening test, pipemidic acid,
  piromidic acid, pradofloxacin, premafloxacin, prulifloxacin,
  rosoxacin, rufloxacin, sarafloxacin, sitafloxacin, sparfloxacin,
  temafloxacin, tilbroquinol, tioxacin, tosufloxacin, and trovafloxacin

- [`rifamycins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  rifabutin, rifampicin, rifampicin/ethambutol/isoniazid,
  rifampicin/isoniazid, rifampicin/pyrazinamide/ethambutol/isoniazid,
  rifampicin/pyrazinamide/isoniazid, rifamycin, and rifapentine

- [`spiropyrimidinetriones()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  zoliflodacin

- [`streptogramins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  pristinamycin and quinupristin/dalfopristin

- [`sulfonamides()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  brodimoprim, sulfadiazine, sulfadiazine/tetroxoprim, sulfadimethoxine,
  sulfadimidine, sulfafurazole, sulfaisodimidine, sulfalene,
  sulfamazone, sulfamerazine, sulfamethizole, sulfamethoxazole,
  sulfamethoxypyridazine, sulfametomidine, sulfametoxydiazine,
  sulfamoxole, sulfanilamide, sulfaperin, sulfaphenazole, sulfapyridine,
  sulfathiazole, and sulfathiourea

- [`tetracyclines()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  cetocycline, chlortetracycline, clomocycline, demeclocycline,
  doxycycline, eravacycline, lymecycline, metacycline, minocycline,
  omadacycline, oxytetracycline, penimepicycline, rolitetracycline,
  sarecycline, tetracycline, tetracycline screening test, and
  tigecycline

- [`trimethoprims()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  brodimoprim, sulfadiazine, sulfadiazine/tetroxoprim,
  sulfadiazine/trimethoprim, sulfadimethoxine, sulfadimidine,
  sulfadimidine/trimethoprim, sulfafurazole, sulfaisodimidine,
  sulfalene, sulfamazone, sulfamerazine, sulfamerazine/trimethoprim,
  sulfamethizole, sulfamethoxazole, sulfamethoxypyridazine,
  sulfametomidine, sulfametoxydiazine, sulfametrole/trimethoprim,
  sulfamoxole, sulfamoxole/trimethoprim, sulfanilamide, sulfaperin,
  sulfaphenazole, sulfapyridine, sulfathiazole, sulfathiourea,
  trimethoprim, and trimethoprim/sulfamethoxazole

- [`ureidopenicillins()`](https://amr-for-r.org/reference/antimicrobial_selectors.md)
  can select:  
  azlocillin, mezlocillin, piperacillin, and piperacillin/tazobactam

## Examples

``` r
x <- custom_mdro_guideline(
  CIP == "R" & age > 60 ~ "Elderly Type A",
  ERY == "R" & age > 60 ~ "Elderly Type B"
)
x
#> A set of custom MDRO rules:
#>   1. If CIP is  R  and age is higher than 60 then: Elderly Type A
#>   2. If ERY is  R  and age is higher than 60 then: Elderly Type B
#>   3. Otherwise: Negative
#> 
#> Unmatched rows will return NA.
#> Results will be of class 'factor', with ordered levels: Negative < Elderly Type A < Elderly Type B

# run the custom rule set (verbose = TRUE will return a logbook instead of the data set):
out <- mdro(example_isolates, guideline = x)
table(out)
#> out
#>       Negative Elderly Type A Elderly Type B 
#>           1070            198            732 

out <- mdro(example_isolates, guideline = x, verbose = TRUE)
head(out)
#>    row_number microorganism           MDRO
#> V1          1          <NA> Elderly Type B
#> V2          2          <NA> Elderly Type B
#> V3          3          <NA>       Negative
#> V4          4          <NA>       Negative
#> V5          5          <NA>       Negative
#> V6          6          <NA>       Negative
#>                                   reason
#> V1 matched rule 2: ERY == "R" & age > 60
#> V2 matched rule 2: ERY == "R" & age > 60
#> V3                      no rules matched
#> V4                      no rules matched
#> V5                      no rules matched
#> V6                      no rules matched
#>                               all_nonsusceptible_columns        guideline
#> V1 PEN, TMP, SXT, LNZ, VAN, TEC, TCY, ERY, CLI, AZM, RIF Custom guideline
#> V2 PEN, TMP, SXT, LNZ, VAN, TEC, TCY, ERY, CLI, AZM, RIF Custom guideline
#> V3                     PEN, FLC, CXM, CAZ, ERY, AZM, COL Custom guideline
#> V4                     PEN, FLC, CXM, CAZ, ERY, AZM, COL Custom guideline
#> V5                PEN, FLC, CXM, CAZ, TMP, ERY, AZM, COL Custom guideline
#> V6           PEN, FLC, CXM, CAZ, TMP, ERY, CLI, AZM, COL Custom guideline

# you can create custom guidelines using selectors (see ?antimicrobial_selectors)
my_guideline <- custom_mdro_guideline(
  AMX == "R" ~ "Custom MDRO 1",
  all(cephalosporins_2nd() == "R") ~ "Custom MDRO 2"
)
my_guideline
#> A set of custom MDRO rules:
#>   1. If AMX is  R  then: Custom MDRO 1
#>   2. If all of cephalosporins_2nd() is  R  then: Custom MDRO 2
#>   3. Otherwise: Negative
#> 
#> Unmatched rows will return NA.
#> Results will be of class 'factor', with ordered levels: Negative < Custom MDRO 1 < Custom MDRO 2

out <- mdro(example_isolates, guideline = my_guideline)
#> ℹ Column 'esbl' is SIR eligible (despite only having empty values), since
#>   it seems to be tazobactam (TAZ)
#> ℹ Column 'mecC' is SIR eligible (despite only having empty values), since
#>   it seems to be mecillinam (MEC)
#> ℹ Column 'vanA' is SIR eligible (despite only having empty values), since
#>   it seems to be lenampicillin (LEN)
#> ℹ Column 'vanB' is SIR eligible (despite only having empty values), since
#>   it seems to be metronidazole (MTR)
#> ℹ For `cephalosporins_2nd()` using columns 'CXM' (cefuroxime) and 'FOX'
#>   (cefoxitin)
#> ℹ Assuming a filter on all 2 cephalosporins_2nd. Wrap around `all()` or
#>   `any()` to prevent this note.
table(out)
#> out
#>      Negative Custom MDRO 1 Custom MDRO 2 
#>          1144           804            52 
```
