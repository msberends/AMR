# Define Custom EUCAST Rules

Define custom EUCAST rules for your organisation or specific analysis
and use the output of this function in
[`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md).

## Usage

``` r
custom_eucast_rules(...)
```

## Arguments

- ...:

  Rules in [formula](https://rdrr.io/r/base/tilde.html) notation, see
  below for instructions, and in *Examples*.

## Value

A [list](https://rdrr.io/r/base/list.html) containing the custom rules

## Details

Some organisations have their own adoption of EUCAST rules. This
function can be used to define custom EUCAST rules to be used in the
[`eucast_rules()`](https://amr-for-r.org/reference/eucast_rules.md)
function.

### Basics

If you are familiar with the
[`case_when()`](https://dplyr.tidyverse.org/reference/case_when.html)
function of the `dplyr` package, you will recognise the input method to
set your own rules. Rules must be set using what R considers to be the
'formula notation'. The rule itself is written *before* the tilde (`~`)
and the consequence of the rule is written *after* the tilde:

    x <- custom_eucast_rules(TZP == "S" ~ aminopenicillins == "S",
                             TZP == "R" ~ aminopenicillins == "R")

These are two custom EUCAST rules: if TZP (piperacillin/tazobactam) is
"S", all aminopenicillins (ampicillin and amoxicillin) must be made "S",
and if TZP is "R", aminopenicillins must be made "R". These rules can
also be printed to the console, so it is immediately clear how they
work:

    x
    #> A set of custom EUCAST rules:
    #>
    #>   1. If TZP is "S" then set to  S :
    #>      amoxicillin (AMX), ampicillin (AMP)
    #>
    #>   2. If TZP is "R" then set to  R :
    #>      amoxicillin (AMX), ampicillin (AMP)

The rules (the part *before* the tilde, in above example `TZP == "S"`
and `TZP == "R"`) must be evaluable in your data set: it should be able
to run as a filter in your data set without errors. This means for the
above example that the column `TZP` must exist. We will create a sample
data set and test the rules set:

    df <- data.frame(mo = c("Escherichia coli", "Klebsiella pneumoniae"),
                     TZP = as.sir("R"),
                     ampi = as.sir("S"),
                     cipro = as.sir("S"))
    df
    #>                      mo TZP ampi cipro
    #> 1      Escherichia coli   R    S     S
    #> 2 Klebsiella pneumoniae   R    S     S

    eucast_rules(df,
                 rules = "custom",
                 custom_rules = x,
                 info = FALSE,
                 overwrite = TRUE)
    #>                      mo TZP ampi cipro
    #> 1      Escherichia coli   R    R     S
    #> 2 Klebsiella pneumoniae   R    R     S

### Using taxonomic properties in rules

There is one exception in columns used for the rules: all column names
of the
[microorganisms](https://amr-for-r.org/reference/microorganisms.md) data
set can also be used, but do not have to exist in the data set. These
column names are: "mo", "fullname", "status", "kingdom", "phylum",
"class", "order", "family", "genus", "species", "subspecies", "rank",
"ref", "oxygen_tolerance", "source", "lpsn", "lpsn_parent",
"lpsn_renamed_to", "mycobank", "mycobank_parent", "mycobank_renamed_to",
"gbif", "gbif_parent", "gbif_renamed_to", "prevalence", and "snomed".
Thus, this next example will work as well, despite the fact that the
`df` data set does not contain a column `genus`:

    y <- custom_eucast_rules(
      TZP == "S" & genus == "Klebsiella" ~ aminopenicillins == "S",
      TZP == "R" & genus == "Klebsiella" ~ aminopenicillins == "R"
    )

    eucast_rules(df,
                 rules = "custom",
                 custom_rules = y,
                 info = FALSE,
                 overwrite = TRUE)
    #>                      mo TZP ampi cipro
    #> 1      Escherichia coli   R    S     S
    #> 2 Klebsiella pneumoniae   R    R     S

### Sharing rules among multiple users

The rules set (the `y` object in this case) could be exported to a
shared file location using
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) if you collaborate
with multiple users. The custom rules set could then be imported using
[`readRDS()`](https://rdrr.io/r/base/readRDS.html).

### Usage of multiple antimicrobials and antimicrobial group names

You can define antimicrobial groups instead of single antimicrobials for
the rule consequence, which is the part *after* the tilde (~). In the
examples above, the antimicrobial group `aminopenicillins` includes both
ampicillin and amoxicillin.

Rules can also be applied to multiple antimicrobials and antimicrobial
groups simultaneously. Use the [`c()`](https://rdrr.io/r/base/c.html)
function to combine multiple antimicrobials. For instance, the following
example sets all aminopenicillins and ureidopenicillins to "R" if column
TZP (piperacillin/tazobactam) is "R":

    x <- custom_eucast_rules(TZP == "R" ~ c(aminopenicillins, ureidopenicillins) == "R")
    x
    #> A set of custom EUCAST rules:
    #>
    #>   1. If TZP is "R" then set to "R":
    #>      amoxicillin (AMX), ampicillin (AMP), azlocillin (AZL), mezlocillin (MEZ), piperacillin (PIP), piperacillin/tazobactam (TZP)

These 35 antimicrobial groups are allowed in the rules
(case-insensitive) and can be used in any combination:

- aminoglycosides  
  (amikacin, amikacin/fosfomycin, apramycin, arbekacin, astromicin,
  bekanamycin, dibekacin, framycetin, gentamicin, gentamicin-high,
  habekacin, hygromycin, isepamicin, kanamycin, kanamycin-high,
  kanamycin/cephalexin, micronomicin, neomycin, netilmicin,
  pentisomicin, plazomicin, propikacin, ribostamycin, sisomicin,
  streptoduocin, streptomycin, streptomycin-high, tobramycin, and
  tobramycin-high)

- aminopenicillins  
  (amoxicillin and ampicillin)

- antifungals  
  (amorolfine, amphotericin B, amphotericin B-high, anidulafungin,
  butoconazole, caspofungin, ciclopirox, clotrimazole, econazole,
  fluconazole, flucytosine, fosfluconazole, griseofulvin, hachimycin,
  ibrexafungerp, isavuconazole, isoconazole, itraconazole, ketoconazole,
  manogepix, micafungin, miconazole, nystatin, oteseconazole, pimaricin,
  posaconazole, rezafungin, ribociclib, sulconazole, terbinafine,
  terconazole, and voriconazole)

- antimycobacterials  
  (4-aminosalicylic acid, calcium aminosalicylate, capreomycin,
  clofazimine, delamanid, enviomycin, ethambutol, ethambutol/isoniazid,
  ethionamide, isoniazid,
  isoniazid/sulfamethoxazole/trimethoprim/pyridoxine, morinamide,
  p-aminosalicylic acid, pretomanid, protionamide, pyrazinamide,
  rifabutin, rifampicin, rifampicin/ethambutol/isoniazid,
  rifampicin/isoniazid, rifampicin/pyrazinamide/ethambutol/isoniazid,
  rifampicin/pyrazinamide/isoniazid, rifamycin, rifapentine, sodium
  aminosalicylate, streptomycin/isoniazid, terizidone, thioacetazone,
  thioacetazone/isoniazid, tiocarlide, and viomycin)

- betalactams  
  (amoxicillin, amoxicillin/clavulanic acid, amoxicillin/sulbactam,
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
  ticarcillin/clavulanic acid, and tigemonam)

- betalactams_with_inhibitor  
  (amoxicillin/clavulanic acid, amoxicillin/sulbactam,
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
  ticarcillin/clavulanic acid)

- carbapenems  
  (biapenem, doripenem, ertapenem, imipenem, imipenem/EDTA,
  imipenem/relebactam, meropenem, meropenem/nacubactam,
  meropenem/vaborbactam, panipenem, razupenem, ritipenem, ritipenem
  acoxil, taniborbactam, and tebipenem)

- cephalosporins  
  (cefacetrile, cefaclor, cefadroxil, cefalexin, cefaloridine,
  cefalotin, cefamandole, cefapirin, cefatrizine, cefazedone, cefazolin,
  cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren
  pivoxil, cefepime, cefepime/amikacin, cefepime/clavulanic acid,
  cefepime/enmetazobactam, cefepime/nacubactam, cefepime/taniborbactam,
  cefepime/tazobactam, cefepime/zidebactam, cefetamet, cefetamet
  pivoxil, cefetecol, cefetrizole, cefiderocol, cefixime, cefmenoxime,
  cefmetazole, cefodizime, cefonicid, cefoperazone,
  cefoperazone/sulbactam, ceforanide, cefoselis, cefotaxime, cefotaxime
  screening test, cefotaxime/clavulanic acid, cefotaxime/sulbactam,
  cefotetan, cefotiam, cefotiam hexetil, cefovecin, cefoxitin, cefoxitin
  screening test, cefozopran, cefpimizole, cefpiramide, cefpirome,
  cefpodoxime, cefpodoxime proxetil, cefpodoxime/clavulanic acid,
  cefprozil, cefquinome, cefroxadine, cefsulodin, cefsumide,
  ceftaroline, ceftaroline/avibactam, ceftazidime,
  ceftazidime/avibactam, ceftazidime/clavulanic acid, cefteram, cefteram
  pivoxil, ceftezole, ceftibuten, ceftiofur, ceftizoxime, ceftizoxime
  alapivoxil, ceftobiprole, ceftobiprole medocaril,
  ceftolozane/tazobactam, ceftriaxone, ceftriaxone/beta-lactamase
  inhibitor, cefuroxime, cefuroxime axetil, cephradine, latamoxef, and
  loracarbef)

- cephalosporins_1st  
  (cefacetrile, cefadroxil, cefalexin, cefaloridine, cefalotin,
  cefapirin, cefatrizine, cefazedone, cefazolin, cefroxadine, ceftezole,
  and cephradine)

- cephalosporins_2nd  
  (cefaclor, cefamandole, cefmetazole, cefonicid, ceforanide, cefotetan,
  cefotiam, cefoxitin, cefoxitin screening test, cefprozil, cefuroxime,
  cefuroxime axetil, and loracarbef)

- cephalosporins_3rd  
  (cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren
  pivoxil, cefetamet, cefetamet pivoxil, cefixime, cefmenoxime,
  cefodizime, cefoperazone, cefoperazone/sulbactam, cefotaxime,
  cefotaxime screening test, cefotaxime/clavulanic acid,
  cefotaxime/sulbactam, cefotiam hexetil, cefovecin, cefpimizole,
  cefpiramide, cefpodoxime, cefpodoxime proxetil, cefpodoxime/clavulanic
  acid, cefsulodin, ceftazidime, ceftazidime/avibactam,
  ceftazidime/clavulanic acid, cefteram, cefteram pivoxil, ceftibuten,
  ceftiofur, ceftizoxime, ceftizoxime alapivoxil, ceftriaxone,
  ceftriaxone/beta-lactamase inhibitor, and latamoxef)

- cephalosporins_4th  
  (cefepime, cefepime/amikacin, cefepime/clavulanic acid,
  cefepime/enmetazobactam, cefepime/nacubactam, cefepime/taniborbactam,
  cefepime/tazobactam, cefepime/zidebactam, cefetecol, cefoselis,
  cefozopran, cefpirome, and cefquinome)

- cephalosporins_5th  
  (ceftaroline, ceftaroline/avibactam, ceftobiprole, ceftobiprole
  medocaril, and ceftolozane/tazobactam)

- cephalosporins_except_caz  
  (cefacetrile, cefaclor, cefadroxil, cefalexin, cefaloridine,
  cefalotin, cefamandole, cefapirin, cefatrizine, cefazedone, cefazolin,
  cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren
  pivoxil, cefepime, cefepime/amikacin, cefepime/clavulanic acid,
  cefepime/enmetazobactam, cefepime/nacubactam, cefepime/taniborbactam,
  cefepime/tazobactam, cefepime/zidebactam, cefetamet, cefetamet
  pivoxil, cefetecol, cefetrizole, cefiderocol, cefixime, cefmenoxime,
  cefmetazole, cefodizime, cefonicid, cefoperazone,
  cefoperazone/sulbactam, ceforanide, cefoselis, cefotaxime, cefotaxime
  screening test, cefotaxime/clavulanic acid, cefotaxime/sulbactam,
  cefotetan, cefotiam, cefotiam hexetil, cefovecin, cefoxitin, cefoxitin
  screening test, cefozopran, cefpimizole, cefpiramide, cefpirome,
  cefpodoxime, cefpodoxime proxetil, cefpodoxime/clavulanic acid,
  cefprozil, cefquinome, cefroxadine, cefsulodin, cefsumide,
  ceftaroline, ceftaroline/avibactam, ceftazidime/avibactam,
  ceftazidime/clavulanic acid, cefteram, cefteram pivoxil, ceftezole,
  ceftibuten, ceftiofur, ceftizoxime, ceftizoxime alapivoxil,
  ceftobiprole, ceftobiprole medocaril, ceftolozane/tazobactam,
  ceftriaxone, ceftriaxone/beta-lactamase inhibitor, cefuroxime,
  cefuroxime axetil, cephradine, latamoxef, and loracarbef)

- fluoroquinolones  
  (besifloxacin, ciprofloxacin, ciprofloxacin/metronidazole,
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
  trovafloxacin)

- glycopeptides  
  (avoparcin, bleomycin, dalbavancin, norvancomycin, oritavancin,
  ramoplanin, teicoplanin, teicoplanin-macromethod, telavancin,
  vancomycin, and vancomycin-macromethod)

- glycopeptides_except_lipo  
  (avoparcin, bleomycin, norvancomycin, ramoplanin, teicoplanin,
  teicoplanin-macromethod, vancomycin, and vancomycin-macromethod)

- isoxazolylpenicillins  
  (cloxacillin, dicloxacillin, flucloxacillin, meticillin, oxacillin,
  and oxacillin screening test)

- lincosamides  
  (clindamycin, lincomycin, and pirlimycin)

- lipoglycopeptides  
  (dalbavancin, oritavancin, and telavancin)

- macrolides  
  (acetylmidecamycin, acetylspiramycin, azithromycin, clarithromycin,
  dirithromycin, erythromycin, flurithromycin, gamithromycin, josamycin,
  kitasamycin, meleumycin, midecamycin, miocamycin, nafithromycin,
  oleandomycin, rokitamycin, roxithromycin, solithromycin, spiramycin,
  telithromycin, tildipirosin, tilmicosin, troleandomycin,
  tulathromycin, tylosin, and tylvalosin)

- monobactams  
  (aztreonam, aztreonam/avibactam, aztreonam/nacubactam, carumonam, and
  tigemonam)

- nitrofurans  
  (furazidin, furazolidone, nifurtoinol, nitrofurantoin, and
  nitrofurazone)

- oxazolidinones  
  (cadazolid, cycloserine, linezolid, tedizolid, and thiacetazone)

- penicillins  
  (amoxicillin, amoxicillin/clavulanic acid, amoxicillin/sulbactam,
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
  acid)

- phenicols  
  (chloramphenicol, florfenicol, and thiamphenicol)

- polymyxins  
  (colistin, polymyxin B, and polymyxin B/polysorbate 80)

- quinolones  
  (besifloxacin, cinoxacin, ciprofloxacin, ciprofloxacin/metronidazole,
  ciprofloxacin/ornidazole, ciprofloxacin/tinidazole, clinafloxacin,
  danofloxacin, delafloxacin, difloxacin, enoxacin, enrofloxacin,
  finafloxacin, fleroxacin, flumequine, garenoxacin, gatifloxacin,
  gemifloxacin, grepafloxacin, lascufloxacin, levofloxacin,
  levofloxacin/ornidazole, levonadifloxacin, lomefloxacin,
  marbofloxacin, metioxate, miloxacin, moxifloxacin, nadifloxacin,
  nalidixic acid, nalidixic acid screening test, nemonoxacin,
  nifuroquine, nitroxoline, norfloxacin, norfloxacin screening test,
  norfloxacin/metronidazole, norfloxacin/tinidazole, ofloxacin,
  ofloxacin/ornidazole, orbifloxacin, oxolinic acid, pazufloxacin,
  pefloxacin, pefloxacin screening test, pipemidic acid, piromidic acid,
  pradofloxacin, premafloxacin, prulifloxacin, rosoxacin, rufloxacin,
  sarafloxacin, sitafloxacin, sparfloxacin, temafloxacin, tilbroquinol,
  tioxacin, tosufloxacin, and trovafloxacin)

- rifamycins  
  (rifabutin, rifampicin, rifampicin/ethambutol/isoniazid,
  rifampicin/isoniazid, rifampicin/pyrazinamide/ethambutol/isoniazid,
  rifampicin/pyrazinamide/isoniazid, rifamycin, and rifapentine)

- streptogramins  
  (pristinamycin and quinupristin/dalfopristin)

- sulfonamides  
  (brodimoprim, sulfadiazine, sulfadiazine/tetroxoprim,
  sulfadimethoxine, sulfadimidine, sulfafurazole, sulfaisodimidine,
  sulfalene, sulfamazone, sulfamerazine, sulfamethizole,
  sulfamethoxazole, sulfamethoxypyridazine, sulfametomidine,
  sulfametoxydiazine, sulfamoxole, sulfanilamide, sulfaperin,
  sulfaphenazole, sulfapyridine, sulfathiazole, and sulfathiourea)

- tetracyclines  
  (cetocycline, chlortetracycline, clomocycline, demeclocycline,
  doxycycline, eravacycline, lymecycline, metacycline, minocycline,
  omadacycline, oxytetracycline, penimepicycline, rolitetracycline,
  sarecycline, tetracycline, tetracycline screening test, and
  tigecycline)

- tetracyclines_except_tgc  
  (cetocycline, chlortetracycline, clomocycline, demeclocycline,
  doxycycline, eravacycline, lymecycline, metacycline, minocycline,
  omadacycline, oxytetracycline, penimepicycline, rolitetracycline,
  sarecycline, tetracycline, and tetracycline screening test)

- trimethoprims  
  (brodimoprim, sulfadiazine, sulfadiazine/tetroxoprim,
  sulfadiazine/trimethoprim, sulfadimethoxine, sulfadimidine,
  sulfadimidine/trimethoprim, sulfafurazole, sulfaisodimidine,
  sulfalene, sulfamazone, sulfamerazine, sulfamerazine/trimethoprim,
  sulfamethizole, sulfamethoxazole, sulfamethoxypyridazine,
  sulfametomidine, sulfametoxydiazine, sulfametrole/trimethoprim,
  sulfamoxole, sulfamoxole/trimethoprim, sulfanilamide, sulfaperin,
  sulfaphenazole, sulfapyridine, sulfathiazole, sulfathiourea,
  trimethoprim, and trimethoprim/sulfamethoxazole)

- ureidopenicillins  
  (azlocillin, mezlocillin, piperacillin, and piperacillin/tazobactam)

## Examples

``` r
x <- custom_eucast_rules(
  AMC == "R" & genus == "Klebsiella" ~ aminopenicillins == "R",
  AMC == "I" & genus == "Klebsiella" ~ aminopenicillins == "I"
)
x
#> A set of custom EUCAST rules:
#> 
#>   1. If AMC is  R  and genus is "Klebsiella" then set to  R :
#>      amoxicillin (AMX), ampicillin (AMP)
#> 
#>   2. If AMC is  I  and genus is "Klebsiella" then set to  I :
#>      amoxicillin (AMX), ampicillin (AMP)

# run the custom rule set (verbose = TRUE will return a logbook instead of the data set):
eucast_rules(example_isolates,
  rules = "custom",
  custom_rules = x,
  info = FALSE,
  overwrite = TRUE,
  verbose = TRUE
)
#> # A tibble: 8 × 9
#>     row col   mo_fullname     old   new   rule  rule_group rule_name rule_source
#>   <int> <chr> <chr>           <ord> <chr> <chr> <chr>      <chr>     <chr>      
#> 1    33 AMP   Klebsiella pne… R     I     "rep… Custom EU… Custom E… Object 'x'…
#> 2    33 AMX   Klebsiella pne… R     I     "rep… Custom EU… Custom E… Object 'x'…
#> 3    34 AMP   Klebsiella pne… R     I     "rep… Custom EU… Custom E… Object 'x'…
#> 4    34 AMX   Klebsiella pne… R     I     "rep… Custom EU… Custom E… Object 'x'…
#> 5   531 AMP   Klebsiella pne… R     I     "rep… Custom EU… Custom E… Object 'x'…
#> 6   531 AMX   Klebsiella pne… R     I     "rep… Custom EU… Custom E… Object 'x'…
#> 7  1485 AMP   Klebsiella oxy… R     I     "rep… Custom EU… Custom E… Object 'x'…
#> 8  1485 AMX   Klebsiella oxy… R     I     "rep… Custom EU… Custom E… Object 'x'…

# combine rule sets
x2 <- c(
  x,
  custom_eucast_rules(TZP == "R" ~ carbapenems == "R")
)
x2
#> A set of custom EUCAST rules:
#> 
#>   1. If AMC is  R  and genus is "Klebsiella" then set to  R :
#>      amoxicillin (AMX), ampicillin (AMP)
#> 
#>   2. If AMC is  I  and genus is "Klebsiella" then set to  I :
#>      amoxicillin (AMX), ampicillin (AMP)
#> 
#>   3. If TZP is  R  then set to  R :
#>      biapenem (BIA), doripenem (DOR), ertapenem (ETP), imipenem (IPM),
#>      imipenem/EDTA (IPE), imipenem/relebactam (IMR), meropenem (MEM),
#>      meropenem/nacubactam (MNC), meropenem/vaborbactam (MEV), panipenem (PAN),
#>      razupenem (RZM), ritipenem (RIT), ritipenem acoxil (RIA), taniborbactam
#>      (TAN), tebipenem (TBP)
```
