# Define Custom EUCAST Rules

Define custom EUCAST rules for your organisation or specific analysis
and use the output of this function in
[`eucast_rules()`](https://amr-for-r.org/release/2.0.0/reference/eucast_rules.md).

## Usage

``` r
custom_eucast_rules(...)
```

## Arguments

- ...:

  rules in [formula](https://rdrr.io/r/base/tilde.html) notation, see
  *Examples*

## Value

A [list](https://rdrr.io/r/base/list.html) containing the custom rules

## Details

Some organisations have their own adoption of EUCAST rules. This
function can be used to define custom EUCAST rules to be used in the
[`eucast_rules()`](https://amr-for-r.org/release/2.0.0/reference/eucast_rules.md)
function.

## How it works

### Basics

If you are familiar with the
[`case_when()`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html)
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

    eucast_rules(df, rules = "custom", custom_rules = x, info = FALSE)
    #>                      mo TZP ampi cipro
    #> 1      Escherichia coli   R    R     S
    #> 2 Klebsiella pneumoniae   R    R     S

### Using taxonomic properties in rules

There is one exception in columns used for the rules: all column names
of the
[microorganisms](https://amr-for-r.org/release/2.0.0/reference/microorganisms.md)
data set can also be used, but do not have to exist in the data set.
These column names are: "mo", "fullname", "status", "kingdom", "phylum",
"class", "order", "family", "genus", "species", "subspecies", "rank",
"ref", "source", "lpsn", "lpsn_parent", "lpsn_renamed_to", "gbif",
"gbif_parent", "gbif_renamed_to", "prevalence", and "snomed". Thus, this
next example will work as well, despite the fact that the `df` data set
does not contain a column `genus`:

    y <- custom_eucast_rules(TZP == "S" & genus == "Klebsiella" ~ aminopenicillins == "S",
                             TZP == "R" & genus == "Klebsiella" ~ aminopenicillins == "R")

    eucast_rules(df, rules = "custom", custom_rules = y, info = FALSE)
    #>                      mo TZP ampi cipro
    #> 1      Escherichia coli   R    S     S
    #> 2 Klebsiella pneumoniae   R    R     S

### Usage of antibiotic group names

It is possible to define antibiotic groups instead of single antibiotics
for the rule consequence, the part *after* the tilde. In above examples,
the antibiotic group `aminopenicillins` is used to include ampicillin
and amoxicillin. The following groups are allowed (case-insensitive).
Within parentheses are the drugs that will be matched when running the
rule.

- "aminoglycosides"  
  (amikacin, amikacin/fosfomycin, amphotericin B-high, apramycin,
  arbekacin, astromicin, bekanamycin, dibekacin, framycetin, gentamicin,
  gentamicin-high, habekacin, hygromycin, isepamicin, kanamycin,
  kanamycin-high, kanamycin/cephalexin, micronomicin, neomycin,
  netilmicin, pentisomicin, plazomicin, propikacin, ribostamycin,
  sisomicin, streptoduocin, streptomycin, streptomycin-high, tobramycin,
  and tobramycin-high)

- "aminopenicillins"  
  (amoxicillin and ampicillin)

- "antifungals"  
  (amphotericin B, anidulafungin, butoconazole, caspofungin, ciclopirox,
  clotrimazole, econazole, fluconazole, flucytosine, fosfluconazole,
  griseofulvin, hachimycin, ibrexafungerp, isavuconazole, isoconazole,
  itraconazole, ketoconazole, manogepix, micafungin, miconazole,
  nystatin, oteseconazole, pimaricin, posaconazole, rezafungin,
  ribociclib, sulconazole, terbinafine, terconazole, and voriconazole)

- "antimycobacterials"  
  (4-aminosalicylic acid, calcium aminosalicylate, capreomycin,
  clofazimine, delamanid, enviomycin, ethambutol, ethambutol/isoniazid,
  ethionamide, isoniazid,
  isoniazid/sulfamethoxazole/trimethoprim/pyridoxine, morinamide,
  p-aminosalicylic acid, pretomanid, protionamide, pyrazinamide,
  rifabutin, rifampicin, rifampicin/ethambutol/isoniazid,
  rifampicin/isoniazid, rifampicin/pyrazinamide/ethambutol/isoniazid,
  rifampicin/pyrazinamide/isoniazid, rifamycin, rifapentine,
  simvastatin/fenofibrate, sodium aminosalicylate,
  streptomycin/isoniazid, terizidone, thioacetazone,
  thioacetazone/isoniazid, tiocarlide, and viomycin)

- "betalactams"  
  (amoxicillin, amoxicillin/clavulanic acid, amoxicillin/sulbactam,
  ampicillin, ampicillin/sulbactam, apalcillin, aspoxicillin, avibactam,
  azidocillin, azlocillin, aztreonam, aztreonam/avibactam,
  aztreonam/nacubactam, bacampicillin, benzathine benzylpenicillin,
  benzathine phenoxymethylpenicillin, benzylpenicillin, biapenem,
  carbenicillin, carindacillin, cefacetrile, cefaclor, cefadroxil,
  cefalexin, cefaloridine, cefalotin, cefamandole, cefapirin,
  cefatrizine, cefazedone, cefazolin, cefcapene, cefcapene pivoxil,
  cefdinir, cefditoren, cefditoren pivoxil, cefepime,
  cefepime/clavulanic acid, cefepime/nacubactam, cefepime/tazobactam,
  cefetamet, cefetamet pivoxil, cefetecol, cefetrizole, cefixime,
  cefmenoxime, cefmetazole, cefodizime, cefonicid, cefoperazone,
  cefoperazone/sulbactam, ceforanide, cefoselis, cefotaxime,
  cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefotetan, cefotiam,
  cefotiam hexetil, cefovecin, cefoxitin, cefoxitin screening,
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
  meticillin, mezlocillin, mezlocillin/sulbactam, nacubactam, nafcillin,
  oxacillin, panipenem, penamecillin, penicillin/novobiocin,
  penicillin/sulbactam, pheneticillin, phenoxymethylpenicillin,
  piperacillin, piperacillin/sulbactam, piperacillin/tazobactam,
  piridicillin, pivampicillin, pivmecillinam, procaine benzylpenicillin,
  propicillin, razupenem, ritipenem, ritipenem acoxil, sarmoxicillin,
  sulbactam, sulbenicillin, sultamicillin, talampicillin, tazobactam,
  tebipenem, temocillin, ticarcillin, and ticarcillin/clavulanic acid)

- "carbapenems"  
  (biapenem, doripenem, ertapenem, imipenem, imipenem/EDTA,
  imipenem/relebactam, meropenem, meropenem/nacubactam,
  meropenem/vaborbactam, panipenem, razupenem, ritipenem, ritipenem
  acoxil, and tebipenem)

- "cephalosporins"  
  (cefacetrile, cefaclor, cefadroxil, cefalexin, cefaloridine,
  cefalotin, cefamandole, cefapirin, cefatrizine, cefazedone, cefazolin,
  cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren
  pivoxil, cefepime, cefepime/clavulanic acid, cefepime/tazobactam,
  cefetamet, cefetamet pivoxil, cefetecol, cefetrizole, cefixime,
  cefmenoxime, cefmetazole, cefodizime, cefonicid, cefoperazone,
  cefoperazone/sulbactam, ceforanide, cefoselis, cefotaxime,
  cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefotetan, cefotiam,
  cefotiam hexetil, cefovecin, cefoxitin, cefoxitin screening,
  cefozopran, cefpimizole, cefpiramide, cefpirome, cefpodoxime,
  cefpodoxime proxetil, cefpodoxime/clavulanic acid, cefprozil,
  cefquinome, cefroxadine, cefsulodin, cefsumide, ceftaroline,
  ceftaroline/avibactam, ceftazidime, ceftazidime/avibactam,
  ceftazidime/clavulanic acid, cefteram, cefteram pivoxil, ceftezole,
  ceftibuten, ceftiofur, ceftizoxime, ceftizoxime alapivoxil,
  ceftobiprole, ceftobiprole medocaril, ceftolozane/tazobactam,
  ceftriaxone, ceftriaxone/beta-lactamase inhibitor, cefuroxime,
  cefuroxime axetil, cephradine, latamoxef, and loracarbef)

- "cephalosporins_1st"  
  (cefacetrile, cefadroxil, cefalexin, cefaloridine, cefalotin,
  cefapirin, cefatrizine, cefazedone, cefazolin, cefroxadine, ceftezole,
  and cephradine)

- "cephalosporins_2nd"  
  (cefaclor, cefamandole, cefmetazole, cefonicid, ceforanide, cefotetan,
  cefotiam, cefoxitin, cefoxitin screening, cefprozil, cefuroxime,
  cefuroxime axetil, and loracarbef)

- "cephalosporins_3rd"  
  (cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren
  pivoxil, cefetamet, cefetamet pivoxil, cefixime, cefmenoxime,
  cefodizime, cefoperazone, cefoperazone/sulbactam, cefotaxime,
  cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefotiam hexetil,
  cefovecin, cefpimizole, cefpiramide, cefpodoxime, cefpodoxime
  proxetil, cefpodoxime/clavulanic acid, cefsulodin, ceftazidime,
  ceftazidime/avibactam, ceftazidime/clavulanic acid, cefteram, cefteram
  pivoxil, ceftibuten, ceftiofur, ceftizoxime, ceftizoxime alapivoxil,
  ceftriaxone, ceftriaxone/beta-lactamase inhibitor, and latamoxef)

- "cephalosporins_4th"  
  (cefepime, cefepime/clavulanic acid, cefepime/tazobactam, cefetecol,
  cefoselis, cefozopran, cefpirome, and cefquinome)

- "cephalosporins_5th"  
  (ceftaroline, ceftaroline/avibactam, ceftobiprole, ceftobiprole
  medocaril, and ceftolozane/tazobactam)

- "cephalosporins_except_caz"  
  (cefacetrile, cefaclor, cefadroxil, cefalexin, cefaloridine,
  cefalotin, cefamandole, cefapirin, cefatrizine, cefazedone, cefazolin,
  cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren
  pivoxil, cefepime, cefepime/clavulanic acid, cefepime/tazobactam,
  cefetamet, cefetamet pivoxil, cefetecol, cefetrizole, cefixime,
  cefmenoxime, cefmetazole, cefodizime, cefonicid, cefoperazone,
  cefoperazone/sulbactam, ceforanide, cefoselis, cefotaxime,
  cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefotetan, cefotiam,
  cefotiam hexetil, cefovecin, cefoxitin, cefoxitin screening,
  cefozopran, cefpimizole, cefpiramide, cefpirome, cefpodoxime,
  cefpodoxime proxetil, cefpodoxime/clavulanic acid, cefprozil,
  cefquinome, cefroxadine, cefsulodin, cefsumide, ceftaroline,
  ceftaroline/avibactam, ceftazidime/avibactam, ceftazidime/clavulanic
  acid, cefteram, cefteram pivoxil, ceftezole, ceftibuten, ceftiofur,
  ceftizoxime, ceftizoxime alapivoxil, ceftobiprole, ceftobiprole
  medocaril, ceftolozane/tazobactam, ceftriaxone,
  ceftriaxone/beta-lactamase inhibitor, cefuroxime, cefuroxime axetil,
  cephradine, latamoxef, and loracarbef)

- "fluoroquinolones"  
  (besifloxacin, ciprofloxacin, clinafloxacin, danofloxacin,
  delafloxacin, difloxacin, enoxacin, enrofloxacin, finafloxacin,
  fleroxacin, garenoxacin, gatifloxacin, gemifloxacin, grepafloxacin,
  lascufloxacin, levofloxacin, levonadifloxacin, lomefloxacin,
  marbofloxacin, metioxate, miloxacin, moxifloxacin, nadifloxacin,
  nifuroquine, norfloxacin, ofloxacin, orbifloxacin, pazufloxacin,
  pefloxacin, pradofloxacin, premafloxacin, prulifloxacin, rufloxacin,
  sarafloxacin, sitafloxacin, sparfloxacin, temafloxacin, tilbroquinol,
  tioxacin, tosufloxacin, and trovafloxacin)

- "glycopeptides"  
  (avoparcin, dalbavancin, norvancomycin, oritavancin, ramoplanin,
  teicoplanin, teicoplanin-macromethod, telavancin, vancomycin, and
  vancomycin-macromethod)

- "glycopeptides_except_lipo"  
  (avoparcin, norvancomycin, ramoplanin, teicoplanin,
  teicoplanin-macromethod, vancomycin, and vancomycin-macromethod)

- "lincosamides"  
  (acetylmidecamycin, acetylspiramycin, clindamycin, gamithromycin,
  kitasamycin, lincomycin, meleumycin, nafithromycin, pirlimycin,
  primycin, solithromycin, tildipirosin, tilmicosin, tulathromycin,
  tylosin, and tylvalosin)

- "lipoglycopeptides"  
  (dalbavancin, oritavancin, and telavancin)

- "macrolides"  
  (acetylmidecamycin, acetylspiramycin, azithromycin, clarithromycin,
  dirithromycin, erythromycin, flurithromycin, gamithromycin, josamycin,
  kitasamycin, meleumycin, midecamycin, miocamycin, nafithromycin,
  oleandomycin, pirlimycin, primycin, rokitamycin, roxithromycin,
  solithromycin, spiramycin, telithromycin, tildipirosin, tilmicosin,
  troleandomycin, tulathromycin, tylosin, and tylvalosin)

- "oxazolidinones"  
  (cadazolid, cycloserine, linezolid, tedizolid, and thiacetazone)

- "penicillins"  
  (amoxicillin, amoxicillin/clavulanic acid, amoxicillin/sulbactam,
  ampicillin, ampicillin/sulbactam, apalcillin, aspoxicillin, avibactam,
  azidocillin, azlocillin, aztreonam, aztreonam/avibactam,
  aztreonam/nacubactam, bacampicillin, benzathine benzylpenicillin,
  benzathine phenoxymethylpenicillin, benzylpenicillin, carbenicillin,
  carindacillin, cefepime/nacubactam, ciclacillin, clometocillin,
  cloxacillin, dicloxacillin, epicillin, flucloxacillin, hetacillin,
  lenampicillin, mecillinam, metampicillin, meticillin, mezlocillin,
  mezlocillin/sulbactam, nacubactam, nafcillin, oxacillin, penamecillin,
  penicillin/novobiocin, penicillin/sulbactam, pheneticillin,
  phenoxymethylpenicillin, piperacillin, piperacillin/sulbactam,
  piperacillin/tazobactam, piridicillin, pivampicillin, pivmecillinam,
  procaine benzylpenicillin, propicillin, sarmoxicillin, sulbactam,
  sulbenicillin, sultamicillin, talampicillin, tazobactam, temocillin,
  ticarcillin, and ticarcillin/clavulanic acid)

- "polymyxins"  
  (colistin, polymyxin B, and polymyxin B/polysorbate 80)

- "quinolones"  
  (besifloxacin, cinoxacin, ciprofloxacin, clinafloxacin, danofloxacin,
  delafloxacin, difloxacin, enoxacin, enrofloxacin, finafloxacin,
  fleroxacin, flumequine, garenoxacin, gatifloxacin, gemifloxacin,
  grepafloxacin, lascufloxacin, levofloxacin, levonadifloxacin,
  lomefloxacin, marbofloxacin, metioxate, miloxacin, moxifloxacin,
  nadifloxacin, nalidixic acid, nemonoxacin, nifuroquine, nitroxoline,
  norfloxacin, ofloxacin, orbifloxacin, oxolinic acid, pazufloxacin,
  pefloxacin, pipemidic acid, piromidic acid, pradofloxacin,
  premafloxacin, prulifloxacin, rosoxacin, rufloxacin, sarafloxacin,
  sitafloxacin, sparfloxacin, temafloxacin, tilbroquinol, tioxacin,
  tosufloxacin, and trovafloxacin)

- "streptogramins"  
  (pristinamycin and quinupristin/dalfopristin)

- "tetracyclines"  
  (cetocycline, chlortetracycline, clomocycline, demeclocycline,
  doxycycline, eravacycline, lymecycline, metacycline, minocycline,
  omadacycline, oxytetracycline, penimepicycline, rolitetracycline,
  sarecycline, tetracycline, and tigecycline)

- "tetracyclines_except_tgc"  
  (cetocycline, chlortetracycline, clomocycline, demeclocycline,
  doxycycline, eravacycline, lymecycline, metacycline, minocycline,
  omadacycline, oxytetracycline, penimepicycline, rolitetracycline,
  sarecycline, and tetracycline)

- "trimethoprims"  
  (brodimoprim, sulfadiazine, sulfadiazine/tetroxoprim,
  sulfadiazine/trimethoprim, sulfadimethoxine, sulfadimidine,
  sulfadimidine/trimethoprim, sulfafurazole, sulfaisodimidine,
  sulfalene, sulfamazone, sulfamerazine, sulfamerazine/trimethoprim,
  sulfamethizole, sulfamethoxazole, sulfamethoxypyridazine,
  sulfametomidine, sulfametoxydiazine, sulfametrole/trimethoprim,
  sulfamoxole, sulfamoxole/trimethoprim, sulfanilamide, sulfaperin,
  sulfaphenazole, sulfapyridine, sulfathiazole, sulfathiourea,
  trimethoprim, and trimethoprim/sulfamethoxazole)

- "ureidopenicillins"  
  (azlocillin, mezlocillin, piperacillin, and piperacillin/tazobactam)

## Examples

``` r
x <- custom_eucast_rules(
  AMC == "R" & genus == "Klebsiella" ~ aminopenicillins == "R",
  AMC == "I" & genus == "Klebsiella" ~ aminopenicillins == "I"
)
x

# run the custom rule set (verbose = TRUE will return a logbook instead of the data set):
eucast_rules(example_isolates,
  rules = "custom",
  custom_rules = x,
  info = FALSE,
  verbose = TRUE
)

# combine rule sets
x2 <- c(
  x,
  custom_eucast_rules(TZP == "R" ~ carbapenems == "R")
)
x2
```
