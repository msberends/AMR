% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/custom_mdro_guideline.R
\name{custom_mdro_guideline}
\alias{custom_mdro_guideline}
\alias{c.custom_mdro_guideline}
\title{Define Custom MDRO Guideline}
\usage{
custom_mdro_guideline(..., as_factor = TRUE)

\method{c}{custom_mdro_guideline}(x, ..., as_factor = NULL)
}
\arguments{
\item{...}{Guideline rules in \link[base:tilde]{formula} notation, see below for instructions, and in \emph{Examples}.}

\item{as_factor}{A \link{logical} to indicate whether the returned value should be an ordered \link{factor} (\code{TRUE}, default), or otherwise a \link{character} vector. For combining rules sets (using \code{\link[=c]{c()}}) this value will be inherited from the first set at default.}

\item{x}{Existing custom MDRO rules}
}
\value{
A \link{list} containing the custom rules
}
\description{
Define custom a MDRO guideline for your organisation or specific analysis and use the output of this function in \code{\link[=mdro]{mdro()}}.
}
\details{
Using a custom MDRO guideline is of importance if you have custom rules to determine MDROs in your hospital, e.g., rules that are dependent on ward, state of contact isolation or other variables in your data.
\subsection{Basics}{

If you are familiar with the \code{\link[dplyr:case_when]{case_when()}} function of the \code{dplyr} package, you will recognise the input method to set your own rules. Rules must be set using what \R considers to be the 'formula notation'. The rule itself is written \emph{before} the tilde (\code{~}) and the consequence of the rule is written \emph{after} the tilde:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{custom <- custom_mdro_guideline(CIP == "R" & age > 60 ~ "Elderly Type A",
                                ERY == "R" & age > 60 ~ "Elderly Type B")
}\if{html}{\out{</div>}}

If a row/an isolate matches the first rule, the value after the first \code{~} (in this case \emph{'Elderly Type A'}) will be set as MDRO value. Otherwise, the second rule will be tried and so on. The number of rules is unlimited.

You can print the rules set in the console for an overview. Colours will help reading it if your console supports colours.

\if{html}{\out{<div class="sourceCode r">}}\preformatted{custom
#> A set of custom MDRO rules:
#>   1. If CIP is R and age is higher than 60 then: Elderly Type A
#>   2. If ERY is R and age is higher than 60 then: Elderly Type B
#>   3. Otherwise: Negative

#> Unmatched rows will return NA.
#> Results will be of class 'factor', with ordered levels: Negative < Elderly Type A < Elderly Type B
}\if{html}{\out{</div>}}

The outcome of the function can be used for the \code{guideline} argument in the \code{\link[=mdro]{mdro()}} function:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{x <- mdro(example_isolates, guideline = custom)
#> Determining MDROs based on custom rules, resulting in factor levels: Negative < Elderly Type A < Elderly Type B.
#> - Custom MDRO rule 1: CIP == "R" & age > 60 (198 rows matched)
#> - Custom MDRO rule 2: ERY == "R" & age > 60 (732 rows matched)
#> => Found 930 custom defined MDROs out of 2000 isolates (46.5\%)

table(x)
#> x
#>       Negative  Elderly Type A  Elderly Type B
#>           1070             198             732
}\if{html}{\out{</div>}}

Rules can also be combined with other custom rules by using \code{\link[=c]{c()}}:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{x <- mdro(example_isolates,
          guideline = c(custom,
                        custom_mdro_guideline(ERY == "R" & age > 50 ~ "Elderly Type C")))
#> Determining MDROs based on custom rules, resulting in factor levels: Negative < Elderly Type A < Elderly Type B < Elderly Type C.
#> - Custom MDRO rule 1: CIP == "R" & age > 60 (198 rows matched)
#> - Custom MDRO rule 2: ERY == "R" & age > 60 (732 rows matched)
#> - Custom MDRO rule 3: ERY == "R" & age > 50 (109 rows matched)
#> => Found 1039 custom defined MDROs out of 2000 isolates (52.0\%)

table(x)
#> x
#>       Negative  Elderly Type A  Elderly Type B  Elderly Type C
#>            961             198             732             109
}\if{html}{\out{</div>}}
}

\subsection{Sharing rules among multiple users}{

The rules set (the \code{custom} object in this case) could be exported to a shared file location using \code{\link[=saveRDS]{saveRDS()}} if you collaborate with multiple users. The custom rules set could then be imported using \code{\link[=readRDS]{readRDS()}}.
}

\subsection{Usage of multiple antimicrobials and antimicrobial group names}{

You can define antimicrobial groups instead of single antimicrobials for the rule itself, which is the part \emph{before} the tilde (~). Use \code{\link[=any]{any()}} or \code{\link[=all]{all()}} to specify the scope of the antimicrobial group:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{custom_mdro_guideline(
  AMX == "R"                       ~ "My MDRO #1",
  any(cephalosporins_2nd() == "R") ~ "My MDRO #2",
  all(glycopeptides() == "R")      ~ "My MDRO #3"
)
}\if{html}{\out{</div>}}

All 35 antimicrobial selectors are supported for use in the rules:
\itemize{
\item \code{\link[=aminoglycosides]{aminoglycosides()}} can select: \cr amikacin, amikacin/fosfomycin, apramycin, arbekacin, astromicin, bekanamycin, dibekacin, framycetin, gentamicin, gentamicin-high, habekacin, hygromycin, isepamicin, kanamycin, kanamycin-high, kanamycin/cephalexin, micronomicin, neomycin, netilmicin, pentisomicin, plazomicin, propikacin, ribostamycin, sisomicin, streptoduocin, streptomycin, streptomycin-high, tobramycin, and tobramycin-high
\item \code{\link[=aminopenicillins]{aminopenicillins()}} can select: \cr amoxicillin and ampicillin
\item \code{\link[=antifungals]{antifungals()}} can select: \cr amorolfine, amphotericin B, amphotericin B-high, anidulafungin, butoconazole, caspofungin, ciclopirox, clotrimazole, econazole, fluconazole, flucytosine, fosfluconazole, griseofulvin, hachimycin, ibrexafungerp, isavuconazole, isoconazole, itraconazole, ketoconazole, manogepix, micafungin, miconazole, nystatin, oteseconazole, pimaricin, posaconazole, rezafungin, ribociclib, sulconazole, terbinafine, terconazole, and voriconazole
\item \code{\link[=antimycobacterials]{antimycobacterials()}} can select: \cr 4-aminosalicylic acid, calcium aminosalicylate, capreomycin, clofazimine, delamanid, enviomycin, ethambutol, ethambutol/isoniazid, ethionamide, isoniazid, isoniazid/sulfamethoxazole/trimethoprim/pyridoxine, morinamide, p-aminosalicylic acid, pretomanid, protionamide, pyrazinamide, rifabutin, rifampicin, rifampicin/ethambutol/isoniazid, rifampicin/isoniazid, rifampicin/pyrazinamide/ethambutol/isoniazid, rifampicin/pyrazinamide/isoniazid, rifamycin, rifapentine, simvastatin/fenofibrate, sodium aminosalicylate, streptomycin/isoniazid, terizidone, thioacetazone, thioacetazone/isoniazid, tiocarlide, and viomycin
\item \code{\link[=betalactams]{betalactams()}} can select: \cr amoxicillin, amoxicillin/clavulanic acid, amoxicillin/sulbactam, ampicillin, ampicillin/sulbactam, apalcillin, aspoxicillin, azidocillin, azlocillin, aztreonam, aztreonam/avibactam, aztreonam/nacubactam, bacampicillin, benzathine benzylpenicillin, benzathine phenoxymethylpenicillin, benzylpenicillin, benzylpenicillin screening test, biapenem, carbenicillin, carindacillin, carumonam, cefacetrile, cefaclor, cefadroxil, cefalexin, cefaloridine, cefalotin, cefamandole, cefapirin, cefatrizine, cefazedone, cefazolin, cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren pivoxil, cefepime, cefepime/amikacin, cefepime/clavulanic acid, cefepime/enmetazobactam, cefepime/nacubactam, cefepime/tazobactam, cefepime/zidebactam, cefetamet, cefetamet pivoxil, cefetecol, cefetrizole, cefiderocol, cefixime, cefmenoxime, cefmetazole, cefodizime, cefonicid, cefoperazone, cefoperazone/sulbactam, ceforanide, cefoselis, cefotaxime, cefotaxime screening test, cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefotetan, cefotiam, cefotiam hexetil, cefovecin, cefoxitin, cefoxitin screening test, cefozopran, cefpimizole, cefpiramide, cefpirome, cefpodoxime, cefpodoxime proxetil, cefpodoxime/clavulanic acid, cefprozil, cefquinome, cefroxadine, cefsulodin, cefsumide, ceftaroline, ceftaroline/avibactam, ceftazidime, ceftazidime/avibactam, ceftazidime/clavulanic acid, cefteram, cefteram pivoxil, ceftezole, ceftibuten, ceftiofur, ceftizoxime, ceftizoxime alapivoxil, ceftobiprole, ceftobiprole medocaril, ceftolozane/tazobactam, ceftriaxone, ceftriaxone/beta-lactamase inhibitor, cefuroxime, cefuroxime axetil, cephradine, ciclacillin, clometocillin, cloxacillin, dicloxacillin, doripenem, epicillin, ertapenem, flucloxacillin, hetacillin, imipenem, imipenem/EDTA, imipenem/relebactam, latamoxef, lenampicillin, loracarbef, mecillinam, meropenem, meropenem/nacubactam, meropenem/vaborbactam, metampicillin, meticillin, mezlocillin, mezlocillin/sulbactam, nafcillin, oxacillin, oxacillin screening test, panipenem, penamecillin, penicillin/novobiocin, penicillin/sulbactam, pheneticillin, phenoxymethylpenicillin, piperacillin, piperacillin/sulbactam, piperacillin/tazobactam, piridicillin, pivampicillin, pivmecillinam, procaine benzylpenicillin, propicillin, razupenem, ritipenem, ritipenem acoxil, sarmoxicillin, sulbenicillin, sultamicillin, talampicillin, tebipenem, temocillin, ticarcillin, ticarcillin/clavulanic acid, and tigemonam
\item \code{\link[=betalactams_with_inhibitor]{betalactams_with_inhibitor()}} can select: \cr amoxicillin/clavulanic acid, amoxicillin/sulbactam, ampicillin/sulbactam, aztreonam/avibactam, aztreonam/nacubactam, cefepime/amikacin, cefepime/clavulanic acid, cefepime/enmetazobactam, cefepime/nacubactam, cefepime/tazobactam, cefepime/zidebactam, cefoperazone/sulbactam, cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefpodoxime/clavulanic acid, ceftaroline/avibactam, ceftazidime/avibactam, ceftazidime/clavulanic acid, ceftolozane/tazobactam, ceftriaxone/beta-lactamase inhibitor, imipenem/relebactam, meropenem/nacubactam, meropenem/vaborbactam, mezlocillin/sulbactam, penicillin/novobiocin, penicillin/sulbactam, piperacillin/sulbactam, piperacillin/tazobactam, and ticarcillin/clavulanic acid
\item \code{\link[=carbapenems]{carbapenems()}} can select: \cr biapenem, doripenem, ertapenem, imipenem, imipenem/EDTA, imipenem/relebactam, meropenem, meropenem/nacubactam, meropenem/vaborbactam, panipenem, razupenem, ritipenem, ritipenem acoxil, and tebipenem
\item \code{\link[=cephalosporins]{cephalosporins()}} can select: \cr cefacetrile, cefaclor, cefadroxil, cefalexin, cefaloridine, cefalotin, cefamandole, cefapirin, cefatrizine, cefazedone, cefazolin, cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren pivoxil, cefepime, cefepime/amikacin, cefepime/clavulanic acid, cefepime/enmetazobactam, cefepime/nacubactam, cefepime/tazobactam, cefepime/zidebactam, cefetamet, cefetamet pivoxil, cefetecol, cefetrizole, cefiderocol, cefixime, cefmenoxime, cefmetazole, cefodizime, cefonicid, cefoperazone, cefoperazone/sulbactam, ceforanide, cefoselis, cefotaxime, cefotaxime screening test, cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefotetan, cefotiam, cefotiam hexetil, cefovecin, cefoxitin, cefoxitin screening test, cefozopran, cefpimizole, cefpiramide, cefpirome, cefpodoxime, cefpodoxime proxetil, cefpodoxime/clavulanic acid, cefprozil, cefquinome, cefroxadine, cefsulodin, cefsumide, ceftaroline, ceftaroline/avibactam, ceftazidime, ceftazidime/avibactam, ceftazidime/clavulanic acid, cefteram, cefteram pivoxil, ceftezole, ceftibuten, ceftiofur, ceftizoxime, ceftizoxime alapivoxil, ceftobiprole, ceftobiprole medocaril, ceftolozane/tazobactam, ceftriaxone, ceftriaxone/beta-lactamase inhibitor, cefuroxime, cefuroxime axetil, cephradine, latamoxef, and loracarbef
\item \code{\link[=cephalosporins_1st]{cephalosporins_1st()}} can select: \cr cefacetrile, cefadroxil, cefalexin, cefaloridine, cefalotin, cefapirin, cefatrizine, cefazedone, cefazolin, cefroxadine, ceftezole, and cephradine
\item \code{\link[=cephalosporins_2nd]{cephalosporins_2nd()}} can select: \cr cefaclor, cefamandole, cefmetazole, cefonicid, ceforanide, cefotetan, cefotiam, cefoxitin, cefoxitin screening test, cefprozil, cefuroxime, cefuroxime axetil, and loracarbef
\item \code{\link[=cephalosporins_3rd]{cephalosporins_3rd()}} can select: \cr cefcapene, cefcapene pivoxil, cefdinir, cefditoren, cefditoren pivoxil, cefetamet, cefetamet pivoxil, cefixime, cefmenoxime, cefodizime, cefoperazone, cefoperazone/sulbactam, cefotaxime, cefotaxime screening test, cefotaxime/clavulanic acid, cefotaxime/sulbactam, cefotiam hexetil, cefovecin, cefpimizole, cefpiramide, cefpodoxime, cefpodoxime proxetil, cefpodoxime/clavulanic acid, cefsulodin, ceftazidime, ceftazidime/avibactam, ceftazidime/clavulanic acid, cefteram, cefteram pivoxil, ceftibuten, ceftiofur, ceftizoxime, ceftizoxime alapivoxil, ceftriaxone, ceftriaxone/beta-lactamase inhibitor, and latamoxef
\item \code{\link[=cephalosporins_4th]{cephalosporins_4th()}} can select: \cr cefepime, cefepime/amikacin, cefepime/clavulanic acid, cefepime/enmetazobactam, cefepime/nacubactam, cefepime/tazobactam, cefepime/zidebactam, cefetecol, cefoselis, cefozopran, cefpirome, and cefquinome
\item \code{\link[=cephalosporins_5th]{cephalosporins_5th()}} can select: \cr ceftaroline, ceftaroline/avibactam, ceftobiprole, ceftobiprole medocaril, and ceftolozane/tazobactam
\item \code{\link[=fluoroquinolones]{fluoroquinolones()}} can select: \cr besifloxacin, ciprofloxacin, ciprofloxacin/metronidazole, ciprofloxacin/ornidazole, ciprofloxacin/tinidazole, clinafloxacin, danofloxacin, delafloxacin, difloxacin, enoxacin, enrofloxacin, finafloxacin, fleroxacin, garenoxacin, gatifloxacin, gemifloxacin, grepafloxacin, lascufloxacin, levofloxacin, levofloxacin/ornidazole, levonadifloxacin, lomefloxacin, marbofloxacin, metioxate, miloxacin, moxifloxacin, nadifloxacin, nemonoxacin, nifuroquine, nitroxoline, norfloxacin, norfloxacin screening test, norfloxacin/metronidazole, norfloxacin/tinidazole, ofloxacin, ofloxacin/ornidazole, orbifloxacin, pazufloxacin, pefloxacin, pefloxacin screening test, pradofloxacin, premafloxacin, prulifloxacin, rufloxacin, sarafloxacin, sitafloxacin, sparfloxacin, temafloxacin, tilbroquinol, tioxacin, tosufloxacin, and trovafloxacin
\item \code{\link[=glycopeptides]{glycopeptides()}} can select: \cr avoparcin, bleomycin, dalbavancin, norvancomycin, oritavancin, ramoplanin, teicoplanin, teicoplanin-macromethod, telavancin, vancomycin, and vancomycin-macromethod
\item \code{\link[=isoxazolylpenicillins]{isoxazolylpenicillins()}} can select: \cr cloxacillin, dicloxacillin, flucloxacillin, meticillin, oxacillin, and oxacillin screening test
\item \code{\link[=lincosamides]{lincosamides()}} can select: \cr clindamycin, lincomycin, and pirlimycin
\item \code{\link[=lipoglycopeptides]{lipoglycopeptides()}} can select: \cr dalbavancin, oritavancin, and telavancin
\item \code{\link[=macrolides]{macrolides()}} can select: \cr acetylmidecamycin, acetylspiramycin, azithromycin, clarithromycin, dirithromycin, erythromycin, flurithromycin, gamithromycin, josamycin, kitasamycin, meleumycin, midecamycin, miocamycin, nafithromycin, oleandomycin, rokitamycin, roxithromycin, solithromycin, spiramycin, telithromycin, tildipirosin, tilmicosin, troleandomycin, tulathromycin, tylosin, and tylvalosin
\item \code{\link[=monobactams]{monobactams()}} can select: \cr aztreonam, aztreonam/avibactam, aztreonam/nacubactam, carumonam, and tigemonam
\item \code{\link[=nitrofurans]{nitrofurans()}} can select: \cr furazidin, furazolidone, nifurtoinol, nitrofurantoin, and nitrofurazone
\item \code{\link[=oxazolidinones]{oxazolidinones()}} can select: \cr cadazolid, cycloserine, linezolid, tedizolid, and thiacetazone
\item \code{\link[=penicillins]{penicillins()}} can select: \cr amoxicillin, amoxicillin/clavulanic acid, amoxicillin/sulbactam, ampicillin, ampicillin/sulbactam, apalcillin, aspoxicillin, azidocillin, azlocillin, bacampicillin, benzathine benzylpenicillin, benzathine phenoxymethylpenicillin, benzylpenicillin, benzylpenicillin screening test, carbenicillin, carindacillin, ciclacillin, clometocillin, cloxacillin, dicloxacillin, epicillin, flucloxacillin, hetacillin, lenampicillin, mecillinam, metampicillin, meticillin, mezlocillin, mezlocillin/sulbactam, nafcillin, oxacillin, oxacillin screening test, penamecillin, penicillin/novobiocin, penicillin/sulbactam, pheneticillin, phenoxymethylpenicillin, piperacillin, piperacillin/sulbactam, piperacillin/tazobactam, piridicillin, pivampicillin, pivmecillinam, procaine benzylpenicillin, propicillin, sarmoxicillin, sulbenicillin, sultamicillin, talampicillin, temocillin, ticarcillin, and ticarcillin/clavulanic acid
\item \code{\link[=phenicols]{phenicols()}} can select: \cr chloramphenicol, florfenicol, and thiamphenicol
\item \code{\link[=polymyxins]{polymyxins()}} can select: \cr colistin, polymyxin B, and polymyxin B/polysorbate 80
\item \code{\link[=quinolones]{quinolones()}} can select: \cr besifloxacin, cinoxacin, ciprofloxacin, ciprofloxacin/metronidazole, ciprofloxacin/ornidazole, ciprofloxacin/tinidazole, clinafloxacin, danofloxacin, delafloxacin, difloxacin, enoxacin, enrofloxacin, finafloxacin, fleroxacin, flumequine, garenoxacin, gatifloxacin, gemifloxacin, grepafloxacin, lascufloxacin, levofloxacin, levofloxacin/ornidazole, levonadifloxacin, lomefloxacin, marbofloxacin, metioxate, miloxacin, moxifloxacin, nadifloxacin, nalidixic acid, nalidixic acid screening test, nemonoxacin, nifuroquine, nitroxoline, norfloxacin, norfloxacin screening test, norfloxacin/metronidazole, norfloxacin/tinidazole, ofloxacin, ofloxacin/ornidazole, orbifloxacin, oxolinic acid, pazufloxacin, pefloxacin, pefloxacin screening test, pipemidic acid, piromidic acid, pradofloxacin, premafloxacin, prulifloxacin, rosoxacin, rufloxacin, sarafloxacin, sitafloxacin, sparfloxacin, temafloxacin, tilbroquinol, tioxacin, tosufloxacin, and trovafloxacin
\item \code{\link[=rifamycins]{rifamycins()}} can select: \cr rifabutin, rifampicin, rifampicin/ethambutol/isoniazid, rifampicin/isoniazid, rifampicin/pyrazinamide/ethambutol/isoniazid, rifampicin/pyrazinamide/isoniazid, rifamycin, and rifapentine
\item \code{\link[=streptogramins]{streptogramins()}} can select: \cr pristinamycin and quinupristin/dalfopristin
\item \code{\link[=sulfonamides]{sulfonamides()}} can select: \cr brodimoprim, sulfadiazine, sulfadiazine/tetroxoprim, sulfadimethoxine, sulfadimidine, sulfafurazole, sulfaisodimidine, sulfalene, sulfamazone, sulfamerazine, sulfamethizole, sulfamethoxazole, sulfamethoxypyridazine, sulfametomidine, sulfametoxydiazine, sulfamoxole, sulfanilamide, sulfaperin, sulfaphenazole, sulfapyridine, sulfathiazole, and sulfathiourea
\item \code{\link[=tetracyclines]{tetracyclines()}} can select: \cr cetocycline, chlortetracycline, clomocycline, demeclocycline, doxycycline, eravacycline, lymecycline, metacycline, minocycline, omadacycline, oxytetracycline, penimepicycline, rolitetracycline, sarecycline, tetracycline, tetracycline screening test, and tigecycline
\item \code{\link[=trimethoprims]{trimethoprims()}} can select: \cr brodimoprim, sulfadiazine, sulfadiazine/tetroxoprim, sulfadiazine/trimethoprim, sulfadimethoxine, sulfadimidine, sulfadimidine/trimethoprim, sulfafurazole, sulfaisodimidine, sulfalene, sulfamazone, sulfamerazine, sulfamerazine/trimethoprim, sulfamethizole, sulfamethoxazole, sulfamethoxypyridazine, sulfametomidine, sulfametoxydiazine, sulfametrole/trimethoprim, sulfamoxole, sulfamoxole/trimethoprim, sulfanilamide, sulfaperin, sulfaphenazole, sulfapyridine, sulfathiazole, sulfathiourea, trimethoprim, and trimethoprim/sulfamethoxazole
\item \code{\link[=ureidopenicillins]{ureidopenicillins()}} can select: \cr azlocillin, mezlocillin, piperacillin, and piperacillin/tazobactam
}
}
}
\examples{
x <- custom_mdro_guideline(
  CIP == "R" & age > 60 ~ "Elderly Type A",
  ERY == "R" & age > 60 ~ "Elderly Type B"
)
x

# run the custom rule set (verbose = TRUE will return a logbook instead of the data set):
out <- mdro(example_isolates, guideline = x)
table(out)

out <- mdro(example_isolates, guideline = x, verbose = TRUE)
head(out)

# you can create custom guidelines using selectors (see ?antimicrobial_selectors)
my_guideline <- custom_mdro_guideline(
  AMX == "R" ~ "Custom MDRO 1",
  all(cephalosporins_2nd() == "R") ~ "Custom MDRO 2"
)
my_guideline

out <- mdro(example_isolates, guideline = my_guideline)
table(out)
}
