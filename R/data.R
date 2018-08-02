# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# AUTHORS                                                              #
# Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)           #
#                                                                      #
# LICENCE                                                              #
# This program is free software; you can redistribute it and/or modify #
# it under the terms of the GNU General Public License version 2.0,    #
# as published by the Free Software Foundation.                        #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
# ==================================================================== #

#' Dataset with 420 antibiotics
#'
#' A dataset containing all antibiotics with a J0 code, with their DDD's. Properties were downloaded from the WHO, see Source.
#' @format A data.frame with 420 observations and 18 variables:
#' \describe{
#'   \item{\code{atc}}{ATC code, like \code{J01CR02}}
#'   \item{\code{molis}}{MOLIS code, like \code{amcl}}
#'   \item{\code{umcg}}{UMCG code, like \code{AMCL}}
#'   \item{\code{abbr}}{Abbreviation as used by many countries, to be used for \code{\link{guess_atc}}}
#'   \item{\code{official}}{Official name by the WHO, like \code{"Amoxicillin and enzyme inhibitor"}}
#'   \item{\code{official_nl}}{Official name in the Netherlands, like \code{"Amoxicilline met enzymremmer"}}
#'   \item{\code{trivial_nl}}{Trivial name in Dutch, like \code{"Amoxicilline/clavulaanzuur"}}
#'   \item{\code{trade_name}}{Trade name as used by many countries, to be used for \code{\link{guess_atc}}}
#'   \item{\code{oral_ddd}}{Defined Daily Dose (DDD), oral treatment}
#'   \item{\code{oral_units}}{Units of \code{ddd_units}}
#'   \item{\code{iv_ddd}}{Defined Daily Dose (DDD), parenteral treatment}
#'   \item{\code{iv_units}}{Units of \code{iv_ddd}}
#'   \item{\code{atc_group1}}{ATC group, like \code{"Macrolides, lincosamides and streptogramins"}}
#'   \item{\code{atc_group2}}{Subgroup of \code{atc_group1}, like \code{"Macrolides"}}
#'   \item{\code{atc_group1_nl}}{ATC group in Dutch, like \code{"Macroliden, lincosamiden en streptograminen"}}
#'   \item{\code{atc_group2_nl}}{Subgroup of \code{atc_group1} in Dutch, like \code{"Macroliden"}}
#'   \item{\code{useful_gramnegative}}{\code{FALSE} if not useful according to EUCAST, \code{NA} otherwise (see Source)}
#'   \item{\code{useful_grampositive}}{\code{FALSE} if not useful according to EUCAST, \code{NA} otherwise (see Source)}
#' }
#' @source - World Health Organization: \url{https://www.whocc.no/atc_ddd_index/} \cr - EUCAST - Expert rules intrinsic exceptional V3.1 \cr - MOLIS (LIS of Certe): \url{https://www.certe.nl} \cr - GLIMS (LIS of UMCG): \url{https://www.umcg.nl}
#' @seealso \code{\link{microorganisms}}
# abbr and trade_name created with:
# https://hs.unr.edu/Documents/dhs/chs/NVPHTC/antibiotic_refeference_guide.pdf
# antibiotics %>%
#   mutate(abbr =
#            case_when(
#              official == 'Amikacin' ~ 'Ak|AN|AMI|AMK',
#              official == 'Amoxicillin' ~ 'AMX|AMOX|AC',
#              official == 'Amoxicillin and beta-lactamase inhibitor' ~ 'AUG|A/C|XL|AML',
#              official == 'Ampicillin' ~ 'AM|AMP',
#              official == 'Ampicillin and beta-lactamase inhibitor' ~ 'A/S|SAM|AMS|AB',
#              official == 'Azithromycin' ~ 'Azi|AZM|AZ',
#              official == 'Azlocillin' ~ 'AZ|AZL',
#              official == 'Aztreonam' ~ 'Azt|ATM|AT|AZM',
#              official == 'Carbenicillin' ~ 'Cb|BAR',
#              official == 'Cefaclor' ~ 'Ccl|CEC|Cfr|FAC|CF',
#              official == 'Cefadroxil' ~ 'CFR|FAD',
#              official == 'Cefazolin' ~ 'Cfz|CZ|FAZ|KZ',
#              official == 'Cefdinir' ~ 'Cdn|CDR|DIN|CD|CFD',
#              official == 'Cefditoren' ~ 'CDN',
#              official == 'Cefepime' ~ 'Cpe|FEP|PM|CPM',
#              official == 'Cefixime' ~ 'Cfe|DCFM|FIX|IX',
#              official == 'Cefoperazone' ~ 'Cfp|CPZ|PER|FOP|CP',
#              official == 'Cefotaxime' ~ 'Cft|CTX|TAX|FOT|CT',
#              official == 'Cefotetan' ~ 'Ctn|CTT|CTE|TANS|CN',
#              official == 'Cefoxitin' ~ 'Cfx|FOX|CX|FX',
#              official == 'Cefpodoxime' ~ 'Cpd|POD|PX',
#              official == 'Cefprozil' ~ 'Cpz|CPR|FP',
#              official == 'Ceftaroline' ~ 'CPT',
#              official == 'Ceftazidime' ~ 'Caz|TAZ|TZ',
#              official == 'Ceftibuten' ~ 'CTB|TIB|CB',
#              official == 'Ceftizoxime' ~ 'Cz|ZOX|CZX|CZ|CTZ|TIZ',
#              official == 'Ceftriaxone' ~ 'Cax|CRO|CTR|FRX|AXO|TX',
#              official == 'Cefuroxime' ~ 'Crm|CXM|CFX|ROX|FUR|XM',
#              official == 'Cephalexin' ~ 'CN|LX|CFL',
#              official == 'Cephalothin' ~ 'Cf',
#              official == 'Chloramphenicol' ~ 'C|CHL|CL',
#              official == 'Ciprofloxacin' ~ 'Cp|CIP|CI',
#              official == 'Clarithromycin' ~ 'Cla|CLR|CLM|CH',
#              official == 'Clindamycin' ~ 'Cd|CC|CM|CLI|DA',
#              official == 'Colistin' ~ 'CL|CS|CT',
#              official == 'Daptomycin' ~ 'Dap',
#              official == 'Doxycycline' ~ 'Dox',
#              official == 'Doripenem' ~ 'DOR|Dor',
#              official == 'Ertapenem' ~ 'Etp',
#              official == 'Erythromycin' ~ 'E|ERY|EM',
#              official == 'Fosfomycin' ~ 'FOS|FF|FO|FM',
#              official == 'Flucloxacillin' ~ 'CLOX',
#              official == 'Gentamicin' ~ 'Gm|CN|GEN',
#              official == 'Imipenem' ~ 'Imp|IPM|IMI|IP',
#              official == 'Kanamycin' ~ 'K|KAN|HLK|KM',
#              official == 'Levofloxacin' ~ 'Lvx|LEV|LEVO|LE',
#              official == 'Linezolid' ~ 'Lzd|LNZ|LZ',
#              official == 'Lomefloxacin' ~ 'Lmf|LOM',
#              official == 'Meropenem' ~ 'Mer|MEM|MERO|MRP|MP',
#              official == 'Metronidazole' ~ 'MNZ',
#              official == 'Mezlocillin' ~ 'Mz|MEZ',
#              official == 'Minocycline' ~ 'Min|MI|MN|MNO|MC|MH',
#              official == 'Moxifloxacin' ~ 'Mox|MXF',
#              official == 'Mupirocin' ~ 'MUP',
#              official == 'Nafcillin' ~ 'Naf|NF',
#              official == 'Nalidixic acid' ~ 'NA|NAL',
#              official == 'Nitrofurantoin' ~ 'Fd|F/M|FT|NIT|NI|F',
#              official == 'Norfloxacin' ~ 'Nxn|NOR|NX',
#              official == 'Ofloxacin' ~ 'Ofl|OFX|OF',
#              official == 'Oxacillin' ~ 'Ox|OXS|OXA',
#              official == 'Benzylpenicillin' ~ 'P|PEN|PV',
#              official == 'Penicillins, combinations with other antibacterials' ~ 'P|PEN|PV',
#              official == 'Piperacillin' ~ 'Pi|PIP|PP',
#              official == 'Piperacillin and beta-lactamase inhibitor' ~ 'PT|TZP|PTZ|P/T|PTc',
#              official == 'Polymyxin B' ~ 'PB',
#              official == 'Quinupristin/dalfopristin' ~ 'Syn|Q/D|QDA|RP',
#              official == 'Rifampin' ~ 'Rif|RA|RI|RD',
#              official == 'Spectinomycin' ~ 'SPT|SPE|SC',
#              official == 'Streptomycin' ~ 'S|STR',
#              official == 'Teicoplanin' ~ 'Tei|TEC|TPN|TP|TPL',
#              official == 'Telavancin' ~ 'TLV',
#              official == 'Telithromcyin' ~ 'Tel',
#              official == 'Tetracycline' ~ 'Te|TET|TC',
#              official == 'Ticarcillin' ~ 'Ti|TIC|TC',
#              official == 'Ticarcillin and beta-lactamase inhibitor' ~ 'Tim|T/C|TCC|TLc',
#              official == 'Tigecycline' ~ 'TGC',
#              official == 'Tobramycin' ~ 'To|NN|TM|TOB',
#              official == 'Trimethoprim' ~ 'T|TMP|TR|W',
#              official == 'Sulfamethoxazole and trimethoprim' ~ 'T/S|SXT|SxT|TS|COT',
#              official == 'Vancomycin' ~ 'Va|VAN',
#              TRUE ~ NA_character_),
#
#          trade_name =
#            case_when(
#              official == 'Amikacin' ~ 'Amikin',
#              official == 'Amoxicillin' ~ 'Amoxil|Dispermox|Larotid|Trimox',
#              official == 'Amoxicillin and beta-lactamase inhibitor' ~ 'Augmentin',
#              official == 'Ampicillin' ~ 'Pfizerpen-A|Principen',
#              official == 'Ampicillin and beta-lactamase inhibitor' ~ 'Unasyn',
#              official == 'Azithromycin' ~ 'Zithromax',
#              official == 'Azlocillin' ~ 'Azlin',
#              official == 'Aztreonam' ~ 'Azactam',
#              official == 'Carbenicillin' ~ 'Geocillin',
#              official == 'Cefaclor' ~ 'Ceclor',
#              official == 'Cefadroxil' ~ 'Duricef',
#              official == 'Cefazolin' ~ 'Ancef',
#              official == 'Cefdinir' ~ 'Omnicef',
#              official == 'Cefditoren' ~ 'Spectracef',
#              official == 'Cefepime' ~ 'Maxipime',
#              official == 'Cefixime' ~ 'Suprax',
#              official == 'Cefoperazone' ~ 'Cefobid',
#              official == 'Cefotaxime' ~ 'Claforan',
#              official == 'Cefotetan' ~ 'Cefotan',
#              official == 'Cefoxitin' ~ 'Mefoxin',
#              official == 'Cefpodoxime' ~ 'Vantin',
#              official == 'Cefprozil' ~ 'Cefzil',
#              official == 'Ceftaroline' ~ 'Teflaro',
#              official == 'Ceftazidime' ~ 'Fortaz|Tazicef|Tazidime',
#              official == 'Ceftibuten' ~ 'Cedax',
#              official == 'Ceftizoxime' ~ 'Cefizox',
#              official == 'Ceftriaxone' ~ 'Rocephin',
#              official == 'Cefuroxime' ~ 'Ceftin|Zinacef',
#              official == 'Cephalexin' ~ 'Keflex|Panixine',
#              official == 'Cephalothin' ~ 'Keflin',
#              official == 'Chloramphenicol' ~ 'Chloromycetin',
#              official == 'Ciprofloxacin' ~ 'Cipro|Ciloxan|Ciproxin',
#              official == 'Clarithromycin' ~ 'Biaxin',
#              official == 'Clindamycin' ~ 'Cleocin|Clinda-Derm|Clindagel|Clindesse|Clindets|Evoclin',
#              official == 'Colistin' ~ 'Coly-Mycin',
#              official == 'Daptomycin' ~ 'Cubicin',
#              official == 'Doxycycline' ~ 'Doryx|Monodox|Vibramycin|Atridox|Oracea|Periostat|Vibra-Tabs',
#              official == 'Doripenem' ~ 'Doribax',
#              official == 'Ertapenem' ~ 'Invanz',
#              official == 'Erythromycin' ~ 'Eryc|EryPed|Erythrocin|E-Base|E-Glades|E-Mycin|E.E.S.|Ery-Tab|Eryderm|Erygel|Erythra-derm|Eryzole|Pediamycin',
#              official == 'Fosfomycin' ~ 'Monurol',
#              official == 'Flucloxacillin' ~ 'Flopen|Floxapen|Fluclox|Sesamol|Softapen|Staphylex',
#              official == 'Gentamicin' ~ 'Garamycin|Genoptic',
#              official == 'Imipenem' ~ 'Primaxin',
#              official == 'Kanamycin' ~ 'Kantrex',
#              official == 'Levofloxacin' ~ 'Levaquin|Quixin',
#              official == 'Linezolid' ~ 'Zyvox',
#              official == 'Lomefloxacin' ~ 'Maxaquin',
#              official == 'Meropenem' ~ 'Merrem',
#              official == 'Metronidazole' ~ 'Flagyl|MetroGel|MetroCream|MetroLotion',
#              official == 'Mezlocillin' ~ 'Mezlin',
#              official == 'Minocycline' ~ 'Arestin|Solodyn',
#              official == 'Moxifloxacin' ~ 'Avelox|Vigamox',
#              official == 'Mupirocin' ~ 'Bactroban|Centany',
#              official == 'Nafcillin' ~ 'Unipen',
#              official == 'Nalidixic acid' ~ 'NegGram',
#              official == 'Nitrofurantoin' ~ 'Furadantin|Macrobid|Macrodantin',
#              official == 'Norfloxacin' ~ 'Noroxin',
#              official == 'Ofloxacin' ~ 'Floxin|Ocuflox|Ophthalmic',
#              official == 'Oxacillin' ~ 'Bactocill',
#              official == 'Benzylpenicillin' ~ 'Permapen|Pfizerpen|Veetids',
#              official == 'Penicillins, combinations with other antibacterials' ~ 'Permapen|Pfizerpen|Veetids',
#              official == 'Piperacillin' ~ 'Pipracil',
#              official == 'Piperacillin and beta-lactamase inhibitor' ~ 'Zosyn',
#              official == 'Polymyxin B' ~ 'Poly-RX',
#              official == 'Quinupristin/dalfopristin' ~ 'Synercid',
#              official == 'Rifampin' ~ 'Rifadin|Rifamate|Rimactane',
#              official == 'Spectinomycin' ~ 'Trobicin',
#              official == 'Streptomycin' ~ 'Streptomycin Sulfate',
#              official == 'Teicoplanin' ~ 'Targocid',
#              official == 'Telavancin' ~ 'Vibativ',
#              official == 'Telithromcyin' ~ 'Ketek',
#              official == 'Tetracycline' ~ 'Sumycin|Bristacycline|Tetrex',
#              official == 'Ticarcillin' ~ 'Ticar',
#              official == 'Ticarcillin and beta-lactamase inhibitor' ~ 'Timentin',
#              official == 'Tigecycline' ~ 'Tygacil',
#              official == 'Tobramycin' ~ 'Tobi|Aktob|Tobre',
#              official == 'Trimethoprim' ~ 'Primsol|Proloprim',
#              official == 'Sulfamethoxazole and trimethoprim' ~ 'Bactrim|Septra|Sulfatrim',
#              official == 'Vancomycin' ~ 'Vancocin|Vancomycin Hydrochloride',
#              TRUE ~ NA_character_)
#   )
# last two columns created with:
# antibiotics %>%
#   mutate(useful_gramnegative =
#            if_else(
#              atc_group1 %like% '(fusidic|glycopeptide|macrolide|lincosamide|daptomycin|linezolid)' |
#                atc_group2 %like% '(fusidic|glycopeptide|macrolide|lincosamide|daptomycin|linezolid)' |
#                official %like% '(fusidic|glycopeptide|macrolide|lincosamide|daptomycin|linezolid)',
#              FALSE,
#              NA
#            ),
#          useful_grampositive =
#            if_else(
#              atc_group1 %like% '(aztreonam|temocillin|polymyxin|colistin|nalidixic)' |
#                atc_group2 %like% '(aztreonam|temocillin|polymyxin|colistin|nalidixic)' |
#                official %like% '(aztreonam|temocillin|polymyxin|colistin|nalidixic)',
#              FALSE,
#              NA
#            )
#   )
"antibiotics"

#' Dataset with ~2500 microorganisms
#'
#' A dataset containing 2456 microorganisms. MO codes of the UMCG can be looked up using \code{\link{microorganisms.umcg}}.
#' @format A data.frame with 2456 observations and 12 variables:
#' \describe{
#'   \item{\code{bactid}}{ID of microorganism}
#'   \item{\code{bactsys}}{Bactsyscode of microorganism}
#'   \item{\code{family}}{Family name of microorganism}
#'   \item{\code{genus}}{Genus name of microorganism, like \code{"Echerichia"}}
#'   \item{\code{species}}{Species name of microorganism, like \code{"coli"}}
#'   \item{\code{subspecies}}{Subspecies name of bio-/serovar of microorganism, like \code{"EHEC"}}
#'   \item{\code{fullname}}{Full name, like \code{"Echerichia coli (EHEC)"}}
#'   \item{\code{type}}{Type of microorganism, like \code{"Bacteria"} and \code{"Fungus/yeast"}}
#'   \item{\code{gramstain}}{Gram of microorganism, like \code{"Negative rods"}}
#'   \item{\code{aerobic}}{Logical whether bacteria is aerobic}
#'   \item{\code{type_nl}}{Type of microorganism in Dutch, like \code{"Bacterie"} and \code{"Schimmel/gist"}}
#'   \item{\code{gramstain_nl}}{Gram of microorganism in Dutch, like \code{"Negatieve staven"}}
#' }
#  source MOLIS (LIS of Certe) - \url{https://www.certe.nl}
#' @seealso \code{\link{guess_bactid}} \code{\link{antibiotics}} \code{\link{microorganisms.umcg}}
"microorganisms"

#' Translation table for UMCG with ~1100 microorganisms
#'
#' A dataset containing all bacteria codes of UMCG MMB. These codes can be joined to data with an ID from \code{\link{microorganisms}$bactid} (using \code{\link{left_join_microorganisms}}). GLIMS codes can also be translated to valid \code{bactid}'s with \code{\link{guess_bactid}}.
#' @format A data.frame with 1090 observations and 2 variables:
#' \describe{
#'   \item{\code{mocode}}{Code of microorganism according to UMCG MMB}
#'   \item{\code{bactid}}{Code of microorganism in \code{\link{microorganisms}}}
#' }
# source MOLIS (LIS of Certe) - \url{https://www.certe.nl} \cr \cr GLIMS (LIS of UMCG) - \url{https://www.umcg.nl}
#' @seealso \code{\link{guess_bactid}} \code{\link{microorganisms}}
"microorganisms.umcg"

#' Dataset with 2000 blood culture isolates of septic patients
#'
#' An anonymised dataset containing 2000 microbial blood culture isolates with their antibiogram of septic patients found in 5 different hospitals in the Netherlands, between 2001 and 2017. This data.frame can be used to practice AMR analysis. For examples, press F1.
#' @format A data.frame with 2000 observations and 49 variables:
#' \describe{
#'   \item{\code{date}}{date of receipt at the laboratory}
#'   \item{\code{hospital_id}}{ID of the hospital}
#'   \item{\code{ward_icu}}{logical to determine if ward is an intensive care unit}
#'   \item{\code{ward_clinical}}{logical to determine if ward is a regular clinical ward}
#'   \item{\code{ward_outpatient}}{logical to determine if ward is an outpatient clinic}
#'   \item{\code{age}}{age of the patient}
#'   \item{\code{sex}}{sex of the patient}
#'   \item{\code{patient_id}}{ID of the patient, first 10 characters of an SHA hash containing irretrievable information}
#'   \item{\code{bactid}}{ID of microorganism, see \code{\link{microorganisms}}}
#'   \item{\code{peni:rifa}}{40 different antibiotics with class \code{rsi} (see \code{\link{as.rsi}}); these column names occur in \code{\link{antibiotics}} data set and can be translated with \code{\link{abname}}}
#' }
# source MOLIS (LIS of Certe) - \url{https://www.certe.nl}
#' @examples
#' # ----------- #
#' # PREPARATION #
#' # ----------- #
#'
#' # Save this example dataset to an object, so we can edit it:
#' my_data <- septic_patients
#'
#' # load the dplyr package to make data science A LOT easier
#' library(dplyr)
#'
#' # Add first isolates to our dataset:
#' my_data <- my_data %>%
#'   mutate(first_isolates = first_isolate(my_data, "date", "patient_id", "bactid"))
#'
#' # -------- #
#' # ANALYSIS #
#' # -------- #
#'
#' # 1. Get the amoxicillin resistance percentages (p)
#' #     and numbers (n) of E. coli, divided by hospital:
#'
#' my_data %>%
#'   filter(bactid == guess_bactid("E. coli"),
#'          first_isolates == TRUE) %>%
#'   group_by(hospital_id) %>%
#'   summarise(n = n_rsi(amox),
#'             p = resistance(amox))
#'
#'
#' # 2. Get the amoxicillin/clavulanic acid resistance
#' #    percentages of E. coli, trend over the years:
#'
#' my_data %>%
#'   filter(bactid == guess_bactid("E. coli"),
#'          first_isolates == TRUE) %>%
#'   group_by(year = format(date, "%Y")) %>%
#'   summarise(n = n_rsi(amcl),
#'             p = resistance(amcl, minimum = 20))
"septic_patients"
