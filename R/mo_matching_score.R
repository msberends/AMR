# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Data Analysis for R                   #
#                                                                      #
# SOURCE                                                               #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2022 Berends MS, Luz CF et al.                              #
# Developed at the University of Groningen, the Netherlands, in        #
# collaboration with non-profit organisations Certe Medical            #
# Diagnostics & Advice, and University Medical Center Groningen.       # 
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                      #
# Visit our website for the full manual and a complete tutorial about  #
# how to conduct AMR data analysis: https://msberends.github.io/AMR/   #
# ==================================================================== #

#' Calculate the Matching Score for Microorganisms
#' 
#' This algorithm is used by [as.mo()] and all the [`mo_*`][mo_property()] functions to determine the most probable match of taxonomic records based on user input. 
#' @inheritSection lifecycle Stable Lifecycle
#' @author Dr Matthijs Berends
#' @param x Any user input value(s)
#' @param n A full taxonomic name, that exists in [`microorganisms$fullname`][microorganisms] or [`microorganisms.old$fullname`][microorganisms.old]
#' @section Matching Score for Microorganisms:
#' With ambiguous user input in [as.mo()] and all the [`mo_*`][mo_property()] functions, the returned results are chosen based on their matching score using [mo_matching_score()]. This matching score \eqn{m}, is calculated as:
#' 
#' \ifelse{latex}{\deqn{m_{(x, n)} = \frac{l_{n} - 0.5 \cdot \min \begin{cases}l_{n} \\ \textrm{lev}(x, n)\end{cases}}{l_{n} \cdot p_{n} \cdot k_{n}}}}{\ifelse{html}{\figure{mo_matching_score.png}{options: width="300" alt="mo matching score"}}{m(x, n) = ( l_n * min(l_n, lev(x, n) ) ) / ( l_n * p_n * k_n )}}
#' 
#' where:
#' 
#' * \ifelse{html}{\out{<i>x</i> is the user input;}}{\eqn{x} is the user input;}
#' * \ifelse{html}{\out{<i>n</i> is a taxonomic name (genus, species, and subspecies);}}{\eqn{n} is a taxonomic name (genus, species, and subspecies);}
#' * \ifelse{html}{\out{<i>l<sub>n</sub></i> is the length of <i>n</i>;}}{l_n is the length of \eqn{n};}
#' * \ifelse{html}{\out{<i>lev</i> is the <a href="https://en.wikipedia.org/wiki/Levenshtein_distance">Levenshtein distance function</a>, which counts any insertion, deletion and substitution as 1 that is needed to change <i>x</i> into <i>n</i>;}}{lev is the Levenshtein distance function, which counts any insertion, deletion and substitution as 1 that is needed to change \eqn{x} into \eqn{n};}
#' * \ifelse{html}{\out{<i>p<sub>n</sub></i> is the human pathogenic prevalence group of <i>n</i>, as described below;}}{p_n is the human pathogenic prevalence group of \eqn{n}, as described below;}
#' * \ifelse{html}{\out{<i>k<sub>n</sub></i> is the taxonomic kingdom of <i>n</i>, set as Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5.}}{l_n is the taxonomic kingdom of \eqn{n}, set as Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5.}
#' 
#' The grouping into human pathogenic prevalence (\eqn{p}) is based on experience from several microbiological laboratories in the Netherlands in conjunction with international reports on pathogen prevalence. **Group 1** (most prevalent microorganisms) consists of all microorganisms where the taxonomic class is Gammaproteobacteria or where the taxonomic genus is *Enterococcus*, *Staphylococcus* or *Streptococcus*. This group consequently contains all common Gram-negative bacteria, such as *Pseudomonas* and *Legionella* and all species within the order Enterobacterales. **Group 2** consists of all microorganisms where the taxonomic phylum is Proteobacteria, Firmicutes, Actinobacteria or Sarcomastigophora, or where the taxonomic genus is *Absidia*, *Acholeplasma*, *Acremonium*, *Actinotignum*, *Aedes*, *Alistipes*, *Alloprevotella*, *Alternaria*, *Anaerosalibacter*, *Ancylostoma*, *Angiostrongylus*, *Anisakis*, *Anopheles*, *Apophysomyces*, *Arachnia*, *Aspergillus*, *Aureobasidium*, *Bacteroides*, *Basidiobolus*, *Beauveria*, *Bergeyella*, *Blastocystis*, *Blastomyces*, *Borrelia*, *Brachyspira*, *Branhamella*, *Butyricimonas*, *Candida*, *Capillaria*, *Capnocytophaga*, *Catabacter*, *Cetobacterium*, *Chaetomium*, *Chlamydia*, *Chlamydophila*, *Chryseobacterium*, *Chrysonilia*, *Cladophialophora*, *Cladosporium*, *Conidiobolus*, *Contracaecum*, *Cordylobia*, *Cryptococcus*, *Curvularia*, *Deinococcus*, *Demodex*, *Dermatobia*, *Diphyllobothrium*, *Dirofilaria*, *Dysgonomonas*, *Echinostoma*, *Elizabethkingia*, *Empedobacter*, *Enterobius*, *Exophiala*, *Exserohilum*, *Fasciola*, *Flavobacterium*, *Fonsecaea*, *Fusarium*, *Fusobacterium*, *Giardia*, *Haloarcula*, *Halobacterium*, *Halococcus*, *Hendersonula*, *Heterophyes*, *Histoplasma*, *Hymenolepis*, *Hypomyces*, *Hysterothylacium*, *Lelliottia*, *Leptosphaeria*, *Leptotrichia*, *Lucilia*, *Lumbricus*, *Malassezia*, *Malbranchea*, *Metagonimus*, *Microsporum*, *Mortierella*, *Mucor*, *Mycocentrospora*, *Mycoplasma*, *Myroides*, *Necator*, *Nectria*, *Ochroconis*, *Odoribacter*, *Oesophagostomum*, *Oidiodendron*, *Opisthorchis*, *Ornithobacterium*, *Parabacteroides*, *Pediculus*, *Pedobacter*, *Phlebotomus*, *Phocaeicola*, *Phocanema*, *Phoma*, *Piedraia*, *Pithomyces*, *Pityrosporum*, *Porphyromonas*, *Prevotella*, *Pseudallescheria*, *Pseudoterranova*, *Pulex*, *Rhizomucor*, *Rhizopus*, *Rhodotorula*, *Riemerella*, *Saccharomyces*, *Sarcoptes*, *Scolecobasidium*, *Scopulariopsis*, *Scytalidium*, *Sphingobacterium*, *Spirometra*, *Spiroplasma*, *Sporobolomyces*, *Stachybotrys*, *Streptobacillus*, *Strongyloides*, *Syngamus*, *Taenia*, *Tannerella*, *Tenacibaculum*, *Terrimonas*, *Toxocara*, *Treponema*, *Trichinella*, *Trichobilharzia*, *Trichoderma*, *Trichomonas*, *Trichophyton*, *Trichosporon*, *Trichostrongylus*, *Trichuris*, *Tritirachium*, *Trombicula*, *Tunga*, *Ureaplasma*, *Victivallis*, *Wautersiella*, *Weeksella*, or *Wuchereria*. **Group 3** consists of all other microorganisms.
#' 
#' All characters in \eqn{x} and \eqn{n} are ignored that are other than A-Z, a-z, 0-9, spaces and parentheses.
#' 
#' All matches are sorted descending on their matching score and for all user input values, the top match will be returned. This will lead to the effect that e.g., `"E. coli"` will return the microbial ID of *Escherichia coli* (\eqn{m = `r round(mo_matching_score("E. coli", "Escherichia coli"), 3)`}, a highly prevalent microorganism found in humans) and not *Entamoeba coli* (\eqn{m = `r round(mo_matching_score("E. coli", "Entamoeba coli"), 3)`}, a less prevalent microorganism in humans), although the latter would alphabetically come first. 
#' @export
#' @inheritSection AMR Reference Data Publicly Available
#' @inheritSection AMR Read more on Our Website!
#' @examples 
#' as.mo("E. coli")
#' mo_uncertainties()
#' 
#' mo_matching_score(x = "E. coli",
#'                   n = c("Escherichia coli", "Entamoeba coli"))
#'                   
#' # old taxonomic names work as well:
#' mo_matching_score(x = "E. blattae", n = "Escherichia blattae")
mo_matching_score <- function(x, n) {
  meet_criteria(x, allow_class = c("character", "data.frame", "list"))
  meet_criteria(n, allow_class = "character")
  
  # no dots and other non-whitespace characters
  x <- gsub("[^a-zA-Z0-9 \\(\\)]+", "", x)
  
  # only keep one space
  x <- gsub(" +", " ", x)
  
  # n is always a taxonomically full name
  if (length(n) == 1) {
    n <- rep(n, length(x))
  }
  if (length(x) == 1) {
    x <- rep(x, length(n))
  }
  
  # length of fullname
  l_n <- nchar(n)
  lev <- double(length = length(x))
  l_n.lev <- double(length = length(x))
  for (i in seq_len(length(x))) {
    # determine Levenshtein distance, but maximise to nchar of n
    lev[i] <- utils::adist(x[i], n[i], ignore.case = FALSE, fixed = TRUE, costs = c(ins = 1, del = 1, sub = 1))
    # minimum of (l_n, Levenshtein distance)
    l_n.lev[i] <- min(l_n[i], as.double(lev[i]))
  }
  # support n to be a previously valid name - replace here with current name to get prevalence and kingdom
  n_to_current <- which(!n %in% MO_lookup$fullname)
  if (length(n_to_current) > 0) {
    n[n_to_current] <- MO.old_lookup$fullname_new[match(n[n_to_current], MO.old_lookup$fullname)]
  }
  # human pathogenic prevalence (1 to 3), see ?as.mo
  p_n <- MO_lookup$prevalence[match(n, MO_lookup$fullname)]
  # kingdom index (Bacteria = 1, Fungi = 2, Protozoa = 3, Archaea = 4, others = 5)
  k_n <- MO_lookup$kingdom_index[match(n, MO_lookup$fullname)]
  
  # matching score (save Levenshtein part as attribute):
  structure((l_n - 0.5 * l_n.lev) / (l_n * p_n * k_n),
            lev = (l_n - 0.5 * l_n.lev) / l_n)
}
