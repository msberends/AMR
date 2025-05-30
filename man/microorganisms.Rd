% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{microorganisms}
\alias{microorganisms}
\title{Data Set with 78 679 Taxonomic Records of Microorganisms}
\format{
A \link[tibble:tibble]{tibble} with 78 679 observations and 26 variables:
\itemize{
\item \code{mo}\cr ID of microorganism as used by this package. \emph{\strong{This is a unique identifier.}}
\item \code{fullname}\cr Full name, like \code{"Escherichia coli"}. For the taxonomic ranks genus, species and subspecies, this is the 'pasted' text of genus, species, and subspecies. For all taxonomic ranks higher than genus, this is the name of the taxon. \emph{\strong{This is a unique identifier.}}
\item \code{status} \cr Status of the taxon, either "accepted", "not validly published", "synonym", or "unknown"
\item \code{kingdom}, \code{phylum}, \code{class}, \code{order}, \code{family}, \code{genus}, \code{species}, \code{subspecies}\cr Taxonomic rank of the microorganism. Note that for fungi, \emph{phylum} is equal to their taxonomic \emph{division}. Also, for fungi, \emph{subkingdom} and \emph{subdivision} were left out since they do not occur in the bacterial taxonomy.
\item \code{rank}\cr Text of the taxonomic rank of the microorganism, such as \code{"species"} or \code{"genus"}
\item \code{ref}\cr Author(s) and year of related scientific publication. This contains only the \emph{first surname} and year of the \emph{latest} authors, e.g. "Wallis \emph{et al.} 2006 \emph{emend.} Smith and Jones 2018" becomes "Smith \emph{et al.}, 2018". This field is directly retrieved from the source specified in the column \code{source}. Moreover, accents were removed to comply with CRAN that only allows ASCII characters.
\item \code{oxygen_tolerance} \cr Oxygen tolerance, either "aerobe", "anaerobe", "anaerobe/microaerophile", "facultative anaerobe", "likely facultative anaerobe", or "microaerophile". These data were retrieved from BacDive (see \emph{Source}). Items that contain "likely" are missing from BacDive and were extrapolated from other species within the same genus to guess the oxygen tolerance. Currently 68.3\% of all ~39 000 bacteria in the data set contain an oxygen tolerance.
\item \code{source}\cr Either "GBIF", "LPSN", "Manually added", "MycoBank", or "manually added" (see \emph{Source})
\item \code{lpsn}\cr Identifier ('Record number') of List of Prokaryotic names with Standing in Nomenclature (LPSN). This will be the first/highest LPSN identifier to keep one identifier per row. For example, \emph{Acetobacter ascendens} has LPSN Record number 7864 and 11011. Only the first is available in the \code{microorganisms} data set. \emph{\strong{This is a unique identifier}}, though available for only ~33 000 records.
\item \code{lpsn_parent}\cr LPSN identifier of the parent taxon
\item \code{lpsn_renamed_to}\cr LPSN identifier of the currently valid taxon
\item \code{mycobank}\cr Identifier ('MycoBank #') of MycoBank. \emph{\strong{This is a unique identifier}}, though available for only ~18 000 records.
\item \code{mycobank_parent}\cr MycoBank identifier of the parent taxon
\item \code{mycobank_renamed_to}\cr MycoBank identifier of the currently valid taxon
\item \code{gbif}\cr Identifier ('taxonID') of Global Biodiversity Information Facility (GBIF). \emph{\strong{This is a unique identifier}}, though available for only ~49 000 records.
\item \code{gbif_parent}\cr GBIF identifier of the parent taxon
\item \code{gbif_renamed_to}\cr GBIF identifier of the currently valid taxon
\item \code{prevalence}\cr Prevalence of the microorganism based on Bartlett \emph{et al.} (2022, \doi{10.1099/mic.0.001269}), see \code{\link[=mo_matching_score]{mo_matching_score()}} for the full explanation
\item \code{snomed}\cr Systematized Nomenclature of Medicine (SNOMED) code of the microorganism, version of July 16th, 2024 (see \emph{Source}). Use \code{\link[=mo_snomed]{mo_snomed()}} to retrieve it quickly, see \code{\link[=mo_property]{mo_property()}}.
}
}
\source{
Taxonomic entries were imported in this order of importance:
\enumerate{
\item List of Prokaryotic names with Standing in Nomenclature (LPSN):\cr\cr
Parte, AC \emph{et al.} (2020). \strong{List of Prokaryotic names with Standing in Nomenclature (LPSN) moves to the DSMZ.} International Journal of Systematic and Evolutionary Microbiology, 70, 5607-5612; \doi{10.1099/ijsem.0.004332}. Accessed from \url{https://lpsn.dsmz.de} on June 24th, 2024.
\item MycoBank:\cr\cr
Vincent, R \emph{et al} (2013). \strong{MycoBank gearing up for new horizons.} IMA Fungus, 4(2), 371-9; \doi{10.5598/imafungus.2013.04.02.16}. Accessed from \url{https://www.mycobank.org} on June 24th, 2024.
\item Global Biodiversity Information Facility (GBIF):\cr\cr
GBIF Secretariat (2023). GBIF Backbone Taxonomy. Checklist dataset \doi{10.15468/39omei}. Accessed from \url{https://www.gbif.org} on June 24th, 2024.
}

Furthermore, these sources were used for additional details:
\itemize{
\item BacDive:\cr\cr
Reimer, LC \emph{et al.} (2022). \strong{\emph{BacDive} in 2022: the knowledge base for standardized bacterial and archaeal data.} Nucleic Acids Res., 50(D1):D741-D74; \doi{10.1093/nar/gkab961}. Accessed from \url{https://bacdive.dsmz.de} on July 16th, 2024.
\item Systematized Nomenclature of Medicine - Clinical Terms (SNOMED-CT):\cr\cr
Public Health Information Network Vocabulary Access and Distribution System (PHIN VADS). US Edition of SNOMED CT from 1 September 2020. Value Set Name 'Microorganism', OID 2.16.840.1.114222.4.11.1009 (v12). Accessed from \url{https://www.cdc.gov/phin/php/phinvads/} on July 16th, 2024.
\item Grimont \emph{et al.} (2007). Antigenic Formulae of the Salmonella Serovars, 9th Edition. WHO Collaborating Centre for Reference and Research on \emph{Salmonella} (WHOCC-SALM).
\item Bartlett \emph{et al.} (2022). \strong{A comprehensive list of bacterial pathogens infecting humans} \emph{Microbiology} 168:001269; \doi{10.1099/mic.0.001269}
}
}
\usage{
microorganisms
}
\description{
A data set containing the full microbial taxonomy (\strong{last updated: June 24th, 2024}) of six kingdoms. This data set is the backbone of this \code{AMR} package. MO codes can be looked up using \code{\link[=as.mo]{as.mo()}} and microorganism properties can be looked up using any of the \code{\link[=mo_property]{mo_*}} functions.

This data set is carefully crafted, yet made 100\% reproducible from public and authoritative taxonomic sources (using \href{https://github.com/msberends/AMR/blob/main/data-raw/_reproduction_scripts/reproduction_of_microorganisms.R}{this script}), namely: \emph{List of Prokaryotic names with Standing in Nomenclature (LPSN)} for bacteria, \emph{MycoBank} for fungi, and \emph{Global Biodiversity Information Facility (GBIF)} for all others taxons.
}
\details{
Please note that entries are only based on LPSN, MycoBank, and GBIF (see below). Since these sources incorporate entries based on (recent) publications in the International Journal of Systematic and Evolutionary Microbiology (IJSEM), it can happen that the year of publication is sometimes later than one might expect.

For example, \emph{Staphylococcus pettenkoferi} was described for the first time in Diagnostic Microbiology and Infectious Disease in 2002 (\doi{10.1016/s0732-8893(02)00399-1}), but it was not until 2007 that a publication in IJSEM followed (\doi{10.1099/ijs.0.64381-0}). Consequently, the \code{AMR} package returns 2007 for \code{mo_year("S. pettenkoferi")}.
}
\section{Included Taxa}{

Included taxonomic data from \href{https://lpsn.dsmz.de}{LPSN}, \href{https://www.mycobank.org}{MycoBank}, and \href{https://www.gbif.org}{GBIF} are:
\itemize{
\item All ~39 000 (sub)species from the kingdoms of Archaea and Bacteria
\item ~28 000 species from the kingdom of Fungi. The kingdom of Fungi is a very large taxon with almost 300,000 different (sub)species, of which most are not microbial (but rather macroscopic, like mushrooms). Because of this, not all fungi fit the scope of this package. Only relevant fungi are covered (such as all species of \emph{Aspergillus}, \emph{Candida}, \emph{Cryptococcus}, \emph{Histoplasma}, \emph{Pneumocystis}, \emph{Saccharomyces} and \emph{Trichophyton}).
\item ~8 100 (sub)species from the kingdom of Protozoa
\item ~1 600 (sub)species from 39 other relevant genera from the kingdom of Animalia (such as \emph{Strongyloides} and \emph{Taenia})
\item All ~22 000 previously accepted names of all included (sub)species (these were taxonomically renamed)
\item The complete taxonomic tree of all included (sub)species: from kingdom to subspecies
\item The identifier of the parent taxons
\item The year and first author of the related scientific publication
}
\subsection{Manual additions}{

For convenience, some entries were added manually:
\itemize{
\item ~1 500 entries of \emph{Salmonella}, such as the city-like serovars and groups A to H
\item 37 species groups (such as the beta-haemolytic \emph{Streptococcus} groups A to K, coagulase-negative \emph{Staphylococcus} (CoNS), \emph{Mycobacterium tuberculosis} complex, etc.), of which the group compositions are stored in the \link{microorganisms.groups} data set
\item 1 entry of \emph{Blastocystis} (\emph{B.  hominis}), although it officially does not exist (Noel \emph{et al.} 2005, PMID 15634993)
\item 1 entry of \emph{Moraxella} (\emph{M. catarrhalis}), which was formally named \emph{Branhamella catarrhalis} (Catlin, 1970) though this change was never accepted within the field of clinical microbiology
\item 8 other 'undefined' entries (unknown, unknown Gram-negatives, unknown Gram-positives, unknown yeast, unknown fungus, and unknown anaerobic Gram-pos/Gram-neg bacteria)
}

The syntax used to transform the original data to a cleansed \R format, can be \href{https://github.com/msberends/AMR/blob/main/data-raw/_reproduction_scripts/reproduction_of_microorganisms.R}{found here}.
}
}

\section{Download Our Reference Data}{

All reference data sets in the AMR package - including information on microorganisms, antimicrobials, and clinical breakpoints - are freely available for download in multiple formats: R, MS Excel, Apache Feather, Apache Parquet, SPSS, and Stata.

For maximum compatibility, we also provide machine-readable, tab-separated plain text files suitable for use in any software, including laboratory information systems.

Visit \href{https://amr-for-r.org/articles/datasets.html}{our website for direct download links}, or explore the actual files in \href{https://github.com/msberends/AMR/tree/main/data-raw/datasets}{our GitHub repository}.
}

\examples{
microorganisms
}
\seealso{
\code{\link[=as.mo]{as.mo()}}, \code{\link[=mo_property]{mo_property()}}, \link{microorganisms.groups}, \link{microorganisms.codes}, \link{intrinsic_resistant}
}
\keyword{datasets}
