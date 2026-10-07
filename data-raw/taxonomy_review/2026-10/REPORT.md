# Taxonomy update of `microorganisms`, October 2026

Unattended rebuild of the `microorganisms` data set with
`data-raw/_reproduction_scripts/reproduction_of_microorganisms.R`, following
`reproduction_of_microorganisms_CLAUDE.md`. All review CSV files referred to below are in this folder.

**Baseline.** On request of Matthijs S. Berends, the previous data set is the one of release **v3.0.1**, not the
development version on `main`, which has known flaws. The script now takes it from the last git tag (see script change
1). The comparisons below are therefore against v3.0.1 (78,679 records).

## 1. Sources

| Source | Version / accessed | Citation |
|---|---|---|
| GBIF/COL | Catalogue of Life **2026-04-18 XR** (COL26.4 XR), issued 2026-04-18, `doi:10.48580/dgxjw`; downloaded as COL Data Package (ColDP), file date 2026-04-30 | Bánki, O. *et al.* (2026). Catalogue of Life (2026-04-18 XR). Catalogue of Life Foundation, Amsterdam, Netherlands. (unchanged, `main` already cited this release) |
| LPSN | `taxonomy.csv`, file date 2026-10-05 (34,524 records) | Freese, HM *et al.* (2026), unchanged |
| MycoBank | `MBList.xlsx`, file date 2026-01-07 | Vincent, R *et al.* (2013), unchanged |
| BacDive | oxygen tolerance and cell shape exports, file date 2026-10-05 | Reimer, LC *et al.* (2022), unchanged |
| Bartlett *et al.* | latest `Tab 6 Full List` from GitHub (replaces the tracked file; 8 spelling corrections, e.g. *Arcanobacterium haemolyticum*) | Bartlett *et al.* (2022), `doi:10.1099/mic.0.001269` |

The LPSN scrape (lineages of all 6,546 genera, 0 failures) ran once; its cache (`data-raw/lpsn_scrape_cache.rds`) is
git-ignored and not committed.

## 2. Summary

**107,253 records** (v3.0.1: 78,679; development version on `main`: 96,982). Per domain, old (v3.0.1) vs new, see
also [CSV](comparison_per_domain_rank_status.csv) and [CSV](new_records_per_domain_and_kind.csv):

| Domain | Total old | Total new | Accepted old | Accepted new | Synonym old | Synonym new | Unknown/other old / new | Change |
|---|---|---|---|---|---|---|---|---|
| Animalia | 1,628 | 2,579 | 1,379 | 1,754 | 244 | 810 | 5 / 15 | +58.4% |
| Archaea | 1,419 | 1,573 | 1,225 | 1,359 | 181 | 201 | 13 / 13 | +10.9% |
| Bacteria | 39,249 | 48,239 | 29,853 | 38,689 | 7,076 | 8,777 | 2,320 / 773 | +22.9% |
| Chromista | 178 | 418 | 157 | 342 | 17 | 68 | 4 / 8 | +134.8% |
| Fungi | 28,137 | 47,412 | 14,385 | 22,886 | 8,746 | 19,294 | 5,006 / 5,232 | +68.5% |
| Protozoa | 8,067 | 7,030 | 6,056 | 5,097 | 1,880 | 1,805 | 131 / 128 | -12.9% |

Explanation of the changes of more than 5%:

- **Bacteria (+22.9%)**: new validly published species in LPSN since 2024 (accepted +8,836), and v3.0.1 'not validly
  published' names now as accepted, synonym or unknown.
- **Fungi (+68.5%)**: mostly synonyms (+10,548): old names of the kept fungal species, from MycoBank and COL, which
  `as.mo()` uses to translate old names (kept on purpose, decision of 5 October 2026). Accepted +8,501: new species in
  relevant genera, current names of species of relevant genera (script change 5), and the microsporidia, now in the
  Fungi (1,394 names of v3.0.1 moved from the Protozoa, script change 19). Lichens are removed (script change 17), except
  the 84 that were released before.
- **Animalia (+58.4%)**, **Chromista (+134.8%)**: helminths, vectors and protists of the relevant genera that COL now
  contains (and the COL records without authors, script change 3), for Chromista also *Balantidium* and *Isospora*;
  small absolute numbers (+951 and +240).
- **Archaea (+10.9%)**: new species in LPSN.
- **Protozoa (-12.9%)**: the microsporidia moved to the Fungi (script change 19); the Protozoa follow the relevance
  rules of the Fungi (script change 16), most protozoa of v3.0.1 are restored as released taxa.

Records by rank: species 83,950, subspecies 8,682, genus 11,286, family 2,093, order 726, class 318, phylum 131, and
37 species groups. MO codes of existing names (same name and domain) that changed: **0**
([CSV](changed_mo_codes_of_existing_names.csv)); names that moved to another domain (and so got another code prefix):
1,481, of which 1,394 microsporidia (Protozoa to Fungi), the others the corrections in the checklist
([CSV](names_that_moved_to_another_domain.csv)). Codes in the package's other data sets missing from `microorganisms`:
**0**.

Run time: the source chunks take about 15 minutes (LPSN scrape 10.5 minutes, ran once), each full rebuild from a
checkpoint 15 to 25 minutes. 18 runs were needed because of the script changes below (logs available on request).

## 3. Script changes

All changes are in `data-raw/_reproduction_scripts/reproduction_of_microorganisms.R` unless stated otherwise.

1. **Baseline from the last release.** `microorganisms_old` is now read from `data/microorganisms.rda` of the last git
   tag (v3.0.1), with `domain` set to `kingdom` for releases before the domain column existed. The development data set
   is kept as `microorganisms_dev`, only to translate the codes in the package's other data sets (`genera_of_mo()`,
   `fix_old_mos()`). Affects all relevance lists, carried-over and manually added records, and all comparisons.
2. **COL Data Package and file names.** COL is now downloaded as ColDP (your updated download notes): `NameUsage.tsv` is
   read with explicit column selection and its columns are translated to the Darwin Core names the script was built on
   (`parentID` is the accepted name of a synonym). Verified: all 7,834,411 rows read without parsing problems. File
   names now match the downloads (`MBList.xlsx`, `bacdive_oxygen_tolerance.csv`, `bacdive_cell_shape.csv`).
3. **COL records without authors were dropped.** `filter(ref != "AmSOD")` removed every record without authors, as
   `NA != "AmSOD"` is `NA`. Now `is.na(ref) | ref != "AmSOD"`: 16,508 records more in `taxonomy_gbif.rds`, e.g.
   *Balamuthia mandrillaris*.
4. **Genus homonyms across domains** (#309):
   - MycoBank records of a genus without any phylum are no authoritative fungal records any more (MycoBank also indexes
     *Plasmodium* Marchiafava & Celli, *Sarcocystis*, *Pseudospirillum*, which made these 'Fungi');
   - new rule 2: the domain of the genus in the last release wins, if an authoritative source supports it (or no domain
     has one), and never on unclassified MycoBank records alone. Without this, the released bacterial genera *Moorella*,
     *Serpula*, *Pirella*, *Chainia*, *Bogoriella*, *Frondicola*, *Pseudospirillum* (now LPSN synonyms, same organisms)
     lost to fungal namesakes and their released codes would have been retired. See
     [CSV](homonyms_decided_by_the_domain_of_the_last_release.csv) (column `released_decides`).
5. **Relevant genera no longer include whole 'current genera'.** The October 2026 step that added the current genera of
   all species of relevant genera made e.g. *Cercospora* and *Colletotrichum* relevant entirely (the prevalence of
   v3.0.1 was set per genus, so every historical *Fusarium* name counted): the data set grew to 222,329 records. Now the
   current names themselves are protected (`relevant_current_species`) with their genus record.
6. **Wrong links from `mo_current()` and `mo_synonyms()`.** Names missing from the development data (e.g. *Brugia*,
   *Ascaris*, *Babesia*) were fuzzy-matched to other organisms, e.g. *Brugia* to the puffball *Lycoperdon*, which made
   420 puffball records 'potentially pathogenic'. New helper `current_names_same_domain()`; `mo_synonyms()` is only
   used for names that exist in the development data.
7. **MycoBank names must be in Latin script** (letters, space, hyphen): removes names such as
   *Cordyceps \*jezoensoides* and *Pyrenopeziza poae* written with Cyrillic letters (integrity rule `invalid_mo_characters`).
8. **Old records outside the scope of their source** get source `manually added` (e.g. the microsporidian families
   Culicosporidae, Golbergiidae, Janacekiidae, Protozoa with source MycoBank in v3.0.1; integrity rule
   `mycobank_outside_fungi`, known defect).
9. **Other organisms with the same genus name from v3.0.1** (the *Graphium* butterflies), which entered via three routes:
   carried-over records, previously 'manually added' records (until 2026 also inferred records) and restored released
   taxa. New helper `names_in_sources()`; (sub)species of which the species only exists in another domain in the
   sources are not carried over; removed homonyms are not re-added; new explicit list `released_homonym_genera`
   (*Graphium* = Fungi): of these genera, released (sub)species are only kept if a current source has them in that
   domain, the others are retired with `NEEDS REVIEW`.
10. **Restoring released taxa**:
    - missing taxa are detected by name *and* rank, and a name at two ranks gets the `{rank}` suffix again (the genus
      *Kapabacteria* and the class 'Kapabacteria {class}'; integrity rule `species_without_genus`);
    - a released code in the winning domain is not retired as a homonym if its name still exists in that domain in the
      sources (e.g. the fungal genus *Trichurus*, as 'Trichurus' in Animalia is a misspelling of *Trichuris*);
    - a name registered in two domains (e.g. *Octospora*, a fungus in v3.0.1 and a microsporidian in v2.1.1) is only
      restored from the most recent release; the older codes are retired with `NEEDS REVIEW`.
11. **Empty families of genera** are filled from the COL genus record in the same order (e.g. *Lacazia*:
    Ajellomycetaceae), and for genera without a family in any source from the new commented list
    `genus_family_override` (*Ascaris* = Ascarididae, *Balamuthia* = Balamuthiidae).
12. ***Brugia*** species added like the existing *Leishmania* block (COL has no *Brugia* species at all):
    *B. malayi*, *B. timori*. Before, `as.mo("Brugia malayi")` gave *Serratia microhaemolytica* (also on `main`).
13. ***Salmonella* species synonyms that are serovars** in this data set (*S. typhi*, *S. paratyphi*, *S. typhimurium*)
    are removed: they made `as.mo("S. typhi")` return *S. enterica* instead of *S.* Typhi (regression against `main`
    and v3.0.1; until 2026 they were removed as duplicate codes).
14. **Oxygen tolerance from BacDive**: a minority of strain records of 'facultative anaerobe' no longer decides for the
    species (only if at least half): *Clostridioides difficile* had 6 of 52 such records and became 'facultative
    anaerobe' (unit test).
15. **Codes in the other data sets** (`fix_old_mos()`): translated with the development data set; two names that no
    longer exist are mapped (*Trichophyton mentagrophytes indotineae* to *T. indotineae*, as in EUCAST and MycoBank;
    *Mycobacteroides stephanolepidis* to *Mycobacterium stephanolepidis*, LPSN); unmatched codes are reviewed per data
    set; in `intrinsic_resistant`, unmatched rows are dropped (non-bacteria that the development version had as
    bacteria, see [CSV](codes_without_a_match_in_amr_intrinsic_resistant_mo.csv)).
16. **Scope decisions of 5 October 2026** (Matthijs S. Berends): names that the COL eXtended Release added
    programmatically (`clb:merged`) are only used if their genus is clinically relevant, they are a protected current
    name, or they were in the last release (40,897 fewer records in `taxonomy_gbif.rds`); Protozoa follow the relevance
    rules of the Fungi instead of being kept entirely (that rule dated from before the XR filled this domain).
17. **No lichens** (decision of 5 October 2026): fungi in the lichen-forming classes (Arthoniomycetes, Candelariomycetes,
    Lecanoromycetes, Lichinomycetes) and the order Verrucariales are removed unless their genus is clinically relevant,
    with the synonyms that point to them: 4,155 records (528 accepted, 3,627 synonyms). They entered as current names of
    kept synonyms. The 84 lichens of v3.0.1 are restored as released taxa. See
    [CSV](lichens_and_their_synonyms_these_will_be_removed.csv).
18. **Synonym cycles decided** (decision of 5 October 2026): in a cycle that the sources cannot settle, the member in a
    genus of the new list `synonym_cycle_current_genera` becomes the current name: *Capillidium* (former
    *Conidiobolus* species) and *Scolecobasidium* (former *Ochroconis* and *Pseudosigmoidea* species).
19. **Microsporidia are Fungi** (decision of 5 October 2026): new helpers `is_microsporidian()` and
    `microsporidia_to_fungi()`, applied to COL, the carried-over, previously manually added and restored records, and
    the released domain of genera. COL and most releases had them in the Protozoa, so their codes change from `P_` to
    `F_` (e.g. *Enterocytozoon bieneusi*); the `P_` codes keep working, translated by name (registry rule 5). Codes of
    v2.x of the same microsporidia are no longer retired, but translated by name as well.
20. ***Enterocytozoon bieneusi*** **and** ***Encephalitozoon cuniculi*** **are the accepted names** (decision of
    5 October 2026): COL had them as synonyms of '*Encephalitozoon bieneusi*' (as did v3.0.1) and '*Nosema cuniculi*'.
    New commented list `accepted_name_override` for such source flaws: the accepted name and its synonym are swapped,
    and other synonyms point to the accepted name.
21. ***Cystoisospora belli*** is added as a synonym of *Isospora belli* (decision of 5 October 2026), as COL does not
    contain this name at all (`as.mo()` returned an unknown).

Other files:

- `data-raw/_pre_commit_checks.R`: *Staphylococcus dromedarii* (abstract: "a novel coagulase-negative
  *Staphylococcus* species", `doi:10.1099/ijsem.0.007230`), *S. parequorum* (Baek *et al.*, J Microbiol 2025) and
  *S. xeri* (Belhout *et al.* 2026, `doi:10.1099/ijsem.0.007253`) added to CoNS; *S. parequorum* and *S. xeri*
  confirmed as coagulase-negative by Matthijs S. Berends.
- `R/mo.R` (on request): `as.mo()` now finds *Salmonella* serovars written with the species, such as
  "Salmonella enterica (subsp. enterica) serovar Typhi", which returned *S. bongori* (also on `main` and in v3.0.1),
  and no longer turns "Paratyphi A/B/C" into the serogroups A/B/C; with tests in `test-mo.R` and a `NEWS.md` bullet.
- `data-raw/_pre_commit_checks.R`: `microorganisms.dta` is written with `strl_threshold = 255` (as
  `clinical_breakpoints.dta` already was), as it otherwise grew to 130 MB, beyond the file size limit of GitHub (now
  78 MiB, verified to read back identically).
- `data-raw/_pre_commit_checks.R`: *Balantidium* and *Isospora* added to `MO_RELEVANT_GENERA`, the genera in which COL
  has *Balantidium coli* (*Balantioides coli*) and *Isospora belli* (*Cystoisospora belli*); before, `as.mo()` returned
  *Campylobacter pinnipediorum caledonicus* for *Balantidium coli* and *Stenotrophomonas beteli* for *Isospora belli*.
- `tests/testthat/test-mo_property.R` and `test-data-microorganisms.R`: see the checklist items on the unit tests.
- `R/aa_globals.R`: accessed dates in `TAXONOMY_VERSION`; `NEWS.md`: the existing taxonomy bullet of this series updated.

## 4. Checklist

### Review blocks of the build (3a)

- [x] **LPSN genera with multiple lineages**: 6 rows, the homonyms *Halalkalibacterium*, *Pusillimonas* and
  *Rhodococcus*; the script removes the Balneolaceae, Oscillospiraceae and Chroococcaceae lineages explicitly, the kept
  families match v3.0.1, see [CSV](lpsn_genera_with_multiple_lineages.csv)
- [x] **LPSN higher taxa that could not be downloaded**: 0 rows
- [x] **Homonyms of relevant genera**: 219 rows (about 100 genera in more than one domain); after script change 4 all
  winners plausible, e.g. *Plasmodium*, *Sarcocystis*, *Babesia* Chromista (COL; *Plasmodium* is set to Protozoa by
  the curated block later on), *Moorella*, *Serpula*, *Chainia*, *Pirella* Bacteria, *Necator*, *Capillaria*, *Hymenolepis* Animalia (override list), *Trichophyton*, *Graphium*, *Cryptococcus*,
  *Mucor* and the microsporidia (e.g. *Enterocytozoon*) Fungi, see [CSV](homonyms_of_relevant_genera_check_best_domain.csv)
- [x] **Homonyms decided by the domain of the last release** (new): 21 rows, of which the release decides where an
  authoritative source supports it (column `released_decides`), see
  [CSV](homonyms_decided_by_the_domain_of_the_last_release.csv); doubtful ones under "Needs a human"
- [x] **Records removed as homonyms of a genus in another domain**: 173 rows, other organisms with the same genus name
  (e.g. the polychaete *Serpula*, fungal *Moorella*, *Necator*, the unclassified MycoBank index entries of *Plasmodium*
  and *Sarcocystis*), see [CSV](records_removed_as_homonyms_of_a_genus_in_another_domain.csv)
- [x] **Genera with multiple families**: 0 rows
- [x] **Not validly published or unknown status, but kept as clinically relevant**: 21 rows, all plausible (e.g.
  *Tropheryma whipplei*, *Neoehrlichia mikurensis*, *Rickettsia felis*, *Mycobacterium lepromatosis*), see
  [CSV](not_validly_published_or_unknown_status_but_kept_as_clinically_relevant.csv)
- [x] **Duplicate full names after adding missing parents**: 8 rows, inherited from v3.0.1 (labyrinthulids in both
  Fungi and Chromista, Microsporea, Myxomycetes), resolved by source priority, none clinically relevant, see
  [CSV](duplicate_full_names_after_adding_missing_parents.csv)
- [x] **Released MO codes that are retired, since their taxon is another organism with the same genus name**: 477
  rows: *Graphium* butterflies 406 (248 of them `NEEDS REVIEW`, script change 9), *Capillaria* fungi 20, sea star
  *Nectria* 18, fungal *Morganella* 8, scale insect *Cryptococcus* 6, fungal *Necator* 5, and 14 others, see
  [CSV](released_mo_codes_that_are_retired_since_their_taxon_is_another_organism_with_the_same_genus_name.csv)
- [x] **Released MO codes that are retired, since a more recent release has their name in another domain** (new): 17
  rows, all `NEEDS REVIEW`, v2.x Protozoa codes of microsporidia whose name v3.0.x has for another fungus (e.g.
  *Caudospora*, *Campanulospora*, *Octospora*); the v2.x codes of the same microsporidia are translated by name, see
  [CSV](released_mo_codes_that_are_retired_since_a_more_recent_release_has_their_name_in_another_domain.csv)
- [x] **Released MO codes that are retired, since their genus is now only known in another domain**: 0 rows
- [x] **Released taxa that were missing and are restored**: 15,643 taxa of v3.0.1 that are not in the new selection (or
  not in the sources anymore) are restored with their old code and status (README rule 4), of which 492 clinically
  relevant (almost all bacterial species absent from current LPSN, e.g. *Neisseria bergeri*,
  *Streptococcus halitosis*), see [CSV](released_taxa_that_were_missing_and_are_restored_by_domain_rank_and_status.csv)
  and [CSV](restored_released_taxa_that_are_clinically_relevant.csv)
- [x] **Orthographic variants linked to their accepted name**: 982 rows, sample of 25 correct (Latin gender
  agreement), one in the wrong direction under "Needs a human", see
  [CSV](orthographic_variants_linked_to_their_accepted_name.csv)
- [x] **Genera without a family that get a family from GBIF or the override list** (new): 204 rows, e.g. *Ascaris*
  Ascarididae, *Balamuthia* Balamuthiidae, *Lacazia* Ajellomycetaceae, see
  [CSV](genera_without_a_family_that_get_a_family_from_gbif_or_the_override_list.csv)
- [x] **Cycles of synonyms**: 10 pairs (*Conidiobolus*/*Capillidium*, *Ochroconis*/*Scolecobasidium*), decided by
  `synonym_cycle_current_genera` (script change 18), see [CSV](cycles_of_synonyms.csv)
- [x] **Synonym genera with accepted species, these will become accepted**: 239 rows, plausible (relevant ones:
  *Debaryozyma*, *Hansenula*, *Ochroconis*, *Pseudallescheria*, *Saprochaete*), see
  [CSV](synonym_genera_with_accepted_species_these_will_become_accepted.csv)
- [x] **Salmonella species synonyms that are serovars** (new): 3 rows removed (script change 13), see
  [CSV](salmonella_species_synonyms_that_are_serovars_in_this_data_set_these_will_be_removed.csv)
- [x] **Lichens and their synonyms, these will be removed** (new): 4,155 records (script change 17), see
  [CSV](lichens_and_their_synonyms_these_will_be_removed.csv)
- [x] **Duplicate full names / duplicate MO codes / MO codes with repeated elements / records without a valid MO
  code**: all 0 rows
- [x] **New taxa** and **Removed taxa**: 29,126 and 492 rows, see [CSV](new_taxa.csv) and [CSV](removed_taxa.csv)
- [x] **Removed taxa that were clinically relevant**: 207 rows, all retired codes of other organisms with the same genus
  name (*Graphium* 164, *Capillaria* 11, *Nectria* 9, *Morganella* 8, *Cryptococcus* 6, *Necator* 5, and 4 others), no
  pathogen lost, see [CSV](removed_taxa_that_were_clinically_relevant_prevalence_2.csv)
- [x] **Genera that moved to another domain**: 254 rows: the microsporidia to Fungi (1,394 records, script change
  19), *Capillaria* and *Necator* to Animalia (override list), the labyrinthulids (*Diplophrys*, *Oblongichytrium*,
  *Sorodiplophrys*) to Chromista, *Actinomycodium* and *Tubercularia* to Fungi, see
  [CSV](genera_that_moved_to_another_domain.csv)
- [x] **Relevant genera with a new or empty family**: 188 rows; the empty families are mostly fungal form genera
  without a family (*incertae sedis*) in MycoBank; exceptions under "Needs a human", see
  [CSV](relevant_genera_with_a_new_or_empty_family.csv)
- [x] **Previously manually added taxa that are not in the new data set**: 13 rows: 11 *Graphium* butterflies
  (intended), *Microsphaera penicillata*, and the species group *Mycobacterium avium-intracellulare complex* (approved
  rename to *M. avium complex*), see [CSV](previously_manually_added_taxa_that_are_not_in_the_new_data_set.csv)
- [x] **Synonyms without a current name**: 3,909 rows (names that `as.mo()` cannot update), see
  [CSV](synonyms_without_a_current_name.csv)
- [x] **Codes without a match in the other data sets** (new): 0 for `clinical_breakpoints`, `example_isolates`,
  `microorganisms.codes`, `microorganisms.groups`; 69 for `intrinsic_resistant` (non-bacteria that the development
  version had as bacteria, e.g. fungal *Bogoriella*, *Microsphaera*, *Morganella*), these rows are dropped, see
  [CSV](codes_without_a_match_in_amr_intrinsic_resistant_mo.csv)

### Sentinel organisms (3b)

- [x] **All 32 sentinels pass**: all named organisms present with the expected domain, family, status and code prefix
  (e.g. *Candida auris* synonym of *Candidozyma auris*, *Trichophyton rubrum* `F_TRCHP_RBRM` in Arthrodermataceae,
  *Necator americanus* and *Capillaria* in Animalia, *Plasmodium falciparum* `P_`, *Ascaris lumbricoides* in
  Ascarididae), *Graphium* only fungal, *Microsporidium* more than 100 species; all codes of the package's data sets
  exist, see [CSV](sentinel_organisms.csv) and the check script [R](checks_sentinels_comparison.R)

### Comparison with v3.0.1 (3c)

- [x] **Counts per domain, rank and status**: see section 2 and [CSV](comparison_per_domain_rank_status.csv)
- [x] **MO codes of existing names that changed**: 0, see [CSV](changed_mo_codes_of_existing_names.csv)
- [x] **No new taxon received a registered code of another taxon**: only `B_MYCBC_AVIM-C`, the approved rename in
  `mo_code_renames.csv`, see [CSV](registered_codes_with_another_name.csv)
- [x] **Prevalence by rank**, old vs new, see [CSV](prevalence_by_rank.csv)
- [x] **Most relevant new and removed genera and species** (prevalence up to 1.25): 7,021 new, 207 removed (see above),
  see [CSV](most_relevant_new_taxa.csv) and [CSV](most_relevant_removed_taxa.csv)

### Integrity rules and MO code registry (3d)

- [x] **Integrity rules and registry**: the build passes all rules of `tests/testthat/helper-microorganisms.R`; 494
  codes retired in `mo_code_retirements.csv` (265 `NEEDS REVIEW`, each listed in the review CSV files above, with
  the reason); no new rows in `mo_code_renames.csv`
- [x] **Known defects**: all rows of `microorganisms_known_defects.csv` are solved, so the file is deleted

### Unit tests (3e)

- [x] **Unit tests**: `devtools::test()` passes completely (after adding the three new staphylococci to CoNS)
- [x] **Changed test**, `test-mo_property.R`: Gram-positive classes *Limnocylindria* (new class in Chloroflexota) and
  *Tissierellia* added
- [x] **Changed test**, `test-mo_property.R`: `mo_pathogenicity(example_isolates$mo)` counts from 1915/62/1/22 to
  1912/71/1/16: *Nakaseomyces glabratus* (6 isolates) from 'Unknown' to 'Potentially pathogenic' (a correction),
  *Aerococcus urinaeequi* (2) and *Paenibacillus durus* (1) from 'Pathogenic' to 'Potentially pathogenic' (not in
  Bartlett *et al.*; on `main` only via synonym links in older LPSN data)
- [x] **Changed test**, `test-data-microorganisms.R`: the registry test used `paste0()` on possibly empty vectors, which
  gives `" ()"` and so a false defect when no code is wrong; now `sprintf()`, as elsewhere in that file
- [x] **Changed test**, `test-data-microorganisms.R`: `F_TRCHP_RBRM` is the current code of *Trichophyton rubrum* again,
  so no warning about an earlier version is expected anymore; the name is still checked

## 5. Needs a human

- [x] **Retired codes with `NEEDS REVIEW`** (265): 248 *Graphium* names that no current source has as a fungus (sample of
  40: all swallowtails), and 17 v2.x Protozoa codes whose name v3.0.x has in Fungi for another organism (e.g.
  *Octospora*, *Caudospora*). Accepted (the v2.x codes of the same microsporidia are now translated by name instead,
  script change 19)
- [x] **Retired as 'another organism' although the same organisms reclassified**: *Rozella* (2) and *Nephridiophaga*
  (1) v2.x Protozoa codes; *Thermus profundus*
  (`B_THERMS_PRFN`, no longer in LPSN, retired because of an Archaea namesake); *Hymenolepis leptocephala*
  (`AN_HYMNL_LPTC`, retired because of a plant namesake). Accepted, low impact
- [x] **Released taxa restored with their status of v3.0.1** (15,643, of which 492 clinically relevant): names such as
  *Neisseria bergeri* and *Streptococcus halitosis* are in no current source but keep status 'accepted' (README rule 4).
  Kept as they are
- [x] **Domain decided by the last release**, doubtful cases: *Copromonas* (a euglenoid, released as Bacteria, no
  authoritative source, stays Bacteria); *Aplanochytrium* (a labyrinthulid, Chromista in COL, stays Fungi as released).
  Kept as they are
- [x] **Microsporidia are Fungi** (decided, script change 19): all microsporidia are now in the Fungi; 1,394 names
  of v3.0.1 change from `P_` to `F_` (e.g. *Encephalitozoon intestinalis* `P_ENCPH_INTS` to `F_ENCPH_INTS`), the `P_`
  codes are translated by name, see [CSV](names_that_moved_to_another_domain.csv). *Enterocytozoon bieneusi* and
  *Encephalitozoon cuniculi* are the accepted names (script change 20)
- [x] **Synonym cycles** (10 pairs): decided, *Capillidium* and *Scolecobasidium* are the current genera (script change
  18), see [CSV](cycles_of_synonyms.csv)
- [x] **Orthographic variant in the wrong direction**: *Wickerhamomyces anomalus anomalus* is linked to
  '*W. anomala anomala*', while *Wickerhamomyces* is masculine (source data). Accepted for now
- [x] **Empty families** of relevant genera that I did not fill as I am not certain: *Sappinia* (v3.0.1: Stenamoebidae)
  and *Fenollaria* (no family in LPSN). Left empty
- [x] **Relevant genera without species** (decided): `MO_RELEVANT_GENERA` listed the newer genera *Balantioides*
  and *Cystoisospora*, while COL has the species in *Balantidium* and *Isospora*, so these pathogens were missing and
  `as.mo()` returned *Campylobacter pinnipediorum caledonicus* for *Balantidium coli* and *Stenotrophomonas beteli* for
  *Isospora belli*. Now *Balantidium* and *Isospora* are in `MO_RELEVANT_GENERA`, and *Cystoisospora belli* (not in COL)
  is added as synonym of *Isospora belli* (script change 21); all four names are found
- [x] **`TAXONOMY_VERSION` dates**: MycoBank is set to the file date 2026-01-07, GBIF/COL to the file date of `COL.zip`
  (2026-04-30); the COL citation is unchanged, as `main` already cited release 2026-04-18 XR. Confirmed
- [x] **Stale file**: `data-raw/taxonomy_lpsn0.rds` dates from June 2026 and is no longer written by the script;
  removed
- [x] **One commit instead of five**: the pre-commit hook of the repository stages all changed files in `data-raw/` and
  `man/`, so the first commit took everything; splitting would need `--no-verify`, which I did not use. As PRs are
  squash-merged, this only affects the review per commit
- [x] **Scope decisions of 5 October 2026** (implemented as script changes 16 to 21): COL XR names only if clinically
  relevant or released before (40,897 fewer COL records), Protozoa with the relevance rules of the Fungi, fungal synonyms
  kept, lichens removed, microsporidia in the Fungi

## 6. Second round (6 and 7 October 2026)

Reason: `synonyms_without_a_current_name.csv` was not empty while all tests passed. LPSN showed that e.g. *Eggerthella
lenta*, *Gordonia amarae*, *Budvicia aquatica* and *Burkholderia pyrrocinia* are correct names (and COL has *Anisakis
simplex* and *Enterobius vermicularis* as accepted), while this build had them as synonyms without a current name.

### Script changes

- [x] **Cause (22)**: a source can have the same name twice (LPSN: the correct name and an illegitimate homotypic
  synonym with the same spelling; COL: an accepted name with a subgenus, e.g. *Enterobius (Enterobius) vermicularis*,
  next to the synonym *Enterobius vermicularis*). The record with the lowest identifier was kept, now the accepted
  one, and links to the dropped record move to the kept one (`one_record_per_name()`). LPSN: 34 dangling links to 0
- [x] **LPSN lookup (23)**: `get_lpsn_and_author()` read only the nomenclatural status, so a validly published
  synonym (e.g. *Eubacterium lentum*) became 'accepted'; it now also reads the taxonomic status and the correct name,
  and `apply_lpsn_result()` links such a synonym. Subspecies URLs fixed. The LPSN cache was renewed for this
- [x] **Synonyms without a current name (24, decision of 7 October)**: prokaryotes are looked up in LPSN (linked,
  accepted, or not validly published and then only kept if protected); all others are removed and their released
  codes retired with the reason. Records with children that are kept become 'unknown'. Outcome: 3,827 records, 3,691
  removed, 123 'unknown', 13 kept; 3,668 codes retired, see [CSV](synonyms_without_a_current_name_and_their_outcome.csv)
- [x] **_M. tuberculosis_ complex protected (25, decision of 7 October)**: also *M. canettii* and *M. orygis* (LPSN:
  preferred names, not validly published) are kept, as accepted
- [x] **Subspecies codes reserved (26)**: registered subspecies codes were not reserved, so a removed subspecies could
  pass its code to another spelling (*Candida melibiosi membranaefaciens*); now as for genera and species
- [x] **Childless genera (27)**: a genus that is the current name of a synonym is kept (*Monotosporella*)
- [x] **Missing parent records (28)**: a referenced family, order, class or phylum is added or moved to the domain of
  its members (*Plasmodiidae* from Chromista to Protozoa), see [CSV](missing_parent_records_that_were_added_or_moved.csv)
- [x] **Higher taxonomy harmonised (29)**: every record takes its higher taxonomy from its parent record, outdated
  parent names are replaced by their current name (5,287 records), see
  [CSV](records_whose_higher_taxonomy_was_harmonised_with_the_record_of_their_parent.csv)
- [x] **Shared identifiers (30)**: records of the same domain and rank with the same identifier become one taxon (the
  name in a current source, else the spelling LPSN knows, else alphabetical: 14 marked NEEDS REVIEW); records of
  another domain or rank only lose the identifier (GBIF 439 was both *Septobasidiales* and the nematode order
  *Strongylida*), see [CSV](records_that_shared_a_source_identifier_one_kept_as_the_current_name.csv)
- [x] **Pointers (31)**: only synonyms have a 'renamed to' identifier (1,734 MycoBank-accepted records kept a COL
  pointer), parent identifiers to the record itself or to a lower rank are removed
- [x] **Current names in data sets (32, decision of 7 October)**: `microorganisms.codes`, `clinical_breakpoints`,
  `example_isolates` and `microorganisms.groups` refer only to current names (`intrinsic_resistant` lists synonyms on
  purpose); the groups script takes group members only from current names (otherwise e.g. *Gardnerella vaginalis*
  entered HACEK as 'Haemophilus vaginalis', and *S. aureus* entered CoPS as '*S. roterodami*'); `MO_CONS`, `MO_COPS`,
  `MO_STREP_ABCG` and `MO_LANCEFIELD` keep a synonym only if its current name is in the list too; 'Teleomorph' removed
  from `MO_RELEVANT_GENERA`
- [x] **Tests**: 28 new blocks in `tests/testthat/test-data.R` (synonyms, pointers, identifiers, parents, hierarchy,
  names, sources, clinically important names and higher taxa with their domain, outdated names, regressions, current
  names in data sets, internal lists, `intrinsic_resistant`, interpretive rules, groups, codes, breakpoints,
  antimicrobials). On the data before this round they fail on exactly these defects; now `devtools::test()`: 1,581
  expectations, 0 failures

### Needs a human

- [x] ***Klebsiella quasivariicola*** and other preferred names that are not validly published: decided, all taxa in the
  genera of `MO_WHO_PRIORITY_GENERA` are protected (script change 33, replaces the hard-coded *M. tuberculosis*
  complex); *K. quasivariicola*, *M. canettii*, *M. orygis*, *E. massiliensis*, *M. liflandii*, *P. aestus*, *P.
  stewarti* and *S. periodonticum* are kept as accepted
- [x] ***M. bovis*, *M. africanum*, *M. caprae*, *M. microti*, *M. pinnipedii***: decided, LPSN is followed, these are
  heterotypic synonyms of *M. tuberculosis* (Riojas et al. 2018); a test enforces this. Considered and rejected: keeping
  heterotypic synonyms in the WHO genera separate (98 names, e.g. *Citrobacter diversus*, *Salmonella enteritidis*,
  which would then lose the breakpoints and rules of their current names)
- [x] ***Giardia duodenalis* and *Giardia lamblia*** are both accepted (COL): kept as they are
- [x] **Protist groups in more than one domain**: decided (script change 34): *Aplanochytrium* in the Chromista (with
  the labyrinthulids), *Dictydiaethalium* in the Protozoa (with the *Myxomycetes*), and a name at two ranks gets the
  rank suffix also between domains ('*Acantharia* {class}' next to the fungal genus *Acantharia*)
- [x] **Spelling chosen alphabetically** for 14 pairs with a shared identifier: accepted (Latin rules would be better,
  but the impact is small)
- [x] ***Ameson michaeli*** (released) is removed and retired, while *A. michaelis* exists: accepted
- [x] ***Macrococcus caseolyticus*** is no longer CoNS (also no EUCAST cefoxitin screening row): accepted
- [ ] **Relevance drifts between builds**: the relevant genera include the genera of `intrinsic_resistant`, which lists
  every code of the previous data set, so that a build partly keeps what the previous one had (569 never released COL
  insects such as *Anisoptera* dropped out after `intrinsic_resistant` was regenerated)
