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

**111,103 records** (v3.0.1: 78,679; development version on `main`: 96,982). Per domain, old (v3.0.1) vs new, see
also [CSV](comparison_per_domain_rank_status.csv) and [CSV](new_records_per_domain_and_kind.csv):

| Domain | Total old | Total new | Accepted old | Accepted new | Synonym old | Synonym new | Unknown/other old / new | Change |
|---|---|---|---|---|---|---|---|---|
| Animalia | 1,628 | 2,579 | 1,379 | 1,754 | 244 | 810 | 5 / 15 | +58.4% |
| Archaea | 1,419 | 1,573 | 1,225 | 1,359 | 181 | 201 | 13 / 13 | +10.9% |
| Bacteria | 39,249 | 48,239 | 29,853 | 38,689 | 7,076 | 8,777 | 2,320 / 773 | +22.9% |
| Chromista | 178 | 300 | 157 | 236 | 17 | 57 | 4 / 7 | +68.5% |
| Fungi | 28,137 | 49,716 | 14,385 | 21,877 | 8,746 | 22,585 | 5,006 / 5,254 | +76.7% |
| Protozoa | 8,067 | 8,694 | 6,056 | 6,458 | 1,880 | 2,098 | 131 / 138 | +7.8% |

Explanation of the changes of more than 5%:

- **Bacteria (+22.9%)**: new validly published species in LPSN since 2024 (accepted +8,836), and v3.0.1 'not validly
  published' names now as accepted, synonym or unknown.
- **Fungi (+76.7%)**: mostly synonyms (+13,839): old names of the kept fungal species, from MycoBank and COL, which
  `as.mo()` uses to translate old names. Accepted +7,492: new species in relevant genera and current names of species of
  relevant genera (script change 5). This also includes about 1,000 lichens, see "Needs a human".
- **Animalia (+58.4%)**, **Chromista (+68.5%)**: helminths, vectors and protists of the relevant genera that COL now
  contains (and the COL records without authors, script change 3); small absolute numbers (+951 and +122).
- **Archaea (+10.9%)**: new species in LPSN.
- Protozoa (+7.8%) is below 5% after the decision to apply the Fungi relevance rules to Protozoa; most protozoa of
  v3.0.1 are restored as released taxa.

Records by rank: species 85,841, subspecies 10,507, genus 11,417, family 2,097, order 728, class 317, phylum 129, and
37 species groups. MO codes of existing names (same name and domain) that changed: **0**
([CSV](changed_mo_codes_of_existing_names.csv)). Codes in the package's other data sets missing from `microorganisms`:
**0**.

Run time: the source chunks take about 15 minutes (LPSN scrape 10.5 minutes, ran once), each full rebuild from a
checkpoint 15 to 25 minutes. 14 runs were needed because of the script changes below (logs available on request).

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
- `tests/testthat/test-mo_property.R` and `test-data-microorganisms.R`: see the checklist items on the unit tests.
- `R/aa_globals.R`: accessed dates in `TAXONOMY_VERSION`; `NEWS.md`: the existing taxonomy bullet of this series updated.

## 4. Checklist

### Review blocks of the build (3a)

- [ ] **LPSN genera with multiple lineages**: 6 rows, the homonyms *Halalkalibacterium*, *Pusillimonas* and
  *Rhodococcus*; the script removes the Balneolaceae, Oscillospiraceae and Chroococcaceae lineages explicitly, the kept
  families match v3.0.1, see [CSV](lpsn_genera_with_multiple_lineages.csv)
- [ ] **LPSN higher taxa that could not be downloaded**: 0 rows
- [ ] **Homonyms of relevant genera**: 220 rows (about 100 genera in more than one domain); after script change 4 all
  winners plausible, e.g. *Plasmodium*, *Sarcocystis*, *Babesia* Chromista (COL; *Plasmodium* is set to Protozoa by
  the curated block later on), *Moorella*, *Serpula*, *Chainia*, *Pirella* Bacteria, *Enterocytozoon* Protozoa,
  *Necator*, *Capillaria*, *Hymenolepis* Animalia (override list), *Trichophyton*, *Graphium*, *Cryptococcus*,
  *Mucor* Fungi, see [CSV](homonyms_of_relevant_genera_check_best_domain.csv)
- [ ] **Homonyms decided by the domain of the last release** (new): 27 rows, of which the release decides where an
  authoritative source supports it (column `released_decides`), see
  [CSV](homonyms_decided_by_the_domain_of_the_last_release.csv); doubtful ones under "Needs a human"
- [ ] **Records removed as homonyms of a genus in another domain**: 179 rows, other organisms with the same genus name
  (e.g. the polychaete *Serpula*, fungal *Moorella*, *Necator*, the unclassified MycoBank index entries of *Plasmodium*
  and *Sarcocystis*), see [CSV](records_removed_as_homonyms_of_a_genus_in_another_domain.csv)
- [ ] **Genera with multiple families**: 0 rows
- [ ] **Not validly published or unknown status, but kept as clinically relevant**: 21 rows, all plausible (e.g.
  *Tropheryma whipplei*, *Neoehrlichia mikurensis*, *Rickettsia felis*, *Mycobacterium lepromatosis*), see
  [CSV](not_validly_published_or_unknown_status_but_kept_as_clinically_relevant.csv)
- [ ] **Duplicate full names after adding missing parents**: 13 rows, inherited from v3.0.1 (labyrinthulids in both
  Fungi and Chromista, Microsporea, Myxomycetes), resolved by source priority, none clinically relevant, see
  [CSV](duplicate_full_names_after_adding_missing_parents.csv)
- [ ] **Released MO codes that are retired, since their taxon is another organism with the same genus name**: 495
  rows: *Graphium* butterflies 406 (248 of them `NEEDS REVIEW`, script change 9), *Capillaria* fungi 20, sea star
  *Nectria* 18, microsporidian *Caudospora* 15, fungal *Morganella* 8, scale insect *Cryptococcus* 6, fungal *Necator* 5,
  and 17 single cases, see
  [CSV](released_mo_codes_that_are_retired_since_their_taxon_is_another_organism_with_the_same_genus_name.csv)
- [ ] **Released MO codes that are retired, since a more recent release has their name in another domain** (new): 27
  rows, all `NEEDS REVIEW`, v2.x Protozoa codes of microsporidia that v3.0.x had in Fungi (e.g. *Chytridiopsis*,
  *Octospora*, *Nolleria pulicis*), see
  [CSV](released_mo_codes_that_are_retired_since_a_more_recent_release_has_their_name_in_another_domain.csv)
- [ ] **Released MO codes that are retired, since their genus is now only known in another domain**: 0 rows
- [ ] **Released taxa that were missing and are restored**: 15,323 taxa of v3.0.1 that are not in the new selection (or
  not in the sources anymore) are restored with their old code and status (README rule 4), of which 484 clinically
  relevant (almost all bacterial species absent from current LPSN, e.g. *Neisseria bergeri*,
  *Streptococcus halitosis*), see [CSV](released_taxa_that_were_missing_and_are_restored_by_domain_rank_and_status.csv)
  and [CSV](restored_released_taxa_that_are_clinically_relevant.csv)
- [ ] **Orthographic variants linked to their accepted name**: 982 rows, sample of 25 correct (Latin gender
  agreement), one in the wrong direction under "Needs a human", see
  [CSV](orthographic_variants_linked_to_their_accepted_name.csv)
- [ ] **Genera without a family that get a family from GBIF or the override list** (new): 199 rows, e.g. *Ascaris*
  Ascarididae, *Balamuthia* Balamuthiidae, *Lacazia* Ajellomycetaceae, see
  [CSV](genera_without_a_family_that_get_a_family_from_gbif_or_the_override_list.csv)
- [ ] **Cycles of synonyms**: 20 rows (10 pairs, *Conidiobolus*/*Capillidium*, *Ochroconis*/*Scolecobasidium*), both
  names made accepted by the script, see [CSV](cycles_of_synonyms.csv) and "Needs a human"
- [ ] **Synonym genera with accepted species, these will become accepted**: 239 rows, plausible (relevant ones:
  *Debaryozyma*, *Hansenula*, *Ochroconis*, *Pseudallescheria*, *Saprochaete*), see
  [CSV](synonym_genera_with_accepted_species_these_will_become_accepted.csv)
- [ ] **Salmonella species synonyms that are serovars** (new): 3 rows removed (script change 13), see
  [CSV](salmonella_species_synonyms_that_are_serovars_in_this_data_set_these_will_be_removed.csv)
- [ ] **Duplicate full names / duplicate MO codes / MO codes with repeated elements / records without a valid MO
  code**: all 0 rows
- [ ] **New taxa** and **Removed taxa**: 32,953 and 500 rows, see [CSV](new_taxa.csv) and [CSV](removed_taxa.csv)
- [ ] **Removed taxa that were clinically relevant**: 207 rows, all retired codes of other organisms with the same genus
  name (*Graphium* 164, *Capillaria* 11, *Nectria* 9, *Morganella* 8, *Cryptococcus* 6, *Necator* 5, and 4 others), no
  pathogen lost, see [CSV](removed_taxa_that_were_clinically_relevant_prevalence_2.csv)
- [ ] **Genera that moved to another domain**: 14 rows: corrections are *Capillaria* and *Necator* to Animalia
  (override list), the labyrinthulids (*Diplophrys*, *Oblongichytrium*, *Sorodiplophrys*) to Chromista, *Actinomycodium*
  and *Tubercularia* to Fungi; debatable is the microsporidian *Pleistosporidium* to Protozoa (as in COL), see "Needs a
  human" and [CSV](genera_that_moved_to_another_domain.csv)
- [ ] **Relevant genera with a new or empty family**: 188 rows; the empty families are mostly fungal form genera
  without a family (*incertae sedis*) in MycoBank; exceptions under "Needs a human", see
  [CSV](relevant_genera_with_a_new_or_empty_family.csv)
- [ ] **Previously manually added taxa that are not in the new data set**: 13 rows: 11 *Graphium* butterflies
  (intended), *Microsphaera penicillata*, and the species group *Mycobacterium avium-intracellulare complex* (approved
  rename to *M. avium complex*), see [CSV](previously_manually_added_taxa_that_are_not_in_the_new_data_set.csv)
- [ ] **Synonyms without a current name**: 3,878 rows (names that `as.mo()` cannot update), see
  [CSV](synonyms_without_a_current_name.csv)
- [ ] **Codes without a match in the other data sets** (new): 0 for `clinical_breakpoints`, `example_isolates`,
  `microorganisms.codes`, `microorganisms.groups`; 69 for `intrinsic_resistant` (non-bacteria that the development
  version had as bacteria, e.g. fungal *Bogoriella*, *Microsphaera*, *Morganella*), these rows are dropped, see
  [CSV](codes_without_a_match_in_amr_intrinsic_resistant_mo.csv)

### Sentinel organisms (3b)

- [ ] **All 32 sentinels pass**: all named organisms present with the expected domain, family, status and code prefix
  (e.g. *Candida auris* synonym of *Candidozyma auris*, *Trichophyton rubrum* `F_TRCHP_RBRM` in Arthrodermataceae,
  *Necator americanus* and *Capillaria* in Animalia, *Plasmodium falciparum* `P_`, *Ascaris lumbricoides* in
  Ascarididae), *Graphium* only fungal, *Microsporidium* more than 100 species; all codes of the package's data sets
  exist, see [CSV](sentinel_organisms.csv) and the check script [R](checks_sentinels_comparison.R)

### Comparison with v3.0.1 (3c)

- [ ] **Counts per domain, rank and status**: see section 2 and [CSV](comparison_per_domain_rank_status.csv)
- [ ] **MO codes of existing names that changed**: 0, see [CSV](changed_mo_codes_of_existing_names.csv)
- [ ] **No new taxon received a registered code of another taxon**: only `B_MYCBC_AVIM-C`, the approved rename in
  `mo_code_renames.csv`, see [CSV](registered_codes_with_another_name.csv)
- [ ] **Prevalence by rank**, old vs new, see [CSV](prevalence_by_rank.csv)
- [ ] **Most relevant new and removed genera and species** (prevalence up to 1.25): 6,951 new, 207 removed (see above),
  see [CSV](most_relevant_new_taxa.csv) and [CSV](most_relevant_removed_taxa.csv)

### Integrity rules and MO code registry (3d)

- [ ] **Integrity rules and registry**: the build passes all rules of `tests/testthat/helper-microorganisms.R`; 522
  codes retired in `mo_code_retirements.csv` (275 `NEEDS REVIEW`, each listed in the review CSV files above, with
  the reason); no new rows in `mo_code_renames.csv`
- [ ] **Known defects**: all rows of `microorganisms_known_defects.csv` are solved, so the file is deleted

### Unit tests (3e)

- [ ] **Unit tests**: `devtools::test()` passes completely (after adding the three new staphylococci to CoNS)
- [ ] **Changed test**, `test-mo_property.R`: Gram-positive classes *Limnocylindria* (new class in Chloroflexota) and
  *Tissierellia* added
- [ ] **Changed test**, `test-mo_property.R`: `mo_pathogenicity(example_isolates$mo)` counts from 1915/62/1/22 to
  1912/71/1/16: *Nakaseomyces glabratus* (6 isolates) from 'Unknown' to 'Potentially pathogenic' (a correction),
  *Aerococcus urinaeequi* (2) and *Paenibacillus durus* (1) from 'Pathogenic' to 'Potentially pathogenic' (not in
  Bartlett *et al.*; on `main` only via synonym links in older LPSN data)
- [ ] **Changed test**, `test-data-microorganisms.R`: the registry test used `paste0()` on possibly empty vectors, which
  gives `" ()"` and so a false defect when no code is wrong; now `sprintf()`, as elsewhere in that file
- [ ] **Changed test**, `test-data-microorganisms.R`: `F_TRCHP_RBRM` is the current code of *Trichophyton rubrum* again,
  so no warning about an earlier version is expected anymore; the name is still checked

## 5. Needs a human

- [ ] **Fungal synonyms and lichens**: Fungi grow by 13,839 synonyms (old names of kept species). Kept as they are,
  following your decision (useful for `as.mo()`, no extra layer of assessment). About 1,000 lichen records (*Lecidea*
  558, *Lecanora* 444) enter as current names of kept synonyms, with their own synonyms; recommend to accept for now
- [ ] **Retired codes with `NEEDS REVIEW`** (275): 248 *Graphium* names that no current source has as a fungus (sample of
  40: all swallowtails), and 27 v2.x Protozoa codes of microsporidia whose name v3.0.x has in Fungi. Recommend to accept;
  for the microsporidia, the alternative is to keep them translated by name (`as.mo()` would then give the Fungi record)
- [ ] **Retired as 'another organism' although the same organisms reclassified**: *Rozella* (2) and *Nephridiophaga* (1)
  v2.x Protozoa codes (now Fungi); *Thermus profundus* (`B_THERMS_PRFN`, no longer in LPSN, retired because of an
  Archaea namesake); *Hymenolepis leptocephala* (`AN_HYMNL_LPTC`, retired because of a plant namesake). Recommend to
  accept, low impact
- [ ] **Released taxa restored with their status of v3.0.1** (15,323, of which 484 clinically relevant): names such as
  *Neisseria bergeri* and *Streptococcus halitosis* are in no current source but keep status 'accepted' (README rule 4).
  Decide whether names absent from all sources should become 'unknown'
- [ ] **Domain decided by the last release**, doubtful cases: *Copromonas* (a euglenoid, released as Bacteria, no
  authoritative source, stays Bacteria); *Aplanochytrium* (a labyrinthulid, Chromista in COL, stays Fungi as released)
- [ ] **Microsporidia**: *Pleistosporidium* moved from Fungi to Protozoa as in COL (codes change from `F_` to `P_`);
  current consensus places microsporidia with the fungi
- [ ] **Synonym cycles** (10 pairs): both names are accepted now; recommend *Capillidium* and *Scolecobasidium* as the
  current genera
- [ ] **Orthographic variant in the wrong direction**: *Wickerhamomyces anomalus anomalus* is linked to
  '*W. anomala anomala*', while *Wickerhamomyces* is masculine (source data)
- [ ] **Empty families** of relevant genera that I did not fill as I am not certain: *Sappinia* (v3.0.1: Stenamoebidae)
  and *Fenollaria* (no family in LPSN)
- [ ] **Relevant genera without species**: *Balantioides coli* is only a synonym in COL, *Cystoisospora* is a synonym
  genus in COL without *C. belli*; consider adding them like *Brugia* (script change 12)
- [ ] **`TAXONOMY_VERSION` dates**: MycoBank is set to the file date 2026-01-07 (probably the server date of
  `MBList.zip`, not the download date), GBIF/COL to the file date of `COL.zip` (2026-04-30); the COL citation is
  unchanged, as `main` already cited release 2026-04-18 XR
- [ ] **Stale file**: `data-raw/taxonomy_lpsn0.rds` dates from June 2026 and is no longer written by the script; left
  untouched
- [ ] **One commit instead of five**: the pre-commit hook of the repository stages all changed files in `data-raw/` and
  `man/`, so the first commit took everything; splitting would need `--no-verify`, which I did not use. As PRs are
  squash-merged, this only affects the review per commit
- [ ] **Scope decisions of 5 October 2026** (yours, implemented as script changes): COL XR names only if clinically
  relevant or released before (40,897 fewer COL records), and Protozoa with the relevance rules of the Fungi
