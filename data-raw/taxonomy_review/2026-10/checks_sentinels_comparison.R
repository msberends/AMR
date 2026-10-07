# Sentinel and comparison checks of the new `microorganisms` data set (steps 3b and 3c)
suppressMessages({
  library(dplyr)
  devtools::load_all(".", quiet = TRUE)
})
review_dir <- "data-raw/taxonomy_review/2026-10"
scratch <- "/tmp/claude-1000/-home-uscloud-AMR/cb1d4493-f2a9-47e9-b397-cf7900d03e5f/scratchpad"

env <- new.env()
load(file.path(scratch, "microorganisms_old.rda"), envir = env)
old <- env$microorganisms
old$domain <- old$kingdom
new <- AMR::microorganisms
old$mo <- as.character(old$mo)
new$mo <- as.character(new$mo)

# 3b. sentinels ----
get_row <- function(name) {
  x <- new[new$fullname == name, , drop = FALSE]
  if (nrow(x) == 0) {
    return(data.frame(name = name, present = FALSE))
  }
  x <- x[order(x$status != "accepted"), ][1, ]
  current <- tryCatch(suppressWarnings(suppressMessages(mo_name(mo_current(name)))), error = function(e) NA_character_)
  data.frame(
    name = name, present = TRUE, mo = x$mo, status = x$status, domain = x$domain, kingdom = x$kingdom,
    family = x$family, rank = x$rank, prevalence = x$prevalence, current_name = current,
    n_records = sum(new$fullname == name)
  )
}
expect <- tibble::tribble(
  ~name, ~exp_domain, ~exp_family, ~exp_status, ~exp_prefix, ~exp_max_prev,
  "Escherichia coli", "Bacteria", "Enterobacteriaceae", "accepted", "B_", NA,
  "Klebsiella pneumoniae", "Bacteria", "Enterobacteriaceae", "accepted", "B_", NA,
  "Staphylococcus aureus", "Bacteria", "Staphylococcaceae", "accepted", NA, NA,
  "Pseudomonas aeruginosa", "Bacteria", "Pseudomonadaceae", NA, NA, NA,
  "Acinetobacter baumannii", "Bacteria", "Moraxellaceae", NA, NA, NA,
  "Enterococcus faecium", "Bacteria", "Enterococcaceae", NA, NA, NA,
  "Streptococcus pneumoniae", "Bacteria", "Streptococcaceae", NA, NA, NA,
  "Mycobacterium tuberculosis", "Bacteria", "Mycobacteriaceae", NA, NA, NA,
  "Clostridioides difficile", "Bacteria", "Peptostreptococcaceae", NA, NA, NA,
  "Tropheryma whipplei", NA, NA, NA, NA, NA,
  "Candida albicans", "Fungi", NA, "accepted", NA, NA,
  "Candida auris", NA, NA, "synonym", NA, NA,
  "Candidozyma auris", "Fungi", NA, "accepted", NA, 1.99,
  "Nakaseomyces glabratus", "Fungi", NA, "accepted", NA, 1.99,
  "Aspergillus fumigatus", "Fungi", "Aspergillaceae", NA, NA, NA,
  "Trichophyton rubrum", "Fungi", "Arthrodermataceae", NA, "F_", NA,
  "Mucor", "Fungi", "Mucoraceae", NA, NA, NA,
  "Blastomyces", "Fungi", NA, "accepted", NA, NA,
  "Blastomyces dermatitidis", "Fungi", NA, "accepted", NA, NA,
  "Cryptococcus neoformans", "Fungi", NA, NA, NA, NA,
  "Pneumocystis jirovecii", "Fungi", NA, NA, NA, NA,
  "Necator americanus", "Animalia", NA, NA, NA, NA,
  "Capillaria", "Animalia", NA, NA, NA, NA,
  "Ascaris lumbricoides", "Animalia", "Ascarididae", NA, NA, NA,
  "Echinococcus granulosus", "Animalia", "Taeniidae", NA, NA, NA,
  "Plasmodium falciparum", "Protozoa", NA, NA, "P_", NA,
  "Toxoplasma gondii", NA, NA, "accepted", NA, NA,
  "Leishmania donovani", "Protozoa", NA, NA, NA, NA,
  "Trypanosoma cruzi", "Protozoa", NA, NA, NA, NA,
  "Giardia", "Protozoa", NA, NA, NA, NA
)
sentinels <- bind_rows(lapply(expect$name, get_row)) %>%
  left_join(expect, by = "name") %>%
  mutate(
    problems = paste0(
      if_else(!present, "missing; ", ""),
      if_else(present & !is.na(exp_domain) & domain != exp_domain, "domain; ", ""),
      if_else(present & !is.na(exp_family) & family != exp_family, "family; ", ""),
      if_else(present & !is.na(exp_status) & status != exp_status, "status; ", ""),
      if_else(present & !is.na(exp_prefix) & !startsWith(mo, exp_prefix), "mo prefix; ", ""),
      if_else(present & !is.na(exp_max_prev) & prevalence > exp_max_prev, "prevalence; ", "")
    ),
    passed = problems == ""
  )
# Candida auris must have a current name in Candidozyma
ca <- sentinels$name == "Candida auris"
if (!isTRUE(sentinels$current_name[ca] %like% "^Candidozyma")) {
  sentinels$problems[ca] <- paste0(sentinels$problems[ca], "current name not in Candidozyma; ")
  sentinels$passed[ca] <- FALSE
}
graphium <- new %>% filter(genus == "Graphium")
graphium_check <- data.frame(
  name = c("Graphium (only fungal species)", "Microsporidium (> 100 species)"),
  present = c(nrow(graphium) > 0, TRUE),
  domain = c(toString(unique(graphium$domain)), "Fungi"),
  n_records = c(nrow(graphium), sum(new$genus == "Microsporidium" & new$rank == "species")),
  problems = c(
    if_else(any(graphium$domain != "Fungi") | "Graphium sarpedon" %in% graphium$fullname, "non-fungal Graphium; ", ""),
    if_else(sum(new$genus == "Microsporidium" & new$rank == "species") <= 100, "too few Microsporidium species; ", "")
  )
)
graphium_check$passed <- graphium_check$problems == ""
sentinels <- bind_rows(sentinels, graphium_check)
write.csv(sentinels, file.path(review_dir, "sentinel_organisms.csv"), row.names = FALSE, na = "")
print(as.data.frame(sentinels %>% select(name, present, mo, status, domain, family, prevalence, current_name, problems)))

# codes in package data sets
used <- c(
  as.character(clinical_breakpoints$mo), as.character(intrinsic_resistant$mo),
  as.character(microorganisms.groups$mo), as.character(microorganisms.groups$mo_group),
  as.character(microorganisms.codes$mo), as.character(example_isolates$mo)
)
cat("\nCodes in package data sets missing from microorganisms:", sum(!unique(used) %in% new$mo), "\n")
print(setdiff(unique(used), new$mo))

# 3c. comparison ----
cmp <- full_join(
  old %>% count(domain, rank, status, name = "n_old"),
  new %>% count(domain, rank, status, name = "n_new"),
  by = c("domain", "rank", "status")
) %>%
  mutate(across(c(n_old, n_new), ~ coalesce(.x, 0L)), diff = n_new - n_old) %>%
  arrange(domain, rank, status)
write.csv(cmp, file.path(review_dir, "comparison_per_domain_rank_status.csv"), row.names = FALSE, na = "")
per_domain <- full_join(
  old %>% count(domain, name = "n_old"),
  new %>% count(domain, name = "n_new"),
  by = "domain"
) %>%
  mutate(across(c(n_old, n_new), ~ coalesce(.x, 0L)), diff = n_new - n_old, pct = round(100 * diff / n_old, 1))
write.csv(per_domain, file.path(review_dir, "comparison_per_domain.csv"), row.names = FALSE, na = "")
cat("\nPer domain:\n")
print(as.data.frame(per_domain))

# MO codes of existing names (same fullname and domain) that changed
changed <- inner_join(
  old %>% select(fullname, domain, mo_old = mo, rank_old = rank),
  new %>% select(fullname, domain, mo_new = mo, rank_new = rank, status_new = status, source_new = source),
  by = c("fullname", "domain"),
  relationship = "many-to-many"
) %>%
  group_by(fullname, domain) %>%
  filter(!any(mo_old == mo_new)) %>%
  ungroup()
write.csv(changed, file.path(review_dir, "changed_mo_codes_of_existing_names.csv"), row.names = FALSE, na = "")
cat("\nChanged MO codes of existing names:", n_distinct(changed$fullname), "\n")
print(head(as.data.frame(changed), 50))

# new taxa must not get a code that was ever used for another taxon
registry <- read.csv("data-raw/microorganisms_files/mo_code_registry.csv", stringsAsFactors = FALSE)
reused <- new %>%
  filter(mo %in% registry$mo) %>%
  left_join(registry %>% select(mo, reg_fullname = fullname), by = "mo") %>%
  filter(fullname != reg_fullname) %>%
  select(mo, fullname, rank, domain, reg_fullname)
write.csv(reused, file.path(review_dir, "registered_codes_with_another_name.csv"), row.names = FALSE, na = "")
cat("\nRegistered codes with another name:", nrow(reused), "\n")
print(head(as.data.frame(reused %>% select(mo, fullname, reg_fullname)), 50))

# prevalence by rank
prev <- bind_rows(
  old %>% mutate(version = "v3.0.1"),
  new %>% mutate(version = "new")
) %>%
  count(version, rank, prevalence) %>%
  tidyr::pivot_wider(names_from = version, values_from = n, values_fill = 0) %>%
  arrange(rank, prevalence)
write.csv(prev, file.path(review_dir, "prevalence_by_rank.csv"), row.names = FALSE, na = "")
cat("\nPrevalence by rank:\n")
print(as.data.frame(prev))

# most relevant new and removed genera and species
new_taxa <- new %>%
  filter(!fullname %in% old$fullname, rank %in% c("genus", "species"), prevalence <= 1.25) %>%
  arrange(prevalence, rank, fullname) %>%
  select(mo, fullname, rank, status, domain, family, source, prevalence)
removed_taxa <- old %>%
  filter(!fullname %in% new$fullname, rank %in% c("genus", "species"), prevalence <= 1.25) %>%
  arrange(prevalence, rank, fullname) %>%
  select(mo, fullname, rank, status, domain, family, source, prevalence)
write.csv(new_taxa, file.path(review_dir, "most_relevant_new_taxa.csv"), row.names = FALSE, na = "")
write.csv(removed_taxa, file.path(review_dir, "most_relevant_removed_taxa.csv"), row.names = FALSE, na = "")
cat("\nNew relevant taxa:", nrow(new_taxa), " removed relevant taxa:", nrow(removed_taxa), "\n")
print(head(as.data.frame(new_taxa), 30))
print(head(as.data.frame(removed_taxa), 30))
