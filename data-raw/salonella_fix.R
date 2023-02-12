snomed2 <- microorganisms %>%
  filter(mo %in% c("B_SLMNL_TYPH", "B_SLMNL_HMRM", "B_SLMNL_PRTY")) %>%
  pull(snomed)

new_typhi <- microorganisms %>%
  filter(mo == "B_SLMNL_THSS") %>%
  slice(c(1, 1, 1)) %>%
  mutate(
    mo = c("B_SLMNL_TYPH", "B_SLMNL_HMRM", "B_SLMNL_PRTY"),
    fullname = c("Salmonella Typhi", "Salmonella Typhimurium", "Salmonella Paratyphi"),
    subspecies = c("Typhi", "Typhimurium", "Paratyphi"),
    snomed = snomed2
  )

new_groupa <- microorganisms %>%
  filter(mo == "B_SLMNL_GRPB") %>%
  mutate(
    mo = "B_SLMNL_GRPA",
    fullname = gsub("roup B", "roup A", fullname),
    species = gsub("roup B", "roup A", species)
  )

microorganisms$mo <- as.character(microorganisms$mo)

microorganisms <- microorganisms %>%
  filter(!mo %in% c("B_SLMNL_TYPH", "B_SLMNL_HMRM", "B_SLMNL_PRTY")) %>%
  bind_rows(new_typhi, new_groupa) %>%
  arrange(fullname)

microorganisms$lpsn_parent[which(microorganisms$genus == "Salmonella" & microorganisms$rank == "species")] <- "516547"
microorganisms$gbif_parent[which(microorganisms$genus == "Salmonella" & microorganisms$rank == "species")] <- "3221815"

class(microorganisms$mo) <- c("mo", "character")
