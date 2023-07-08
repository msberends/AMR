g <- certetoolbox::import_clipboard()
g <- g %>% as_tibble()
g
for (i in seq_len(nrow(g))) {
  if (is.na(g$taxonomy_from[i])) next
  mo <- as.mo(g$taxonomy_from[i])
  g$phylum[i] <- mo_phylum(mo)
  g$class[i] <- mo_class(mo)
  g$order[i] <- mo_order(mo)
  g$family[i] <- mo_family(mo)
  g$genus[i] <- mo_genus(mo)
  g$species[i] <- paste(mo_species(mo), "complex")
  g$lpsn_parent[i] <- mo_property(mo, property = "lpsn_parent")
  g$gbif_parent[i] <- mo_property(mo, property = "gbif_parent")
  g$oxygen_tolerance[i] <- mo_property(mo, property = "oxygen_tolerance")
  g$prevalence[i] <- mo_property(mo, property = "prevalence")
}
certetoolbox::export_clipboard(g)
