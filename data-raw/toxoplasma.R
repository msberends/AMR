microorganisms <- microorganisms |>
  bind_rows(
    # Toxoplasma
    data.frame(
      mo = "P_TXPL_GOND", # species
      fullname = "Toxoplasma gondii",
      kingdom = "(unknown kingdom)",
      phylum = "Apicomplexa",
      class = "Conoidasida",
      order = "Eucoccidiorida",
      family = "Sarcocystidae",
      genus = "Toxoplasma",
      species = "gondii",
      subspecies = "",
      rank = "species",
      ref = "Nicolle et al., 1908",
      species_id = NA_real_,
      source = "manually added",
      prevalence = 2,
      stringsAsFactors = FALSE
    ),
    data.frame(
      mo = "P_TXPL", # genus
      fullname = "Toxoplasma",
      kingdom = "(unknown kingdom)",
      phylum = "Apicomplexa",
      class = "Conoidasida",
      order = "Eucoccidiorida",
      family = "Sarcocystidae",
      genus = "Toxoplasma",
      species = "",
      subspecies = "",
      rank = "genus",
      ref = "Nicolle et al., 1909",
      species_id = NA_real_,
      source = "manually added",
      prevalence = 2,
      stringsAsFactors = FALSE
    ),
    data.frame(
      mo = "[FAM]_SRCCYSTD", # family
      fullname = "Sarcocystidae",
      kingdom = "(unknown kingdom)",
      phylum = "Apicomplexa",
      class = "Conoidasida",
      order = "Eucoccidiorida",
      family = "Sarcocystidae",
      genus = "",
      species = "",
      subspecies = "",
      rank = "family",
      ref = "Poche, 1913",
      species_id = NA_real_,
      source = "manually added",
      prevalence = 2,
      stringsAsFactors = FALSE
    ),
    data.frame(
      mo = "[ORD]_EUCCCDRD", # order
      fullname = "Eucoccidiorida",
      kingdom = "(unknown kingdom)",
      phylum = "Apicomplexa",
      class = "Conoidasida",
      order = "Eucoccidiorida",
      family = "",
      genus = "",
      species = "",
      subspecies = "",
      rank = "order",
      ref = "Leger et al., 1910",
      species_id = NA_real_,
      source = "manually added",
      prevalence = 2,
      stringsAsFactors = FALSE
    ),
  ) |>
  arrange(fullname)
