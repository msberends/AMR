# ==================================================================== #
# TITLE:                                                               #
# AMR: An R Package for Working with Antimicrobial Resistance Data     #
#                                                                      #
# SOURCE CODE:                                                         #
# https://github.com/msberends/AMR                                     #
#                                                                      #
# PLEASE CITE THIS SOFTWARE AS:                                        #
# Berends MS, Luz CF, Friedrich AW, et al. (2022).                     #
# AMR: An R Package for Working with Antimicrobial Resistance Data.    #
# Journal of Statistical Software, 104(3), 1-31.                       #
# https://doi.org/10.18637/jss.v104.i03                                #
#                                                                      #
# Developed at the University of Groningen and the University Medical  #
# Center Groningen in the Netherlands, in collaboration with many      #
# colleagues from around the world, see our website.                   #
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
# how to conduct AMR data analysis: https://amr-for-r.org              #
# ==================================================================== #

test_that("test-data.R", {
  skip_on_cran()

  # IDs should always be unique
  expect_identical(nrow(microorganisms), length(unique(microorganisms$mo)))
  expect_identical(class(microorganisms$mo), c("mo", "character"))
  expect_identical(nrow(antimicrobials), length(unique(AMR::antimicrobials$ab)))
  expect_identical(class(AMR::antimicrobials$ab), c("ab", "character"))
  expect_identical(
    nrow(antimicrobials[!is.na(antimicrobials$cid), ]),
    length(unique(AMR::antimicrobials$cid[!is.na(antimicrobials$cid)]))
  )

  # check cross table reference
  expect_true(all(microorganisms.codes$mo %in% microorganisms$mo))
  expect_true(all(example_isolates$mo %in% microorganisms$mo))
  expect_true(all(microorganisms.groups$mo %in% microorganisms$mo))
  expect_true(all(microorganisms.groups$mo_group %in% microorganisms$mo))
  expect_true(all(clinical_breakpoints$mo %in% microorganisms$mo))
  expect_true(all(clinical_breakpoints$ab %in% AMR::antimicrobials$ab))
  expect_true(all(intrinsic_resistant$mo %in% microorganisms$mo))
  expect_true(all(intrinsic_resistant$ab %in% AMR::antimicrobials$ab))
  expect_false(anyNA(microorganisms.codes$code))
  expect_false(anyNA(microorganisms.codes$mo))
  expect_true(all(dosage$ab %in% AMR::antimicrobials$ab))
  expect_true(all(dosage$name %in% AMR::antimicrobials$name))
  interpretive_abx <- AMR:::INTERPRETIVE_RULES_DF$and_these_antibiotics
  interpretive_abx <- unique(unlist(strsplit(interpretive_abx[!is.na(interpretive_abx)], ", +")))
  expect_true(all(interpretive_abx %in% AMR::antimicrobials$ab),
    info = paste0(
      "Missing in `antimicrobials` data set: ",
      toString(interpretive_abx[which(!interpretive_abx %in% AMR::antimicrobials$ab)])
    )
  )

  # check valid disks/MICs
  # rows without any breakpoint are EUCAST's blocking rows ("no breakpoint" for that organism), which must state their reason
  is_blocking <- is.na(clinical_breakpoints$breakpoint_S) & is.na(clinical_breakpoints$breakpoint_R)
  expect_true(all(grepl("^\\[No breakpoint\\] ", clinical_breakpoints$note[is_blocking])))
  expect_true(all(clinical_breakpoints$guideline[is_blocking] %like% "EUCAST"))
  expect_false(anyNA(as.mic(clinical_breakpoints[which(clinical_breakpoints$method == "MIC" & clinical_breakpoints$ref_tbl != "ECOFF" & !is_blocking), "breakpoint_S", drop = TRUE])))
  expect_true(anyNA(as.mic(clinical_breakpoints[which(clinical_breakpoints$method == "MIC" & clinical_breakpoints$ref_tbl != "ECOFF" & !is_blocking), "breakpoint_R", drop = TRUE])))
  expect_false(anyNA(as.disk(clinical_breakpoints[which(clinical_breakpoints$method == "DISK" & clinical_breakpoints$ref_tbl != "ECOFF" & !is_blocking), "breakpoint_S", drop = TRUE])))
  expect_true(anyNA(as.disk(clinical_breakpoints[which(clinical_breakpoints$method == "DISK" & clinical_breakpoints$ref_tbl != "ECOFF" & !is_blocking), "breakpoint_R", drop = TRUE])))

  # antibiotic names must always be coercible to their original AB code
  expect_identical(as.ab(AMR::antimicrobials$name), AMR::antimicrobials$ab)

  if (AMR:::pkg_is_available("tibble")) {
    # there should be no diacritics (i.e. non ASCII) characters in the datasets (CRAN policy)
    datasets <- data(package = "AMR", envir = asNamespace("AMR"))$results[, "Item", drop = TRUE]
    datasets <- datasets[datasets != "antibiotics"]
    for (i in seq_len(length(datasets))) {
      dataset <- get(datasets[i], envir = asNamespace("AMR"))
      expect_identical(AMR:::dataset_UTF8_to_ASCII(dataset), dataset, info = datasets[i])
    }
  }

  AMR:::add_MO_lookup_to_AMR_env()
  df <- AMR:::AMR_env$MO_lookup
  expect_true(all(c(
    "mo", "fullname", "status", "kingdom", "phylum", "class", "order",
    "family", "genus", "species", "subspecies", "rank", "ref", "source",
    "lpsn", "lpsn_parent", "lpsn_renamed_to", "gbif", "gbif_parent", "gbif_renamed_to", "prevalence",
    "snomed", "domain_index", "fullname_lower", "full_first", "species_first"
  ) %in% colnames(df)))

  expect_inherits(AMR:::MO_CONS, "mo")

  # (only taxa that are staphylococci: e.g. 'Staphylococcus rhodochrous' is a synonym of Rhodococcus rhodochrous)
  # (a synonym is categorised by its current name, e.g. Staphylococcus roterodami is S. aureus)
  current_mo <- AMR:::synonym_mo_to_accepted_mo(as.character(microorganisms$mo), fill_in_accepted = TRUE)
  current_genus <- microorganisms$genus[match(current_mo, as.character(microorganisms$mo))]
  current_species <- microorganisms$species[match(current_mo, as.character(microorganisms$mo))]
  uncategorised <- subset(
    microorganisms,
    genus == "Staphylococcus" &
      current_genus %in% "Staphylococcus" &
      !species %in% c("", "aureus") &
      !current_species %in% "aureus" &
      !mo %in% c(AMR:::MO_CONS, AMR:::MO_COPS) &
      !current_mo %in% c(AMR:::MO_CONS, AMR:::MO_COPS)
  )
  expect_true(NROW(uncategorised) == 0,
    info = ifelse(NROW(uncategorised) == 0,
      "All staphylococcal species categorised as CoNS/CoPS.",
      paste0(
        "Staphylococcal species not categorised as CoNS/CoPS: S. ",
        uncategorised$species, " (", uncategorised$mo, ")",
        collapse = "\n"
      )
    )
  )

  # THIS WILL CHECK NON-ASCII STRINGS IN ALL FILES:

  # check_non_ascii <- function() {
  #   purrr::map_df(
  #     .id = "file",
  #     # list common text files
  #     .x = fs::dir_ls(
  #       recurse = TRUE,
  #       type = "file",
  #       # ignore images, compressed
  #       regexp = "\\.(png|ico|rda|ai|tar.gz|zip|xlsx|csv|pdf|psd)$",
  #       invert = TRUE
  #     ),
  #     .f = function(path) {
  #       x <- readLines(path, warn = FALSE)
  #       # from tools::showNonASCII()
  #       asc <- iconv(x, "latin1", "ASCII")
  #       ind <- is.na(asc) | asc != x
  #       # make data frame
  #       if (any(ind)) {
  #         tibble::tibble(
  #           row = which(ind),
  #           line = iconv(x[ind], "latin1", "ASCII", sub = "byte")
  #         )
  #       } else {
  #         tibble::tibble()
  #       }
  #     }
  #   )
  # }
  # x <- check_non_ascii() %>%
  #   filter(file %unlike% "^(data-raw|docs|git_)")
})

test_that("taxonomic name columns contain no NA (empty string is used instead)", {
  for (col in c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies")) {
    expect_false(anyNA(microorganisms[[col]]), info = col)
  }
})

# ==================================================================== #
# Strict integrity of the taxonomy and of every data set that refers to it. Each test names the
# records that break it. These tests exist because in October 2026 a taxonomy build passed all
# earlier tests while Eggerthella lenta, Gordonia amarae and Anisakis simplex were synonyms
# without a current name (a source contained the same name twice), and while data sets referred
# to outdated codes.
# ==================================================================== #

# describes failing records for test output, e.g. "B_EGGRT_LENT (Eggerthella lenta)"
describe_mo <- function(mo, n = 25) {
  mo <- as.character(mo)
  nm <- microorganisms$fullname[match(mo, as.character(microorganisms$mo))]
  out <- paste0(mo, " (", nm, ")")
  paste0(paste(utils::head(out, n), collapse = ", "), ifelse(length(out) > n, paste0(", ... (", length(out), " in total)"), ""))
}
mo_df <- function() {
  df <- as.data.frame(microorganisms, stringsAsFactors = FALSE)
  df$mo <- as.character(df$mo)
  df
}
# the row of the current name of each record (its own row if it is not a synonym), with the priority of
# synonym_mo_to_accepted_mo(): LPSN > MycoBank > GBIF
current_row <- function(df) {
  target <- rep(NA_integer_, nrow(df))
  for (id in c("gbif", "mycobank", "lpsn")) {
    t <- match(df[[paste0(id, "_renamed_to")]], df[[id]], incomparables = NA)
    target[!is.na(t)] <- t[!is.na(t)]
  }
  target[df$status != "synonym"] <- NA_integer_
  out <- ifelse(df$status != "synonym", seq_len(nrow(df)), target)
  for (step in seq_len(10)) {
    next_row <- ifelse(!is.na(out) & df$status[out] %in% "synonym", target[out], out)
    out <- next_row
  }
  out
}

test_that("taxonomy: every synonym leads to a current name", {
  skip_on_cran()
  df <- mo_df()
  cur <- current_row(df)
  orphans <- df$mo[df$status == "synonym" & (is.na(cur) | df$status[cur] %in% "synonym")]
  expect_true(length(orphans) == 0, info = paste("Synonyms without a current name:", describe_mo(orphans)))
})

test_that("taxonomy: synonym_mo_to_accepted_mo() never returns a synonym, NA or the input itself", {
  skip_on_cran()
  df <- mo_df()
  syn <- df$mo[df$status == "synonym"]
  res <- AMR:::synonym_mo_to_accepted_mo(syn)
  expect_true(!anyNA(res), info = paste("No current name:", describe_mo(syn[is.na(res)])))
  expect_true(all(res != syn, na.rm = TRUE), info = paste("Renamed to itself:", describe_mo(syn[which(res == syn)])))
  res_status <- df$status[match(res, df$mo)]
  expect_true(!any(res_status %in% "synonym"), info = paste("Current name is a synonym:", describe_mo(syn[res_status %in% "synonym"])))
})

test_that("taxonomy: synonym chains are short and free of cycles", {
  skip_on_cran()
  df <- mo_df()
  target <- rep(NA_integer_, nrow(df))
  for (id in c("gbif", "mycobank", "lpsn")) {
    t <- match(df[[paste0(id, "_renamed_to")]], df[[id]], incomparables = NA)
    target[!is.na(t)] <- t[!is.na(t)]
  }
  target[df$status != "synonym"] <- NA_integer_
  pos <- ifelse(df$status == "synonym", target, NA_integer_)
  for (step in seq_len(5)) {
    pos <- ifelse(!is.na(pos) & df$status[pos] %in% "synonym", target[pos], pos)
  }
  too_long <- df$mo[!is.na(pos) & df$status[pos] %in% "synonym"]
  expect_true(length(too_long) == 0, info = paste("Chains longer than 5 steps or cycles:", describe_mo(too_long)))
})

test_that("taxonomy: records that are not synonyms do not point to a current name", {
  skip_on_cran()
  df <- mo_df()
  points <- rep(FALSE, nrow(df))
  for (id in c("lpsn", "mycobank", "gbif")) {
    points <- points | (!is.na(df[[paste0(id, "_renamed_to")]]) & df[[paste0(id, "_renamed_to")]] %in% df[[id]])
  }
  wrong <- df$mo[df$status != "synonym" & points]
  expect_true(length(wrong) == 0, info = paste("Accepted or unknown records with a 'renamed to' target:", describe_mo(wrong)))
})

test_that("taxonomy: no record points to itself", {
  skip_on_cran()
  df <- mo_df()
  for (id in c("lpsn", "mycobank", "gbif")) {
    self <- df$mo[which(df[[paste0(id, "_renamed_to")]] == df[[id]])]
    expect_true(length(self) == 0, info = paste(id, "renamed to itself:", describe_mo(self)))
  }
})

test_that("taxonomy: the identifier of a source belongs to one record only", {
  skip_on_cran()
  df <- mo_df()
  for (id in c("lpsn", "mycobank", "gbif")) {
    x <- df[[id]]
    dup <- df$mo[!is.na(x) & (duplicated(x) | duplicated(x, fromLast = TRUE))]
    expect_true(length(dup) == 0, info = paste(id, "identifier used by more than one record:", describe_mo(dup)))
  }
})

test_that("taxonomy: parent identifiers refer to existing records of a higher rank", {
  skip_on_cran()
  df <- mo_df()
  rank_order <- c("domain" = 0, "kingdom" = 1, "phylum" = 2, "class" = 3, "order" = 4, "family" = 5, "genus" = 6, "species" = 7, "subspecies" = 8)
  for (id in c("lpsn", "mycobank", "gbif")) {
    parent <- df[[paste0(id, "_parent")]]
    has <- !is.na(parent) & df$status == "accepted"
    p <- match(parent, df[[id]], incomparables = NA)
    missing_parent <- df$mo[has & is.na(p)]
    expect_true(length(missing_parent) == 0, info = paste(id, "parent not in the data set:", describe_mo(missing_parent)))
    lower <- df$mo[has & !is.na(p) & df$rank %in% names(rank_order) & df$rank[p] %in% names(rank_order) &
      rank_order[df$rank[p]] >= rank_order[df$rank]]
    expect_true(length(lower) == 0, info = paste(id, "parent not of a higher rank:", describe_mo(lower)))
  }
})

test_that("taxonomy: the taxonomic fields match the rank", {
  skip_on_cran()
  df <- mo_df()
  is_taxon <- df$fullname %unlike% "unknown"
  wrong <- df$mo[is_taxon & (
    (df$rank == "genus" & (df$genus == "" | df$species != "" | df$subspecies != "")) |
      (df$rank == "species" & (df$genus == "" | df$species == "" | df$subspecies != "")) |
      (df$rank == "subspecies" & (df$genus == "" | df$species == "" | df$subspecies == "")) |
      (df$rank == "family" & (df$family == "" | df$genus != "")) |
      (df$rank == "order" & (df$order == "" | df$family != "")) |
      (df$rank == "class" & (df$class == "" | df$order != "")) |
      (df$rank == "phylum" & (df$phylum == "" | df$class != ""))
  )]
  expect_true(length(wrong) == 0, info = paste("Fields do not match the rank:", describe_mo(wrong)))
})

test_that("taxonomy: the full name follows from the taxonomic fields", {
  skip_on_cran()
  df <- mo_df()
  is_taxon <- df$fullname %unlike% "unknown"
  low <- df$rank %in% c("genus", "species", "subspecies") & is_taxon
  # Salmonella serovars are named without their species, e.g. 'Salmonella Typhi' (genus Salmonella, species enterica)
  serovar <- df$genus == "Salmonella" & df$subspecies %like_case% "^[A-Z]"
  # (curated group labels, such as 'Milleri Group Streptococcus (MGS)', are named on purpose)
  low <- low & df$fullname %unlike_case% " Group "
  expected <- trimws(gsub(" +", " ", ifelse(serovar, paste(df$genus, df$subspecies), paste(df$genus, df$species, df$subspecies))))
  wrong <- df$mo[low & mo_name_without_suffix(df$fullname) != expected]
  expect_true(length(wrong) == 0, info = paste("Full name does not match genus/species/subspecies:", describe_mo(wrong)))
  for (r in c("family", "order", "class", "phylum")) {
    wrong <- df$mo[df$rank == r & is_taxon & mo_name_without_suffix(df$fullname) != df[[r]]]
    expect_true(length(wrong) == 0, info = paste("Full name of a", r, "does not match the", r, "field:", describe_mo(wrong)))
  }
})

test_that("taxonomy: every species and subspecies has a genus with the same higher taxonomy", {
  skip_on_cran()
  df <- mo_df()
  genera <- df[df$rank == "genus" & df$fullname %unlike% "unknown", , drop = FALSE]
  low <- which(df$rank %in% c("species", "subspecies") & df$fullname %unlike% "unknown")
  g <- match(paste(df$domain[low], df$genus[low]), paste(genera$domain, genera$genus))
  without <- df$mo[low[is.na(g)]]
  expect_true(length(without) == 0, info = paste("No genus record:", describe_mo(without)))
  for (r in c("kingdom", "phylum", "class", "order", "family")) {
    differs <- df$mo[low[!is.na(g) & df[[r]][low] != genera[[r]][g]]]
    expect_true(length(differs) == 0, info = paste("Other", r, "than its genus:", describe_mo(differs)))
  }
})

test_that("taxonomy: every subspecies has a species", {
  skip_on_cran()
  df <- mo_df()
  sub <- which(df$rank == "subspecies" & df$fullname %unlike% "unknown")
  sp <- paste(df$domain, df$genus, df$species)[df$rank == "species"]
  without <- df$mo[sub[!paste(df$domain[sub], df$genus[sub], df$species[sub]) %in% sp]]
  expect_true(length(without) == 0, info = paste("No species record:", describe_mo(without)))
})

test_that("taxonomy: every genus has its family, and every family its order, with the same higher taxonomy", {
  skip_on_cran()
  df <- mo_df()
  # PENDING A DECISION (October 2026): protist groups that the sources place in more than one domain, so that their
  # higher taxa exist in another domain than (some of) their members; a policy for their domain is needed first
  # (the labyrinthulids, the radiolarian Acantharia and the plasmodial slime moulds)
  pending_classes <- c("Labyrinthulea", "Acantharia", "Myxomycetes")
  pending <- c(
    pending_classes, "Thraustochytrida", "Thraustochytriaceae", "Thraustochytriidae", "Amphifilidae",
    "Diplophryidae", "Oblongichytriidae", "Sorodiplophryidae", "Dictydiaethaliaceae", "Arthracanthida"
  )
  df <- df[!df$fullname %in% pending & !df$class %in% pending_classes & !df$order %in% pending & !df$family %in% pending, , drop = FALSE]
  pairs <- list(c("genus", "family"), c("family", "order"), c("order", "class"), c("class", "phylum"))
  for (pr in pairs) {
    child <- which(df$rank == pr[1] & df[[pr[2]]] != "" & df$fullname %unlike% "unknown")
    parents <- df[df$rank == pr[2], , drop = FALSE]
    p <- match(paste(df$domain[child], df[[pr[2]]][child]), paste(parents$domain, parents[[pr[2]]]))
    without <- df$mo[child[is.na(p)]]
    expect_true(length(without) == 0, info = paste("A", pr[1], "without a", pr[2], "record:", describe_mo(without)))
    higher <- c("kingdom", "phylum", "class", "order")
    higher <- higher[seq_len(which(c("kingdom", "phylum", "class", "order", "family") == pr[2]) - 1)]
    for (r in higher) {
      differs <- df$mo[child[!is.na(p) & df[[r]][child] != parents[[r]][p]]]
      expect_true(length(differs) == 0, info = paste("A", pr[1], "with another", r, "than its", pr[2], ":", describe_mo(differs)))
    }
  }
})

test_that("taxonomy: references are clean", {
  skip_on_cran()
  df <- mo_df()
  ref <- df$ref[!is.na(df$ref) & df$ref != ""]
  mo <- df$mo[!is.na(df$ref) & df$ref != ""]
  unbalanced <- mo[lengths(regmatches(ref, gregexpr("(", ref, fixed = TRUE))) != lengths(regmatches(ref, gregexpr(")", ref, fixed = TRUE)))]
  expect_true(length(unbalanced) == 0, info = paste("Unbalanced parentheses in ref:", describe_mo(unbalanced)))
  quoted <- mo[ref %like_case% "^[\"' ]|[\"' ]$"]
  expect_true(length(quoted) == 0, info = paste("Quotes or spaces around ref:", describe_mo(quoted)))
})

test_that("taxonomy: LPSN records are prokaryotes with an LPSN identifier", {
  skip_on_cran()
  df <- mo_df()
  wrong <- df$mo[df$source == "LPSN" & (is.na(df$lpsn) | !df$domain %in% c("Bacteria", "Archaea"))]
  expect_true(length(wrong) == 0, info = describe_mo(wrong))
  wrong <- df$mo[df$source == "MycoBank" & is.na(df$mycobank)]
  expect_true(length(wrong) == 0, info = describe_mo(wrong))
  wrong <- df$mo[df$source == "GBIF" & is.na(df$gbif)]
  expect_true(length(wrong) == 0, info = describe_mo(wrong))
})

test_that("taxonomy: clinically important names are current names", {
  skip_on_cran()
  # each must exist with status 'accepted'; this list includes names that were wrongly a synonym in a build
  must_be_current <- c(
    # bacteria
    "Escherichia coli", "Klebsiella pneumoniae", "Klebsiella oxytoca", "Klebsiella aerogenes", "Enterobacter cloacae",
    "Citrobacter freundii", "Serratia marcescens", "Proteus mirabilis", "Morganella morganii", "Salmonella enterica",
    "Shigella sonnei", "Yersinia enterocolitica", "Pseudomonas aeruginosa", "Acinetobacter baumannii",
    "Stenotrophomonas maltophilia", "Burkholderia cepacia", "Burkholderia pyrrocinia", "Haemophilus influenzae",
    "Moraxella catarrhalis", "Neisseria gonorrhoeae", "Neisseria meningitidis", "Legionella pneumophila",
    "Campylobacter jejuni", "Helicobacter pylori", "Bordetella pertussis", "Brucella melitensis", "Brucella anthropi",
    "Budvicia aquatica", "Bacteroides fragilis", "Staphylococcus aureus", "Staphylococcus epidermidis",
    "Staphylococcus lugdunensis", "Streptococcus pneumoniae", "Streptococcus pyogenes", "Streptococcus agalactiae",
    "Enterococcus faecalis", "Enterococcus faecium", "Listeria monocytogenes", "Clostridioides difficile",
    "Clostridium perfringens", "Thomasclavelia ramosa", "Eggerthella lenta", "Cutibacterium acnes",
    "Corynebacterium diphtheriae", "Nocardia farcinica", "Gordonia amarae", "Lacticaseibacillus rhamnosus",
    "Mycobacterium tuberculosis", "Mycobacterium canettii", "Mycobacterium orygis", "Mycobacterium avium",
    "Mycobacterium abscessus", "Mycobacterium leprae",
    "Tropheryma whipplei", "Treponema pallidum", "Chlamydia trachomatis", "Mycoplasmoides pneumoniae",
    # fungi
    "Candida albicans", "Candida parapsilosis", "Candida tropicalis", "Nakaseomyces glabratus", "Candidozyma auris",
    "Pichia kudriavzevii", "Cryptococcus neoformans", "Aspergillus fumigatus", "Aspergillus flavus", "Aspergillus niger",
    "Pneumocystis jirovecii", "Trichophyton rubrum", "Histoplasma capsulatum",
    # parasites
    "Plasmodium falciparum", "Plasmodium vivax", "Toxoplasma gondii", "Giardia duodenalis", "Anisakis simplex",
    "Enterobius vermicularis", "Ascaris lumbricoides", "Echinococcus granulosus",
    # genera and higher taxa (in a build of October 2026, the nematode order Strongylida became a synonym of a
    # fungal order with the same old GBIF identifier)
    "Staphylococcus", "Streptococcus", "Enterococcus", "Klebsiella", "Pseudomonas", "Acinetobacter", "Mycobacterium",
    "Candida", "Aspergillus", "Plasmodium", "Ancylostoma", "Necator", "Strongyloides", "Schistosoma", "Taenia",
    "Enterobacterales", "Enterobacteriaceae", "Staphylococcaceae", "Mycobacteriaceae", "Strongylida",
    "Eurotiales", "Saccharomycetes", "Apicomplexa", "Nematoda", "Platyhelminthes"
  )
  # higher taxa are often of 'unknown' status in the sources, they must not be a synonym
  higher <- must_be_current[!grepl(" ", must_be_current)]
  status <- microorganisms$status[match(must_be_current, microorganisms$fullname)]
  missing <- must_be_current[is.na(status)]
  not_current <- must_be_current[!is.na(status) & (status == "synonym" | (status != "accepted" & !must_be_current %in% higher))]
  expect_true(length(missing) == 0, info = paste("Missing:", paste(missing, collapse = ", ")))
  expect_true(length(not_current) == 0, info = paste("Not 'accepted':", paste(paste0(not_current, " (", status[match(not_current, must_be_current)], ")"), collapse = ", ")))
  domains <- c(
    "Strongylida" = "Animalia", "Nematoda" = "Animalia", "Platyhelminthes" = "Animalia",
    "Ancylostoma" = "Animalia", "Necator" = "Animalia", "Eurotiales" = "Fungi", "Saccharomycetes" = "Fungi",
    "Apicomplexa" = "Protozoa", "Plasmodium" = "Protozoa", "Enterobacterales" = "Bacteria", "Mycobacterium" = "Bacteria"
  )
  actual <- microorganisms$domain[match(names(domains), microorganisms$fullname)]
  wrong <- names(domains)[!actual %in% domains]
  expect_true(length(wrong) == 0, info = paste("In the wrong domain:", paste(paste0(wrong, " (", actual[match(wrong, names(domains))], ")"), collapse = ", ")))
  # these must exist, but LPSN has them as synonyms of M. tuberculosis (Riojas et al. 2018), so they can be either
  must_exist <- c("Mycobacterium bovis", "Mycobacterium africanum", "Mycobacterium caprae", "Mycobacterium microti", "Mycobacterium pinnipedii")
  missing <- must_exist[!must_exist %in% microorganisms$fullname]
  expect_true(length(missing) == 0, info = paste("Missing:", paste(missing, collapse = ", ")))
})

test_that("taxonomy: outdated names lead to their current names", {
  skip_on_cran()
  renamed <- c(
    "Clostridium difficile" = "Clostridioides difficile",
    "Enterobacter aerogenes" = "Klebsiella aerogenes",
    "Propionibacterium acnes" = "Cutibacterium acnes",
    "Eubacterium lentum" = "Eggerthella lenta",
    "Clostridium ramosum" = "Thomasclavelia ramosa",
    "Erysipelatoclostridium ramosum" = "Thomasclavelia ramosa",
    "Ochrobactrum anthropi" = "Brucella anthropi",
    "Lactobacillus rhamnosus" = "Lacticaseibacillus rhamnosus",
    "Candida glabrata" = "Nakaseomyces glabratus",
    "Candida krusei" = "Pichia kudriavzevii",
    "Mycobacteroides abscessus" = "Mycobacterium abscessus",
    "Mycoplasma pneumoniae" = "Mycoplasmoides pneumoniae",
    "Nocardia amarae" = "Gordonia amarae",
    "Pseudomonas pyrrocinia" = "Burkholderia pyrrocinia",
    "Giardia intestinalis" = "Giardia duodenalis",
    "Bacillus coli" = "Escherichia coli"
  )
  current <- suppressWarnings(suppressMessages(mo_name(names(renamed), language = NULL, keep_synonyms = FALSE, info = FALSE)))
  wrong <- names(renamed)[current != unname(renamed) | is.na(current)]
  expect_true(length(wrong) == 0,
    info = paste(sprintf("%s -> %s (expected %s)", wrong, current[match(wrong, names(renamed))], renamed[wrong]), collapse = "; ")
  )
})

test_that("taxonomy: a name that a source has twice was taken from its correct record", {
  skip_on_cran()
  # regression: LPSN has e.g. the correct name Eggerthella lenta Wade et al. 1999 (775884) and the illegitimate
  # homotypic synonym Eggerthella lenta Kageyama et al. 1999 (7150); the latter was taken until October 2026
  expected <- c(
    "Eggerthella lenta" = "775884", "Gordonia amarae" = "776579", "Budvicia aquatica" = "774260",
    "Burkholderia pyrrocinia" = "774298"
  )
  actual <- microorganisms$lpsn[match(names(expected), microorganisms$fullname)]
  expect_identical(unname(actual), unname(expected))
  expect_identical(microorganisms$gbif[match("Anisakis simplex", microorganisms$fullname)], "KZNMP")
})

test_that("data sets refer only to current names", {
  skip_on_cran()
  status_of <- function(mo) microorganisms$status[match(as.character(mo), as.character(microorganisms$mo))]
  check <- list(
    "example_isolates$mo" = example_isolates$mo,
    "clinical_breakpoints$mo" = clinical_breakpoints$mo,
    "microorganisms.codes$mo" = microorganisms.codes$mo,
    "microorganisms.groups$mo" = microorganisms.groups$mo,
    "microorganisms.groups$mo_group" = microorganisms.groups$mo_group
  )
  for (nm in names(check)) {
    mo <- unique(as.character(check[[nm]]))
    mo <- mo[!is.na(mo) & mo != "UNKNOWN"]
    st <- status_of(mo)
    expect_true(!anyNA(st), info = paste(nm, "contains codes that do not exist:", paste(utils::head(mo[is.na(st)], 25), collapse = ", ")))
    expect_true(!any(st %in% "synonym"), info = paste(nm, "contains synonyms:", describe_mo(mo[st %in% "synonym"])))
  }
})

test_that("internal organism lists contain existing codes and are complete for synonyms", {
  skip_on_cran()
  # the genus lists may contain outdated genera on purpose (e.g. Wangiella), so that these are recognised as well
  genera <- microorganisms$fullname[microorganisms$rank == "genus"]
  for (nm in c("MO_RELEVANT_GENERA", "MO_WHO_PRIORITY_GENERA")) {
    x <- get(nm, envir = asNamespace("AMR"))
    wrong <- x[!x %in% genera]
    expect_true(length(wrong) == 0, info = paste(nm, "contains names that are no genus in the data:", paste(wrong, collapse = ", ")))
  }
  # the code lists contain synonyms on purpose (so that outdated codes are classified as well), but then also
  # their current names
  for (nm in c("MO_CONS", "MO_COPS", "MO_STREP_ABCG", "MO_LANCEFIELD")) {
    x <- as.character(get(nm, envir = asNamespace("AMR")))
    absent <- x[!x %in% as.character(microorganisms$mo)]
    expect_true(length(absent) == 0, info = paste(nm, "contains codes that do not exist:", paste(absent, collapse = ", ")))
    cur <- AMR:::synonym_mo_to_accepted_mo(x)
    incomplete <- unique(cur[!is.na(cur) & !cur %in% x])
    expect_true(length(incomplete) == 0, info = paste(nm, "lacks the current names of its synonyms:", describe_mo(incomplete)))
  }
})

test_that("CoNS and CoPS lists exist, do not overlap and are staphylococci", {
  skip_on_cran()
  expect_true(all(AMR:::MO_CONS %in% microorganisms$mo), info = describe_mo(AMR:::MO_CONS[!AMR:::MO_CONS %in% microorganisms$mo]))
  expect_true(all(AMR:::MO_COPS %in% microorganisms$mo), info = describe_mo(AMR:::MO_COPS[!AMR:::MO_COPS %in% microorganisms$mo]))
  expect_length(intersect(AMR:::MO_CONS, AMR:::MO_COPS), 0)
  cur <- AMR:::synonym_mo_to_accepted_mo(as.character(c(AMR:::MO_CONS, AMR:::MO_COPS)), fill_in_accepted = TRUE)
  genus <- microorganisms$genus[match(cur, microorganisms$mo)]
  # (Staphylococcus caseolyticus, now Macrococcus caseolyticus, is in the CoNS species list of
  # create_species_cons_cops() in data-raw/_pre_commit_checks.R)
  wrong <- c(AMR:::MO_CONS, AMR:::MO_COPS)[!genus %in% c("Staphylococcus", "Mammaliicoccus") &
    !microorganisms$fullname[match(cur, microorganisms$mo)] %in% "Macrococcus caseolyticus"]
  expect_true(length(wrong) == 0, info = paste("Not a staphylococcus:", describe_mo(wrong)))
})

test_that("intrinsic_resistant: synonyms have exactly the resistances of their current name", {
  skip_on_cran()
  ir <- data.frame(mo = as.character(intrinsic_resistant$mo), ab = as.character(intrinsic_resistant$ab), stringsAsFactors = FALSE)
  expect_false(any(duplicated(ir)), info = "intrinsic_resistant contains duplicated rows")
  syn <- unique(ir$mo[microorganisms$status[match(ir$mo, microorganisms$mo)] %in% "synonym"])
  cur <- AMR:::synonym_mo_to_accepted_mo(syn)
  key <- function(m) vapply(split(ir$ab, ir$mo)[m], function(x) paste(sort(x), collapse = ","), character(1))
  differ <- syn[!is.na(cur) & key(syn) != vapply(cur, function(m) if (m %in% ir$mo) key(m) else "", character(1))]
  expect_true(length(differ) == 0, info = paste("Synonyms with other resistances than their current name:", describe_mo(differ)))
  # every current name of a synonym in the table is in the table as well
  absent <- unique(cur[!is.na(cur) & !cur %in% ir$mo])
  expect_true(length(absent) == 0, info = paste("Current names missing while their synonym is listed:", describe_mo(absent)))
})

test_that("interpretive rules name organisms that exist", {
  skip_on_cran()
  r <- AMR:::INTERPRETIVE_RULES_DF
  r <- r[r$like.is.one_of %in% c("is", "one_of") & r$if_mo_property %in% c("genus", "genus_species", "fullname"), , drop = FALSE]
  nms <- unique(trimws(unlist(strsplit(r$this_value, ",", fixed = TRUE))))
  known <- tolower(nms) %in% tolower(c(microorganisms$fullname, microorganisms.groups$mo_group_name))
  # 'Staphylococcus coagulase-negative' is matched as a group name in interpretive_rules()
  unknown <- setdiff(nms[!known], "Staphylococcus coagulase-negative")
  expect_true(length(unknown) == 0, info = paste("Rules name organisms that are not in the data:", paste(unknown, collapse = ", ")))
})

test_that("microorganisms.groups: consistent groups of current names", {
  skip_on_cran()
  g <- microorganisms.groups
  expect_false(any(duplicated(paste(g$mo_group, g$mo))), info = "duplicated group members")
  per_group <- tapply(g$mo_group_name, as.character(g$mo_group), function(x) length(unique(x)))
  expect_true(all(per_group == 1), info = paste("Groups with more than one name:", paste(names(per_group)[per_group != 1], collapse = ", ")))
  # (WHONET also defines a group by genus, e.g. Chilomastix)
  group_rank <- microorganisms$rank[match(g$mo_group, microorganisms$mo)]
  expect_true(all(group_rank %in% c("species group", "genus")), info = paste(unique(g$mo_group_name[!group_rank %in% c("species group", "genus")]), collapse = ", "))
  expect_identical(g$mo_name, microorganisms$fullname[match(g$mo, microorganisms$mo)])
  expect_identical(g$mo_group_name, microorganisms$fullname[match(g$mo_group, microorganisms$mo)])
})

test_that("microorganisms.codes: unique codes for existing current names", {
  skip_on_cran()
  dup <- microorganisms.codes$code[duplicated(toupper(microorganisms.codes$code))]
  expect_true(length(dup) == 0, info = paste("Duplicated codes:", paste(utils::head(dup, 25), collapse = ", ")))
  expect_false(anyNA(microorganisms.codes$mo))
})

test_that("clinical_breakpoints: consistent and unambiguous", {
  skip_on_cran()
  cb <- clinical_breakpoints
  expect_true(all(cb$guideline %like_case% "^(EUCAST|CLSI) [0-9]{4}$"), info = paste(unique(cb$guideline[cb$guideline %unlike_case% "^(EUCAST|CLSI) [0-9]{4}$"]), collapse = ", "))
  expect_true(all(cb$type %in% c("human", "animal", "ECOFF")))
  expect_true(all(cb$method %in% c("MIC", "DISK")))
  expect_false(any(duplicated(cb)), info = "fully duplicated rows")
  key <- paste(cb$guideline, cb$type, cb$host, cb$method, cb$site, cb$mo, cb$ab, cb$uti)
  dup <- unique(key[duplicated(key)])
  expect_true(length(dup) == 0, info = paste("More than one breakpoint for the same guideline, type, host, method, site, organism, agent and UTI:", paste(utils::head(dup, 15), collapse = " | ")))
  has_both <- !is.na(cb$breakpoint_S) & !is.na(cb$breakpoint_R)
  mic_wrong <- which(has_both & cb$method == "MIC" & cb$breakpoint_S > cb$breakpoint_R)
  disk_wrong <- which(has_both & cb$method == "DISK" & cb$breakpoint_S < cb$breakpoint_R)
  expect_length(mic_wrong, 0)
  expect_length(disk_wrong, 0)
  no_dose <- which(cb$method == "DISK" & !is.na(cb$breakpoint_S) & cb$type != "ECOFF" & cb$disk_dose %in% c("", NA))
  expect_true(length(no_dose) == 0, info = paste("Disk breakpoints without a disk content:", paste(utils::head(key[no_dose], 10), collapse = " | ")))
  # rank_index must follow the rank of the organism (1 = subspecies, 2 = species, 2.5 = species group, 3 = genus, ...)
  rank <- microorganisms$rank[match(cb$mo, microorganisms$mo)]
  expected <- c("subspecies" = 1, "species" = 2, "species group" = 2.5, "genus" = 3, "family" = 4, "order" = 5)[rank]
  wrong <- which(!is.na(expected) & floor(cb$rank_index * 2) / 2 != expected & cb$rank_index != expected)
  expect_true(length(wrong) == 0, info = paste("rank_index does not match the rank:", describe_mo(unique(cb$mo[wrong]))))
})

test_that("antimicrobials: unique and well-formed", {
  skip_on_cran()
  ab <- AMR::antimicrobials
  expect_false(any(duplicated(ab$ab)))
  expect_false(any(duplicated(tolower(ab$name))), info = paste(ab$name[duplicated(tolower(ab$name))], collapse = ", "))
  expect_false(anyNA(ab$name))
  expect_false(anyNA(ab$group))
  atc <- unlist(ab$atc)
  atc <- atc[!is.na(atc) & atc != ""]
  # (ATCvet codes start with Q, e.g. QJ01CA04)
  expect_true(all(atc %like_case% "^Q?[A-Z][0-9]{2}[A-Z]{2}[0-9]{2}$"), info = paste(atc[atc %unlike_case% "^Q?[A-Z][0-9]{2}[A-Z]{2}[0-9]{2}$"], collapse = ", "))
  loinc <- unlist(ab$loinc)
  loinc <- loinc[!is.na(loinc) & loinc != ""]
  expect_true(all(loinc %like_case% "^[0-9]+-[0-9]$"), info = paste(utils::head(loinc[loinc %unlike_case% "^[0-9]+-[0-9]$"], 20), collapse = ", "))
})
