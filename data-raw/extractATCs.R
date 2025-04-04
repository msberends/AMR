library(rvest)
library(dplyr)
library(tidyr)
library(AMR)
# we need J01, J02 and J04 (J03 does not exist)

url <- "https://atcddd.fhi.no/atc_ddd_index/?code={code}&showdescription=no"
complete_vector <- character(0)
for (Jxx in c("J01", "J02", "J04")) {
  site <- gsub("{code}", Jxx, url, fixed = TRUE)
  tbl <- site |>
    read_html() |>
    html_element("#content") |>
    html_text() |>
    strsplit(Jxx, fixed = TRUE)
  out <- paste0(Jxx, tbl[[1]])
  out <- out[out %like% paste0(Jxx, "[A-Z0-9]")]
  out <- gsub("(.*?)\n.*", "\\1", out)
  complete_vector <- c(complete_vector, out)
}

complete_vector.bak <- complete_vector
next_codes <- gsub("(.*?) .*", "\\1", complete_vector)
for (Jxxx in next_codes) {
  message("Checking ", Jxxx)
  site <- gsub("{code}", Jxxx, url, fixed = TRUE)
  tbl <- site |>
    read_html() |>
    html_element("#content") |>
    html_text() |>
    strsplit(Jxxx, fixed = TRUE)
  out <- paste0(Jxxx, tbl[[1]])
  out <- out[out %like% paste0(Jxxx, "[A-Z0-9]")]
  out <- gsub("(.*?)\n.*", "\\1", out)
  complete_vector <- c(complete_vector, out)
}

next_codes <- gsub("(.*?) .*", "\\1", complete_vector[complete_vector %like% "[A-Z][0-9][0-9][A-Z][A-Z]"])
complete_tbl <- NULL
for (Jxxxxx in next_codes) {
  message("Checking ", Jxxxxx)
  site <- gsub("{code}", Jxxxxx, url, fixed = TRUE)
  tbl <- site |>
    read_html() |>
    html_element("table") |>
    html_table(header = TRUE) |>
    mutate_if(is.character, function(x) ifelse(x == "", NA_character_, x)) |>
    fill(`ATC code`, Name, .direction = "down")
  if (is.null(complete_tbl)) {
    complete_tbl <- tbl
  } else {
    complete_tbl <- bind_rows(complete_tbl, tbl)
  }
}

new_ab <- complete_tbl |>
  select(atc = `ATC code`, name = Name, ddd = DDD, ddd_units = U, route = Adm.R) |>
  filter(route %in% c("O", "P", NA)) |>
  mutate(route = case_when(
    route == "O" ~ "oral",
    route == "P" ~ "iv",
    TRUE ~ NA_character_
  )) |>
  pivot_wider(names_from = route, values_from = c(ddd, ddd_units), values_fn = first) |>
  select(atc, name,
    oral_ddd = ddd_oral, oral_units = ddd_units_oral,
    iv_ddd = ddd_iv, iv_units = ddd_units_iv
  ) |>
  mutate(name = paste0(substr(toupper(name), 1, 1), substr(name, 2, 999))) |>
  mutate(name = gsub(" and ", "/", name)) |>
  filter(name %unlike% "^Combinations",
         name %unlike% "/beta[-]lactamase inhibitor",
         name %unlike% "combinations") |>
  arrange(name)

new_atcs <- new_ab |>
  mutate(name_gen = generalise_antibiotic_name(name)) |>
  filter(name_gen %in% AMR_env$AB_lookup$generalised_name) |>
  mutate(ab = as.ab(name_gen))
new_atcs$atc_group1 <- complete_vector
for (i in seq_len(nrow(new_atcs))) {
  # fill in DDDs and ATC code
  antimicrobials$oral_ddd[which(antimicrobials$ab == new_atcs$ab[i])] <- new_atcs$oral_ddd[i]
  antimicrobials$oral_units[which(antimicrobials$ab == new_atcs$ab[i])] <- new_atcs$oral_units[i]
  antimicrobials$iv_ddd[which(antimicrobials$ab == new_atcs$ab[i])] <- new_atcs$iv_ddd[i]
  antimicrobials$iv_units[which(antimicrobials$ab == new_atcs$ab[i])] <- new_atcs$iv_units[i]
  antimicrobials$atc[which(antimicrobials$ab == new_atcs$ab[i])] <- list(c(antimicrobials$atc[which(antimicrobials$ab == new_atcs$ab[i])][[1]], new_atcs$atc[i]))
}


# check these - any new?
new_ab |>
  filter(!name %in% antimicrobials$name | !atc %in% unlist(antimicrobials$atc)) |>
  mutate(
    name_old = ab_name(atc, language = NULL),
    new = !atc %in% unlist(antimicrobials$atc)
  ) |>
  View()


atc_ref <- antimicrobials |>
  select(ab, group, atc_group1, atc_group2) |>
  mutate(atc = ab_atc(ab, only_first = TRUE)) |>
  filter(!is.na(atc)) |>
  mutate(atc = substr(atc, 1, 5)) |>
  select(atc, everything(), -ab) |>
  distinct() |>
  arrange(atc)
atc_ref <- atc_ref |>
  arrange(atc, group, atc_group1, atc_group2) |>
  distinct(atc, .keep_all = TRUE)

antimicrobials |>
  bind_rows(new |>
    select(-group, -atc_group1, -atc_group2) |>
    left_join(atc_ref, by = c("atc_part" = "atc")) %>%
    mutate(
      atc = as.list(atc),
      loinc = as.list(character(nrow(.))),
      abbreviations = loinc,
      synonyms = loinc
    ) |>
    select(-atc_part)) |>
  View()
