df <- example_isolates |>
  filter_first_isolate(method = "e", episode_days = 14) |>
  mutate(mo = ifelse(mo_genus(mo) == "Klebsiella", as.mo("Klebsiella"), mo)) |>
  top_n_microorganisms(10)

out_new <- df |> antibiogram(c("TZP","TZP+GEN","TZP+TOB"), wisca = TRUE, syndromic_group = "ward")
out_nonwisca <- df |> antibiogram(c("TZP","TZP+GEN","TZP+TOB"),
                                  syndromic_group = "ward",
                                  mo_transform = function(x) "",
                                  digits = 1,
                                  minimum = 10,
                                  formatting_type = 14) |>
  as_tibble() |>
  select(-Pathogen)

# parameters_amr.R#L110: no filter on ward, so pts are only in 1 ward, depending on order of data
# parameters_amr.R: number of first isolates are determined on the whole data set, while Klebsiella is aggregated afterwards (=duplicates on genus level)

source("~/Downloads/estimate_definition_amr.R")




