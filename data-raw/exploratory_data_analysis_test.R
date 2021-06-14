library(dplyr)
example_isolates %>% 
  select(mo, where(is.rsi)) %>% 
  tidyr::pivot_longer(cols = where(is.rsi)) %>% 
  # remove intrisic R
  filter(!paste(mo, name) %in% AMR:::INTRINSIC_R) %>% 
  mutate(name = as.ab(name),
         value = ifelse(value == "R", 1, 0),
         class = ab_group(name)) %>% 
  group_by(mo, class) %>% 
  summarise(n = n(),
            res = mean(value, na.rm = TRUE)) %>% 
  filter(n > 30, !is.na(res))



df <- example_isolates
search_mo <- "B_ESCHR_COLI"
intrinsic_res <- INTRINSIC_R[INTRINSIC_R %like% search_mo]
intrinsic_res <- gsub(".* (.*)", "\\1", intrinsic_res)

x <- df %>%
  select(mo, where(is.rsi)) %>% 
  filter(mo == search_mo) %>%
  # at least 30 results available
  select(function(x) sum(!is.na(x)) >= 30) %>% 
  # remove intrisic R
  select(!matches(paste(intrinsic_res, collapse = "|")))
