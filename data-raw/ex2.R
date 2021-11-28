ex2 <- example_isolates
for (extra_id in seq_len(50)) {
  ex2 <- ex2 %>% 
    bind_rows(example_isolates %>% mutate(patient_id = paste0(patient_id, extra_id)))
}
# randomly clear antibibiograms of 2%
clr <- sort(sample(x = seq_len(nrow(ex2)),
                   size = nrow(ex2) * 0.02))
for (row in which(is.rsi(ex2))) {
  ex2[clr, row] <- NA_rsi_
}


