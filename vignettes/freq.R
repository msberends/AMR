## ----setup, include = FALSE, results = 'markup'--------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)
library(dplyr)
library(AMR)

## ---- echo = TRUE, results = 'hide'--------------------------------------
# # just using base R
freq(septic_patients$sex)

# # using base R to select the variable and pass it on with a pipe
septic_patients$sex %>% freq()

# # do it all with pipes, using the `select` function of the dplyr package
septic_patients %>%
  select(sex) %>%
  freq()

## ---- echo = TRUE--------------------------------------------------------
freq(septic_patients$sex)

## ---- echo = TRUE, results = 'hide'--------------------------------------
my_patients <- septic_patients %>% 
  left_join_microorganisms()

## ---- echo = TRUE--------------------------------------------------------
colnames(microorganisms)

## ---- echo = TRUE--------------------------------------------------------
dim(septic_patients)
dim(my_patients)

## ---- echo = TRUE--------------------------------------------------------
my_patients %>%
  select(genus, species) %>%
  freq()

## ---- echo = TRUE--------------------------------------------------------
# # get age distribution of unique patients
septic_patients %>% 
  distinct(patient_id, .keep_all = TRUE) %>% 
  select(age) %>% 
  freq(nmax = 5)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  select(hospital_id) %>% 
  freq()

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  select(hospital_id) %>% 
  freq(sort.count = TRUE)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  select(amox) %>% 
  freq()

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  select(date) %>% 
  freq(nmax = 5)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  select(amox) %>% 
  freq(na.rm = FALSE)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  select(hospital_id) %>% 
  freq(markdown = TRUE)

## ---- echo = TRUE--------------------------------------------------------
my_df <- septic_patients %>%
  select(hospital_id) %>% 
  freq(as.data.frame = TRUE)

my_df

class(my_df)

## ---- echo = FALSE-------------------------------------------------------
# this will print "2018" in 2018, and "2018-yyyy" after 2018.
yrs <- c(2018:format(Sys.Date(), "%Y"))
yrs <- c(min(yrs), max(yrs))
yrs <- paste(unique(yrs), collapse = "-")

