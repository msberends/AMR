## ----setup, include = FALSE, results = 'markup'--------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)
library(dplyr)
library(AMR)

## ---- echo = TRUE, results = 'hide'--------------------------------------
# just using base R
freq(septic_patients$sex)

# using base R to select the variable and pass it on with a pipe from the dplyr package
septic_patients$sex %>% freq()

# do it all with pipes, using the `select` function from the dplyr package
septic_patients %>%
  select(sex) %>%
  freq()

# or the preferred way: using a pipe to pass the variable on to the freq function
septic_patients %>% freq(sex) # this also shows 'age' in the title


## ---- echo = TRUE--------------------------------------------------------
freq(septic_patients$sex)

## ---- echo = TRUE, results = 'hide'--------------------------------------
my_patients <- septic_patients %>% left_join_microorganisms()

## ---- echo = TRUE--------------------------------------------------------
colnames(microorganisms)

## ---- echo = TRUE--------------------------------------------------------
dim(septic_patients)
dim(my_patients)

## ---- echo = TRUE--------------------------------------------------------
my_patients %>% freq(genus, species)

## ---- echo = TRUE--------------------------------------------------------
# # get age distribution of unique patients
septic_patients %>% 
  distinct(patient_id, .keep_all = TRUE) %>% 
  freq(age, nmax = 5)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  freq(hospital_id)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  freq(hospital_id, sort.count = TRUE)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  select(amox) %>% 
  freq()

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  select(date) %>% 
  freq(nmax = 5)

## ---- echo = TRUE--------------------------------------------------------
my_df <- septic_patients %>% freq(age)
class(my_df)

## ---- echo = TRUE--------------------------------------------------------
dim(my_df)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  freq(amox, na.rm = FALSE)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  freq(hospital_id, row.names = FALSE)

## ---- echo = TRUE--------------------------------------------------------
septic_patients %>%
  freq(hospital_id, markdown = TRUE)

## ---- echo = FALSE-------------------------------------------------------
# this will print "2018" in 2018, and "2018-yyyy" after 2018.
yrs <- c(2018:format(Sys.Date(), "%Y"))
yrs <- c(min(yrs), max(yrs))
yrs <- paste(unique(yrs), collapse = "-")

