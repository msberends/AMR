# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2018-2020 Berends MS, Luz CF et al.                              #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# We created this package for both routine data analysis and academic  #
# research and it was publicly released in the hope that it will be    #
# useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
# Visit our website for more info: https://msberends.gitlab.io/AMR.    #
# ==================================================================== #

# Read and format data ----------------------------------------------------

library(tidyverse)
library(maps)
library(httr)

GET_df <- function(ip) {
  ip <- paste0("https://ipinfo.io/", ip, "?token=", ipinfo_token)
  result <- ip %>% GET()
  stop_for_status(result)
  result %>%
    content(type = "text", encoding = "UTF-8") %>%
    jsonlite::fromJSON(flatten = TRUE) %>% 
    as_tibble()
}

# get website analytics
source("data-raw/country_analysis_url_token.R")
url_json <- paste0(country_analysis_url, 
                   "/index.php?&module=API&token_auth=", 
                   country_analysis_token,
                   "&method=Live.getLastVisitsDetails&idSite=3&language=en&expanded=1&date=2018-01-01,2028-01-01&period=range&filter_limit=-1&format=JSON&segment=&translateColumnNames=1")
data_json <- jsonlite::read_json(url_json)
data <- tibble(
  timestamp_server = as.POSIXct(sapply(data_json, function(x) x$serverTimestamp), origin = "1970-01-01"),
  ipaddress = sapply(data_json, function(x) x$visitIp))
rm(data_json)

# add country data based on IP address and ipinfo.io API
unique_ip <- unique(data$ipaddress)
ip_tbl <- GET_df(unique_ip[1])
p <- AMR:::progress_estimated(n = length(unique_ip) - 1, min_time = 0)
for (i in 2:length(unique_ip)) {
  p$tick()$print()
  ip_tbl <- ip_tbl %>% 
    bind_rows(GET_df(unique_ip[i]))
}

ip_tbl.bak <- ip_tbl

# add long and lat
ip_tbl <- ip_tbl %>% 
  separate(loc, into = c("y", "x"), sep = ",", remove = FALSE, convert = TRUE)

# Plot world map ----------------------------------------------------------

countries_geometry <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE)) %>% 
  mutate(countries_code = countrycode::countrycode(ID, 
                                                   origin = 'country.name', 
                                                   destination = 'iso2c', 
                                                   custom_match = c("Ascension Island" = "GB", # Great Britain
                                                                    "Azores" = "PT", # Portugal
                                                                    "Barbuda" = "GB", # Great Britain
                                                                    "Bonaire" = "BQ", # Bonaire, Saint Eustatius and Saba
                                                                    "Canary Islands" = "ES", # Spain
                                                                    "Chagos Archipelago" = "MU", # Mauritius
                                                                    "Grenadines" = "VC", # Saint Vincent and the Grenadines
                                                                    "Heard Island" = "AU", # Australia
                                                                    "Kosovo" = "XK",
                                                                    "Madeira Islands" = "PT", # Portugal
                                                                    "Micronesia" = "FM",
                                                                    "Saba" = "BQ", # Bonaire, Saint Eustatius and Saba
                                                                    "Saint Martin" = "MF",
                                                                    "Siachen Glacier" = "IN", # India
                                                                    "Sint Eustatius" = "BQ" # Bonaire, Saint Eustatius and Saba
                                                   )),
         included = as.integer(countries_code %in% ip_tbl$country),
         not_antarctica = as.integer(ID != "Antarctica"),
         countries_name = ifelse(included == 1, as.character(ID), NA))

# how many?
countries_geometry %>% filter(included == 1) %>% nrow()

countries_plot <- ggplot(countries_geometry) +
    geom_sf(aes(fill = included, colour = not_antarctica),
            size = 0.25, 
            show.legend = FALSE) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank()) +
    scale_fill_gradient(low = "white", high = "#CAD6EA", ) +
    # this makes the border Antarctica turn white (invisible):
    scale_colour_gradient(low = "white", high = "#81899B")

countries_plot_mini <- countries_plot
countries_plot_mini$data <- countries_plot_mini$data %>% filter(ID != "Antarctica")
countries_plot_mini <- countries_plot_mini + scale_colour_gradient(low = "#81899B", high = "#81899B")
countries_plot_big <- countries_plot +
  labs(title = tools::toTitleCase("Countries where the AMR package for R was downloaded from"),
       subtitle = paste0("Between March 2018 (first release) and ", format(Sys.Date(), "%B %Y"), "." #,
                         #"The dots denote visitors on our website https://gitlab.io/msberends/AMR."
                         )) +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)) +
  geom_text(aes(x = -170,
                y = -75,
                label = stringr::str_wrap(paste0("Countries (n = ", 
                                                 length(countries_name[!is.na(countries_name)]), "): ", 
                                                 paste(sort(countries_name[!is.na(countries_name)]), collapse = ", ")),
                                          200)),
            hjust = 0,
             size = 4) # +
  # # points of visitors
  # geom_point(data = ip_tbl,
  #            aes(x = x, y = y), 
  #            size = 1,
  #            colour = "#81899B")

# main website page
ggsave("pkgdown/logos/countries.png",
       width = 6, 
       height = 2.5, 
       units = "in", 
       dpi = 100, 
       plot = countries_plot_mini, 
       scale = 1)
# when clicked - a high res enlargement
ggsave("pkgdown/logos/countries_large.png",
       width = 11,
       height = 6, 
       units = "in", 
       dpi = 300, 
       plot = countries_plot_big,
       scale = 1.5)


# Gibberish ---------------------------------------------------------------

data %>%
  left_join(ip_tbl, by = c("ipaddress" = "ip")) %>% 
  group_by(country = countrycode::countrycode(country, 
                                              origin = 'iso2c', 
                                              destination = 'country.name')) %>%
  summarise(first = min(timestamp_server)) %>%
  arrange(desc(first)) %>% 
  mutate(frame = case_when(first <= as.POSIXct("2019-06-30") ~ "Q1-Q2 2019",
                           first <= as.POSIXct("2019-12-31") ~ "Q3-Q4 2019", 
                           TRUE ~ "Q1-Q2 2020")) %>% 
  View()
# 
# p1 <- data %>%
#   group_by(country) %>%
#   summarise(first = min(timestamp_server)) %>%
#   arrange(first) %>%
#   mutate(n = row_number()) %>%
#   ggplot(aes(x = first, y = n)) + 
#   geom_line() +
#   geom_point(aes(x = max(first), y = max(n)), size = 3) + 
#   scale_x_datetime(date_breaks = "2 months", date_labels = "%B %Y") +
#   labs(x = NULL, y = "Number of countries")
# 
# package_releases <- read_html("https://cran.r-project.org/src/contrib/Archive/AMR/") %>%
#   rvest::html_table() %>%
#   .[[1]] %>%
#   as_tibble(.name_repair = "unique") %>%
#   filter(`Last modified` != "") %>%
#   transmute(version = gsub("[^0-9.]", "", 
#                            gsub(".tar.gz", "", Name)), 
#             datetime = as.POSIXct(`Last modified`)) %>% 
#   # add current
#   bind_rows(tibble(version = as.character(packageVersion("AMR")),
#                    datetime = as.POSIXct(packageDate("AMR")))) %>% 
#   # remove the ones not plottable
#   filter(datetime > min(p1$data$first))
# 
# p1 + geom_linerange(data = package_releases, aes(x = datetime, ymin = 0, ymax = 80), colour = "red", inherit.aes = FALSE)
# 
