
# Read and format data ----------------------------------------------------

library(tidyverse)
library(maps)

# get website analytics
source("data-raw/country_analysis_url_token.R")
url_json <- paste0(country_analysis_url, 
                   "/index.php?&module=API&token_auth=", 
                   country_analysis_token,
                   "&method=Live.getLastVisitsDetails&idSite=3&language=en&expanded=1&date=2018-01-01,2028-01-01&period=range&filter_limit=-1&format=JSON&segment=&translateColumnNames=1")

data_json <- jsonlite::read_json(url_json)
data <- tibble(
  timestamp_server = as.POSIXct(sapply(data_json, function(x) x$serverTimestamp), origin = "1970-01-01"),
  country = sapply(data_json, function(x) x$country))


# Plot world map ----------------------------------------------------------

countries_name <- sort(unique(data$country))
countries_name <- countries_name[countries_name != "Unknown"]
countries_iso <- countrycode::countrycode(countries_name, 'country.name', 'iso3c')

world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE)) %>% 
  mutate(countries_code = countrycode::countrycode(ID, 'country.name', 'iso3c'),
         included = as.integer(countries_code %in% countries_iso)) %>% 
  mutate(not_antarctica = as.integer(ID != "Antarctica"))

(ggplot(world1) +
    geom_sf(aes(fill = included, colour = not_antarctica), size = 0.25) +
    theme_minimal() +
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank()) +
    scale_fill_gradient(low = "white", high = "#CAD6EA") +
    # this makes the border Antarctica turn white (invisible):
    scale_colour_gradient(low = "white", high = "#81899B") +
    geom_text(aes(x = -170,
                  y = -70,
                  label = stringr::str_wrap(paste0("Accented countries (n = ", 
                                                   length(countries_name), "): ", 
                                                   paste(countries_name, collapse = ", ")),
                                            225)),
              hjust = 0,
              size = 4)) %>% 
  ggsave("pkgdown/logos/countries.png", dpi = 300, plot = ., scale = 1.5)


# Gibberish ---------------------------------------------------------------

p1 <- data %>%
  group_by(country) %>%
  summarise(first = min(timestamp_server)) %>%
  arrange(first) %>%
  mutate(n = row_number()) %>%
  ggplot(aes(x = first, y = n)) + 
  geom_line() +
  geom_point(aes(x = max(first), y = max(n)), size = 3) + 
  scale_x_datetime(date_breaks = "2 months", date_labels = "%B %Y") +
  labs(x = NULL, y = "Number of countries")

package_releases <- read_html("https://cran.r-project.org/src/contrib/Archive/AMR/") %>%
  rvest::html_table() %>%
  .[[1]] %>%
  as_tibble(.name_repair = "unique") %>%
  filter(`Last modified` != "") %>%
  transmute(version = gsub("[^0-9.]", "", 
                           gsub(".tar.gz", "", Name)), 
            datetime = as.POSIXct(`Last modified`)) %>% 
  # add current
  bind_rows(tibble(version = as.character(packageVersion("AMR")),
                   datetime = as.POSIXct(packageDate("AMR")))) %>% 
  # remove the ones not plottable
  filter(datetime > min(p1$data$first))

p1 + geom_linerange(data = package_releases, aes(x = datetime, ymin = 0, ymax = 80), colour = "red", inherit.aes = FALSE)

