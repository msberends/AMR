# get complete filenames of all R files in the GitHub repository of nathaneastwood/poorman

commit <- "7d76d77f8f7bc663bf30fb5a161abb49801afa17"

files <- xml2::read_html(paste0("https://github.com/nathaneastwood/poorman/tree/", commit, "/R")) %>% 
  rvest::html_nodes("table") %>%
  rvest::html_table()
files <- files[[1]][,"Name"]

# remove files with only pkg specific code
files <- files[!files %in% c("zzz.R", "init.R")]
files <- paste0("https://raw.githubusercontent.com/nathaneastwood/poorman/", commit, "/R/",
                files[grepl("[.]R$", files)])

# add our prepend file, containing info about the source of the data
files <- c("data-raw/poorman_prepend.R", files)

# read all contents to a character vector
contents <- character(0)
sapply(files, function(file) {
  contents <<- c(contents, readLines(file))
  invisible()
})

# remove lines starting with "#'" and NULL and write to file
contents <- contents[!grepl("^(#'|NULL|\"_PACKAGE)", contents)]

# now make it independent on UseMethod, since we will not export these functions
contents <- gsub('UseMethod[(]"(.*?)"[)]',
                 'if ("grouped_data" %in% class(.data)) {|||    \\1.grouped_data(.data, ...)|||  } else {|||    \\1.default(.data, ...)|||  }', 
                 paste(contents, collapse = "|||"),
                 perl = TRUE) %>% 
  # add commit to intro part
  gsub("{commit}", commit, ., fixed = TRUE) %>% 
  strsplit(split = "|||", fixed = TRUE) %>% 
  unlist()

writeLines(contents, "R/aa_helper_functions_dplyr.R")
