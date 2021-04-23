# get complete filenames of all R files in the GitHub repository of nathaneastwood/poorman

commit <- "52eb6947e0b4430cd588976ed8820013eddf955f"

files <- xml2::read_html(paste0("https://github.com/nathaneastwood/poorman/tree/", commit, "/R")) %>%
  rvest::html_nodes("a") %>%
  rvest::html_attr("href")

# get full URLs of all raw R files
files <- sort(paste0("https://raw.githubusercontent.com", gsub("blob/", "", files[files %like% "/R/.*.R$"])))
# remove files with only pkg specific code
files <- files[files %unlike% "(zzz|init)[.]R$"]
# also, there's a lot of functions we don't use
files <- files[files %unlike% "(slice|glimpse|recode|replace_na|coalesce)[.]R$"]

# add our prepend file, containing info about the source of the data
intro <- readLines("data-raw/poorman_prepend.R")
# copyright info:
copyright <- paste0("# ", readLines("https://raw.githubusercontent.com/nathaneastwood/poorman/master/LICENSE"))

# read all contents to a character vector
contents <- character(0)
sapply(files, function(file) {
  message("reading ", basename(file))
  contents <<- c(contents, readLines(file))
  invisible()
})
contents <- c(intro,
              copyright,
              "",
              contents)

# remove lines starting with "#'" and NULL and write to file
contents <- contents[!grepl("^(#'|NULL|\"_PACKAGE)", contents)]

# now make it independent on UseMethod, since we will not export these functions

contents <- gsub('UseMethod[(]"(.*?)"[)]',
                 'if ("grouped_data" %in% class(.data)) {|||    \\1.grouped_data(.data, ...)|||  } else {|||    \\1.default(.data, ...)|||  }', 
                 paste(contents, collapse = "|||"),
                 perl = TRUE) %>%
  # add commit to intro part
  gsub("{commit}", commit, ., fixed = TRUE) %>%
  # add date to intro part
  gsub("{date}", format(Sys.Date(), "%e %B %Y"), ., fixed = TRUE) %>%
  strsplit(split = "|||", fixed = TRUE) %>%
  unlist() %>% 
  # add "pm_" as prefix to all functions
  gsub("^([a-z_.]+) <- function", "pm_\\1 <- function", .)

# now get all those pm_* functions to replace all untransformed function name calls as well
new_pm_names <- sort(gsub("pm_(.*?) <-.*", "\\1", contents[grepl("^pm_", contents)]))
for (i in seq_len(length(new_pm_names))) {
  contents <- gsub(paste0("([^a-z._])", new_pm_names[i], "([^a-z._])"), paste0("\\1pm_", new_pm_names[i], "\\2"), contents)
  # starting with a space or a straight bracket or an opening parenthesis, ending with nothing or a non-character or a closing parenthesis
  contents <- gsub(paste0("( |\\[|\\()", new_pm_names[i], "($|[^a-z]|\\))"), paste0("\\1pm_", new_pm_names[i], "\\2"), contents)
}

# replace %>% with %pm>% 
contents <- gsub("%>%", "%pm>%", contents, fixed = TRUE)
# fix for new lines, since n() also existed
contents <- gsub("\\pm_n", "\\n", contents, fixed = TRUE)
# prefix other functions also with "pm_"
contents <- gsub("^([a-z_]+)(\\$|)", "pm_\\1\\2", contents)
# prefix environments
contents <- gsub("eval_env", "pm_eval_env", contents, fixed = TRUE)
contents <- gsub("select_env", "pm_select_env", contents, fixed = TRUE)
contents <- gsub("context", "pm_context", contents, fixed = TRUE)
# now some items are overprefixed
contents <- gsub("(pm_)+", "pm_", contents)
# special case for pm_distinct(), we need '.keep_all' to work
contents <- gsub("pm_distinct <- function(.data, ..., .keep_all = FALSE)", "pm_distinct <- function(.data, ...)", contents, fixed = TRUE)

# who needs US spelling?
contents <- contents[!grepl("summarize", contents)]

writeLines(contents, "R/aa_helper_pm_functions.R")

# after this, comment out:
# pm_left_join() since we use a faster version
# pm_group_split() since we don't use it and it relies on R 3.5.0 for the use of ...length(), which is hard to support with C++ code
