# get complete filenames of all R files in the GitHub repository of nathaneastwood/poorman

library(magrittr)
commit <- "3cc0a9920b1eb559dd166f548561244189586b3a"

files <- xml2::read_html(paste0("https://github.com/nathaneastwood/poorman/tree/", commit, "/R")) %>%
  rvest::html_nodes("a") %>%
  rvest::html_attr("href")
files <- files[files %like% "/blob/.*R$"]

# get full URLs of all raw R files
files <- sort(paste0("https://raw.githubusercontent.com", gsub("blob/", "", files[files %like% "/R/.*.R$"])))
# remove files with only pkg specific code
files <- files[files %unlike% "(zzz|init)[.]R$"]
# also, there's a lot of functions we don't use
files <- files[files %unlike% "/(between|coalesce|cumulative|fill|glimpse|gluestick|group_cols|na_if|near|nest_by|poorman-package|print|recode|reconstruct|replace_na|replace_with|rownames|slice|union_all|unite|window_rank|with_groups)[.]R$"]

# add our prepend file, containing info about the source of the data
intro <- readLines("data-raw/poorman_prepend.R") %>% 
  # add commit to intro part
  gsub("{commit}", commit, ., fixed = TRUE) %>%
  # add date to intro part
  gsub("{date}", trimws(format(Sys.Date(), "%e %B %Y")), ., fixed = TRUE)
# copyright info:
copyright <- paste0("# ", readLines(paste0("https://raw.githubusercontent.com/nathaneastwood/poorman/", commit, "/LICENSE")))

# read all contents to a character vector
contents <- character(0)
sapply(files, function(file) {
  message("reading ", basename(file))
  contents <<- c(contents, readLines(file))
  invisible()
})

# remove lines starting with "#'" and NULL and write to file
contents <- contents[!grepl("^(#'|NULL|\"_PACKAGE)", contents)]
contents.bak <- contents

# grouped attributes same as dplyr
contents <- gsub("grouped_data", "grouped_df", contents, fixed = TRUE)
# now make it independent on UseMethod, since we will not export these functions
has_usemethods <- gsub("^([a-z_]+).*", "\\1", contents[which(contents %like% "usemethod") - 1])
for (use in has_usemethods) {
  relevant_row <- which(contents %like% paste0("^", use, " <- function")) + 1
  function_call <- trimws(gsub(".*function(.*)\\{.*", "\\1", contents[relevant_row - 1]))
  function_call1 <- trimws(gsub("[()]", "", strsplit(function_call, ",")[[1]][1]))
  if (any(contents %like% paste0(use, ".grouped_df"))) {
    # this function will have methods for data.frame and grouped_df
    contents[relevant_row] <- paste0("  if (\"grouped_df\" %in% class(", function_call1, ")) ", use, ".grouped_df", function_call, " else ", use, ".data.frame", function_call)
  } else {
    # this function will only have data.frame as method
    contents[relevant_row] <- paste0("  ", use, ".data.frame", function_call)
  }
  # add pm_ prefix
  contents[relevant_row - 1] <- paste0("pm_", contents[relevant_row - 1])
  
}
# correct for NextMethod
contents <- gsub("NextMethod\\(\"(.*)\"\\)", "\\1.data.frame(...)", contents)
# correct for 'default' method
contents <- gsub(".default <-", ".data.frame <-", contents, fixed = TRUE)
contents <- gsub("pm_group_by_drop.data.frame", "pm_group_by_drop", contents, fixed = TRUE)
# now get all those pm_* functions to replace all untransformed function name calls as well
new_pm_names <- sort(gsub("pm_(.*?) <-.*", "\\1", contents[grepl("^pm_", contents)]))
for (i in seq_len(length(new_pm_names))) {
  contents <- gsub(paste0("([^a-z._])", new_pm_names[i], "([^a-z._])"), paste0("\\1pm_", new_pm_names[i], "\\2"), contents)
  # starting with a space or a straight bracket or an opening parenthesis, ending with nothing or a non-character or a closing parenthesis
  contents <- gsub(paste0("( |\\[|\\()", new_pm_names[i], "($|[^a-z]|\\))"), paste0("\\1pm_", new_pm_names[i], "\\2"), contents)
}
# replace %>% with %pm>%
contents[which(contents %like% "^\\|\\|") - 1] <- paste0(contents[which(contents %like% "^\\|\\|") - 1], " ||")
contents[which(contents %like% "^\\|\\|")] <- gsub("^\\|\\|", "", contents[which(contents %like% "^\\|\\|")])
contents <- gsub("%>%", "%pm>%", contents, fixed = TRUE)
# fix for new lines, since n() also existed
contents <- gsub("\\pm_n", "\\n", contents, fixed = TRUE)
# prefix other functions also with "pm_"
contents <- gsub("^([a-z_]+)(\\$|)", "pm_\\1\\2", contents)
# prefix environmental objects and functions
contents <- gsub("(eval_env|select_env|select_context|context|dotdotdot|as_symbols|insert_dot|deparse_|groups_set|apply_grouped_function|split_into_groups|calculate_groups|has_groups|eval_select_pos|select_positions|eval_expr|eval_call|add_group_columns|find_used|is_nested|setup_|select_|group_)", "pm_\\1", contents)
# now some items are overprefixed
contents <- gsub("(pm_)+", "pm_", contents)
contents <- gsub("pm_if (\"grouped_df", "if (\"grouped_df", contents, fixed = TRUE)
# remove comments and empty lines
contents <- gsub("#.*", "", contents)
contents <- contents[trimws(contents) != ""]
# fix for their relocate()
contents <- gsub("if (!missing(.before))", "if (!missing(.before) && !is.null(.before))", contents, fixed = TRUE)
contents <- gsub("if (!missing(.after))", "if (!missing(.after) && !is.null(.after))", contents, fixed = TRUE)
contents[which(contents %like% "reshape\\($") + 1] <- gsub("data", "as.data.frame(data, stringsAsFactors = FALSE)", contents[which(contents %like% "reshape\\($") + 1])
contents <- gsub('pm_relocate(.data = long, values_to, .after = -1)', 'pm_relocate(.data = long, "value", .after = -1)', contents, fixed = TRUE)

# who needs US spelling?
contents <- contents[contents %unlike% "summarize"]

# add intro
contents <- c(
  intro,
  copyright,
  "",
  contents
)

writeLines(contents, "R/aa_helper_pm_functions.R")

# note: pm_left_join() will be overwritten by aaa_helper_functions.R, which contains a faster implementation
