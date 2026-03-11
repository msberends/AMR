# Download and Unpack an AMR Course Repository

Downloads and unpacks a GitHub repository containing course materials,
using `usethis::use_course()`. This is a convenience wrapper intended
for use in educational settings, such as workshops or tutorials
associated with the AMR package.

## Usage

``` r
amr_course(github_repo, branch = "main", ...)
```

## Arguments

- github_repo:

  A character string specifying the GitHub repository with username and
  repo name, e.g. `"https://github.com/username/repo"`.

- branch:

  A character string specifying the branch to download. Defaults to
  `"main"`.

- ...:

  Additional arguments passed on to `usethis::use_course()`.

## Value

Called for its side effect. `usethis::use_course()` will prompt the user
to choose a destination and open the extracted project. Returns
invisibly whatever `usethis::use_course()` returns.

## Details

This function constructs a ZIP archive URL from the provided
`github_repo` and `branch`, then delegates to `usethis::use_course()` to
handle the download and extraction.

The function is designed for interactive use in course or workshop
settings and is not intended for use in non-interactive or automated
pipelines.

## See also

`usethis::use_course()`

## Examples

``` r
if (FALSE) { # \dontrun{

# Let this run by users, e.g., webinar participants
amr_course("https://github.com/my_user_name/our_AMR_course")
} # }
```
