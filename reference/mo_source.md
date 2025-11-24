# User-Defined Reference Data Set for Microorganisms

These functions can be used to predefine your own reference to be used
in [`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and
consequently all
[`mo_*`](https://amr-for-r.org/reference/mo_property.md) functions (such
as [`mo_genus()`](https://amr-for-r.org/reference/mo_property.md) and
[`mo_gramstain()`](https://amr-for-r.org/reference/mo_property.md)).

This is **the fastest way** to have your organisation (or analysis)
specific codes picked up and translated by this package, since you don't
have to bother about it again after setting it up once.

## Usage

``` r
set_mo_source(path, destination = getOption("AMR_mo_source",
  "~/mo_source.rds"))

get_mo_source(destination = getOption("AMR_mo_source", "~/mo_source.rds"))
```

## Arguments

- path:

  Location of your reference file, this can be any text file (comma-,
  tab- or pipe-separated) or an Excel file (see *Details*). Can also be
  `""`, `NULL` or `FALSE` to delete the reference file.

- destination:

  Destination of the compressed data file - the default is the user's
  home directory.

## Details

The reference file can be a text file separated with commas (CSV) or
tabs or pipes, an Excel file (either 'xls' or 'xlsx' format) or an R
object file (extension '.rds'). To use an Excel file, you will need to
have the `readxl` package installed.

`set_mo_source()` will check the file for validity: it must be a
[data.frame](https://rdrr.io/r/base/data.frame.html), must have a column
named `"mo"` which contains values from
[`microorganisms$mo`](https://amr-for-r.org/reference/microorganisms.md)
or
[`microorganisms$fullname`](https://amr-for-r.org/reference/microorganisms.md)
and must have a reference column with your own defined values. If all
tests pass, `set_mo_source()` will read the file into R and will ask to
export it to `"~/mo_source.rds"`. The CRAN policy disallows packages to
write to the file system, although '*exceptions may be allowed in
interactive sessions if the package obtains confirmation from the
user*'. For this reason, this function only works in interactive
sessions so that the user can **specifically confirm and allow** that
this file will be created. The destination of this file can be set with
the `destination` argument and defaults to the user's home directory. It
can also be set with the package option
[`AMR_mo_source`](https://amr-for-r.org/reference/AMR-options.md), e.g.
`options(AMR_mo_source = "my/location/file.rds")`.

The created compressed data file `"mo_source.rds"` will be used at
default for MO determination (function
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md) and consequently
all `mo_*` functions like
[`mo_genus()`](https://amr-for-r.org/reference/mo_property.md) and
[`mo_gramstain()`](https://amr-for-r.org/reference/mo_property.md)). The
location and timestamp of the original file will be saved as an
[attribute](https://rdrr.io/r/base/attributes.html) to the compressed
data file.

The function `get_mo_source()` will return the data set by reading
`"mo_source.rds"` with
[`readRDS()`](https://rdrr.io/r/base/readRDS.html). If the original file
has changed (by checking the location and timestamp of the original
file), it will call `set_mo_source()` to update the data file
automatically if used in an interactive session.

Reading an Excel file (`.xlsx`) with only one row has a size of 8-9 kB.
The compressed file created with `set_mo_source()` will then have a size
of 0.1 kB and can be read by `get_mo_source()` in only a couple of
microseconds (millionths of a second).

## How to Setup

Imagine this data on a sheet of an Excel file. The first column contains
the organisation specific codes, the second column contains valid
taxonomic names:

      |         A          |            B          |
    --|--------------------|-----------------------|
    1 | Organisation XYZ   | mo                    |
    2 | lab_mo_ecoli       | Escherichia coli      |
    3 | lab_mo_kpneumoniae | Klebsiella pneumoniae |
    4 |                    |                       |

We save it as `"/Users/me/Documents/ourcodes.xlsx"`. Now we have to set
it as a source:

    set_mo_source("/Users/me/Documents/ourcodes.xlsx")
    #> NOTE: Created mo_source file '/Users/me/mo_source.rds' (0.3 kB) from
    #>       '/Users/me/Documents/ourcodes.xlsx' (9 kB), columns
    #>       "Organisation XYZ" and "mo"

It has now created a file `"~/mo_source.rds"` with the contents of our
Excel file. Only the first column with foreign values and the 'mo'
column will be kept when creating the RDS file.

And now we can use it in our functions:

    as.mo("lab_mo_ecoli")
    #> Class 'mo'
    #> [1] B_ESCHR_COLI

    mo_genus("lab_mo_kpneumoniae")
    #> [1] "Klebsiella"

    # other input values still work too
    as.mo(c("Escherichia coli", "E. coli", "lab_mo_ecoli"))
    #> NOTE: Translation to one microorganism was guessed with uncertainty.
    #>       Use mo_uncertainties() to review it.
    #> Class 'mo'
    #> [1] B_ESCHR_COLI B_ESCHR_COLI B_ESCHR_COLI

If we edit the Excel file by, let's say, adding row 4 like this:

      |         A          |            B          |
    --|--------------------|-----------------------|
    1 | Organisation XYZ   | mo                    |
    2 | lab_mo_ecoli       | Escherichia coli      |
    3 | lab_mo_kpneumoniae | Klebsiella pneumoniae |
    4 | lab_Staph_aureus   | Staphylococcus aureus |
    5 |                    |                       |

...any new usage of an MO function in this package will update your data
file:

    as.mo("lab_mo_ecoli")
    #> NOTE: Updated mo_source file '/Users/me/mo_source.rds' (0.3 kB) from
    #>       '/Users/me/Documents/ourcodes.xlsx' (9 kB), columns
    #>        "Organisation XYZ" and "mo"
    #> Class 'mo'
    #> [1] B_ESCHR_COLI

    mo_genus("lab_Staph_aureus")
    #> [1] "Staphylococcus"

To delete the reference data file, just use `""`, `NULL` or `FALSE` as
input for `set_mo_source()`:

    set_mo_source(NULL)
    #> Removed mo_source file '/Users/me/mo_source.rds'

If the original file (in the previous case an Excel file) is moved or
deleted, the `mo_source.rds` file will be removed upon the next use of
[`as.mo()`](https://amr-for-r.org/reference/as.mo.md) or any
[`mo_*`](https://amr-for-r.org/reference/mo_property.md) function.
