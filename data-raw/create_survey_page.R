
license_text <- readLines("docs/LICENSE-text.html")
license_text <- paste(license_text, collapse = "|||")
license_text <- gsub("licen(s|c)e", "Survey", license_text, ignore.case = TRUE)
license_text <- gsub("<pre>.*</pre>", '<p>If no form is showing below, please <a href="https://forms.office.com/Pages/ResponsePage.aspx?id=-SJRM_TUZ02i_M1twg3ecDlnO1BBtdxGi-GnYu7DKfdUNTFLQ0xVSUlRVVlXTlVTNjZFMjdRUVpCSy4u" target="_blank">click here to open it</a>.</p><iframe width="100%" height= "500px" src= "https://forms.office.com/Pages/ResponsePage.aspx?id=-SJRM_TUZ02i_M1twg3ecDlnO1BBtdxGi-GnYu7DKfdUNTFLQ0xVSUlRVVlXTlVTNjZFMjdRUVpCSy4u&embed=true" frameborder= "0" marginwidth= "0" marginheight= "0" style= "border: none; max-width:100%; max-height:100vh" allowfullscreen webkitallowfullscreen mozallowfullscreen msallowfullscreen> </iframe>', license_text)
writeLines(unlist(strsplit(license_text, "|||", fixed = TRUE)), "docs/survey.html")
