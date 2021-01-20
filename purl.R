dir.create("docs/scripts")

files <- list.files()

for(i in 1:length(files))
{
  if(endsWith(files[i], '.Rmd'))
  {
    R_file <- stringr::str_c("docs/scripts/", substr(files[i], 1, nchar(files[i])-2))
    knitr::purl(files[i], R_file)
  }
}
