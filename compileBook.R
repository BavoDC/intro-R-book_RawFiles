################################################
##               Compilation                  ##
################################################

# Set correct working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Compile the book (directory: docs)
bookdown::render_book('index.Rmd', 'bookdown::gitbook')
if(!file.exists('docs/.nojekyll'))
  file.create('docs/.nojekyll')

# Create R scripts (directory: docs/scripts)
source("purl.R")


################################################
##               Instructions                 ##
################################################

# Write chapters in separate .Rmd files, file name format:
#
#   'XX-YYYYYY.RMD'
#
# Where:
#
#   XX : The chapter number
#   YY : Chosen name for the Rmd file -- Has no effect on the compiled result
#

# The RMD file for Each Chapter starts with:
#
#   '# XXXXX {#YYYYY}'
# 
# Where:
#
#   XX : The name of the chapter - This is the name shown in the compiled book
#   YY : Used for referencing the chapter in other parts of the book


################################################
##               Directories                  ##
################################################

# _book : The compiled book
# _bookdown_files : ???
# bib   : bibtex references
# data  : all data files used in the project. Only use relative paths no absolute paths
# docs  : documentation --> Contains the R scripts for each chapter
# images : all images used in the book
# rds   : simulated data, created during compilation

