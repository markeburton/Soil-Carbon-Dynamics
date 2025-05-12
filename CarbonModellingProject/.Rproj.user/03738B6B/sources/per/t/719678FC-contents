
#### Load Libraries ####
library(ggplot2)
library(readxl)
library(openxlsx)
library(tidyverse)
library(inTextSummaryTable)
library(Hmisc)

#library(drcarlate)

#### Directories

#Fill in Project Name to be used in output files names. 
Project_initials <- "XXX"

# UPDATE TO LINK TO THE FOLDER STRUCTURE ON YOUR MACHINE.
parent_dir <- "FILE PATH TO THE FOLDER STRUCTURE GOES HERE/CarbonModellingProject/"

# Folder where all data file structure is located.
source_dir <- paste0(parent_dir, "InputFiles/")
# Folder where all the map outputs will live.
out_dir <- paste0(parent_dir, "OutputFiles/")
#Today's date for automating the naming convention of the files.
date <- paste(format(Sys.Date(), format="%Y"), format(Sys.Date(), format="%m"), format(Sys.Date(), format="%d"), sep="")
date

#Turn off scientific notation in tables if desired
options(scipen = 9999)

