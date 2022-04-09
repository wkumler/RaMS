## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
options(tidyverse.quiet = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "80%",
  fig.align = 'center',
  fig.height = 3,
  fig.width = 6.5
)

## ----listfiles----------------------------------------------------------------
files_to_convert <- list.files(
  system.file("extdata", package = "RaMS"), full.names = TRUE, pattern = "mzML"
)[2:4]

## ----maketmzmls---------------------------------------------------------------
library(RaMS)
# Create a folder to hold the new files
dir.create("tmzMLs")

# Convert a single file
file_to_create <- "tmzMLs/LB12HL_AB.tmzML"
tmzmlMaker(input_filename = files_to_convert[1], output_filename = file_to_create)

# Convert multiple files
files_to_create <- paste0("tmzMLs/", basename(files_to_convert))
# Make sure they end in .tmzML!
files_to_create <- gsub(x = files_to_create, "\\.mzML.*", ".tmzML")
# Loop over each file input/output pair
created_files <- mapply(tmzmlMaker, files_to_convert, files_to_create)

## ----grabmsdata---------------------------------------------------------------
msdata <- grabMSdata(created_files, verbosity=0)

## ----subset$------------------------------------------------------------------
ms_data_table <- msdata$MS1[mz%between%pmppm(152.05723, 5)]

## ----ggplot, fig.width=8, warning=FALSE---------------------------------------
library(ggplot2)
ggplot(ms_data_table) + geom_line(aes(x=rt, y=int, color=filename)) + xlim(8, 9.5)

## ----dplyr--------------------------------------------------------------------
library(dplyr)
ms_data_table %>%
  filter(rt%between%c(8.4, 8.9)) %>%
  group_by(filename) %>%
  summarise(area=sum(int))

## ----print msdata_connection--------------------------------------------------
print(msdata)

## ----str msdata_connection----------------------------------------------------
str(msdata)

## ----expected errors, error=TRUE----------------------------------------------
# Cannot order data by mass
msdata$MS1[-mz]
# Cannot request the top few rows
msdata[1:10,]
# Cannot subset by retention time
msdata$MS1[rt %between% c(8, 9.5)]

## ----cleanup------------------------------------------------------------------
unlink("tmzMLs", recursive = TRUE)

