## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
options(tidyverse.quiet = TRUE)
data.table::setDTthreads(2)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "80%",
  fig.align = 'center',
  fig.height = 3,
  fig.width = 6.5
)

## ----minifyfile, results='hide'-----------------------------------------------
library(RaMS)
msdata_files <- list.files(
  system.file("extdata", package = "RaMS"), full.names = TRUE, pattern = "mzML"
)[1:4]

initial_filename <- msdata_files[1]
output_filename <- gsub(x=paste0("minified_", basename(initial_filename)), "\\.gz", "")

masses_of_interest <- c(118.0865, 138.0555)
minifyMSdata(files = initial_filename, output_files = output_filename, 
             mz_include = masses_of_interest, ppm = 10, warn = FALSE)

## ----opendata, results='hide'-------------------------------------------------
init_msdata <- grabMSdata(initial_filename)
msdata <- grabMSdata(output_filename)

## ----showmini-----------------------------------------------------------------
knitr::kable(head(msdata$MS1, 3))
knitr::kable(head(msdata$MS2, 3))

## ----TICBPCdiffs--------------------------------------------------------------
par(mfrow=c(2, 1), mar=c(2.1, 2.1, 1.1, 0.1))
plot(init_msdata$BPC$rt, init_msdata$BPC$int, type = "l", main = "Initial BPC")
plot(msdata$BPC$rt, msdata$BPC$int, type = "l", main = "New BPC")

## ----cleanmini, include=FALSE-------------------------------------------------
layout(1)
unlink(output_filename)

## ----vectorizedmini-----------------------------------------------------------
dir.create("mini_mzMLs/")
output_files <- paste0("mini_mzMLs/", basename(msdata_files))
output_files <- gsub(x=output_files, "\\.gz", "")

minifyMSdata(files = msdata_files, output_files = output_files, verbosity = 0,
             mz_include = masses_of_interest, ppm = 10, warn = FALSE)

## ----ggplotmini, fig.height=2-------------------------------------------------
mini_msdata <- grabMSdata(output_files, verbosity = 0)

library(ggplot2)
ggplot(mini_msdata$BPC) + geom_line(aes(x=rt, y=int, color=filename)) + theme_bw()

## ----cleanminimulti, include=FALSE--------------------------------------------
unlink("mini_mzMLs", recursive = TRUE)

## ----mzs to include, eval=FALSE-----------------------------------------------
#  raw_stans <- read.csv(paste0("https://raw.githubusercontent.com/",
#                               "IngallsLabUW/Ingalls_Standards/",
#                               "b098927ea0089b6e7a31e1758e7c7eaad5408535/",
#                               "Ingalls_Lab_Standards_NEW.csv"))
#  
#  mzs_to_include <- as.numeric(unique(raw_stans[raw_stans$Fraction1=="HILICPos",]$m.z))
#  # Include glycine betaine isotopes for README demo
#  mzs_to_include <- c(mzs_to_include, 119.0899, 119.0835)

## ----download the mzXMLs, eval=FALSE------------------------------------------
#  if(!dir.exists("vignettes/data"))dir.create("vignettes/data")
#  base_url <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS703/"
#  chosen_files <- paste0(base_url, "170223_Smp_LB12HL_", c("AB", "CD", "EF"), "_pos.mzXML")
#  new_names <- gsub(x=basename(chosen_files), "170223_Smp_", "")
#  
#  mapply(download.file, chosen_files, paste0("vignettes/data/", new_names),
#         mode = "wb", method = "libcurl")

## ----get DDA, eval=FALSE------------------------------------------------------
#  file.copy(from = paste0("Z:/1_QEdata/2016/2016_Katherine_1335_LightB12_",
#                          "Experiment/170223_KRH_Rerun_1335_LightB12_Exp_HILIC/",
#                          "positive/",
#                          "170223_Poo_AllCyanoAqExtracts_DDApos_2.mzXML"),
#            to = "vignettes/data/DDApos_2.mzXML", overwrite = TRUE)

## ----minify, warning=FALSE, eval=FALSE----------------------------------------
#  library(RaMS)
#  
#  if(!dir.exists("inst/extdata"))dir.create("inst/extdata", recursive = TRUE)
#  init_files <- list.files("vignettes/data/", full.names = TRUE)
#  out_files <- paste0("inst/extdata/", basename(init_files))
#  minifyMSdata(files = init_files, output_files = out_files, warn = FALSE,
#               mz_include = mzs_to_include, ppm = 20)

## ----msconvertchunk, eval=FALSE-----------------------------------------------
#  system("msconvert inst/extdata/*.mzXML -o inst/extdata/temp --noindex")
#  system("msconvert --mzXML inst/extdata/*.mzXML -o inst/extdata/temp --noindex")
#  system('msconvert inst/extdata/temp/*.mzML --filter \"scanTime [240,900]\" -o inst/extdata -g')
#  system('msconvert inst/extdata/temp/*.mzXML --mzXML --filter \"scanTime [240,900]\" -o inst/extdata -g')

## ----renameremove, eval=FALSE-------------------------------------------------
#  init_files <- list.files("inst/extdata", full.names = TRUE)
#  new_names <- paste0("inst/extdata/", gsub(x=init_files, ".*(Smp_|Extracts_)", ""))
#  file.rename(init_files, new_names)
#  
#  unlink("inst/extdata/temp", recursive = TRUE)
#  file.remove(list.files("inst/extdata", pattern = "mzXML$", full.names = TRUE))
#  file.remove(paste0("inst/extdata/", c("LB12HL_CD.mzXML.gz", "LB12HL_EF.mzXML.gz")))

## ----checkquality, eval=FALSE-------------------------------------------------
#  MSnbase::readMSData(list.files("inst/extdata", full.names = TRUE)[1], msLevel. = 1)
#  RaMS::grabMSdata(new_names[1])

## ----vigdata cleanup, eval=FALSE----------------------------------------------
#  unlink("vignettes/data", recursive = TRUE)

