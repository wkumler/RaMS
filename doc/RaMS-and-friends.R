## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
library(RaMS)

# Locate an MS file
single_file <- system.file("extdata", "LB12HL_AB.mzML.gz", package = "RaMS")

# Grab the MS data
msdata <- grabMSdata(single_file, grab_what = "everything")

# Write out MS1 data to .csv file
write.csv(x = msdata$MS1, file = "MS1_data.csv")

# Clean up afterward
file.remove("MS1_data.csv")

## -----------------------------------------------------------------------------
library(openxlsx)

# Locate an MS2 file
MS2_file <- system.file("extdata", "DDApos_2.mzML.gz", package = "RaMS")

# Grab the MS1 and MS2 data
msdata <- grabMSdata(MS2_file, grab_what=c("MS1", "MS2"))

# Write out MS data to Excel file
# openxlsx writes each object in a list to a unique sheet
# Produces one sheet for MS1 and one for MS2
write.xlsx(msdata, file = "MS2_data.xlsx")

# Clean up afterward
file.remove("MS2_data.xlsx")

## -----------------------------------------------------------------------------
library(DBI)
# Get data from multiple files to show off
mzml_files <- system.file(c("extdata/LB12HL_AB.mzML.gz", 
                            "extdata/LB12HL_CD.mzML.gz"), 
                          package = "RaMS")
msdata <- grabMSdata(mzml_files)

# Create the sqlite database and connect to it
MSdb <- dbConnect(RSQLite::SQLite(), "MSdata.sqlite")

# Export MS1 and MS2 data to sqlite tables
dbWriteTable(MSdb, "MS1", msdata$MS1)
dbWriteTable(MSdb, "MS2", msdata$MS2)
dbListTables(MSdb)

# Perform a simple query to ensure data was exported correctly
dbGetQuery(MSdb, 'SELECT * FROM MS1 LIMIT 3')

# Perform EIC extraction in SQL rather than in R
EIC_query <- 'SELECT * FROM MS1 WHERE mz BETWEEN :lower_bound AND :upper_bound'
query_params <- list(lower_bound=118.086, upper_bound=118.087)
EIC <- dbGetQuery(MSdb, EIC_query, params = query_params)

# Append with additional files
extra_file <- system.file("extdata", "LB12HL_EF.mzML.gz", package = "RaMS")
extra_msdata <- grabMSdata(extra_file, grab_what = "everything")
dbGetQuery(MSdb, 'SELECT COUNT(*) FROM MS1') # Initial number of rows
dbAppendTable(MSdb, "MS1", extra_msdata$MS1)
# Confirm three different files exist in DB
dbGetQuery(MSdb, 'SELECT DISTINCT filename FROM MS1')
# Confirm new rows have been added
dbGetQuery(MSdb, 'SELECT COUNT(*) FROM MS1')

# Disconnect after export
dbDisconnect(MSdb)

# Clean up afterward
unlink("MSdata.sqlite")

## -----------------------------------------------------------------------------
# Locate a couple MS files
data_dir <- system.file("extdata", package = "RaMS")
file_paths <- list.files(data_dir, pattern = "HL.*mzML", full.names = TRUE)

msdata <- grabMSdata(files = file_paths, grab_what = "BPC")$BPC

