library(testthat)
library(data.table)
library(RaMS)

mzML_filenames <- list.files(system.file("extdata", package = "RaMS"),
                             pattern = "mzML", full.names = TRUE)
mzXML_filenames <- list.files(system.file("extdata", package = "RaMS"),
                              pattern = "mzXML", full.names = TRUE)

mzML_data <- grabMSdata(files = mzML_filenames, grab_what = "everything", verbosity = "none")
mzML_EICs <- grabMSdata(files = mzML_filenames, grab_what = c("EIC", "EIC_MS2"),
                        verbosity = "none", mz = 118.0865, ppm=5)
mzML_trimmed <- grabMSdata(files = mzML_filenames, grab_what = "everything",
                           verbosity = "none", rtrange = c(45, 46))

test_check("RaMS")
