library(testthat)
library(data.table)
library(RaMS)

mzML_filenames <- list.files(system.file("extdata", package = "RaMS"),
                             pattern = "mzML", full.names = TRUE)
mzXML_filenames <- list.files(system.file("extdata", package = "RaMS"),
                              pattern = "mzXML", full.names = TRUE)

mzML_data <- grabMSdata(files = mzML_filenames, grab_what = "everything", verbosity = "minimal")
mzML_EICs <- grabMSdata(files = mzML_filenames, grab_what = c("EIC", "EIC_MS2"),
                        verbosity = "minimal", mz = 118.0865, ppm=5)
mzML_trimmed <- grabMSdata(files = mzML_filenames, grab_what = "everything",
                           verbosity = "minimal", rtrange = c(5, 6))

test_check("RaMS")
