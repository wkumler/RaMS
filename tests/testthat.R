library(testthat)
library(RaMS)

mzML_filenames <- list.files(system.file("extdata", package = "RaMS"),
                             pattern = "mzML", full.names = TRUE)
mzXML_filenames <- list.files(system.file("extdata", package = "RaMS"),
                              pattern = "mzXML", full.names = TRUE)

mzML_data <- grabMSdata(files = mzML_filenames, grab_what = c("MS1", "MS2", "TIC", "BPC"), verbosity = "none")
mzML_EICs <- grabMSdata(files = mzML_filenames, grab_what = c("EIC", "EIC_MS2"),
                        verbosity = "none", mz = 118.0865, ppm=5)

mzXML_data <- grabMSdata(files = mzXML_filenames, grab_what = c("MS1", "MS2", "TIC", "BPC"), verbosity = "none")
mzXML_EICs <- grabMSdata(files = mzXML_filenames, c("EIC", "EIC_MS2"),
                         verbosity = "none", mz = 118.0865, ppm=5)

test_check("RaMS")
