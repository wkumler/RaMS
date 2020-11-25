library(testthat)
library(MSnbase)
library(data.table)

mzML_filename <- system.file("proteomics", "MS3TMT11.mzML", package = "msdata")
mzXML_filename <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                              package = "msdata")

mzML_data <- grabMSdata(files = mzML_filename, grab_what = "everything", verbosity = "none")
mzML_EICs <- grabMSdata(files = mzML_filename, grab_what = c("EIC", "EIC_MS2"),
                        verbosity = "none", mz = 118.0865, ppm=5)
mzML_trimmed <- grabMSdata(files = mzML_filename, grab_what = "everything",
                           verbosity = "none", rtrange = c(45, 46))

test_check("RaMS")
