library(RaMS)
library(data.table)

mzML_filenames <- list.files(system.file("extdata", package = "RaMS"),
                             pattern = "mzML", full.names = TRUE)
mzXML_filenames <- list.files(system.file("extdata", package = "RaMS"),
                              pattern = "mzXML", full.names = TRUE)

mzML_everything <- grabMSdata(mzML_filenames[1:2], grab_what = "everything", verbosity=2)

mzML_data <- grabMSdata(mzML_filenames[1:2], grab_what = c("BPC", "TIC", "MS1", "MS2", "metadata"))

mzXML_data <- grabMSdata(mzXML_filenames[1:2], grab_what = "everything", verbosity=2)

mzML_multi_data <- grabMSdata(mzML_filenames[2:4], grab_what = c("TIC", "MS1", "metadata"))

EIC_mz <- 118.0865
EIC_ppm <- 5

mzML_EIC_data <- grabMSdata(mzML_filenames[1:2], grab_what = c("EIC", "EIC_MS2"),
                            mz = EIC_mz, ppm = EIC_ppm, verbosity = 2)
mzXML_EIC_data <- grabMSdata(mzXML_filenames[1:2], grab_what = c("EIC", "EIC_MS2"),
                             mz = EIC_mz, ppm = EIC_ppm, verbosity = 2)

test_rtrange <- c(5, 6)
mzML_trimmed <- grabMSdata(files = mzML_filenames[1:2],
                           grab_what = c("MS1", "MS2", "TIC", "BPC"),
                           rtrange = test_rtrange)
mzXML_trimmed <- grabMSdata(files = mzXML_filenames[1],
                           grab_what = c("MS1", "MS2", "TIC", "BPC"),
                           rtrange = test_rtrange)


test_zerorange <- c(5, 5)
mzML_zeroed <- grabMSdata(files = mzML_filenames[1:2],
                          grab_what = c("MS1", "MS2", "TIC", "BPC"),
                          rtrange = test_zerorange)


