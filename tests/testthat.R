library(testthat)
library(RaMS)

mzML_filename <- system.file("extdata",
                             "180205_Poo_TruePoo_Full2.mzML",
                             package = "RaMS")

test_check("RaMS")
