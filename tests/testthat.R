library(testthat)
library(RaMS)
library(mzR)

mzML_filename <- system.file("extdata",
                             "190715_Poo_TruePooFK180310_Full2.mzML",
                             package = "RaMS")
mzXML_filename <- system.file("extdata",
                              "190715_Poo_TruePooFK180310_Full2.mzXML",
                              package = "RaMS")

test_check("RaMS")
