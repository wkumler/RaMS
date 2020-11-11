library(testthat)
library(RaMS)
library(MSnbase)

mzML_filename <- system.file("proteomics", "MS3TMT11.mzML", package = "msdata")
mzXML_filename <- system.file("proteomics", "threonine_i2_e35_pH_tree.mzXML",
                              package = "msdata")
mzML_MS2_filename <- system.file(
  "proteomics", "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz",
  package = "msdata"
)

test_check("RaMS")
