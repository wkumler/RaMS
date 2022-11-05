
msdata_MRM_chroms <- grabMSdata(files = list.files(system.file("extdata", package = "RaMS"),
                              pattern = "wk_chrom.mzML.gz", full.names = TRUE))

test_that("chroms are grabbed when using chroms", {
  expect_identical(mzML_everything, mzML_data)
})

test_that("warning thrown when requesting chroms from mzXML", {
  expect_warning(grabMSdata(mzXML_filenames[1]))
})
