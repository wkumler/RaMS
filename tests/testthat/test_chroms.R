
msdata_MRM_chroms <- grabMSdata(files = list.files(system.file("extdata", package = "RaMS"),
                              pattern = "wk_chrom.mzML.gz", full.names = TRUE),
                              grab_what="chroms")

test_that("chroms are grabbed when using chroms", {
  expect_identical(names(msdata_MRM_chroms), "chroms")
  expect_identical(names(msdata_MRM_chroms$chroms),
                   c("chrom_type", "chrom_index", "target_mz", "product_mz",
                     "rt", "int", "filename"))
})

msdata_MRM_everything <- grabMSdata(files = list.files(system.file("extdata", package = "RaMS"),
                                                       pattern = "wk_chrom.mzML.gz", full.names = TRUE),
                                    grab_what=c("everything", "chroms"))
test_that("everything works with chroms", {
  expect_identical(names(msdata_MRM_everything), c("MS1", "MS2", "BPC", "TIC", "chroms", "metadata"))
})

test_that("warning thrown when requesting chroms from mzXML", {
  expect_warning(grabMSdata(mzXML_filenames[2], grab_what = "chroms"),
                 regexp = "grab_what = 'chroms' not available for mzXML documents")
})
