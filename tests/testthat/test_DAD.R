
ext_filepath <- system.file("extdata", package = "RaMS")
DAD_filepath <- list.files(ext_filepath, full.names = TRUE,
                           pattern = "uv_test_mini.mzML")
msdata_DAD <- grabMSdata(files = DAD_filepath, grab_what = c("everything", "DAD"))
msdata_DAD_only <- grabMSdata(files = DAD_filepath, grab_what = "DAD")

test_that("DAD grab works", {
  expect_gt(nrow(msdata_DAD$MS1), 0)
  expect_equal(nrow(msdata_DAD$TIC), 5)
  expect_equal(nrow(msdata_DAD$BPC), 5)
  expect_gt(nrow(msdata_DAD$DAD), 0)
})

test_that("DAD data shaped correctly", {
  expect_identical(names(msdata_DAD$MS1), c("rt", "mz", "int", "filename"))
  expect_identical(names(msdata_DAD$DAD), c("rt", "lambda", "int", "filename"))
  expect_gte(min(msdata_DAD$DAD$lambda), 200)
  expect_lte(max(msdata_DAD$DAD$lambda), 600)
})

test_that("DAD only data still works", {
  expect_length(msdata_DAD_only, 1)
  expect_identical(names(msdata_DAD_only), "DAD")
  expect_identical(names(msdata_DAD_only$DAD), c("rt", "lambda", "int", "filename"))
})

test_that("warning thrown when using mzXML", {
  expect_warning(grabMSdata(mzXML_filenames[2], grab_what = "DAD"),
                 regexp = "grab_what = 'DAD' not available for mzXML documents")
})

test_that("both polarities detected", {
  expect_length(msdata_DAD$metadata$polarity, 2)
})
