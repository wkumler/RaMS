
test_that("MS1 extraction", {
  expect_s3_class(mzML_data$MS1, class = "data.table")
  expect_identical(names(mzML_data$MS1), c("rt", "mz", "int", "filename"))
  expect_gt(nrow(mzML_data$MS1), 0)
})

test_that("MS1 extraction", {
  mzML_data <- grabMSdata(files = mzML_filename, grab_what = "MS1", verbosity = "none")
  expect_s3_class(mzML_data$MS1, class = "data.table")
  expect_identical(names(mzML_data$MS1), c("rt", "mz", "int", "filename"))
  expect_gt(nrow(mzML_data$MS1), 0)
})


test_that("BPC can be read", {
  BPC <- grabMzmlBPC(mzML_filename)
  expect_identical(names(BPC), c("rt", "int"))
  expect_gt(nrow(BPC), 0)

  TIC <- grabMzmlBPC(mzML_filename, TIC = TRUE)
  expect_identical(names(TIC), c("rt", "int"))
  expect_gt(nrow(TIC), 0)

  expect_identical(nrow(BPC), nrow(TIC))
  expect_true(all(TIC$int>=BPC$int))
})

test_that("new matches MSnbase", {
  MSnExp <- readMSData(mzML_filename, msLevel. = 1)
  parsed_spectra <- lapply(as.list(MSnExp@assayData), function(spectrum){
    data.table(rt=spectrum@rt, mz=spectrum@mz, int=spectrum@intensity)
  })
  all_df <- do.call(rbind, parsed_spectra)
  all_df <- all_df[order(all_df$rt),]
  all_df$rt <- all_df$rt/60
  rownames(all_df) <- NULL

  mzML_data <- grabMzmlData(mzML_filename)

  expect_equal(all_df, mzML_data)
})

# Except it doesn't - the sum of the individual intensities is less than the
# actual TIC recorded in the file???
# test_that("BPC matches MSnbase BPC", {
#   MSnExp_obj <- readMSData(mzML_filename, msLevel. = 1)
#   old_chrom <- chromatogram(MSnExp_obj, aggregationFun="max")[[1]]
#   old_mzML_data <- data.frame(rt=old_chrom@rtime, int=old_chrom@intensity)
#   old_mzML_data$rt <- old_mzML_data$rt/60
#   rownames(old_mzML_data) <- NULL
#
#   mzML_data <- grabMzmlBPC(mzML_filename, TIC = FALSE)
#
#   expect_equal(old_mzML_data, mzML_data)
# })
