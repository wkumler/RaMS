test_that("file can be read", {
  mzML_data <- grabMzmlMS2(mzML_filename)
  expect_s3_class(mzML_data, class = "data.frame")
  expect_identical(names(mzML_data), c("rt", "premz", "fragmz", "int"))
  expect_gt(nrow(mzML_data), 0)
})

test_that("new matches MSnbase", {
  MSnExp <- readMSData(mzML_filename, msLevel. = 2)
  parsed_spectra <- lapply(as.list(MSnExp@assayData), function(spectrum){
    data.frame(rt=spectrum@rt, premz=spectrum@precursorMz,
               fragmz=spectrum@mz, int=spectrum@intensity)
  })
  all_df <- do.call(rbind, parsed_spectra)
  all_df <- all_df[order(all_df$rt, all_df$premz),]
  all_df$rt <- all_df$rt/60
  rownames(all_df) <- NULL

  mzML_data <- grabMzmlMS2(mzML_filename)

  expect_equal(all_df, mzML_data)
})
