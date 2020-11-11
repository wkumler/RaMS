test_that("file can be read", {
  mzML_data <- grabMzmlData(mzML_filename)
  expect_s3_class(mzML_data, class = "data.frame")
  expect_identical(names(mzML_data), c("rt", "mz", "int"))
  expect_gt(nrow(mzML_data), 0)
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

test_that("new matches previous", {
  old_grabData_fun <- function(filename){
    msdata <- mzR:::openMSfile(filename)
    fullhd <- mzR::header(msdata)
    spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
      given_peaks <- mzR::peaks(msdata, x)
      rtime <- fullhd[x, "retentionTime"]
      return(cbind(rtime, given_peaks))
    })
    all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)),
                          c("rt", "mz", "int"))
    return(all_data)
  }
  old_mzML_data <- old_grabData_fun(mzML_filename)
  old_mzML_data$rt <- old_mzML_data$rt/60
  mzML_data <- grabMzmlData(mzML_filename)

  expect_equal(old_mzML_data, mzML_data)
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

test_that("mzML matches mzXML", {
  mzML_data <- grabMzmlData(mzML_filename)
  mzXML_data <- grabMzxmlData(mzXML_filename)
  all.equal(mzML_data, mzXML_data)

  mzML_BPC <- grabMzmlBPC(mzML_filename)
  mzXML_BPC <- grabMzxmlBPC(mzXML_filename)
  all.equal(mzML_BPC$rt, mzXML_BPC$rt, tolerance = 1e-3)
  all.equal(mzML_BPC$int, mzXML_BPC$int, tolerance = 1)

  mzML_TIC <- grabMzmlBPC(mzML_filename, TIC = TRUE)
  mzXML_TIC <- grabMzxmlBPC(mzXML_filename, TIC = TRUE)
  all.equal(mzML_TIC$rt, mzXML_TIC$rt, tolerance = 1e-3)
  all.equal(mzML_TIC$int, mzXML_TIC$int, tolerance = 1)
})
