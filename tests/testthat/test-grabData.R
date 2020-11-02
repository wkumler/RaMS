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
