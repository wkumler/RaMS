
# mzML tests ----
filename <- mzML_filenames[2]
output_dir <- tempdir()
output_filename <- paste0(output_dir, "\\mini_mzML.mzML")


# Test blacklist
exclude_mzs <- c(118.0865, 138.0555)
minifyMzml(filename, output_filename, mz_blacklist=exclude_mzs, ppm=5, warn = FALSE)
mini_msdata <- grabMSdata(output_filename)


test_that("minified file is smaller", {
  expect_lt(nrow(mini_msdata$MS1),
            nrow(mzML_data$MS1[filename==basename(mzML_filenames[2])]))
})
test_that("blacklisted file has no data within bounds", {
  for(i in exclude_mzs){
    expect_equal(0, nrow(mini_msdata$MS1[mz%between%pmppm(i, ppm = 5)]))
  }
})
# # Don't want to deal with MSnbase dependencies, but these are good to run personally
# test_that("blacklist file opens in MSnbase", {
#   library(xcms)
#   raw_data <- readMSData(output_filename, pdata = NULL, msLevel. = 1)
#   bpis <- chromatogram(raw_data, aggregationFun = "max")
#   expect_s4_class(bpis, "MChromatograms")
#   plot(bpis)
# })




# Test whitelist
include_mzs <- c(118.0865, 138.0555)
minifyMzml(filename, output_filename, mz_whitelist=include_mzs, ppm=5, warn = FALSE)
mini_msdata <- grabMSdata(output_filename)

test_that("minified file is smaller", {
  expect_lt(nrow(mini_msdata$MS1),
            nrow(mzML_data$MS1[filename==basename(filename)]))
})
test_that("whitelisted file has data only within bounds", {
  expect_lte(max(mini_msdata$MS1$mz), max(pmppm(max(exclude_mzs))))
  expect_gte(min(mini_msdata$MS1$mz), min(pmppm(min(exclude_mzs))))

  mini_list <- lapply(include_mzs, function(mz_i){
    mini_msdata$MS1[mz%between%pmppm(mz_i, ppm = 5)]
  })
  eic_mini <- do.call(rbind, mini_list)
  eic_mini <- eic_mini[order(eic_mini$rt),]
  expect_equal(eic_mini, mini_msdata$MS1)

  full_list <- lapply(include_mzs, function(mz_i){
    eic_start <- mzML_data$MS1[mz%between%pmppm(mz_i, ppm = 5)]
    chosen_filename <- basename(filename)
    eic_start[filename==chosen_filename]
  })
  eic_full <- do.call(rbind, full_list)
  eic_full <- eic_full[order(eic_full$rt),]

  expect_equal(eic_full$rt, mini_msdata$MS1$rt)
  expect_equal(eic_full$mz, mini_msdata$MS1$mz)
  expect_equal(eic_full$int, mini_msdata$MS1$int)
})
# # Don't want to deal with MSnbase dependencies, but these are good to run personally
# test_that("whitelist file opens in MSnbase", {
#   library(xcms)
#   raw_data <- readMSData(output_filename, pdata = NULL, msLevel. = 1)
#   bpis <- chromatogram(raw_data, aggregationFun = "max")
#   expect_s4_class(bpis, "MChromatograms")
# })



unlink(output_filename)
unlink(output_dir)



# mzXML tests ----
filename <- mzXML_filenames[2]
output_dir <- tempdir()
output_filename <- paste0(output_dir, "\\mini_mzXML.mzXML")


# Test blacklist
exclude_mzs <- c(118.0865, 138.0555)
minifyMzxml(filename, output_filename, mz_blacklist=exclude_mzs, ppm=5)
mini_msdata <- grabMSdata(output_filename)


test_that("minified file is smaller", {
  expect_lt(nrow(mini_msdata$MS1),
            nrow(mzXML_data$MS1[filename==basename(mzXML_filenames[2])]))
})
test_that("blacklisted file has no data within bounds", {
  for(i in exclude_mzs){
    expect_equal(0, nrow(mini_msdata$MS1[mz%between%pmppm(i, ppm = 5)]))
  }
})
# # Don't want to deal with MSnbase dependencies, but these are good to run personally
# test_that("blacklist file opens in MSnbase", {
#   library(xcms)
#   raw_data <- MSnbase::readMSData(output_filename, pdata = NULL, msLevel. = 1)
#   bpis <- xcms::chromatogram(raw_data, aggregationFun = "max")
#   expect_s4_class(bpis, "MChromatograms")
# })




# Test whitelist
include_mzs <- c(118.0865, 138.0555)
minifyMzxml(filename, output_filename, mz_whitelist=include_mzs, ppm=5)
mini_msdata <- grabMSdata(output_filename)

test_that("minified file is smaller", {
  expect_lt(nrow(mini_msdata$MS1),
            nrow(mzXML_data$MS1[filename==basename(filename)]))
})
test_that("whitelisted file has data only within bounds", {
  expect_lte(max(mini_msdata$MS1$mz), max(pmppm(max(exclude_mzs))))
  expect_gte(min(mini_msdata$MS1$mz), min(pmppm(min(exclude_mzs))))

  mini_list <- lapply(include_mzs, function(mz_i){
    mini_msdata$MS1[mz%between%pmppm(mz_i, ppm = 5)]
  })
  eic_mini <- do.call(rbind, mini_list)
  eic_mini <- eic_mini[order(eic_mini$rt),]
  expect_equal(eic_mini, mini_msdata$MS1)

  full_list <- lapply(include_mzs, function(mz_i){
    eic_start <- mzXML_data$MS1[mz%between%pmppm(mz_i, ppm = 5)]
    chosen_filename <- basename(filename)
    eic_start[filename==chosen_filename]
  })
  eic_full <- do.call(rbind, full_list)
  eic_full <- eic_full[order(eic_full$rt),]

  expect_equal(eic_full$rt, mini_msdata$MS1$rt)
  expect_equal(eic_full$mz, mini_msdata$MS1$mz)
  expect_equal(eic_full$int, mini_msdata$MS1$int)
})
# # Don't want to deal with MSnbase dependencies, but these are good to run personally
# test_that("whitelist file opens in MSnbase", {
#   raw_data <- MSnbase::readMSData(output_filename, pdata = NULL, msLevel. = 1)
#   bpis <- xcms::chromatogram(raw_data, aggregationFun = "max")
#   expect_s4_class(bpis, "MChromatograms")
# })



unlink(output_filename)
unlink(output_dir)
