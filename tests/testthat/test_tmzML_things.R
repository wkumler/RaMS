
output_dir <- paste0(tempdir(), "/temp_tmzMLs")
dir.create(output_dir)

# Writing tmzMLs ----
output_filename <- paste(output_dir, basename(mzML_filenames[1]), sep = "/")
tmzml_filename <- gsub(x = output_filename, "\\.mzML.*", ".tmzML")
output_filenames <- paste(output_dir, basename(mzML_filenames[1:4]), sep = "/")
tmzml_filenames <- gsub(x = output_filenames, "\\.mzML.*", ".tmzML")

mapply(tmzmlMaker, mzML_filenames[1:4], tmzml_filenames[1:4], verbosity=0)

test_that("tmzML conversion warns if no spectra found", {
  expect_warning(expect_error(tmzmlMaker(mzML_filenames[6])))
})

# test_that("tmzML conversion works for mzMLs", {
#   expect_identical(tmzmlMaker(
#     input_filename = mzML_filenames[1],
#     output_filename = tmzml_filename,
#     verbosity = 0), "temp_tmzMLs/DDApos_2.tmzML")
# })
#
# test_that("tmzML conversion works for mzXMLs", {
#   expect_identical(tmzmlMaker(
#     input_filename = mzML_filenames[1],
#     output_filename = tmzml_filename,
#     verbosity = 0), "temp_tmzMLs/DDApos_2.tmzML")
# })

test_that("tmzML warns if filename doesn't end in .tmzML", {
  expect_warning(tmzmlMaker(
    input_filename = mzML_filenames[1],
    output_filename = output_filename,
    verbosity = 0))
})

test_that("tmzmlMaker stops with multiple files", {
  expect_error(tmzmlMaker(mzML_filenames, tmzml_filenames))
})


# Reading tmzMLs ----
msdata <- grabMSdata(tmzml_filename)
test_that("tmzML can be read", {
  expect_s3_class(msdata, "msdata_connection")
  expect_type(msdata, "list")
})

test_that("tmzML prints expectedly (mostly for coverage)", {
  expect_message(print(msdata))
})

test_that("Requesting nonexistent dt throws error", {
  expect_error(msdata$files)
  expect_error(msdata$grab_what)
  expect_error(msdata$blah)
})

test_that("Requesting connection returns expected", {
  expect_type(msdata$connection, "list")
  expect_named(msdata$connection, c("files", "grab_what", "verbosity"))
})

test_that("Verbosity flags work", {
  expect_silent(grabMSdata(tmzml_filename)$MS1[mz%between%pmppm(118.0865)])
  expect_silent(grabMSdata(tmzml_filename, verbosity = 0)$MS1[mz%between%pmppm(118.0865)])
  expect_output(grabMSdata(tmzml_filename, verbosity = 1)$MS1[mz%between%pmppm(118.0865)])
  expect_output(grabMSdata(tmzml_filename, verbosity = 2)$MS1[mz%between%pmppm(118.0865)])
  expect_silent(grabMSdata(tmzml_filenames[1:4], verbosity = 0)$MS1[mz%between%pmppm(118.0865)])
  expect_output(grabMSdata(tmzml_filenames[1:4])$MS1[mz%between%pmppm(118.0865)])
})

test_that("tmzML files not found throw error", {
  expect_error(grabMSdata("blah.tmzML"))
  expect_error(grabMSdata(c("blah.tmzML", "blah2.tmzML")))
})

test_that("mz(X)ML and tmzML mixing throws error", {
  expect_error(grabMSdata(c(mzML_filenames, tmzml_filename)))
})

test_that("Additional args throw error with tmzMLs", {
  expect_warning(grabMSdata(tmzml_filename, mz=118.0865))
  expect_warning(grabMSdata(tmzml_filename, ppm=5))
  expect_warning(grabMSdata(tmzml_filename, rtrange = c(0, 10)))
  expect_warning(grabMSdata(tmzml_filename, prefilter = 1))
})

test_that("MS2 data also works", {
  expect_length(grabMSdata(tmzml_filename)$MS2[premz%between%pmppm(118.0865)], 6)
})


test_that("Same MS1 data after transposing", {
  init_data <- grabMSdata(mzML_filenames[1], verbosity = 0)$MS1[mz%between%pmppm(118.0865)]
  trans_data <- grabMSdata(tmzml_filenames[1], verbosity = 0)$MS1[mz%between%pmppm(118.0865)]
  expect_identical(init_data$rt, trans_data$rt)
  expect_identical(init_data$mz, trans_data$mz)
  expect_identical(init_data$int, trans_data$int)
})

test_that("Same MS2 data after transposing", {
  init_data <- grabMSdata(mzML_filenames[1], verbosity = 0)$MS2[premz%between%pmppm(118.0865)]
  trans_data <- grabMSdata(tmzml_filenames[1], verbosity = 0)$MS2[premz%between%pmppm(118.0865)]
  expect_identical(init_data$rt, trans_data$rt)
  expect_identical(init_data$premz, trans_data$premz)
  expect_identical(init_data$fragmz, trans_data$fragmz)
  expect_equal(init_data$voltage, trans_data$voltage) # because one's integer and one's double
  expect_identical(init_data$int, trans_data$int)
  # Can't just test equality of dt as a whole bc columns are in different orders
})

test_that("Same MS1 data after transposing, mass #2", {
  init_data <- grabMSdata(mzML_filenames[1], verbosity = 0)$MS1[mz%between%pmppm(90.0555)]
  trans_data <- grabMSdata(tmzml_filenames[1], verbosity = 0)$MS1[mz%between%pmppm(90.0555)]
  expect_identical(init_data$rt, trans_data$rt)
  expect_identical(init_data$mz, trans_data$mz)
  expect_identical(init_data$int, trans_data$int)
})

test_that("Same MS2 data after transposing, mass #2", {
  init_data <- grabMSdata(mzML_filenames[1], verbosity = 0)$MS2[premz%between%pmppm(90.0555)]
  trans_data <- grabMSdata(tmzml_filenames[1], verbosity = 0)$MS2[premz%between%pmppm(90.0555)]
  expect_identical(init_data$rt, trans_data$rt)
  expect_identical(init_data$premz, trans_data$premz)
  expect_identical(init_data$fragmz, trans_data$fragmz)
  expect_equal(init_data$voltage, trans_data$voltage) # because one's integer and one's double
  expect_identical(init_data$int, trans_data$int)
  # Can't just test equality of dt as a whole bc columns are in different orders
})


unlink(output_dir, recursive = TRUE, force = TRUE)
