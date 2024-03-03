test_that("everything grabs everything", {
  expect_identical(mzML_everything, mzML_data)
})

test_that("everything got grabbed", {
  expect_setequal(names(mzML_data), c("MS1", "MS2", "metadata", "TIC", "BPC"))
})

test_that("mzML MS1 data looks normal", {
  # Check that there's data there at all
  expect_gt(nrow(mzML_data$MS1), 0)

  # Check that only expected columns exist
  expect_setequal(names(mzML_data$MS1), c("rt", "mz", "int", "filename"))

  # Check on column types
  expect_type(mzML_data$MS1$rt, "double")
  expect_type(mzML_data$MS1$mz, "double")
  expect_type(mzML_data$MS1$int, "double")
  expect_type(mzML_data$MS1$filename, "character")

  # Check that the expected files were loaded
  expect_length(unique(mzML_data$MS1$filename), 5)
  expect_identical(unique(mzML_data$MS1$filename), basename(mzML_filenames[1:5]))
})

test_that("mzML MS2 data looks normal", {
  # Check that there's data there at all
  expect_gt(nrow(mzML_data$MS2), 0)

  # Check that only expected columns exist
  expect_setequal(names(mzML_data$MS2), c("rt", "fragmz", "premz", "voltage",
                                          "int", "filename"))

  # Check on column types
  expect_type(mzML_data$MS2$rt, "double")
  expect_type(mzML_data$MS2$fragmz, "double")
  expect_type(mzML_data$MS2$premz, "double")
  expect_type(mzML_data$MS2$int, "double")
  expect_type(mzML_data$MS2$voltage, "integer")
  expect_setequal(unique(mzML_data$MS2$voltage), c(40L, NA_integer_))
  expect_type(mzML_data$MS2$filename, "character")

  # Check that the expected files were loaded
  expect_length(unique(mzML_data$MS2$filename), 2)
  expect_identical(unique(mzML_data$MS2$filename), basename(mzML_filenames[c(1,5)]))
})

test_that("mzML MS3 data looks normal", {
  ms3data <- grabMSdata(mzML_filenames[1], grab_what = "MS3")

  expect_gt(nrow(ms3data$MS3), 0)

  # Check that only expected columns exist
  expect_setequal(names(ms3data$MS3), c("rt", "fragmz", "prepremz", "premz",
                                        "voltage", "int", "filename"))

  # Check on column types
  expect_type(ms3data$MS3$rt, "double")
  expect_type(ms3data$MS3$fragmz, "double")
  expect_type(ms3data$MS3$premz, "double")
  expect_type(ms3data$MS3$int, "double")
  expect_type(ms3data$MS3$voltage, "integer")
  expect_identical(unique(ms3data$MS3$voltage), 60L)
  expect_type(ms3data$MS3$filename, "character")

  # Check that the expected files were loaded
  expect_length(unique(ms3data$MS3$filename), 1)
})

test_that("mzML BPC/TIC looks normal", {
  # Check that there's data there at all
  expect_gt(nrow(mzML_data$BPC), 0)

  # Check that only expected columns exist
  expect_setequal(names(mzML_data$BPC), c("rt", "int", "filename"))

  # Check column types
  expect_type(mzML_data$BPC$rt, "double")
  expect_type(mzML_data$BPC$int, "double")
  expect_type(mzML_data$BPC$filename, "character")



  # Check that there's data there at all
  expect_gt(nrow(mzML_data$TIC), 0)

  # Check that only expected columns exist
  expect_setequal(names(mzML_data$TIC), c("rt", "int", "filename"))

  # Check column types
  expect_type(mzML_data$TIC$rt, "double")
  expect_type(mzML_data$TIC$int, "double")
  expect_type(mzML_data$TIC$filename, "character")



  # Check that TIC <= BPC
  expect_true(all(mzML_data$BPC$int<=mzML_data$TIC$int))
})

test_that("multifile works", {
  expect_length(unique(mzML_data$MS1$filename), 5)
  expect_length(unique(mzML_data$TIC$filename), 5)
})

test_that("incl_polarity works on mzXMLs", {
  polarity_file <- grabMSdata(mzXML_filenames[3], incl_polarity = TRUE)
  expect_contains(names(polarity_file$MS1), "polarity")
  expect_contains(names(polarity_file$MS2), "polarity")
  expect_contains(names(polarity_file$BPC), "polarity")
  expect_contains(names(polarity_file$TIC), "polarity")

  expect_setequal(unique(polarity_file$MS1$polarity), c(1, -1))
  expect_setequal(unique(polarity_file$MS2$polarity), c(1, -1))
  expect_setequal(unique(polarity_file$BPC$polarity), c(1, -1))
  expect_setequal(unique(polarity_file$TIC$polarity), c(1, -1))

  ms3_file <- grabMSdata(mzXML_filenames[1], incl_polarity = TRUE,
                         grab_what = c("everything", "MS3"))
  expect_contains(names(ms3_file$MS3), "polarity")
  expect_equal(unique(ms3_file$MS3$polarity), 1)
})

test_that("incl_polarity works on mzMLs", {
  polarity_file <- grabMSdata(mzML_filenames[5], incl_polarity = TRUE)

  expect_contains(names(polarity_file$MS1), "polarity")
  expect_contains(names(polarity_file$MS2), "polarity")
  expect_contains(names(polarity_file$BPC), "polarity")
  expect_contains(names(polarity_file$TIC), "polarity")

  expect_setequal(unique(polarity_file$MS1$polarity), c(1, -1))
  expect_setequal(unique(polarity_file$MS2$polarity), c(1, -1))
  expect_setequal(unique(polarity_file$BPC$polarity), c(1, -1))
  expect_setequal(unique(polarity_file$TIC$polarity), c(1, -1))

  ms3_file <- grabMSdata(mzML_filenames[1], incl_polarity = TRUE,
                         grab_what = c("everything", "MS3"))
  expect_contains(names(ms3_file$MS3), "polarity")
  expect_equal(unique(ms3_file$MS3$polarity), 1)
})
