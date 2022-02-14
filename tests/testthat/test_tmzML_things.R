output_folder <- "temp_tmzMLs"
dir.create(output_folder)

# Writing tmzMLs ----
output_filename <- paste(output_folder, basename(mzML_filenames[1]), sep = "/")
tmzml_filename <- gsub(x = output_filename, "\\.mzML.*", ".tmzML")

test_that("tmzML conversion works for mzMLs", {
  expect_identical(RaMS:::tmzmlMaker(
    input_filename = mzML_filenames[1],
    output_filename = tmzml_filename,
    verbosity = 0), "temp_tmzMLs/DDApos_2.tmzML")
})

test_that("tmzML conversion works for mzXMLs", {
  expect_identical(RaMS:::tmzmlMaker(
    input_filename = mzML_filenames[1],
    output_filename = tmzml_filename,
    verbosity = 0), "temp_tmzMLs/DDApos_2.tmzML")
})




# Reading tmzMLs ----

unlink(output_folder, recursive = TRUE, force = TRUE)
