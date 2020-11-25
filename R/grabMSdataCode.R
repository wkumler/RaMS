# Welcome to RaMS!

# grabMSdata ----

#' Grab mass-spectrometry data from file(s)
#'
#' The main `RaMS` function. This function accepts a list of the files that will
#' be read into R's working memory and returns a list of `data.table`s
#' containing the requested information. What information is requested is
#' determined by the `grab_what` argument, which can include MS1, MS2, BPC, TIC,
#' or metadata information. This function serves as a wrapper around both
#' `grabMzmlData` and `grabMzxmlData` and handles multiple files, but those two
#' have also been exposed to the user in case super-simple handling is desired.
#'
#' @param files A character vector of filenames to read into R's memory. Both
#'   absolute and relative paths are acceptable.
#' @param grab_what What data should be read from the file? Options include
#'   "MS1" for data only from the first spectrometer, "MS2" for fragmentation
#'   data, "BPC" for rapid access to the base peak chromatogram, and "TIC" for
#'   rapid access to the total ion chromatogram. These options can be combined
#'   (i.e. `grab_data=c("MS1", "MS2", "BPC")`) or this argument can be set to
#'   "everything" to extract all of the above. Option "EIC" is useful when
#'   working with files whose total size exceeds working memory - it first
#'   extracts all relevant MS1 and MS2 data, then discards data outside of the
#'   mass range(s) calculated from the provided mz and ppm.
#' @param verbosity Three levels of processing output to the R console: "very",
#'   which provides information about each file as it's read in; "kinda", for a
#'   progress bar but no individual file information; and "none" for no output
#'   of any kind.
#' @param mz A vector of the mass-to-charge ratio for compounds of interest.
#'   Only used when combined with `grab_what = "EIC"` (see above). Multiple
#'   masses can be provided.
#' @param ppm A single number corresponding to the mass accuracy (in parts per
#'   million) of the instrument on which the data was collected. Only used when
#'   combined with `grab_what = "EIC"` (see above).
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#'
#' @return A list of `data.table`s, each named after the arguments requested in
#'   grab_what. $MS1 contains MS1 information, MS2 contains fragmentation info,
#'   etc. MS1 data has four columns: retention time (rt), mass-to-charge (mz),
#'   intensity (int), and filename. MS2 data has six: retention time (rt),
#'   precursor m/z (premz), fragment m/z (fragmz), fragment intensity (int),
#'   collision energy (voltage), and filename. Data requested that does not
#'   exist in the provided files (such as MS2 data requested from MS1-only
#'   files) will return an empty (length zero) data.table. The data.tables
#'   extracted from each of the individual files are collected into one large
#'   table using data.table's `rbindlist`.
#'
#' @export
#'
#' @examples
grabMSdata <- function(files, grab_what=c("MS1", "MS2"), verbosity="very",
                       mz=NULL, ppm=NULL, rtrange=NULL){
  # Check file quality
  checkFiles(files)

  # Define outer control loop so multiple files can be read in simultaneously
  all_file_data <- list()
  if(verbosity=="very"|verbosity=="kinda"){
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  }
  for(i in seq_along(files)){
    filename <- files[i]

    verbose <- ifelse(verbosity=="very", TRUE, FALSE)
    if(grepl("mzML", filename, ignore.case = TRUE)){
      out_data <- grabMzmlData(filename, grab_what, verbose, mz, ppm, rtrange)
    } else if(grepl("mzXML", filename, ignore.case = TRUE)){
      out_data <- grabMzxmlData(filename, grab_what, verbose, mz, ppm, rtrange)
    } else {
      message(paste("Unable to determine file type for", filename))
      stopQuietly()
    }
    all_file_data[[i]] <- out_data
    if(verbosity=="very"|verbosity=="kinda"){
      setTxtProgressBar(pb, i)
    }
  }
  if(verbosity=="very"|verbosity=="kinda"){
    close(pb)
  }

  # Bind all the similar pieces together (e.g. stack MS1, MS2 from different files)
  all_file_data_output <- mapply(rbindlist, all_file_data)
  checkOutputQuality(all_file_data_output, grab_what)
}

# checkFiles ----

#' Check that filenames are acceptable.
#'
#' @param files The vector of files provided to grabMSdata.
#'
#' @return NULL (invisibly). The goal of this function is its side effects, i.e.
#'   throwing errors and providing info when the files are not found.
#'
#' @examples
checkFiles <- function(files){
  # Check that files were actually provided
  if(!length(files))stop("No files provided")
  # Check that files exist and stop with messages if not
  # I could be nicer here and just read in the files that exist but I think not
  files_exist <- file.exists(files)
  if(any(!files_exist)){
    message("Error: Unable to find all files:")
    if(sum(!files_exist)>6){
      # Avoid spamming output if lots of files missing
      message(paste("Couldn't find file", head(files[!files_exist]), "\n"))
      message(paste("... and", sum(!files_exist)-6, "others"))
    } else {
      message(paste("Couldn't find file", files[!files_exist]), "\n")
    }
    stopQuietly()
  }
  return(invisible(NULL))
}

# checkOutputQuality ----

#' Check that the output data is properly formatted.
#'
#' This function checks that data produced by repeated calls to the
#' `grabMzmlData()` and `grabMzxmlData()` functions is formatted properly before
#' it's provided to the user. It checks that all of the requested data has been
#' obtained and warns if data is found to be empty, misnamed, or has columns of
#' the wrong type.
#'
#' @param output_data The collected data resulting from repeated calls to
#'   `grabMzmlData()`, after being bound together.
#' @param grab_what The names of the data requested by the user.
#'
#' @return NULL (invisibly). The goal of this function is its side effects, i.e.
#'   throwing errors and providing info when the files are not found.
#'
#' @examples
checkOutputQuality <- function(output_data, grab_what){
  missing_data <- !grab_what%in%names(output_data)
  if(any(missing_data)){
    message(paste("Not all data collected; missing",
                  paste(grab_what[missing_data], collapse = ", ")))
    # Maybe convert above message into full stop?
    stopQuietly()
  }

  missing_data <- !as.logical(sapply(output_data, length))
  if(any(missing_data)){
    message(paste("Some bits seem to be empty:",
                  paste(names(output_data)[missing_data], collapse = ", ")))
  }

  zero_rows <- !as.logical(sapply(output_data, nrow))
  if(any(missing_data)){
    message(paste("Some data seem to be empty:",
                  paste(names(output_data)[missing_data], collapse = ", ")))
  }

  weird_classes <- !sapply(output_data, class)=="data.frame" # Change to data.table later
  if(any(weird_classes)){
    message(paste("Some data aren't data frames:",
                  paste(names(output_data)[weird_classes], collapse = ", ")))
  }

  misnamed <- mapply(function(nms, dt){
    if(nms=="BPC"|nms=="TIC"){
      proper_names <- c("rt", "int", "filename")
    } else if(nms=="MS1"){
      proper_names <- c("rt", "mz", "int", "filename")
    } else if(nms=="MS2"){
      proper_names <- c("rt", "premz", "fragmz", "int", "voltage", "filename")
    } else if(nms=="metadata"){
      stop("Update the proper names for metadata in checkQualityOutput")
    } else {
      message("Unexpected data subset name in checkQualityOutput, may be malformed")
      proper_names <- nms
    }
    !all(proper_names%in%names(dt))&all(names(dt)%in%proper_names)
  }, nms=names(output_data), dt=output_data)
  if(any(misnamed)){
    message(paste("Some data seem to have incorrect names:",
                  paste(names(output_data)[misnamed], collapse = ", ")))
  }

  wrong_coltype <- sapply(output_data, function(dt){
    col_classes <- sapply(dt, class)
    proper_class <- c(rt="numeric", mz="numeric", int="integer",
                      filename="character", premz="numeric", fragmz="numeric",
                      voltage="numeric")
    which(proper_class[names(col_classes)] != col_classes)
  })
  if(any(sapply(wrong_coltype, length))){
    message("Some data seem to have incorrect column types:")
    invisible(mapply(function(incorrect_coltypes, incorrect_names){
      if(length(incorrect_coltypes)){
        message(paste(incorrect_names, names(incorrect_coltypes), "\n"))
      }
    }, wrong_coltype, names(wrong_coltype)))
  }
  return(invisible(NULL))
}

# Other functions ----
stopQuietly <- function(){
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
