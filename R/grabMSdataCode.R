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
#' Retention times are reported in minutes, and will be converted automatically
#' if they are encoded in seconds.
#'
#' @param files A character vector of filenames to read into R's memory. Both
#'   absolute and relative paths are acceptable.
#' @param grab_what What data should be read from the file? Options include
#'   "MS1" for data only from the first spectrometer, "MS2" for fragmentation
#'   data, "BPC" for rapid access to the base peak chromatogram, and "TIC" for
#'   rapid access to the total ion chromatogram. These options can be combined
#'   (i.e. `grab_data=c("MS1", "MS2", "BPC")`) or this argument can be set to
#'   "everything" to extract all of the above. Options "EIC" and "EIC_MS2" are
#'   useful when working with files whose total size exceeds working memory -
#'   they first extracts all relevant MS1 and MS2 data, then discard data
#'   outside of the mass range(s) calculated from the provided mz and ppm.
#' @param verbosity Three levels of processing output to the R console are
#'   available, with increasing verbosity corresponding to higher integers. A
#'   verbosity of zero means that no output will be produced, useful when
#'   wrapping within larger functions. A verbosity of 1 will produce a progress
#'   bar using base R's txtProgressBar function. A verbosity of 2 or higher will
#'   produce timing output for each individual file read in. The default, NULL,
#'   will select between 1 and 2 depending on the number of files being read: if
#'   a single file, verbosity is set to 2; if multiple files, verbosity is set
#'   to 1.
#' @param mz A vector of the mass-to-charge ratio for compounds of interest.
#'   Only used when combined with `grab_what = "EIC"` (see above). Multiple
#'   masses can be provided.
#' @param ppm A single number corresponding to the mass accuracy (in parts per
#'   million) of the instrument on which the data was collected. Only used when
#'   combined with `grab_what = "EIC"` (see above).
#' @param rtrange Only available when parsing mzML files. A vector of length 2
#'   containing an upper and lower bound on retention times of interest.
#'   Providing a range here can speed up load times (although not enormously, as
#'   the entire file must still be read) and reduce the final object's size.
#' @param prefilter A single number corresponding to the minimum intensity of
#'   interest. Data points with intensities below this threshold will be
#'   silently dropped, which can dramatically reduce the size of the final
#'   object.
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
#' library(RaMS)
#' # Extract MS1 data from a couple files
#' sample_dir <- system.file("extdata", package = "RaMS")
#' sample_files <- list.files(sample_dir, full.names=TRUE)
#' multifile_data <- grabMSdata(sample_files[2:4], grab_what="MS1")
#'
#' # "Stream" data from the internet (i.e. Metabolights)
#' \dontrun{
#' access_url <- "https://www.ebi.ac.uk/metabolights/MTBLS703/files"
#'
#' # URL below obtained by right-clicking site download button and copying
#' # link address
#' sample_url <- paste0("https://www.ebi.ac.uk/metabolights/ws/studies/",
#'                      "MTBLS703/download/acefcd61-a634-4f35-9c3c-c572",
#'                      "ade5acf3?file=161024_Smp_LB12HL_AB_pos.mzXML")
#' file_data <- grabMSdata(sample_url, grab_what="everything", verbosity=2)
#' }
grabMSdata <- function(files, grab_what="everything", verbosity=NULL,
                       mz=NULL, ppm=NULL, rtrange=NULL, prefilter=-1){
  # Check that files were provided
  if(!length(files)>0)stop("No files provided")

  # Add sanity check for EIC extraction
  if(!is.null(mz) & !any(c("EIC", "EIC_MS2")%in%grab_what)){
    warning(paste0('Argument "mz" should be used with grab_what = "EIC" or',
            '"EIC_MS2" and will be ignored in the current call'))
  }
  if(!is.null(ppm) & !any(c("EIC", "EIC_MS2")%in%grab_what)){
    warning(paste0('Argument "mz" should be used with grab_what = "EIC" or',
                   '"EIC_MS2" and will be ignored in the current call'))
  }

  # Define outer control loop so multiple files can be read in simultaneously
  all_file_data <- list()
  if(is.null(verbosity)){
    verbosity <- ifelse(length(files)==1, 2, 1)
  }
  if(verbosity>0){
    if(length(files)>=2){
      pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    }
    start_time <- Sys.time()
  }
  for(i in seq_along(files)){
    filename <- files[i]

    if(grepl("mzML", basename(filename), ignore.case = TRUE)){
      out_data <- grabMzmlData(filename = filename, grab_what = grab_what,
                               verbosity = verbosity, mz = mz, ppm = ppm,
                               rtrange = rtrange, prefilter = prefilter)
    } else if(grepl("mzXML", basename(filename), ignore.case = TRUE)){
      out_data <- grabMzxmlData(filename = filename, grab_what = grab_what,
                                verbosity = verbosity, mz = mz, ppm = ppm,
                                rtrange = rtrange, prefilter = prefilter)
    } else {
      stop(paste("Unable to determine file type for", filename))
    }
    out_data_filenamed <- mapply(function(dt_i, fname_i){
      if(nrow(dt_i)){
        dt_i$filename <- fname_i
      } else {
        dt_i$filename <- character()
      }
      dt_i
    }, dt_i=out_data, fname_i=basename(filename), SIMPLIFY = FALSE)
    all_file_data[[i]] <- out_data_filenamed
    names(all_file_data)[[i]] <- basename(filename)
    if(verbosity>0 & length(files)>2){
      setTxtProgressBar(pb, i)
    }
  }
  if(verbosity>0){
    if(length(files)>2){
      close(pb)
    }
  }

  # Bind all the similar pieces together (e.g. stack MS1 from different files)
  all_file_data_output <- Reduce(function(x,y) Map(rbind, x, y, fill=TRUE), all_file_data)
  invisible(checkOutputQuality(
    output_data = all_file_data_output, grab_what = grab_what
  ))

  if(verbosity>0){
    time_total <- round(difftime(Sys.time(), start_time), digits = 2)
    cat("Total time:", time_total, units(time_total), "\n")
  }

  all_file_data_output
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

checkOutputQuality <- function(output_data, grab_what){
  if("everything"%in%grab_what){
    grab_what <- c("MS1", "MS2", "BPC", "TIC")
  }
  missing_data <- !grab_what%in%names(output_data)
  if(any(missing_data)){
    stop(paste("Not all data collected; missing",
                  paste(grab_what[missing_data], collapse = ", ")))
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

  misnamed <- mapply(function(nms, dt){
    if(nms=="BPC"|nms=="TIC"){
      proper_names <- c("rt", "int", "filename")
    } else if(nms=="MS1"){
      proper_names <- c("rt", "mz", "int", "filename")
    } else if(nms=="MS2"){
      proper_names <- c("rt", "premz", "fragmz", "int", "voltage", "filename")
    } else if (nms=="EIC"){
      proper_names <- c("rt", "mz", "int", "filename")
    } else if (nms=="EIC_MS2"){
      proper_names <- c("rt", "premz", "fragmz", "int", "voltage", "filename")
    } else if(nms=="metadata"){
      return(FALSE) # Because we don't know what names metadata may have
    } else {
      message(paste0("Unexpected data subset name in ",
                     "checkQualityOutput, may be malformed"))
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
    proper_class <- c(rt="numeric", mz="numeric", int="numeric",
                      filename="character", premz="numeric", fragmz="numeric",
                      voltage="integer")
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
checkFileType <- function(xml_data, node_to_check){
  # Check for mzML node
  # Length works because external pointer has length 2
  if(!length(xml2::xml_find_first(xml_data, paste0("//d1:", node_to_check)))){
    stop(paste0("No ", node_to_check, " node found in this file"))
  }
}

checkProvidedMzPpm <- function(mz, ppm){
  if(is.null(mz)){
    stop("Please provide an m/z value when using grab_what = EIC")
  }
  if(class(mz)!="numeric"&&class(mz)!="integer"){
    stop("Please provide a numeric m/z value")
  }
  if(any(is.na(mz))){
    stop("It looks like you have an NA among your m/z input")
  }
  if(any(mz<0)){
    stop("m/z must be positive")
  }

  if(is.null(ppm)){
    stop("Please provide a ppm value when using grab_what = EIC")
  }
  if(class(ppm)!="numeric"&&class(ppm)!="integer"){
    stop("Please provide a numeric ppm value")
  }
  if(ppm<0){
    stop("ppm must be positive")
  }
}


#' Plus/minus parts per million
#'
#' It shouldn't be hard to translate a point mass into a mass window bounded by
#' spectrometer accuracy.
#'
#' @param mass A length-1 numeric representing the mass of
#' interest for which a mass range is desired.
#' @param ppm The parts-per-million accuracy of the mass spectrometer on which
#' the data was collected.
#'
#' @return A length-2 numeric representing the mass range requested
#'
#' @export
#'
#' @examples
#' pmppm(100, 5)
#' pmppm(1000000, 5)
#' pmppm(118.0865, 2.5)
#' pmppm(892.535313, 10)
pmppm <- function(mass, ppm=4)c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))


timeReport <- function(last_time, text=NULL){
  time_total <- round(difftime(Sys.time(), last_time), digits = 2)
  cat(time_total, units(time_total), "\n")
  cat(text)
  Sys.time()
}

checkRTrange <- function(rtrange){
  if(!is.null(rtrange)){
    if("matrix"%in%class(rtrange)){
      rtrange <- as.vector(rtrange)
    }
    if(length(rtrange)!=2){
      stop("Please provide an rtrange of length 2")
    }
    if(class(rtrange)!="numeric"&&class(rtrange)!="integer"){
      stop("Please provide a numeric rtrange")
    }
  }
  rtrange
}


# Import area ----

#' @import utils
#' @import xml2
#' @import data.table
#' @importFrom base64enc base64decode
NULL
