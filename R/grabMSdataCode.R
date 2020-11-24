
grabMSdata <- function(files, grab_what=c("MS1", "MS2"),
                       verbosity=c("very", "kinda", "none"),
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
  setTxtProgressBar(pb, i){
    close(pb)
  }

  # Bind all the similar pieces together (e.g. stack MS1, MS2 from different files)
  all_file_data_output <- mapply(rbindlist, all_file_data)
  checkOutputQuality(all_file_data_output, grab_what)
}

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
}

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
}

stopQuietly <- function(){
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
