# TO-DO:

# minifyMSdata ----
#' Shrink MS data by including only data points near masses of interest
#'
#' MS files can be annoyingly large if only a few masses are of interest. This large size makes it
#' difficult to share them online for debugging purposes and often means that untargeted algorithms
#' spend a lot of time picking peaks in data that's irrelevant. minifyMSdata is a function designed to
#' "minify" MS files by extracting only those data points that are within the ppm error of an m/z value
#' of interest, and returns the file essentially otherwise unchanged.
#'
#' @param files The name of a single file to be minified, usually produced by Proteowizard's `msconvert`
#' or something similar.
#' @param output_files The name of the file to be written out.
#' @param mz_exclude A vector of m/z values that should be excluded from the minified file. This argument
#' must be used with the `ppm` argument and should not be used with mz_include. For each mass provided, an
#' m/z window of +/- `ppm` is calculated, and all data points within that window are removed.
#' @param mz_include A vector of m/z values that should be included in the minified file. This argument
#' must be used with the `ppm` argument and should not be used with mz_exclude. For each mass provided, an
#' m/z window of +/- `ppm` is calculated, and all data points within that window are kept.
#' @param ppm The parts-per-million error of the instrument used to collect the original file.
#' @param warn Boolean. Should the function warn the user when removing an index from an mzML file?
#' @param prefilter A single number corresponding to the minimum intensity of
#'   interest in the MS1 data. Data points with intensities below this threshold
#'   will be silently dropped, which can dramatically reduce the size of the
#'   final object. Currently only works with MS1 data, but could be expanded
#'   easily to handle more.
#' @param verbosity A single number with a sensible default behavior. If larger
#' than 2, will render a progress bar as files are processed.
#'
#' @return Invisibly, the name of the new files.
#' @export
#'
#' @examples
#' \dontrun{
#' library(RaMS)
#' # Extract data corresponding to only valine and homarine
#' # m/z = 118.0865 and 138.0555, respectively
#' filename <- system.file("extdata", "LB12HL_AB.mzML.gz", package = "RaMS")
#' output_filename <- "mini_LB12HL_AB.mzML"
#' include_mzs <- c(118.0865, 138.0555)
#' minifyMSdata(filename, output_filename, mz_include=include_mzs, ppm=5)
#' init_data <- grabMSdata(filename)
#' mini_data <- grabMSdata(output_filename)
#' qplotMS1data(rbind(init_data$BPC, mini_data$BPC), color_col = "filename")
#' unlink(output_filename)
#'
#' # Exclude data corresponding to valine and homarine
#' filename <- system.file("extdata", "LB12HL_AB.mzML.gz", package = "RaMS")
#' output_filename <- "mini_LB12HL_AB.mzML"
#' exclude_mzs <- c(118.0865, 138.0555)
#' minifyMSdata(filename, output_filename, mz_exclude=exclude_mzs, ppm=5)
#' init_data <- grabMSdata(filename)
#' mini_data <- grabMSdata(output_filename)
#' qplotMS1data(rbind(init_data$BPC, mini_data$BPC), color_col = "filename")
#' unlink(output_filename)
#' }
minifyMSdata <- function(files, output_files=NULL, mz_exclude=NULL,
                         mz_include=NULL, ppm=NULL, warn=TRUE,
                         prefilter=-1, verbosity=NULL){
  # Check that files were provided
  if(!length(files)>0)stop("No files provided")

  # Check that all files exist, stop if not
  files_exist <- sapply(files, file.exists)
  if(!all(files_exist)){
    if(sum(!files_exist)>1){
      message("Files not found:")
      message(paste0(files[!files_exist][1:3], collapse = "\n"))
      message("...and ", sum(!files_exist)-3, " more")
      stop("Stopping!")
    } else {
      message(paste("File not found:", files[!files_exist]))
      stop("Stopping!")
    }
  }

  # Check that ppm is not NULL
  if(is.null(ppm)){
    stop("`ppm` must be provided by the user, but seems to still be NULL.")
  }

  # Warn if originals are going to be overwritten
  if(is.null(output_files)){
    output_files <- files
  }
  if(all(sort(files)==sort(output_files))){
    warning(paste("Output files would overwrite originals:",
                  "adding '_mini' to each output filename"))
    output_files <- gsub("\\.(?=mz[X]?ML)", replacement = "_mini\\.",
                             output_files, perl = TRUE)
  }

  # Stop if input not same length as output
  if(length(files)!=length(output_files)){
    stop("`output_files` must be the same length as `files`")
  }

  # Check that output files are same type as input files
  mzmls <- grepl("mzml", basename(files), ignore.case = TRUE)
  mzxmls <- grepl("mzxml", basename(files), ignore.case = TRUE)
  out_mzmls <- grepl("mzml", basename(output_files), ignore.case = TRUE)
  out_mzxmls <- grepl("mzxml", basename(output_files), ignore.case = TRUE)
  if(!all(mzmls==out_mzmls)|!all(mzxmls==out_mzxmls)){
    stop("minifyMSdata cannot convert from mzML to mzXML or vice-versa")
  }

  # Provide some output if requested
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
    if(grepl("mzML", basename(files[i]), ignore.case = TRUE)){
      minifyMzml(files[i], output_filename = output_files[i],
                 mz_exclude = mz_exclude, mz_include = mz_include,
                 ppm = ppm, warn = warn, prefilter = prefilter)
    } else if(grepl("mzXML", basename(files[i]), ignore.case = TRUE)){
      minifyMzxml(files[i], output_filename = output_files[i],
                  mz_exclude = mz_exclude, mz_include = mz_include,
                  ppm = ppm, prefilter = prefilter)
    } else {
      stop(paste("Unable to determine file type for", files[i]))
    }
    if(verbosity>0 & length(files)>=2){
      setTxtProgressBar(pb, i)
    }
  }
  if(verbosity>0){
    if(length(files)>=2){
      close(pb)
    }
    time_total <- round(difftime(Sys.time(), start_time), digits = 2)
    cat("Total time:", time_total, units(time_total), "\n")
  }
  invisible(output_files)
}




# minifyMzml ----

#' Shrink mzML files by including only data points near masses of interest
#'
#' mzML files can be annoyingly large if only a few masses are of interest. This large size makes it
#' difficult to share them online for debugging purposes and often means that untargeted algorithms
#' spend a lot of time picking peaks in data that's irrelevant. minifyMzml is a function designed to
#' "minify" mzML files by extracting only those data points that are within a ppm error of an m/z value
#' of interest, and returns the file essentially otherwise unchanged. This function currently works
#' only on MS1 data, but is reasonably expandable if demand becomes evident.
#'
#' @param filename The name of a single file to be minified, usually produced by Proteowizard's `msconvert`
#' or something similar.
#' @param output_filename The name of the file to be written out.
#' @param mz_exclude A vector of m/z values that should be excluded from the minified file. This argument
#' must be used with the `ppm` argument and should not be used with mz_include. For each mass provided, an
#' m/z window of +/- `ppm` is calculated, and all data points within that window are removed.
#' @param mz_include A vector of m/z values that should be included in the minified file. This argument
#' must be used with the `ppm` argument and should not be used with mz_exclude. For each mass provided, an
#' m/z window of +/- `ppm` is calculated, and all data points within that window are kept.
#' @param ppm The parts-per-million error of the instrument used to collect the original file.
#' @param warn Boolean. Should the function warn the user when removing an index from an mzML file?
#' @param prefilter A single number corresponding to the minimum intensity of
#'   interest in the MS1 data. Data points with intensities below this threshold
#'   will be silently dropped, which can dramatically reduce the size of the
#'   final object. Currently only works with MS1 data, but could be expanded
#'   easily to handle more.
#'
#' @return Invisibly, the name of the new file.
#'
#' @examples
#' \dontrun{
#' library(RaMS)
#' # Extract data corresponding to only valine and homarine
#' # m/z = 118.0865 and 138.0555, respectively
#' filename <- system.file("extdata", "LB12HL_AB.mzML.gz", package = "RaMS")
#' output_filename <- "mini_LB12HL_AB.mzML"
#' include_mzs <- c(118.0865, 138.0555)
#' minifyMzml(filename, output_filename, mz_include=include_mzs, ppm=5)
#' unlink(output_filename)
#'
#' # Exclude data corresponding to valine and homarine
#' filename <- system.file("extdata", "LB12HL_AB.mzML.gz", package = "RaMS")
#' output_filename <- "mini_LB12HL_AB.mzML"
#' exclude_mzs <- c(118.0865, 138.0555)
#' minifyMzml(filename, output_filename, mz_exclude=exclude_mzs, ppm=5)
#' unlink(output_filename)
#' }
minifyMzml <- function(filename, output_filename, ppm,
                       mz_exclude=NULL, mz_include=NULL,
                       warn=TRUE, prefilter=-1){
  if(missing(ppm)){
    stop("Must provide ppm value!")
  }
  if(is.null(mz_exclude) & is.null(mz_include)){
    stop("Either `mz_include` or `mz_exclude` must not be NULL")
  }
  xml_data <- xml2::read_xml(filename)

  checkFileType(xml_data, "mzML")
  prefilter <- checkProvidedPrefilter(prefilter)
  file_metadata <- grabMzmlEncodingData(xml_data)

  # Check for indexed mzML and drop index if present, with warning
  if(xml2::xml_name(xml_data)=="indexedmzML"){
    if(warn){
      warning(paste0("mzML file ", basename(filename), " contains an index. ",
                     "I don't know how to recompile indices so it's ",
                     "getting dropped for the minified file."))
    }
    xml_data <- xml2::xml_new_root(xml2::xml_find_first(xml_data, "//d1:mzML"))
  }
  # Check for chromatogramList and drop if present, with warning
  chromlist_node <- xml_find_all(xml_data, "//d1:chromatogramList")
  if(length(chromlist_node)>0){
    if(warn){
      warning(paste0("mzML file ", basename(filename),
                     " has a compiled chromatogram. ",
                     "This will be outdated after minification so it's ",
                     "getting dropped."))
    }
    xml2::xml_remove(chromlist_node)
  }

  # MS1 things
  ### Find MS1 intensity and m/z nodes
  ms1_xpath <- paste0('//d1:spectrum[d1:cvParam[@name="ms level" and ',
                      '@value="1"]][d1:cvParam[@name="base peak intensity"]]')
  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)
  mz_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[1]/d1:binary'
  ms1_mz_nodes <- xml2::xml_find_all(ms1_nodes, mz_xpath)
  int_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[2]/d1:binary'
  ms1_int_nodes <- xml2::xml_find_all(ms1_nodes, int_xpath)

  ### Convert MS1 nodes into data.tables
  ms1_minified <- mapply(function(ms1_mz_node, ms1_int_node){
    mzs <- getEncoded(xml2::xml_text(ms1_mz_node),
                      compression_type = file_metadata$compression,
                      bin_precision = file_metadata$mz_precision,
                      endi_enc = file_metadata$endi_enc)
    ints <- getEncoded(xml2::xml_text(ms1_int_node),
                       compression_type = file_metadata$compression,
                       bin_precision = file_metadata$int_precision,
                       endi_enc = file_metadata$endi_enc)
    if(!is.null(mz_include)){
      whitelist_data <- lapply(mz_include, function(mz_i){
        range_i <- pmppm(mz_i, ppm)
        mz_idxs <- mzs>min(range_i)&mzs<max(range_i)
        ints <- ints[mz_idxs]
        mzs <- mzs[mz_idxs]
        cbind(mzs=mzs, ints=ints)
      })
      if(all(sapply(whitelist_data, length)==0)){
        output_mat <- cbind(mzs=numeric(0), ints=numeric(0))
      } else {
        output_mat <- do.call(what = "rbind", whitelist_data)
      }
    } else if(!is.null(mz_exclude)){
      iterated_mzs <- mzs
      iterated_ints <- ints
      for(mz_i in unique(mz_exclude)){
        range_i <- pmppm(mz_i, ppm)
        iterate_idxs <- iterated_mzs>min(range_i)&iterated_mzs<max(range_i)
        iterated_mzs <- iterated_mzs[!iterate_idxs]
        iterated_ints <- iterated_ints[!iterate_idxs]
      }
      output_mat <- cbind(mzs=iterated_mzs, ints=iterated_ints)
    } else {
      print("This should never show up!")
    }
    subfilter_idxs <- output_mat[,"ints"]<prefilter
    output_mat <- output_mat[!subfilter_idxs, ,drop=FALSE]
    if(nrow(output_mat)==0){
      recoded_mzs <- ""
      recoded_ints <- ""
      bpmz <- 0
      bpint <- 0
      ticur <- 0
      minmz <- 0
      maxmz <- 0

      arraylength <- 0
    } else {
      recoded_mzs <- giveEncoding(output_mat[,"mzs"],
                                  compression_type = file_metadata$compression,
                                  bin_precision = file_metadata$mz_precision,
                                  endi_enc = file_metadata$endi_enc)
      recoded_ints <- giveEncoding(output_mat[,"ints"],
                                   compression_type = file_metadata$compression,
                                   bin_precision = file_metadata$int_precision,
                                   endi_enc = file_metadata$endi_enc)
      bpmz <- unname(output_mat[which.max(output_mat[,"ints"]), "mzs"])
      bpint <- max(output_mat[,"ints"], na.rm = TRUE)
      ticur <- sum(output_mat[,"ints"], na.rm = TRUE)
      minmz <- min(output_mat[,"mzs"], na.rm = TRUE)
      maxmz <- max(output_mat[,"mzs"], na.rm = TRUE)

      arraylength <- nrow(output_mat)
    }
    data.frame(
      mzs=recoded_mzs, ints=recoded_ints,
      bpmz=bpmz, bpint=bpint,
      ticur=ticur, minmz=minmz, maxmz=maxmz,
      mz_enclength=nchar(recoded_mzs),
      int_enclength=nchar(recoded_ints),
      arraylength=arraylength
    )
  }, ms1_mz_nodes, ms1_int_nodes, SIMPLIFY = FALSE)
  ms1_minified <- do.call(what = "rbind", ms1_minified)

  ### Replace with new data
  xml2::xml_text(ms1_mz_nodes) <- ms1_minified$mzs
  xml2::xml_text(ms1_int_nodes) <- ms1_minified$ints
  xml2::xml_attr(ms1_nodes, "defaultArrayLength") <- ms1_minified$arraylength
  mz_enclength_xpath <- "d1:binaryDataArrayList/d1:binaryDataArray[d1:cvParam[@name='m/z array']]"
  mz_enclength_nodes <- xml2::xml_find_all(ms1_nodes, mz_enclength_xpath)
  xml2::xml_attr(mz_enclength_nodes, "encodedLength") <- ms1_minified$mz_enclength
  int_enclength_xpath <- "d1:binaryDataArrayList/d1:binaryDataArray[d1:cvParam[@name='intensity array']]"
  int_enclength_nodes <- xml2::xml_find_all(ms1_nodes, int_enclength_xpath)
  xml2::xml_attr(int_enclength_nodes, "encodedLength") <- ms1_minified$int_enclength

  ### Update MS1 scan metadata
  param_xpath <- c("base peak m/z", "base peak intensity", "total ion current",
                   "lowest observed m/z", "highest observed m/z")
  param_cols <- c("bpmz", "bpint", "ticur", "minmz", "maxmz")
  mapply(function(xpath, val){
    full_xpath <- paste0("d1:cvParam[@name='", xpath, "']")
    node <- xml2::xml_find_all(ms1_nodes, full_xpath)
    xml2::xml_attr(node, "value") <- ms1_minified[[val]]
  }, param_xpath, param_cols)


  # MS2 things
  ### Drop all nodes with a precursor m/z outside of bounds
  ms2_xpath <- paste0('//d1:spectrum[d1:cvParam[@name="ms level" and ',
                      '@value="2"]][d1:cvParam[@name="base peak intensity"]]')
  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)
  pre_xpath <- paste0('d1:precursorList/d1:precursor/d1:selectedIonList/',
                      'd1:selectedIon/d1:cvParam[@name="selected ion m/z"]')
  ms2_pre_nodes <- xml2::xml_find_all(ms2_nodes, pre_xpath)
  ms2_pre_vals <- as.numeric(xml2::xml_attr(ms2_pre_nodes, "value"))
  if(!is.null(mz_include)){
    ms2_subset <- unlist(lapply(mz_include, function(premz_i){
      mzrange <- pmppm(premz_i, ppm)
      which(data.table::between(ms2_pre_vals, mzrange[1], mzrange[2]))
    }))
    to_remove <- setdiff(seq_along(ms2_pre_vals), ms2_subset)
  } else {
    to_remove <- unlist(lapply(mz_exclude, function(premz_i){
      mzrange <- pmppm(premz_i, ppm)
      which(data.table::between(ms2_pre_vals, mzrange[1], mzrange[2]))
    }))
  }
  xml2::xml_remove(ms2_nodes[to_remove])
  ### Re-number spectra
  spectrum_nodes <- xml2::xml_find_all(xml_data, "//d1:spectrum")
  speclist_node <- xml2::xml_find_first(xml_data, "//d1:spectrumList")
  xml2::xml_attr(speclist_node, "count") <- length(spectrum_nodes)
  xml2::xml_attr(spectrum_nodes, "index") <- seq_along(spectrum_nodes)-1


  # MS3 things (same as above but with every other precursor)
  ### Drop all nodes with a precursor m/z outside of bounds
  ms3_xpath <- paste0('//d1:spectrum[d1:cvParam[@name="ms level" and ',
                      '@value="3"]][d1:cvParam[@name="base peak intensity"]]')
  ms3_nodes <- xml2::xml_find_all(xml_data, ms3_xpath)
  pre_xpath <- paste0('d1:precursorList/d1:precursor/d1:selectedIonList/',
                      'd1:selectedIon/d1:cvParam[@name="selected ion m/z"]')
  ms3_pre_nodes <- xml2::xml_find_all(ms3_nodes, pre_xpath)
  ms3_pre_vals <- as.numeric(xml2::xml_attr(ms3_pre_nodes, "value"))
  ms3_pre_mat <- matrix(ms3_pre_vals, ncol = 2, byrow = TRUE)
  if(!is.null(mz_include)){
    ms3_subset <- unlist(lapply(mz_include, function(premz_i){
      mzrange <- pmppm(premz_i, ppm)
      which(between(ms3_pre_mat[,1], mzrange[1], mzrange[2]) |
              between(ms3_pre_mat[,2], mzrange[1], mzrange[2]))
    }))
    to_remove <- setdiff(seq_len(nrow(ms3_pre_mat)), ms3_subset)
  } else {
    to_remove <- unlist(lapply(mz_exclude, function(premz_i){
      mzrange <- pmppm(premz_i, ppm)
      which(between(ms3_pre_mat[,1], mzrange[1], mzrange[2]) |
              between(ms3_pre_mat[,2], mzrange[1], mzrange[2]))
    }))
  }
  xml2::xml_remove(ms3_nodes[to_remove])
  ### Re-renumber spectra
  spectrum_nodes <- xml2::xml_find_all(xml_data, "//d1:spectrum")
  speclist_node <- xml2::xml_find_first(xml_data, "//d1:spectrumList")
  xml2::xml_attr(speclist_node, "count") <- length(spectrum_nodes)
  xml2::xml_attr(spectrum_nodes, "index") <- seq_along(spectrum_nodes)-1



  # Add note that RaMS was used to shrink the file
  proclist_node <- xml2::xml_find_all(xml_data, "//d1:dataProcessingList")
  xml2::xml_add_child(proclist_node, "dataProcessing", id="RaMS_R_package")
  proc_node <- xml2::xml_find_all(proclist_node, '//dataProcessing[@id="RaMS_R_package"]')
  xml2::xml_add_child(proc_node, "processingMethod", order=0, softwareRef="RaMS")
  meth_node <- xml2::xml_find_all(proclist_node, '//processingMethod[@order="0"]')
  if(!is.null(mz_include)){
    xml2::xml_add_child(meth_node, "userParam", cvRef="MS", accession="MS:1009000",
                        name="Minification by m/z whitelist",
                        value=paste0(mz_include, collapse = "; "))
  } else {
    xml2::xml_add_child(meth_node, "userParam", cvRef="MS", accession="MS:1009001",
                        name="Minification by m/z blacklist",
                        value=paste0(mz_exclude, collapse = "; "))
  }
  xml2::xml_add_child(meth_node, "userParam", cvRef="MS", accession="MS:1009002",
                      name="ppm error for minification", value=ppm)
  process_count <- as.numeric(xml2::xml_attr(proclist_node, "count"))+1
  xml2::xml_attr(proclist_node, "count") <- process_count
  softlist_node <- xml2::xml_find_all(xml_data, "//d1:softwareList")
  xml2::xml_add_child(softlist_node, "software", id="RaMS",
                      version=as.character(packageVersion("RaMS")))
  soft_node <- xml2::xml_find_all(softlist_node, 'software[@id="RaMS"]')
  xml2::xml_add_child(soft_node, "userParam", name="RaMS R package")
  software_count <- as.numeric(xml2::xml_attr(softlist_node, "count"))+1
  xml2::xml_attr(softlist_node, "count") <- software_count


  # And write out the new version
  xml2::write_xml(x = xml_data, file = output_filename)
  return(invisible(output_filename))
}



# minifyMzxml ----
#' Shrink mzXML files by including only data points near masses of interest
#'
#' mzXML files can be annoyingly large if only a few masses are of interest. This large size makes it
#' difficult to share them online for debugging purposes and often means that untargeted algorithms
#' spend a lot of time picking peaks in data that's irrelevant. minifyMzxml is a function designed to
#' "minify" mzXML files by extracting only those data points that are within a ppm error of an m/z value
#' of interest, and returns the file essentially otherwise unchanged. This function currently works
#' only on MS1 data, but is reasonably expandable if demand becomes evident.
#'
#' @param filename The name of a single file to be minified, usually produced by Proteowizard's `msconvert`
#' or something similar.
#' @param output_filename The name of the file to be written out.
#' @param mz_exclude A vector of m/z values that should be excluded from the minified file. This argument
#' must be used with the `ppm` argument and should not be used with mz_include. For each mass provided, an
#' m/z window of +/- `ppm` is calculated, and all data points within that window are removed.
#' @param mz_include A vector of m/z values that should be included in the minified file. This argument
#' must be used with the `ppm` argument and should not be used with mz_exclude. For each mass provided, an
#' m/z window of +/- `ppm` is calculated, and all data points within that window are kept.
#' @param ppm The parts-per-million error of the instrument used to collect the original file.
#' @param warn Boolean. Should the function warn the user when removing an index from an mzML file?
#' @param prefilter A single number corresponding to the minimum intensity of
#'   interest in the MS1 data. Data points with intensities below this threshold
#'   will be silently dropped, which can dramatically reduce the size of the
#'   final object. Currently only works with MS1 data, but could be expanded
#'   easily to handle more.
#'
#' @return Invisibly, the name of the new file.
#'
#' @examples
#' \dontrun{
#' library(RaMS)
#' # Extract data corresponding to only valine and homarine
#' # m/z = 118.0865 and 138.0555, respectively
#' filename <- system.file("extdata", "LB12HL_AB.mzXML.gz", package = "RaMS")
#' output_filename <- "mini_LB12HL_AB.mzXML"
#' include_mzs <- c(118.0865, 138.0555)
#' minifyMzxml(filename, output_filename, mz_include=include_mzs, ppm=5)
#' unlink(output_filename)
#'
#' # Exclude data corresponding to valine and homarine
#' filename <- system.file("extdata", "LB12HL_AB.mzXML.gz", package = "RaMS")
#' output_filename <- "mini_LB12HL_AB.mzXML"
#' exclude_mzs <- c(118.0865, 138.0555)
#' minifyMzxml(filename, output_filename, mz_exclude=exclude_mzs, ppm=5)
#' unlink(output_filename)
#' }
minifyMzxml <- function(filename, output_filename, ppm, mz_exclude=NULL,
                        mz_include=NULL, prefilter=-1, warn=TRUE){
  xml_data <- xml2::read_xml(filename)

  checkFileType(xml_data, "mzXML")
  prefilter <- checkProvidedPrefilter(prefilter)
  file_metadata <- grabMzxmlEncodingData(xml_data)

  # Find MS1 intensity and m/z nodes
  ms1_xpath <- '//d1:scan[@msLevel="1" and @peaksCount>0]'
  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)

  mz_int_vals <- grabMzxmlSpectraMzInt(ms1_nodes, file_metadata)

  # Convert MS1 nodes into data.tables
  ms1_minified <- lapply(mz_int_vals, function(mz_int_val){
    mzs <- mz_int_val[,1]
    ints <- mz_int_val[,2]

    if(!is.null(mz_include)){
      whitelist_data <- lapply(mz_include, function(mz_i){
        range_i <- pmppm(mz_i, ppm)
        mz_idxs <- mzs>min(range_i)&mzs<max(range_i)
        ints <- ints[mz_idxs]
        mzs <- mzs[mz_idxs]
        cbind(mzs=mzs, ints=ints)
      })
      if(all(sapply(whitelist_data, length)==0)){
        output_mat <- cbind(mzs=numeric(0), ints=numeric(0))
      } else {
        output_mat <- do.call(what = "rbind", whitelist_data)
      }
    } else if(!is.null(mz_exclude)){
      iterated_mzs <- mzs
      iterated_ints <- ints
      for(mz_i in unique(mz_exclude)){
        range_i <- pmppm(mz_i, ppm)
        iterate_idxs <- iterated_mzs>min(range_i)&iterated_mzs<max(range_i)
        iterated_mzs <- iterated_mzs[!iterate_idxs]
        iterated_ints <- iterated_ints[!iterate_idxs]
      }
      output_mat <- cbind(mzs=iterated_mzs, ints=iterated_ints)
      if(any(is.na(output_mat))){
        browser("Found you an NA:")
      }
    } else {
      stop("Either `mz_include` or `mz_exclude` must not be NULL")
    }
    subfilter_idxs <- output_mat[,"ints"]<prefilter
    output_mat <- output_mat[!subfilter_idxs, ,drop=FALSE]
    if(nrow(output_mat)==0){
      recoded_mzints <- ""
      bpmz <- 0
      bpint <- 0
      ticur <- 0
      minmz <- 0
      maxmz <- 0

      arraylength <- 0
    } else {
      interped_mzints <- as.numeric(t(output_mat[,c("mzs", "ints")]))
      recoded_mzints <- giveEncoding(interped_mzints,
                                     compression_type = file_metadata$compression,
                                     bin_precision = file_metadata$precision,
                                     endi_enc = file_metadata$endi_enc)
      bpmz <- unname(output_mat[which.max(output_mat[,"ints"]), "mzs"])
      bpint <- max(output_mat[,"ints"], na.rm = TRUE)
      ticur <- sum(output_mat[,"ints"], na.rm = TRUE)
      minmz <- min(output_mat[,"mzs"], na.rm = TRUE)
      maxmz <- max(output_mat[,"mzs"], na.rm = TRUE)

      arraylength <- nrow(output_mat)
    }
    data.frame(mzints=recoded_mzints,
               bpmz=bpmz, bpint=bpint,
               ticur=ticur, minmz=minmz, maxmz=maxmz,
               arraylength=arraylength)
  })
  ms1_minified <- do.call(what = "rbind", ms1_minified)

  peak_nodes <- xml2::xml_find_first(ms1_nodes, "d1:peaks")
  xml2::xml_text(peak_nodes) <- ms1_minified$mzints

  xml2::xml_attr(ms1_nodes, "peaksCount") <- ms1_minified$arraylength
  xml2::xml_attr(ms1_nodes, "lowMz") <- ms1_minified$minmz
  xml2::xml_attr(ms1_nodes, "highMz") <- ms1_minified$maxmz
  xml2::xml_attr(ms1_nodes, "basePeakMz") <- ms1_minified$bpmz
  xml2::xml_attr(ms1_nodes, "basePeakIntensity") <- ms1_minified$bpint
  xml2::xml_attr(ms1_nodes, "totIonCurrent") <- ms1_minified$ticur


  # MS2 things
  ms2_xpath <- '//d1:scan[@msLevel="2"]'
  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)
  ms2_pre_nodes <- xml2::xml_find_all(ms2_nodes, "d1:precursorMz")
  ms2_pre_vals <- as.numeric(xml2::xml_text(ms2_pre_nodes))
  if(!is.null(mz_include)){
    ms2_subset <- unlist(lapply(mz_include, function(premz_i){
      mzrange <- pmppm(premz_i, ppm)
      which(data.table::between(ms2_pre_vals, mzrange[1], mzrange[2]))
    }))
    to_remove <- setdiff(seq_along(ms2_pre_vals), ms2_subset)
  } else {
    to_remove <- unlist(lapply(mz_exclude, function(premz_i){
      mzrange <- pmppm(premz_i, ppm)
      which(data.table::between(ms2_pre_vals, mzrange[1], mzrange[2]))
    }))
  }
  xml2::xml_remove(ms2_nodes[to_remove])

  # MS3 things
  ms3_xpath <- '//d1:scan[@msLevel="3"]'
  ms3_nodes <- xml2::xml_find_all(xml_data, ms3_xpath)
  ms3_pre_nodes <- xml2::xml_find_all(ms3_nodes, "d1:precursorMz")
  ms3_pre_vals <- as.numeric(xml2::xml_text(ms3_pre_nodes))
  ms3_pre_mat <- matrix(ms3_pre_vals, ncol = 2, byrow = TRUE)
  if(!is.null(mz_include)){
    ms3_subset <- unlist(lapply(mz_include, function(premz_i){
      mzrange <- pmppm(premz_i, ppm)
      which(between(ms3_pre_mat[,1], mzrange[1], mzrange[2]) |
              between(ms3_pre_mat[,2], mzrange[1], mzrange[2]))
    }))
    to_remove <- setdiff(seq_along(ms3_pre_vals), ms3_subset)
  } else {
    to_remove <- unlist(lapply(mz_exclude, function(premz_i){
      mzrange <- pmppm(premz_i, ppm)
      which(between(ms3_pre_mat[,1], mzrange[1], mzrange[2]) |
              between(ms3_pre_mat[,2], mzrange[1], mzrange[2]))
    }))
  }
  xml2::xml_remove(ms3_nodes[to_remove])



  # Add note that RaMS was used to shrink the file
  proclist_node <- xml2::xml_find_all(xml_data, "//d1:dataProcessing")
  xml2::xml_add_sibling(proclist_node, "dataProcessing")
  proc_node <- xml2::xml_find_first(xml_data, "//dataProcessing")
  xml2::xml_add_child(proc_node, "software", type="minification", name="RaMS",
                      version=as.character(packageVersion("RaMS")))

  # Remove index because the bytes will be off
  index_node <- xml2::xml_find_first(xml_data, xpath = "/d1:mzXML/d1:indexs")
  if(is.na(xml2::xml_type(index_node))){
    if(warn){
      warning(paste0("mzXML file ", basename(filename), " contains an index. ",
                     "I don't know how to recompile indices so it's ",
                     "getting dropped for the minified file."))
    }
    xml2::xml_remove(index_node)
  }
  offset_node <- xml2::xml_find_first(xml_data, xpath = "/d1:mzXML/d1:indexOffset")
  xml2::xml_remove(offset_node)

  # And write out the new version
  xml2::write_xml(x = xml_data, file = output_filename)
  return(invisible(output_filename))
}


# Helper functions ----
#' Convert from compressed binary to R numeric vector
#'
#' @param mzint_nodes The XML nodes containing the compressed binary string
#' @param compression_type Compression type to be used by memDecompress
#' @param bin_precision The bit (?) precision used by readBin
#' @param endi_enc The byte order (?) of the string. For mzML this is always
#'   "little" but mzXML can also be "big"
#'
#' @return A numeric vector of m/z or intensity values
getEncoded <- function(mzint_nodes, compression_type, bin_precision, endi_enc){
  decoded_mzs <- base64enc::base64decode(mzint_nodes)
  decomp_mzs <- memDecompress(decoded_mzs, type = compression_type)
  readBin(decomp_mzs, what = "double", n=length(decomp_mzs)/bin_precision,
          size = bin_precision, endian = endi_enc)
}

#' Convert from R numeric vector to compressed binary
#'
#' @param mzint_vals A numeric vector of m/z or intensity values
#' @param compression_type Compression type to be used by memCompress
#' @param bin_precision The bit (?) precision used by writeBin
#' @param endi_enc The byte order (?) of the string. For mzML this is always
#'   "little" but mzXML can also be "big"
#'
#' @return A single base64-encoded string of compressed binary values
giveEncoding <- function(mzint_vals, compression_type, bin_precision, endi_enc){
  comp_ints <- writeBin(mzint_vals, raw(0), size = bin_precision,
                        endian = endi_enc)
  new_raw_ints <- memCompress(comp_ints, type=compression_type)
  base64enc::base64encode(new_raw_ints)
}
