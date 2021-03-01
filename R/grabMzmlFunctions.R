# Welcome to RaMS!

# grabMzmlData ----

#' Get mass-spectrometry data from an mzML file
#'
#' This function handles the mzML side of things, reading in files that are
#' written in the mzML format. Much of the code is similar to the mzXML format,
#' but the xpath handles are different and the mz/int array is encoded as two
#' separate entries rather than simultaneously. This function has been exposed
#' to the user in case per-file optimization (such as peakpicking or additional
#' filtering) is desired before the full data object is returned.
#'
#' @param filename A single filename to read into R's memory. Both absolute and
#'   relative paths are acceptable.
#' @param grab_what What data should be read from the file? Options include
#'   "MS1" for data only from the first spectrometer, "MS2" for fragmentation
#'   data, "BPC" for rapid access to the base peak chromatogram, and "TIC" for
#'   rapid access to the total ion chromatogram. These options can be combined
#'   (i.e. `grab_data=c("MS1", "MS2", "BPC")`) or this argument can be set to
#'   "everything" to extract all of the above. Option "EIC" is useful when
#'   working with files whose total size exceeds working memory - it first
#'   extracts all relevant MS1 and MS2 data, then discards data outside of the
#'   mass range(s) calculated from the provided mz and ppm.
#' @param verbose Boolean. If TRUE, R will print information about timing to the
#'   console as the file is read in.
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
#'   grab_what. $MS1 contains MS1 information, $MS2 contains fragmentation info,
#'   etc. MS1 data has three columns: retention time (rt), mass-to-charge (mz),
#'   and intensity (int). MS2 data has five: retention time (rt), precursor m/z
#'   (premz), fragment m/z (fragmz), fragment intensity (int), and collision
#'   energy (voltage). Data requested that does not exist in the provided files
#'   (such as MS2 data requested from MS1-only files) will return an empty
#'   (length zero) data.table.
#'
#' @export
#'
#' @examples
#' sample_file <- system.file("extdata", "LB12HL_AB.mzML.gz", package = "RaMS")
#' file_data <- grabMzmlData(sample_file, grab_what="MS1")
#' # Extract MS1 data and a base peak chromatogram
#' file_data <- grabMzmlData(sample_file, grab_what=c("MS1", "BPC"))
#' # Extract data from a retention time subset
#' file_data <- grabMzmlData(sample_file, grab_what=c("MS1", "BPC"),
#'                           rtrange=c(5, 7))
#' # Extract EIC for a specific mass
#' file_data <- grabMzmlData(sample_file, grab_what="EIC", mz=118.0865, ppm=5)
#' # Extract EIC for several masses simultaneously
#' file_data <- grabMzmlData(sample_file, grab_what="EIC", ppm=5,
#'                           mz=c(118.0865, 146.118104, 189.123918))
#'
#' # Extract MS2 data
#' sample_file <- system.file("extdata", "DDApos_2.mzML.gz", package = "RaMS")
#' MS2_data <- grabMzmlData(sample_file, grab_what="MS2")
grabMzmlData <- function(filename, grab_what, verbose=FALSE,
                         mz=NULL, ppm=NULL, rtrange=NULL){
  if(verbose){
    cat(paste0("\nReading file ", basename(filename), "... "))
    last_time <- Sys.time()
  }
  xml_data <- xml2::read_xml(filename)

  checkFileType(xml_data, "mzML")
  rtrange <- checkRTrange(rtrange)
  file_metadata <- grabMzmlEncodingData(xml_data)

  output_data <- list()

  if("everything"%in%grab_what){
    if(length(setdiff(grab_what, "everything"))&&verbose){
      message(paste("Heads-up: grab_what = `everything` includes",
                    "MS1, MS2, BPC, and TIC data"))
      message("Ignoring additional grabs")
    }
    grab_what <- c("MS1", "MS2", "BPC", "TIC", "metadata")
  }

  if("MS1"%in%grab_what){
    if(verbose)last_time <- timeReport(last_time, text = "Reading MS1 data...")
    output_data$MS1 <- grabMzmlMS1(xml_data = xml_data, rtrange = rtrange,
                                   file_metadata = file_metadata)
  }

  if("MS2"%in%grab_what){
    if(verbose)last_time <- timeReport(last_time, text = "Reading MS2 data...")
    output_data$MS2 <- grabMzmlMS2(xml_data = xml_data, rtrange = rtrange,
                                   file_metadata = file_metadata)
  }

  if("BPC"%in%grab_what){
    if(verbose)last_time <- timeReport(last_time, text = "Reading BPC...")
    output_data$BPC <- grabMzmlBPC(xml_data = xml_data, rtrange = rtrange)
  }

  if("TIC"%in%grab_what){
    if(verbose)last_time <- timeReport(last_time, text = "Reading TIC...")
    output_data$TIC <- grabMzmlBPC(xml_data = xml_data, rtrange = rtrange,
                                   TIC = TRUE)
  }

  if("EIC"%in%grab_what){
    checkProvidedMzPpm(mz, ppm)
    if(verbose){
      last_time <- timeReport(last_time, text = "Extracting EIC...")
    }
    if(!"MS1"%in%grab_what){
      init_dt <- grabMzmlMS1(xml_data = xml_data, rtrange = rtrange,
                             file_metadata = file_metadata)
    } else {
      init_dt <- output_data$MS1
      if(!nrow(init_dt))stop("Something weird - can't find MS1 data to subset")
    }
    EIC_list <- lapply(unique(mz), function(mass){
      init_dt[mz%between%pmppm(mass = mass, ppm = ppm)]
    })
    output_data$EIC <- rbindlist(EIC_list)
  }

  if("EIC_MS2"%in%grab_what){
    checkProvidedMzPpm(mz, ppm)
    if(verbose){
      last_time <- timeReport(last_time, text = "Extracting EIC MS2...")
    }
    if(!"MS2"%in%grab_what){
      init_dt <- grabMzmlMS2(xml_data = xml_data, rtrange = rtrange,
                             file_metadata = file_metadata)
    } else {
      init_dt <- output_data$MS2
    }
    premz <- NULL #To prevent R CMD check "notes"
    EIC_MS2_list <- lapply(unique(mz), function(mass){
      init_dt[premz%between%pmppm(mass = mass, ppm = ppm)]
    })
    output_data$EIC_MS2 <- rbindlist(EIC_MS2_list)
  }

  if("metadata"%in%grab_what){
    if(verbose){
      last_time <- timeReport(last_time, text = "Reading file metadata...")
    }
    output_data$metadata <- grabMzmlMetadata(xml_data = xml_data)
  }

  if(verbose){
    cat(Sys.time()-last_time, "s\n")
  }

  output_data
}

# Get mzML specifics (functions of xml_data) ----

#' Helper function to extract mzML file metadata
#'
#' @param xml_data mzML data as parsed by xml2
#'
#' @return A list of values corresponding to various pieces of metadata
#' for each file
grabMzmlMetadata <- function(xml_data){
  source_node <- xml2::xml_find_first(xml_data, xpath = "//d1:sourceFile")
  if(length(source_node)>0){
    source_file <- xml2::xml_attr(source_node, "name")
  } else {
    source_file <- "None found"
  }

  inst_xpath <- "//d1:referenceableParamGroup/d1:cvParam"
  inst_nodes <- xml2::xml_find_first(xml_data, xpath = inst_xpath)
  if(length(inst_nodes)>0){
    inst_val <- xml2::xml_attr(inst_nodes, "name")
  } else {
    inst_val <- "None found"
  }

  config_xpath <- "//d1:componentList/child::node()"
  config_nodes <- xml2::xml_find_all(xml_data, xpath = config_xpath)
  if(length(inst_nodes)>0){
    config_types <- xml2::xml_name(config_nodes)
    config_order <- xml2::xml_attr(config_nodes, "order")
    config_name_nodes <- xml2::xml_find_first(config_nodes, "d1:cvParam")
    config_names <- xml2::xml_attr(config_name_nodes, "name")
  } else {
    config_types <- "None found"
    config_order <- "None found"
    config_names <- "None found"
  }


  time_node <- xml2::xml_find_first(xml_data, xpath = "//d1:run")
  time_val <- xml2::xml_attr(time_node, "startTimeStamp")
  if(!is.na(time_val)){
    time_stamp <- as.POSIXct(strptime(time_val, "%Y-%m-%dT%H:%M:%SZ"))
  } else {
    time_stamp <- as.POSIXct(NA)
  }

  mslevel_xpath <- '//d1:spectrum/d1:cvParam[@name="ms level"]'
  mslevel_nodes <- xml2::xml_find_all(xml_data, xpath = mslevel_xpath)
  if(length(mslevel_nodes)>0){
    mslevels <- paste0("MS", unique(xml2::xml_attr(mslevel_nodes, "value")),
                      collapse = ", ")
  } else {
    mslevels <- "None found"
  }

  polarity_xpath <- '//d1:spectrum/d1:cvParam[@accession="MS:1000130"]'
  polarity_nodes <- xml2::xml_find_all(xml_data, polarity_xpath)
  if(length(polarity_nodes)>0){
    polarities <- unique(
      gsub(" scan", "", xml2::xml_attr(polarity_nodes, "name"))
    )
  } else {
    polarities <- "None found"
  }

  metadata <- data.table(
    source_file=source_file,
    inst_data=inst_val,
    config_data=list(data.frame(
      order=config_order,
      type=config_types,
      name=config_names
    )),
    timestamp = time_stamp,
    mslevels=paste0(mslevels, collapse = ", "),
    polarity=polarities
  )
}

#' Helper function to extract mzML file encoding data
#'
#' @param xml_data mzML data as parsed by xml2
#'
#' @return A list of values used by other parsing functions, currently
#' compression, mz_precision, int_precision
grabMzmlEncodingData <- function(xml_data){
  init_node <- xml2::xml_find_first(xml_data, xpath = "//d1:spectrum")
  compr_xpath <- paste0('//d1:cvParam[@accession="MS:1000574"]|',
                        '//d1:cvParam[@accession="MS:1000576"]')
  compr_node <- xml2::xml_find_first(init_node, compr_xpath)
  compr_type <- xml2::xml_attr(compr_node, "name")
  compr <- switch(compr_type,
                  `zlib`="gzip",
                  `zlib compression`="gzip",
                  `no compression`="none",
                  `none`="none")

  mz_precision_xpath <- '//d1:cvParam[@accession="MS:1000523"]'
  mz_bit_node <- xml2::xml_find_first(init_node, mz_precision_xpath)
  mz_bit_type <- xml2::xml_attr(mz_bit_node, "name")
  mz_precision <- sub(mz_bit_type, pattern = "-bit float", replacement = "")
  mz_precision <- as.numeric(mz_precision)/8

  int_bit_xpath <- '//d1:cvParam[@accession="MS:1000521"]'
  int_bit_node <- xml2::xml_find_first(init_node, int_bit_xpath)
  int_bit_type <- xml2::xml_attr(int_bit_node, "name")
  int_precision <- sub(int_bit_type, pattern = "-bit float", replacement = "")
  int_precision <- as.numeric(int_precision)/8

  if(is.na(int_precision))int_precision <- mz_precision
  if(is.na(mz_precision))mz_precision <- int_precision

  list(compression=compr, mz_precision=mz_precision,
       int_precision=int_precision)
}


#' Extract the MS1 data from an mzML nodeset
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzML file.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information.
#'
#' @return A `data.table` with columns for retention time (rt), m/z (mz), and
#'   intensity (int).
grabMzmlMS1 <- function(xml_data, rtrange, file_metadata){
  ms1_xpath <- paste0('//d1:spectrum[d1:cvParam[@name="ms level" and ',
                      '@value="1"]][d1:cvParam[@name="base peak intensity"]]')
  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)
  if(!is.null(rtrange)){
    ms1_nodes <- shrinkRTrange(ms1_nodes, rtrange)
  }
  if(!length(ms1_nodes)){
    return(data.table(rt=numeric(), mz=numeric(), int=numeric()))
  }


  rt_vals <- grabSpectraRt(ms1_nodes)
  mz_vals <- grabSpectraMz(ms1_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms1_nodes, file_metadata)

  data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
             mz=unlist(mz_vals), int=as.numeric(unlist(int_vals)))

}


#' Extract the MS2 data from an mzML nodeset
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzML file.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information.
#'
#' @return A `data.table` with columns for retention time (rt),  precursor m/z
#'   (mz), fragment m/z (fragmz), collision energy (voltage), and intensity
#'   (int).
grabMzmlMS2 <- function(xml_data, rtrange, file_metadata){
  ms2_xpath <- paste0('//d1:spectrum[d1:cvParam[@name="ms level" and ',
                      '@value="2"]][d1:cvParam[@name="base peak intensity"]]')

  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)
  if(!is.null(rtrange)){
    ms2_nodes <- shrinkRTrange(ms2_nodes, rtrange)
  }
  if(!length(ms2_nodes)){
    return(data.table(rt=numeric(), premz=numeric(), fragmz=numeric(),
                      int=numeric(), voltage=integer()))
  }

  rt_vals <- grabSpectraRt(ms2_nodes)
  premz_vals <- grabSpectraPremz(ms2_nodes)
  voltage <- grabSpectraVoltage(ms2_nodes)
  mz_vals <- grabSpectraMz(ms2_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms2_nodes, file_metadata)

  data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
             premz=rep(premz_vals, sapply(mz_vals, length)),
             fragmz=unlist(mz_vals), int=as.numeric(unlist(int_vals)),
             voltage=rep(voltage, sapply(mz_vals, length)))
}


#' Grab the BPC or TIC from a file
#'
#' The base peak intensity and total ion current are actually written into the
#' mzML files and aren't encoded, making retrieval of BPC and TIC information
#' blazingly fast if parsed correctly.
#'
#' @param xml_data An `xml2` nodeset, usually created by applying `read_xml` to
#'   an mzML file.
#' @param rtrange A vector of length 2 containing an upper and lower bound on
#'   retention times of interest. Providing a range here can speed up load times
#'   (although not enormously, as the entire file must still be read) and reduce
#'   the final object's size.
#' @param TIC Boolean. If TRUE, the TIC is extracted rather than the BPC.
#'
#' @return A `data.table` with columns for retention time (rt), and intensity
#'   (int).
grabMzmlBPC <- function(xml_data, rtrange, TIC=FALSE){
  ms1_xpath <- paste0('//d1:spectrum[d1:cvParam[@name="ms level" and ',
                      '@value="1"]][d1:cvParam[@name="base peak intensity"]]')

  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)
  if(!is.null(rtrange)){
    ms1_nodes <- shrinkRTrange(ms1_nodes, rtrange)
  }

  rt_vals <- grabSpectraRt(ms1_nodes)

  int_xpath <- ifelse(TIC, "total ion current", "base peak intensity")
  int_xpath_full <- paste0('d1:cvParam[@name="', int_xpath, '"]')
  int_nodes <- xml2::xml_find_all(ms1_nodes, xpath = int_xpath_full)
  int_vals <- as.numeric(xml2::xml_attr(int_nodes, "value"))
  return(data.table(rt=rt_vals, int=int_vals))
}



# Get spectrum things (functions of xml_nodes) ----

#' Extract the retention time from the spectra of an mzML nodeset
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#'
#' @return A numeric vector of retention times, one for each scan
grabSpectraRt <- function(xml_nodes){
  rt_xpath <- 'd1:scanList/d1:scan/d1:cvParam[@name="scan start time"]'
  rt_nodes <- xml2::xml_find_all(xml_nodes, rt_xpath)
  rt_vals <- as.numeric(xml2::xml_attr(rt_nodes, "value"))
  if(max(rt_vals)<150){
    # Guess RT is in minutes if the run is less than 150 units long
    # A 2.5 minute run is unheard of, and a 2.5 hour run is unheard of
    rt_vals <- rt_vals*60
  }
  rt_vals
}


#' Extract the precursor mass from the spectra of an mzML nodeset
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#'
#' @return A numeric vector of precursor masses, one for each scan
grabSpectraPremz <- function(xml_nodes){
  premz_xpath <- paste0('d1:precursorList/d1:precursor/d1:selectedIonList',
                        '/d1:selectedIon/d1:cvParam[@name="selected ion m/z"]')
  premz_nodes <- xml2::xml_find_all(xml_nodes, premz_xpath)
  as.numeric(xml2::xml_attr(premz_nodes, "value"))
}


#' Extract the collison energies from the spectra of an mzML nodeset
#'
#' Although the collision energy is typically fixed per file, it's equally fast
#' (afaik) to just grab them all individually here. Also, I'm worried about
#' these rumors of "ramped" collision energies
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#'
#' @return A numeric vector of collision energies, one for each scan.
grabSpectraVoltage <- function(xml_nodes){
  volt_xpath <- paste0('d1:precursorList/d1:precursor/d1:activation',
                       '/d1:cvParam[@name="collision energy"]')
  volt_nodes <- xml2::xml_find_all(xml_nodes, volt_xpath)
  as.integer(xml2::xml_attr(volt_nodes, "value"))
}


#' Extract the mass-to-charge data from the spectra of an mzML nodeset
#'
#' The mz and intensity information of mzML files are encoded as binary arrays,
#' sometimes compressed via gzip or zlib or numpress. This code finds all the
#' m/z binary arrays and converts them back to the original measurements. See
#' https://github.com/ProteoWizard/pwiz/issues/1301
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information. Here, the compression and
#'   mz precision information is relevant.
#'
#' @return A numeric vector of masses, many for each scan.
grabSpectraMz <- function(xml_nodes, file_metadata){
  mz_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[1]/d1:binary'
  mz_vals <- xml2::xml_text(xml2::xml_find_all(xml_nodes, mz_xpath))
  lapply(mz_vals, function(binary){
    if(!nchar(binary))return(numeric(0))
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/file_metadata$mz_precision,
                            size = file_metadata$mz_precision)
  })
}


#' Extract the intensity information from the spectra of an mzML nodeset
#'
#' The mz and intensity information of mzML files are encoded as binary arrays,
#' sometimes compressed via gzip or zlib or numpress. This code finds all the
#' intensity binary arrays and converts them back to the original measurements.
#' See https://github.com/ProteoWizard/pwiz/issues/1301
#'
#' @param xml_nodes An xml_nodeset object corresponding to the spectra collected
#'   by the mass spectrometer, usually produced by applying `xml_find_all` to an
#'   MS1 or MS2 nodeset.
#' @param file_metadata Information about the file used to decode the binary
#'   arrays containing m/z and intensity information. Here, the compression and
#'   int precision information is relevant.
#'
#' @return A numeric vector of intensities, many for each scan.
grabSpectraInt <- function(xml_nodes, file_metadata){
  int_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[2]/d1:binary'
  int_vals <- xml2::xml_text(xml2::xml_find_all(xml_nodes, int_xpath))
  int_vals <- lapply(int_vals, function(binary){
    if(!nchar(binary))return(numeric(0))
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/file_metadata$int_precision,
                            size = file_metadata$int_precision)
  })
}



# Other helper functions ----
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

shrinkRTrange <- function(xml_nodes, rtrange){
  # Xpath is magic. Basically we walk down from each spectrum node to find its
  # "scan start time" node and can filter directly on that.
  rtrange_xpath <- paste0("d1:scanList/d1:scan/d1:cvParam[",
                          '@name="scan start time"',
                          " and @value>=", min(rtrange),
                          " and @value<=", max(rtrange), "]",
                          "/parent::d1:scan/parent::d1:scanList/",
                          "parent::d1:spectrum")
  xml2::xml_find_all(xml_nodes, rtrange_xpath)
}
