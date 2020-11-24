filename <- "180205_Poo_TruePoo_Full2.mzML"
filename <- "180205_Poo_TruePooPos_dda1.mzML"

grabMzmlData <- function(filename, grab_what, verbose=FALSE,
                         mz=NULL, ppm=NULL, rtrange=NULL){
  if(verbose){
    start_time <- Sys.time()
    last_time <- Sys.time()
    cat("Reading file... ")
  }
  xml_data <- read_xml(filename)

  checkMzmlType(xml_data)
  rtrange <- checkRTrange(rtrange)
  file_metadata <- grabMzmlEncodingData(xml_data)

  output_data <- list()

  if("everything"%in%grab_what){
    if(length(setdiff("everything", grab_what))){
      message("grab_what = `everything` includes MS1, MS2, BPC, and TIC data")
      message("Ignoring additional grabs")
    }
    grab_what <- c("MS1", "MS2", "BPC", "TIC")
  }

  if("MS1"%in%grab_what){
    if(verbose){
      cat(Sys.time()-last_time, "s\n")
      last_time <- Sys.time()
      cat("Reading MS1 data... ")
    }
    output_data$MS1 <- grabMzmlMS1(xml_data, rtrange, file_metadata)
  }

  if("MS2"%in%grab_what){
    if(verbose){
      cat(Sys.time()-last_time, "s\n")
      last_time <- Sys.time()
      cat("Reading MS2 data... ")
    }
    output_data$MS2 <- grabMzmlMS2(xml_data, rtrange, file_metadata)
  }

  if("BPC"%in%grab_what){
    if(verbose){
      cat(Sys.time()-last_time, "s\n")
      last_time <- Sys.time()
      cat("Reading BPC... ")
    }
    output_data$BPC <- grabMzmlBPC(xml_data, rtrange)
  }

  if("TIC"%in%grab_what){
    if(verbose){
      cat(Sys.time()-last_time, "s\n")
      last_time <- Sys.time()
      cat("Reading TIC... ")
    }
    output_data$BPC <- grabMzmlBPC(xml_data, rtrange, TIC = TRUE)
  }

  if("EIC"%in%grab_what){
    checkProvidedMzPpm(mz, ppm)

  }
  if(verbose){
    cat(Sys.time()-last_time, "s\n")
    cat("Total time:", Sys.time()-start_time, "\n")
  }

  output_data
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

checkProvidedMzPpm <- function(mz, ppm){
  if(is.null(mz)){
    stop("Please provide an m/z value when using grab_what = EIC")
  }
  if(class(mz)!="numeric"&&class(mz)!="integer"){
    stop("Please provide a numeric m/z value")
  }
  if(mz<0){
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

checkMzmlType <- function(xml_data){
  # Check for mzML node
  # Length works because external pointer has length 2
  if(!length(xml_find_first(xml_data, "//d1:mzML"))){
    stop("No mzML node found in this file")
  }

  # Check for absence of mzXML node
  if(length(xml_find_first(xml_data, "//d1:mzXML"))){
    stop("This file contains an mzXML node and shouldn't")
  }
}

#' Helper function to extract mzML file metadata
#'
#' @param xml_data mzML data as parsed by xml2
#'
#' @return A list of values used by other parsing functions, currently
#' compression, mz_precision, int_precision
#'
#' @export
grabMzmlEncodingData <- function(xml_data){
  init_node <- xml_find_first(xml_data, xpath = "//d1:spectrum")
  compr_xpath <- paste0('//d1:cvParam[@accession="MS:1000574"]|',
                        '//d1:cvParam[@accession="MS:1000576"]')
  compr_type <- xml2::xml_attr(xml2::xml_find_first(init_node, compr_xpath), "name")
  compr <- switch(compr_type,
                  `zlib compression`="gzip",
                  `no compression`="none",
                  `none`="none")

  mz_precision_xpath <- '//d1:cvParam[@accession="MS:1000523"]'
  mz_bit_type <- xml_attr(xml_find_first(init_node, mz_precision_xpath), "name")
  mz_precision <- sub(mz_bit_type, pattern = "-bit float", replacement = "")
  mz_precision <- as.numeric(mz_precision)/8

  int_bit_xpath <- '//d1:cvParam[@accession="MS:1000521"]'
  int_bit_type <- xml2::xml_attr(xml2::xml_find_first(init_node, int_bit_xpath), "name")
  int_precision <- sub(int_bit_type, pattern = "-bit float", replacement = "")
  int_precision <- as.numeric(int_precision)/8

  list(compression=compr, mz_precision=mz_precision, int_precision=int_precision)
}

grabSpectraRt <- function(xml_nodes){
  rt_xpath <- 'd1:scanList/d1:scan/d1:cvParam[@name="scan start time"]'
  rt_nodes <- xml2::xml_find_all(xml_nodes, rt_xpath)
  as.numeric(xml2::xml_attr(rt_nodes, "value"))
}

grabSpectraPremz <- function(xml_nodes){
  premz_xpath <- paste0('d1:precursorList/d1:precursor/d1:selectedIonList',
                        '/d1:selectedIon/d1:cvParam[@name="selected ion m/z"]')
  premz_nodes <- xml2::xml_find_all(xml_nodes, premz_xpath)
  as.numeric(xml2::xml_attr(premz_nodes, "value"))
}

grabSpectraVoltage <- function(xml_nodes){
  volt_xpath <- paste0('d1:precursorList/d1:precursor/d1:activation',
                       '/d1:cvParam[@name="collision energy"]')
  volt_nodes <- xml2::xml_find_all(xml_nodes, volt_xpath)
  as.numeric(xml2::xml_attr(volt_nodes, "value"))
}

grabSpectraMz <- function(xml_nodes, file_metadata){
  mz_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[1]/d1:binary'
  mz_vals <- xml2::xml_text(xml2::xml_find_all(xml_nodes, mz_xpath))
  lapply(mz_vals, function(binary){
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/file_metadata$mz_precision,
                            size = file_metadata$mz_precision)
  })
}

grabSpectraInt <- function(xml_nodes, file_metadata){
  int_xpath <- 'd1:binaryDataArrayList/d1:binaryDataArray[2]/d1:binary'
  int_vals <- xml2::xml_text(xml2::xml_find_all(xml_nodes, int_xpath))
  int_vals <- lapply(int_vals, function(binary){
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/file_metadata$int_precision,
                            size = file_metadata$int_precision)
  })
}

grabMzmlMS1 <- function(xml_data, rtrange, file_metadata){
  ms1_xpath <- '//d1:cvParam[@name="ms level" and @value="1"]/parent::d1:spectrum'
  ms1_nodes <- xml2::xml_find_all(xml_data, ms1_xpath)
  if(!is.null(rtrange)){
    ms1_nodes <- shrinkRTrange(ms1_nodes, rtrange)
  }

  rt_vals <- grabSpectraRt(ms1_nodes)
  mz_vals <- grabSpectraMz(ms1_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms1_nodes, file_metadata)

  data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
             mz=unlist(mz_vals), int=unlist(int_vals))

}

grabMzmlMS2 <- function(xml_data, rtrange, file_metadata){
  ms2_xpath <- '//d1:cvParam[@name="ms level" and @value="2"]/parent::d1:spectrum'
  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)
  if(!length(ms2_nodes)){
    return(data.table(rt=numeric(0), premz=numeric(0), fragmz=numeric(0),
                      int=numeric(0), voltages=numeric(0)))
  }
  if(!is.null(rtrange)){
    ms2_nodes <- shrinkRTrange(ms2_nodes, rtrange)
  }

  rt_vals <- grabSpectraRt(ms2_nodes)
  premz_vals <- grabSpectraPremz(ms2_nodes)
  voltages <- grabSpectraVoltage(ms2_nodes)
  mz_vals <- grabSpectraMz(ms2_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms2_nodes, file_metadata)

  data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
             premz=rep(premz_vals, sapply(mz_vals, length)),
             fragmz=unlist(mz_vals), int=unlist(int_vals),
             voltages=rep(premz_vals, sapply(mz_vals, length)))
}

grabMzmlBPC <- function(xml_data, rtrange, TIC=FALSE){
  ms1_xpath <- '//d1:cvParam[@name="ms level"][@value="1"]/parent::d1:spectrum'
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

shrinkRTrange <- function(xml_nodes, rtrange){
  rtrange_xpath <- paste0("d1:scanList/d1:scan/d1:cvParam[",
                          '@name="scan start time"',
                          " and @value>=", min(rtrange),
                          " and @value<=", max(rtrange), "]",
                          "/parent::d1:scan/parent::d1:scanList/parent::d1:spectrum")
  xml2::xml_find_all(xml_nodes, rtrange_xpath)
}
