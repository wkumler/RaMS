
#' @import xml2
#' @import data.table
#' @importFrom base64enc base64decode
NULL

#' Grab the BPC or TIC for a given mzML file
#'
#' @param filename The name of the mzML file to be read.
#' @param TIC Should the BPC or TIC be read? If TIC=TRUE,
#' the total ion current (TIC)
#' for each retention time is read in. Otherwise, the BPC
#' (base peak chromatogram) is read for each retention time, corresponding
#' to the *maximum* intensity for each scan.
#'
#' @return A data.frame object with columns for retention time (rt) in minutes and
#' and intensity (int) corresponding to the TIC or BPC for a given file.
#'
#' @export
#'
#' @examples
#' mzML_filename <- system.file("proteomics", "MS3TMT11.mzML", package = "msdata")
#' grabMzmlBPC(mzML_filename)
#' grabMzmlBPC(mzML_filename, TIC=TRUE)
grabMzmlBPC <- function(filename, TIC=FALSE){
  xml_data <- xml2::read_xml(filename)

  ms1_nodes <- xml2::xml_find_all(
    xml_data, '//d1:cvParam[@name="MS1 spectrum"]/parent::d1:spectrum'
  )

  rt_vals <- grabSpectraRt(ms1_nodes)

  int_xpath <- ifelse(TIC, "total ion current", "base peak intensity")
  int_xpath_full <- paste0('d1:cvParam[@name="', int_xpath, '"]')
  int_nodes <- xml2::xml_find_all(ms1_nodes, xpath = int_xpath_full)
  int_vals <- as.numeric(xml2::xml_attr(int_nodes, "value"))
  return(data.table(rt=rt_vals, int=int_vals))
}



#' Read an mzML file into a data.frame
#'
#' @details This function reads an mzML file into R's working memory. mzML files
#' are fundamentally XML documents, which allows rapid access to the data by
#' parsing the XML. The R package `xml2::` is used for this purpose here.
#' Retention time information can be read directly, while *m/z* and intensity
#' information must be decoded from binary.
#'
#' @param filename The name of the mzML file to be read.
#'
#' @return A data.frame object with columns for retention time (rt) in minutes,
#' m/z (mz), and intensity (int).
#'
#' @export
#'
#' @examples
#' mzML_filename <- system.file("proteomics", "MS3TMT11.mzML", package = "msdata")
#' grabMzmlData(mzML_filename)
grabMzmlData <- function(filename){
  xml_data <- xml2::read_xml(filename)

  file_metadata <- grabMzmlMetadata(xml_data)
  ms1_nodes <- xml2::xml_find_all(
    xml_data, '//d1:cvParam[@name="ms level"][@value="1"]/parent::d1:spectrum'
  )

  rt_vals <- grabSpectraRt(ms1_nodes)
  mz_vals <- grabSpectraMz(ms1_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms1_nodes, file_metadata)

  data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
             mz=unlist(mz_vals), int=unlist(int_vals))
}



#' Read an mzML file's MS2 data into a data.frame
#'
#' @details This function reads an mzML file's MS2 data into R's working memory. mzML files
#' are fundamentally XML documents, which allows rapid access to the data by
#' parsing the XML. The R package `xml2::` is used for this purpose here.
#' Retention time and precursor mass information can be read directly, while *m/z* and intensity
#' information must be decoded from binary.
#'
#' @param filename The name of the mzML file to be read.
#'
#' @return A data.frame object with columns for retention time (rt) in minutes,
#' precursor mass (premz), fragment m/z (fragmz), and intensity (int).
#'
#' @export
#'
#' @examples
#' sf <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz"
#' mzML_MS2_filename <- system.file("proteomics", sf, package = "msdata")
#' grabMzmlMS2(mzML_MS2_filename)
grabMzmlMS2 <- function(filename){
  xml_data <- xml2::read_xml(filename)

  file_metadata <- grabMzmlMetadata(xml_data)
  ms2_xpath <- '//d1:cvParam[@name="ms level"][@value="2"]/parent::d1:spectrum'
  ms2_nodes <- xml2::xml_find_all(xml_data, ms2_xpath)

  rt_vals <- grabSpectraRt(ms2_nodes)
  premz_vals <- grabSpectraPremz(ms2_nodes)
  mz_vals <- grabSpectraMz(ms2_nodes, file_metadata)
  int_vals <- grabSpectraInt(ms2_nodes, file_metadata)

  data.table(rt=rep(rt_vals, sapply(mz_vals, length)),
             premz=rep(premz_vals, sapply(mz_vals, length)),
             fragmz=unlist(mz_vals), int=unlist(int_vals))
}


#' Simple wrapper around grabMzmlData to extract a specific mass across multiple files
#'
#' @details This function reads an mzML file's *m/z*, retention time, and intensity
#' data into R using `grabMzmlData` and extracts a mass slice according to the user-
#' provided *m/z* and ppm. It is vectorized via `sapply` across all the files in
#' the `filenames` vector.
#'
#' @param mass The *m/z* of the compound of interest
#' @param ppm The instrument's parts-per-million accuracy
#' @param filenames A vector of complete paths to the mzML files to be read, i.e.
#' those produced by system.file() or list.files(full.names=TRUE)
#'
#' @return A data.frame object with columns for retention time (rt) in minutes,
#' *m/z*, intensity (int), and filename.
#'
#' @export
#'
#' @examples
#' msdata_1 <- system.file("proteomics", "MS3TMT10_01022016_32917-33481.mzML.gz", package = "msdata")
#' msdata_2 <- system.file("proteomics", "MS3TMT11.mzML", package = "msdata")
#' grabMzmlEIC(mass=444.22, ppm=50, filenames=c(msdata_1, msdata_2))
grabMzmlEIC <- function(mass, ppm, filenames){
  if(any(!grepl("\\.mzML|\\.mzML.gz$", filenames))){
    stop("This function is currently only has support for mzML files")
  }
  raw_eics <- sapply(filenames, function(filename){
    all_data <- grabMzmlData(filename)
    all_data[mz%between%pmppm(mass, ppm = ppm)]
  }, simplify = FALSE)
  clean_filenames <- sub(pattern = "\\.mzML.*$", replacement = "", basename(filenames))
  filenamed <- mapply(FUN = cbind, raw_eics, clean_filenames, SIMPLIFY = FALSE)
  out_eic <- do.call(what=rbind, filenamed)
  names(out_eic) <- c("rt", "mz", "int", "filename")
  return(out_eic)
}



#' Helper function to extract mzML file metadata
#'
#' @param xml_data mzML data as parsed by xml2
#'
#' @return A list of values used by other parsing functions, currently
#' compression, mz_precision, int_precision
#'
#' @export
grabMzmlMetadata <- function(xml_data){
  compr_xpath <- paste0('//d1:cvParam[@accession="MS:1000574"]|',
                        '//d1:cvParam[@accession="MS:1000576"]')
  compr_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, compr_xpath), "name")
  compr <- switch(compr_type,
                  `zlib compression`="gzip",
                  `no compression`="none",
                  `none`="none")

  mz_precision_xpath <- '//d1:cvParam[@accession="MS:1000523"]'
  mz_bit_type <- xml_attr(xml_find_first(xml_data, mz_precision_xpath), "name")
  mz_precision <- sub(mz_bit_type, pattern = "-bit float", replacement = "")
  mz_precision <- as.numeric(mz_precision)/8

  int_bit_xpath <- '//d1:cvParam[@accession="MS:1000521"]'
  int_bit_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, int_bit_xpath), "name")
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
    decoded_binary <- base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/file_metadata$int_precision,
                            size = file_metadata$int_precision)
  })
}
