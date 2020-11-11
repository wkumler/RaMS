
#' @import xml2
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
#' mzML_filename <- system.file("extdata", "190715_Poo_TruePooFK180310_Full2.mzML", package = "RaMS")
#' grabMzmlBPC(mzML_filename)
#' grabMzmlBPC(mzML_filename, TIC=TRUE)
grabMzmlBPC <- function(filename, TIC=FALSE){
  mz_xml <- xml2::read_xml(filename)

  rt_nodes <- xml2::xml_find_all(mz_xml, '//d1:cvParam[@name="scan start time"]')
  rt_vals <- as.numeric(xml2::xml_attr(rt_nodes, "value"))

  int_xpath <- ifelse(TIC, "total ion current", "base peak intensity")
  int_xpath_full <- paste0('//d1:cvParam[@name="', int_xpath, '"]')
  int_nodes <- xml2::xml_find_all(mz_xml, xpath = int_xpath_full)
  int_vals <- as.numeric(xml2::xml_attr(int_nodes, "value"))
  return(data.frame(rt=rt_vals, int=int_vals))
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
#' mzML_filename <- system.file("extdata", "190715_Poo_TruePooFK180310_Full2.mzML", package = "RaMS")
#' grabMzmlData(mzML_filename)
grabMzmlData <- function(filename){
  xml_data <- xml2::read_xml(filename)

  file_metadata <- grabMzmlMetadata(xml_data)
  ms1_nodes <- xml2::xml_find_all(
    xml_data, '//d1:cvParam[@name="MS1 spectrum"]/parent::d1:spectrum'
  )

  rt_vals <- grabSpectraRt(ms1_nodes)

  mz_vals <- grabSpectraMz(ms1_nodes, file_metadata)

  int_vals <- grabSpectraInt(ms1_nodes, file_metadata)

  data.frame(rt=rep(rt_vals, sapply(mz_vals, length)),
             mz=unlist(mz_vals), int=unlist(int_vals))
}



#' Read an mzML file's MSn data into a data.frame
#'
#' @details This function reads an mzML file's MSn data into R's working memory. mzML files
#' are fundamentally XML documents, which allows rapid access to the data by
#' parsing the XML. The R package `xml2::` is used for this purpose here.
#' Retention time information can be read directly, while *m/z* and intensity
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
#' mzML_MS2 <- system.file("extdata", "190715_Poo_TruePooFK180310_DDApos50.mzML", package = "RaMS")
#' grabMzmlMS2(mzML_MS2)
grabMzmlMS2 <- function(filename){
  xml_data <- xml2::read_xml(filename)

  file_metadata <- grabMzmlMetadata(xml_data)
  msn_xpath <- '//d1:cvParam[@name="MSn spectrum"]/parent::d1:spectrum'
  msn_nodes <- xml2::xml_find_all(xml_data, msn_xpath)

  rt_vals <- grabSpectraRt(msn_nodes)
  premz_vals <- grabSpectraPremz(msn_nodes)
  mz_vals <- grabSpectraMz(msn_nodes, file_metadata)
  int_vals <- grabSpectraInt(msn_nodes, file_metadata)

  data.frame(rt=rep(rt_vals, sapply(mz_vals, length)),
             premz=rep(premz_vals, sapply(mz_vals, length)),
             fragmz=unlist(mz_vals), int=unlist(int_vals))
}



#' Helper function to extract mzML file metadata
#'
#' @param xml_data mzML data parsed by xml2
#'
#' @return A list of values used by other parsing functions, such as
#' compression, mz_precision, int_precision
#'
#' @export
grabMzmlMetadata <- function(xml_data){
  compr_xpath <- '//d1:cvParam[@accession="MS:1000574"]'
  compr_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, compr_xpath), "name")
  compr <- switch(compr_type,
                  `zlib compression`="gzip",
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
  premz_xpath <- 'd1:scanList/d1:scan/d1:userParam'
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
