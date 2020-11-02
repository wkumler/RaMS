#' @import xml2
#' @importFrom base64enc base64decode


#' Grab the BPC or TIC for a given file
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
#' mzML_filename <- system.file("extdata", "180205_Poo_TruePoo_Full2.mzML", package = "RaMS")
#' grabSingleFileBPC(mzML_filename)
#' grabSingleFileBPC(mzML_filename, TIC=TRUE)
grabSingleFileBPC <- function(filename, TIC=FALSE){
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
#' @import xml2
#' @importFrom base64enc base64decode
#'
#' @export
#'
#' @examples
#' mzML_filename <- system.file("extdata", "180205_Poo_TruePoo_Full2.mzML", package = "RaMS")
#' grabMzmlData(mzML_filename)
grabMzmlData <- function(filename){
  xml_data <- xml2::read_xml(filename)

  rt_xpath <- '//d1:cvParam[@name="scan start time"]'
  rt_nodes <- xml2::xml_find_all(xml_data, rt_xpath)
  rt_vals <- as.numeric(xml2::xml_attr(rt_nodes, "value"))

  compr_xpath <- '//d1:cvParam[@accession="MS:1000574"]'

  compr_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, compr_xpath), "name")
  compr <- switch(compr_type,
                  `zlib compression`="gzip",
                  `none`="none")

  mz_bit_xpath <- '//d1:cvParam[@accession="MS:1000523"]'
  mz_bit_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, mz_bit_xpath), "name")
  mz_precision <- sub(mz_bit_type, pattern = "-bit float", replacement = "")
  mz_precision <- as.numeric(mz_precision)/8

  mz_xpath <- paste0('//d1:spectrum/d1:binaryDataArrayList',
                     '/d1:binaryDataArray[1]/d1:binary')
  mz_vals <- xml2::xml_text(xml2::xml_find_all(xml_data, mz_xpath))
  mz_vals <- lapply(mz_vals, function(binary){
      decoded_binary <- base64enc::base64decode(binary)
      raw_binary <- as.raw(decoded_binary)
      decomp_binary <- memDecompress(raw_binary, type = compr)
      final_binary <- readBin(decomp_binary, what = "double",
                              n=length(decomp_binary)/mz_precision,
                              size = mz_precision)
    })

  int_bit_xpath <- '//d1:cvParam[@accession="MS:1000521"]'
  int_bit_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, int_bit_xpath), "name")
  int_precision <- sub(int_bit_type, pattern = "-bit float", replacement = "")
  int_precision <- as.numeric(int_precision)/8

  int_xpath <- paste0('//d1:spectrum/d1:binaryDataArrayList',
                      '/d1:binaryDataArray[2]/d1:binary')
  int_vals <- xml2::xml_text(xml2::xml_find_all(xml_data, int_xpath))
  int_vals <- lapply(int_vals, function(binary){
    decoded_binary <- base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = compr)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/int_precision,
                            size = int_precision)
  })

  data.frame(rt=rep(rt_vals, sapply(mz_vals, length)),
             mz=unlist(mz_vals), int=unlist(int_vals))
}


filename <- "G:\\Shared drives\\Ingalls Lab\\Collaborative_Projects\\MESO-SCOPE\\Falkor\\HILIC_pos\\190715_Poo_TruePooFK180310_Full2.mzXML"
grabMzXmlData <- function(filename){
  xml_data <- xml2::read_xml(filename)

  rt_xpath <- '//d1:scan'
  rt_nodes <- xml2::xml_find_all(xml_data, rt_xpath)
  rt_attrs <- xml2::xml_attr(rt_nodes, "retentionTime")
  rt_vals <- as.numeric(gsub("PT|S", "", rt_attrs))

  compr_xpath <- '//d1:peaks'
  compr_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, compr_xpath),
                               "compressionType")
  compr <- switch(compr_type,
                  `zlib compression`="gzip",
                  `none`="none")

  bit_xpath <- '//d1:peaks'
  bit_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, bit_xpath), "precision")
  precision <- as.numeric(bit_type)/8

  xpath <- '//d1:peaks'
  vals <- xml2::xml_text(xml2::xml_find_all(xml_data, xpath))
  vals <- lapply(vals, function(binary){
    decoded_binary <- base64enc::base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = compr)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)*2/precision,
                            size = precision)
  })

  int_bit_xpath <- '//d1:cvParam[@accession="MS:1000521"]'
  int_bit_type <- xml2::xml_attr(xml2::xml_find_first(xml_data, int_bit_xpath), "name")
  int_precision <- sub(int_bit_type, pattern = "-bit float", replacement = "")
  int_precision <- as.numeric(int_precision)/8

  int_xpath <- paste0('//d1:spectrum/d1:binaryDataArrayList',
                      '/d1:binaryDataArray[2]/d1:binary')
  int_vals <- xml2::xml_text(xml2::xml_find_all(xml_data, int_xpath))
  int_vals <- lapply(int_vals, function(binary){
    decoded_binary <- base64decode(binary)
    raw_binary <- as.raw(decoded_binary)
    decomp_binary <- memDecompress(raw_binary, type = compr)
    final_binary <- readBin(decomp_binary, what = "double",
                            n=length(decomp_binary)/int_precision,
                            size = int_precision)
  })

  data.frame(rt=rep(rt_vals, sapply(mz_vals, length)),
             mz=unlist(mz_vals), int=unlist(int_vals))
}
