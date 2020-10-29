# Requires xml2, mzR, data.table
library(data.table)

grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- as.data.table(do.call(rbind, spectra_list))
  names(all_data) <- c("rt", "mz", "int")
  return(all_data)
}

grabBPC <- function(filename, TIC=FALSE){
  mz_xml <- xml2::read_xml(filename)

  rt_nodes <- xml2::xml_find_all(mz_xml, '//d1:cvParam[@name="scan start time"]')
  rt_vals <- as.numeric(xml2::xml_attr(rt_nodes, "value"))

  int_xpath <- ifelse(TIC, "total ion current", "base peak intensity")
  int_xpath_full <- paste0('//d1:cvParam[@name="', int_xpath, '"]')
  int_nodes <- xml2::xml_find_all(mz_xml, xpath = int_xpath_full)
  int_vals <- xml2::xml_attr(int_nodes, "value")
  return(data.table)
}

