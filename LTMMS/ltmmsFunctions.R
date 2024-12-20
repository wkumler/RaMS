
library(XML)
library(DBI)
library(duckdb)
library(data.table)

startElemParser <- function(name, attrs, sax_env){
  if(name == "cvParam"){
    if(attrs["name"] == "ms level"){
      sax_env$scan_ms_level <- as.numeric(attrs["value"])
    }
    if(attrs["name"] == "scan start time"){
      sax_env$scan_rt <- as.numeric(attrs["value"])
    }
    if(attrs["name"] == "isolation window target m/z"){
      sax_env$scan_premz <- as.numeric(attrs["value"])
    }
    if(attrs["accession"] %in% c("MS:1000523", "MS:1000521")){
      sax_env$bin_type <- ifelse(attrs["accession"] == "MS:1000523", "m/z", "int")
      bin_prec <- sub("-bit float", "", attrs["name"])
      sax_env$bin_prec <- as.numeric(bin_prec) / 8
    }
    if(attrs["accession"] %in% c("MS:1000574", "MS:1000576")){
      sax_env$bin_compr <- switch(attrs["name"], `zlib` = "gzip", `zlib compression` = "gzip", `no compression` = "none", `none` = "none")
    }
    if(attrs["name"]=="collision energy"){
      sax_env$scan_voltage <- as.numeric(attrs["value"])
    }
  }
  if(name == "binary"){
    sax_env$parsing_binary <- TRUE
    sax_env$binary_bits <- character()
  }
  if(name == "spectrum"){
    sax_env$scan_idx <- as.numeric(attrs["index"])
    # print(paste0("Starting scan #", current_scan))
    if(current_scan %% 1000 == 0){
      print(paste0("Starting scan #", sax_env$scan_idx))
    }
  }
}
textElemParser <- function(content, sax_env){
  if(sax_env$parsing_binary){
    sax_env$binary_bits <- c(sax_env$binary_bits, content)
  }
}
endElemParser <- function(name, attrs, sax_env){
  if(name == "binary"){
    sax_env$parsing_binary <- FALSE
    bin_vals <- paste(sax_env$binary_bits, collapse = "")
    bin_vals <- base64enc::base64decode(bin_vals)
    if(length(bin_vals)==0){
      return(NULL)
    }
    bin_vals <- memDecompress(as.raw(bin_vals), type = sax_env$bin_compr)
    bin_vals <- readBin(bin_vals, what = "double", size = sax_env$bin_prec, n = length(bin_vals) / sax_env$bin_prec)
    if(sax_env$bin_type == "m/z"){
      sax_env$scan_mzs <- bin_vals
    } else {
      sax_env$scan_ints <- bin_vals
    }
  }
  if(name == "spectrum"){
    if(sax_env$scan_ms_level == 1){
      scan_data <- data.table(scan_idx = sax_env$scan_idx, rt = sax_env$scan_rt,
                              mz = sax_env$scan_mzs, int = sax_env$scan_ints)
      sax_env$MS1_scan_data <- c(sax_env$MS1_scan_data, list(scan_data))
    }
    if(sax_env$scan_ms_level == 2){
      scan_data <- data.table(scan_idx = sax_env$scan_idx, rt = sax_env$scan_rt,
                              premz = sax_env$scan_premz, fragmz = sax_env$scan_mzs,
                              int = sax_env$scan_ints,
                              voltage = sax_env$scan_voltage)
      sax_env$MS2_scan_data <- c(sax_env$MS2_scan_data, list(scan_data))
    }
    if((length(sax_env$MS1_scan_data) + length(sax_env$MS2_scan_data)) > 10000){
      print("Writing data to duckdb")
      new_MS1_data <- rbindlist(sax_env$MS1_scan_data)
      new_MS1_data[,filename:=sax_env$filename]
      new_MS2_data <- rbindlist(sax_env$MS2_scan_data)
      new_MS2_data[,filename:=sax_env$filename]
      dbWriteTable(sax_env$msduck, "MS1", new_MS1_data, append = TRUE)
      dbWriteTable(sax_env$msduck, "MS2", new_MS2_data, append = TRUE)
      sax_env$MS1_scan_data <- list()
      sax_env$MS2_scan_data <- list()
    }
  }
  if(name == "mzML"){
    print("Writing data to duckdb")
    new_MS1_data <- rbindlist(sax_env$MS1_scan_data)
    new_MS1_data[,filename:=sax_env$filename]
    new_MS2_data <- rbindlist(sax_env$MS2_scan_data)
    new_MS2_data[,filename:=sax_env$filename]
    dbWriteTable(sax_env$msduck, "MS1", new_MS1_data, append = TRUE)
    dbWriteTable(sax_env$msduck, "MS2", new_MS2_data, append = TRUE)
    sax_env$MS1_scan_data <- list()
    sax_env$MS2_scan_data <- list()
  }
}
convertMzmlToDuckdb <- function(ms_files, output_duckdb_name, overwrite_ok = TRUE){
  sax_env <- new.env()

  sax_env$msduck <- dbConnect(duckdb::duckdb(), output_duckdb_name)
  empty_MS1 <- data.table(filename=character(), rt=numeric(), mz=numeric(), int=numeric())
  dbWriteTable(sax_env$msduck, "MS1", empty_MS1, overwrite=overwrite_ok)
  empty_MS2 <- data.table(filename=character(), rt=numeric(), premz=numeric(),
                          fragmz=numeric(), int=numeric(), voltage=numeric())
  dbWriteTable(sax_env$msduck, "MS2", empty_MS2, overwrite=overwrite_ok)
  on.exit(dbDisconnect(sax_env$msduck))

  sax_env$MS1_scan_data <- list()
  sax_env$MS2_scan_data <- list()
  sax_env$parsing_binary <- FALSE

  sax_handlers <- list(
    startElement = function(name, attrs) startElemParser(name, attrs, sax_env),
    text = function(content) textElemParser(content, sax_env),
    endElement = function(name, attrs) endElemParser(name, attrs, sax_env)
  )

  v <- sapply(ms_files, function(ms_file_i){
    sax_env$filename <- basename(ms_file_i)
    xmlEventParse(ms_file_i, handlers = sax_handlers)
  })
  invisible(v)
}



# convertMzmlToDuckdb(ms_files = "LTMMS/HILIC_PT_Fullstrength_FullMS_ddMS2_CE40_Pos.mzML", output_duckdb_name = "LTMMS/mini_demo.duckdb")
#
# msduck <- dbConnect(duckdb::duckdb(), "LTMMS/mini_demo.duckdb")
# ms1_data <- as.data.table(dbGetQuery(msduck, "SELECT * FROM MS1"))
# ms2_data <- as.data.table(dbGetQuery(msduck, "SELECT * FROM MS2"))
#
# print(ms1_data)
# print(ms2_data)





# library(tidyverse)
# library(patchwork)
# library(RaMS)
# msduck <- dbConnect(duckdb::duckdb(), "LTMMS/mini_demo.duckdb")
# ms1_data <- as.data.table(dbGetQuery(msduck, "SELECT * FROM MS1 WHERE mz > 132.1011 AND mz < 132.1021"))
# ms2_data <- as.data.table(dbGetQuery(msduck, "SELECT * FROM MS2 WHERE premz > 132.1011 AND premz < 132.1021"))
# ggplot(ms1_data) +
#   geom_point(aes(x=rt, y=mz, color=log10(int))) +
#   scale_color_viridis_c()
# ms1_chrom <- qplotMS1data(ms1_data[rt%between%c(7, 11)]) + ggtitle("Leucine/Ile/TMAP")
# ms2_data[rt%between%c(7, 11)][fragmz<200] %>%
#   arrange(desc(int)) %>%
#   mutate(mz_group=mz_group(fragmz, ppm=20, max_groups = 20)) %>%
#   ggplot() +
#   geom_point(aes(x=fragmz, y=log10(int), color=factor(mz_group)))
# ms2_chroms <- ms2_data[rt%between%c(7, 11)] %>%
#   arrange(desc(int)) %>%
#   mutate(mz_group=mz_group(fragmz, ppm=20)) %>%
#   filter(mz_group<20) %>%
#   mutate(mz_group=round(mean(fragmz), 5), .by=mz_group) %>%
#   mutate(mz_group=fct_inorder(as.character(mz_group))) %>%
#   summarise(int=max(int), .by = c(rt, mz_group)) %>%
#   ggplot() +
#   geom_line(aes(x=rt, y=int)) +
#   facet_wrap(~mz_group, scales = "free_y")
# ms1_chrom/ms2_chroms
