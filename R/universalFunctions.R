
checkFileType <- function(xml_data, node_to_check){
  # Check for mzML node
  # Length works because external pointer has length 2
  if(!length(xml_find_first(xml_data, paste0("//d1:", node_to_check)))){
    stop(paste0("No ", node_to_check, " node found in this file"))
  }
}
