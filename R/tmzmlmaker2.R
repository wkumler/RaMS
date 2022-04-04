library(data.table)
library(RaMS)
library(tidyverse)


binwidth <- 1
msdata_init <- grabMSdata("inst/extdata/LB12HL_AB.mzML.gz")
mz_splits_all <- floor(msdata_init$MS1$mz/binwidth)
mz_splits <- sort(unique(mz_splits_all))
data_dict <- paste(mz_splits, seq_along(mz_splits)+2, sep = ", ", collapse = ", ")
split_MS1_mzs <- split(as.data.frame(msdata_init$MS1), mz_splits_all)
outlines <- lapply(split_MS1_mzs, function(x){
  paste(as.matrix(x[,c("rt", "mz", "int")]), sep = ", ", collapse = ", ")
})
doc_out <- paste("tmzML document", data_dict, paste(outlines, collapse = "\n"), sep = "\n")
cat(doc_out, file = "tmzml2_test.tmzML")
msdata_init$MS1[mz%between%pmppm(118.0865)] %>%
  ggplot() +
  geom_line(aes(x=rt, y=int))


mz_i <- pmppm(118.0865)
tmzml_data <- scan("tmzml2_test.tmzML", what = "character", sep = "\n", nmax = 2, quiet = TRUE)
split_data <- as.numeric(strsplit(x = tmzml_data[2], split=", ")[[1]])
mass_idxs <- data.table(matrix(split_data, ncol=2, byrow = TRUE))
needed_lines <- mass_idxs[max(which(V1<min(mz_i))):max(which(V1<max(mz_i)))]$V2
msdata_line <- scan("tmzml2_test.tmzML", what = "character", quiet = TRUE, sep = "\n",
                    skip = needed_lines-1, nlines = length(needed_lines))
msdata_mat <- matrix(as.numeric(unlist(strsplit(msdata_line, split = ", "))), ncol = 3)
msdata <- setNames(as.data.table(msdata_mat), c("rt", "mz", "int"))
ggplot(msdata) +
  geom_line(aes(x=rt, y=int))
