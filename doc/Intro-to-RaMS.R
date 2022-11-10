## ---- include = FALSE---------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
options(tidyverse.quiet = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "80%",
  fig.align = 'center',
  fig.height = 3,
  fig.width = 6.5
)

## ----findfiles, message=FALSE-------------------------------------------------
library(RaMS)

# Locate the file directory
msdata_dir <- system.file("extdata", package = "RaMS")

# Identify the files of interest
data_files <- list.files(msdata_dir, pattern = "mzML", full.names = TRUE)[1:4]

# Check that the files identified are the ones expected
basename(data_files)

## ----loadTIC------------------------------------------------------------------
single_file <- data_files[2]

msdata <- grabMSdata(single_file, grab_what = "TIC")

knitr::kable(head(msdata$TIC, 3))

## ----headerTIC----------------------------------------------------------------
plot(msdata$TIC$rt, msdata$TIC$int, type = "l")

## ----loadBPC, fig.height=3----------------------------------------------------
msdata <- grabMSdata(single_file, grab_what = "BPC")

## ----plotBPC, warning=FALSE---------------------------------------------------
library(tidyverse)
ggplot(msdata$BPC) + geom_line(aes(x=rt, y=int))

## ----loadmultiBPC-------------------------------------------------------------
msdata <- grabMSdata(data_files[2:4], grab_what = "BPC")

ggplot(msdata$BPC) + geom_line(aes(x=rt, y=int, color=filename))

## ----ggplotdemo, dpi=144------------------------------------------------------
ggplot(msdata$BPC) + 
  geom_polygon(aes(x=rt, y=int, color=filename), lwd=1, fill="#FFFFFF44") +
  theme(legend.position = c(0.8, 0.7), plot.title = element_text(face="bold"),
        axis.title = element_text(size = 15)) +
  scale_y_continuous(labels = c(0, "250M", "500M"), breaks = c(0, 2.5e8, 5e8)) +
  scale_colour_manual(values = c("#2596be", "#6c25be", "#bea925")) +
  labs(x="Retention time (minutes)", y="Intensity", 
       title = "Base peak chromatogram", color="Files:") +
  coord_cartesian(xlim = c(7.50, 9), ylim = c(0, 5e8))

## ----metadatademo-------------------------------------------------------------
# Since the minification process strips some metadata, I use the 
# less-minified DDA file here
grabMSdata(files = data_files[1], grab_what = "metadata")

## ----MS1demo------------------------------------------------------------------
msdata <- grabMSdata(data_files[2:4], grab_what = "MS1")

knitr::kable(head(msdata$MS1, 3))

## ----adenineplot, warning=FALSE, message=FALSE--------------------------------
library(data.table)

adenine_mz <- 136.06232

adenine_data <- msdata$MS1[mz%between%pmppm(adenine_mz, ppm=5)]

ggplot(adenine_data) + geom_line(aes(x=rt, y=int, color=filename)) + lims(x=c(4, 9))

## ----multicmpdplot------------------------------------------------------------
mzs_of_interest <- c(adenine=136.06232, valine=118.0865, homarine=138.055503)

mass_data <- imap_dfr(mzs_of_interest, function(mz_i, name){
  cbind(msdata$MS1[mz%between%pmppm(mz_i, ppm=10)], name)
})

ggplot(mass_data) + 
  geom_line(aes(x=rt, y=int, color=filename)) + 
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  lims(x=c(4, 9))

## ----MS2data------------------------------------------------------------------
msdata <- grabMSdata(data_files[1], grab_what = "everything")
knitr::kable(head(msdata$MS2, 3))

## ----MS2demo------------------------------------------------------------------
homarine_mz <- 137.047678+1.007276

homarine_MS2 <- msdata$MS2[premz%between%pmppm(homarine_mz, 5)]
homarine_MS2$int <- homarine_MS2$int/max(homarine_MS2$int)

ggplot(homarine_MS2) +
  geom_point(aes(x=fragmz, y=int)) +
  geom_segment(aes(x=fragmz, xend=fragmz, y=int, yend=0)) +
  scale_y_continuous(breaks = c(0, .5, 1), labels = c("0%", "50%", "100%")) +
  labs(x="Fragment m/z", y="")

## ----plotly, eval=FALSE-------------------------------------------------------
#  ## Not run to save space in the vignette:
#  library(plotly)
#  msdata <- grabMSdata(data_files, grab_what = c("MS1", "MS2", "BPC"))
#  
#  compound_MS1 <- msdata$MS1 %>%
#    filter(mz%between%pmppm(homarine_mz, 10)) %>%
#    filter(!str_detect(filename, "DDA"))
#  
#  compound_MS2 <- msdata$MS2[premz%between%pmppm(homarine_mz, 10)] %>%
#    group_by(rt) %>%
#    arrange(desc(int)) %>%
#    summarise(frags=paste(
#      paste(round(fragmz, digits = 3), round(int), sep = ": "), collapse = "\n"),
#      .groups="drop"
#    ) %>%
#    mutate(int=approx(x = compound_MS1$rt, y=compound_MS1$int, xout = rt)$y)
#  plot_ly(compound_MS1) %>%
#    add_trace(type="scatter", mode="lines", x=~rt, y=~int, color=~filename,
#              hoverinfo="none") %>%
#    add_trace(type="scatter", mode="markers", x=~rt, y=~int,
#              text=~frags, hoverinfo="text", showlegend=FALSE,
#              marker=list(color="black"), data = compound_MS2) %>%
#    layout(annotations=list(x=min(compound_MS2$rt), y=median(compound_MS2$int)*10,
#                            text="Mouse over to see\nMSMS fragments"))

## ----homarine_fragsearch------------------------------------------------------
homarine_frag_mz <- 94.06567
msdata$MS2[fragmz%between%pmppm(homarine_frag_mz, ppm = 5)] %>%
  head() %>% knitr::kable()

## ----homarine_neutsearch------------------------------------------------------
homarine_neutral_loss <- homarine_mz - homarine_frag_mz

msdata$MS2 <- mutate(msdata$MS2, neutral_loss=premz-fragmz) %>%
  select("rt", "premz", "fragmz", "neutral_loss", "int", "voltage", "filename")
msdata$MS2[neutral_loss%between%pmppm(homarine_neutral_loss, ppm = 5)] %>%
  arrange(desc(int)) %>% head() %>% knitr::kable()

## ---- fig.height=3------------------------------------------------------------
chrom_file <- system.file("extdata", "wk_chrom.mzML.gz", package = "RaMS")
msdata_chroms <- grabMSdata(chrom_file, verbosity = 0, grab_what = "chroms")
given_chrom <- msdata_chroms$chroms[chrom_type=="SRM iletter1"]
ptitle <- with(given_chrom, paste0(
  unique(chrom_type), ": Target m/z = ", unique(target_mz), "; Product m/z = ", 
  unique(product_mz)
))
plot(given_chrom$rt, given_chrom$int, type="l", main=ptitle)

## ----EICdemo------------------------------------------------------------------
all_data <- grabMSdata(data_files, grab_what = c("MS1", "MS2"))

mzs_of_interest <- c(adenine=136.06232, valine=118.0865)
small_data <- grabMSdata(data_files, grab_what = c("EIC", "EIC_MS2"),
                         mz=mzs_of_interest, ppm = 5)

all_data$MS1 %>%
  mutate(type="All data") %>%
  rbind(small_data$EIC %>% mutate(type="Extracted data only")) %>%
  filter(!str_detect(filename, "DDA")) %>%
  filter(rt%between%c(5, 15)) %>%
  group_by(rt, filename, type) %>%
  summarise(TIC=sum(int), .groups="drop") %>%
  ggplot() +
  geom_line(aes(x=rt, y=TIC, color=filename)) +
  facet_wrap(~type, ncol = 1)

# Size reduction factor:
as.numeric(object.size(all_data)/object.size(small_data))

## ----rtrangedemo--------------------------------------------------------------
small_data <- grabMSdata(data_files, grab_what = c("MS1", "MS2"), rtrange = c(6, 8))

all_data$MS1 %>%
  mutate(type="All data") %>%
  rbind(small_data$MS1 %>% mutate(type="Extracted data only")) %>%
  filter(!str_detect(filename, "DDA")) %>%
  group_by(rt, filename, type) %>%
  summarise(TIC=sum(int), .groups="drop") %>%
  ggplot() +
  geom_line(aes(x=rt, y=TIC, color=filename)) +
  facet_wrap(~type, ncol = 1)

# Size reduction factor:
as.numeric(object.size(all_data)/object.size(small_data))

## ----verbosedemo--------------------------------------------------------------
all_data <- grabMSdata(data_files, grab_what = c("MS1", "MS2"), verbosity = 2)

## ---- eval=FALSE--------------------------------------------------------------
#  ## Not run:
#  library(parallel)
#  cl <- makeCluster(getOption("cl.cores", detectCores()-1))
#  output_data <- parLapply(data_files, grabMzmlData, grab_what="everything", cl = cl)
#  
#  library(foreach)
#  library(doParallel)
#  registerDoParallel(detectCores()-1)
#  output_data <- foreach (i=data_files) %dopar% {
#    RaMS::grabMzmlData(i, grab_what="everything")
#  }
#  stopImplicitCluster()

## ---- eval=FALSE--------------------------------------------------------------
#  ## Not run:
#  data_nodes <- xml2::xml_find_all(mzML_nodes, xpath="//d1:precursorMz")
#  raw_data <- xml2::xml_attr(data_nodes, "value")

## ---- eval=FALSE--------------------------------------------------------------
#  ## Not run:
#  decoded_binary <- base64enc::base64decode(binary)
#  raw_binary <- as.raw(decoded_binary)
#  decomp_binary <- memDecompress(raw_binary, type = file_metadata$compression)
#  final_binary <- readBin(decomp_binary, what = "double",
#                          n=length(decomp_binary)/file_metadata$mz_precision,
#                          size = file_metadata$mz_precision)
#  
#  # See https://github.com/ProteoWizard/pwiz/issues/1301

