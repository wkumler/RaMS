## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "80%",
  fig.align = 'center',
  fig.height = 3,
  dpi = 36
)

## ----findfiles, message=FALSE-------------------------------------------------
library(RaMS)
library(data.table)

# Locate the file directory
msdata_dir <- system.file("extdata", package = "RaMS")

# Identify the files of interest
data_files <- list.files(msdata_dir, pattern = "Full.*mzML", full.names = TRUE)

# Check that the files identified are the ones expected
basename(data_files)

## ----loadTIC------------------------------------------------------------------
single_file <- data_files[1]

file_data <- grabMSdata(single_file, grab_what = "TIC")

knitr::kable(head(file_data$TIC, 3))

## ----headerTIC----------------------------------------------------------------
par(mar=c(4.1, 4.1, 0.1, 0.1))
plot(file_data$TIC$rt, file_data$TIC$int, type = "l")

## ----loadBPC------------------------------------------------------------------
file_data <- grabMSdata(single_file, grab_what = "BPC")

## ----plotBPC------------------------------------------------------------------
library(tidyverse)
ggplot(file_data$BPC) + geom_line(aes(x=rt, y=int))

## ----loadmultiBPC-------------------------------------------------------------
file_data <- grabMSdata(data_files, grab_what = "BPC", verbosity = "minimal")

ggplot(file_data$BPC) + geom_line(aes(x=rt, y=int, color=filename))

## ----ggplotRaMS, dev.args=list(png  = list(type = "cairo"))-------------------
ggplot(file_data$BPC) + 
  geom_line(aes(x=rt, y=int, color=filename), lwd=1.2) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text = element_text(color = "#AA0000"),
        axis.title = element_text(family = "serif"),
        plot.title = element_text(face = "bold")) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x="Retention time (min)", y="Intensity", 
       title = "My chromatogram", color="File names:")

## -----------------------------------------------------------------------------
# Since the minification process strips a lot of useful data, I use DDA here
metadata_files <- list.files(msdata_dir, pattern = "DDA", full.names = TRUE)
grabMSdata(metadata_files, grab_what = "metadata")

## -----------------------------------------------------------------------------
file_data <- grabMSdata(data_files, grab_what = "MS1")

knitr::kable(head(file_data$MS1, 3))

## -----------------------------------------------------------------------------
adenine_mz <- 136.06232

adenine_data <- file_data$MS1[mz%between%pmppm(adenine_mz, ppm=5)]

ggplot(adenine_data) + geom_line(aes(x=rt, y=int, color=filename))

## -----------------------------------------------------------------------------
masses_of_interest <- c(adenine=136.06232, valine=118.0865, homarine=138.055503)

mass_data <- imap(masses_of_interest, function(mz_i, name){
  cbind(file_data$MS1[mz%between%pmppm(mz_i, ppm=5)], name)
}) %>% rbindlist()

ggplot(mass_data) + 
  geom_line(aes(x=rt, y=int, color=filename)) + 
  facet_wrap(~name, ncol = 1, scales = "free_y")

## -----------------------------------------------------------------------------
DDA_file <- list.files(msdata_dir, pattern = "DDA.*mzML", full.names = TRUE)
DDA_data <- grabMSdata(DDA_file, grab_what = c("MS2"))
knitr::kable(head(DDA_data$MS2, 3))

## -----------------------------------------------------------------------------
betaine_mass <- 118.0865

betaine_MS2 <- DDA_data$MS2[premz%between%pmppm(betaine_mass, 5)]
betaine_MS2$int <- betaine_MS2$int/max(betaine_MS2$int)*100

ggplot(betaine_MS2) +
  geom_point(aes(x=fragmz, y=int)) +
  geom_segment(aes(x=fragmz, xend=fragmz, y=int, yend=0)) +
  labs(x="Fragment m/z", y="Relative intensity (%)")

## ----plotly, eval=FALSE-------------------------------------------------------
#  ## Not run to save space in the vignette:
#  library(plotly)
#  data_files <- list.files(msdata_dir, pattern = "mzML", full.names = TRUE)
#  file_data <- grabMSdata(data_files, grab_what = c("MS1", "MS2", "BPC"))
#  
#  clean_MS2 <- file_data$MS2 %>%
#    filter(premz%between%pmppm(betaine_mass)) %>%
#    group_by(rt) %>%
#    arrange(desc(int)) %>%
#    summarise(frags=paste(
#      paste(round(fragmz, digits = 3), round(int), sep = ": "), collapse = "\n"),
#      .groups="drop"
#    )
#  file_data$MS1 %>%
#    filter(mz%between%pmppm(betaine_mass)) %>%
#    filter(!str_detect(filename, "DDA")) %>%
#    plot_ly() %>%
#    add_trace(type="scatter", mode="lines", x=~rt, y=~int, color=~filename,
#              hoverinfo="none") %>%
#    add_trace(type="scatter", mode="markers", x=~rt, y=0,
#              text=~frags, hoverinfo="text", showlegend=FALSE,
#              marker=list(color="black"), data = clean_MS2) %>%
#    layout(annotations=list(x=min(clean_MS2$rt), y=0,
#                            text="Mouse over to see\nMSMS fragments"))

## -----------------------------------------------------------------------------
betaine_frag_mz <- 58.0660
knitr::kable(head(DDA_data$MS2[fragmz%between%pmppm(betaine_frag_mz, ppm = 5)]))

## -----------------------------------------------------------------------------
betaine_mass <- 118.0865
betaine_neutral_loss <- betaine_mass - betaine_frag_mz

file_data$MS2 <- mutate(DDA_data$MS2, neutral_loss=premz-fragmz) %>%
  select("rt", "premz", "fragmz", "neutral_loss", "int", "voltage", "filename")
file_data$MS2[neutral_loss%between%pmppm(betaine_neutral_loss, ppm = 5)] %>%
  head() %>% knitr::kable()

## -----------------------------------------------------------------------------
data_files <- list.files(msdata_dir, pattern = "mzML", full.names = TRUE)
all_data <- grabMSdata(data_files, grab_what = c("MS1", "MS2"))

masses_of_interest <- c(adenine=136.06232, valine=118.0865, homarine=138.055503)
small_data <- grabMSdata(data_files, grab_what = c("EIC", "EIC_MS2"),
                         mz=masses_of_interest, ppm = 5)

all_data$MS1 %>%
  mutate(type="All data") %>%
  rbind(small_data$EIC %>% mutate(type="Extracted data only")) %>%
  filter(!str_detect(filename, "DDA")) %>%
  group_by(rt, filename, type) %>%
  summarise(TIC=sum(int), .groups="drop") %>%
  ggplot() +
  geom_line(aes(x=rt, y=TIC, color=filename)) +
  facet_wrap(~type, ncol = 1) +
  theme(legend.position = "top")

# Size reduction factor:
as.numeric(object.size(all_data)/object.size(small_data))

## -----------------------------------------------------------------------------
small_data <- grabMSdata(data_files, grab_what = c("MS1", "MS2"), rtrange = c(5, 6))

all_data$MS1 %>%
  mutate(type="All data") %>%
  rbind(small_data$MS1 %>% mutate(type="Extracted data only")) %>%
  filter(!str_detect(filename, "DDA")) %>%
  group_by(rt, filename, type) %>%
  summarise(TIC=sum(int), .groups="drop") %>%
  ggplot() +
  geom_line(aes(x=rt, y=TIC, color=filename)) +
  facet_wrap(~type, ncol = 1) +
  theme(legend.position = "top")

# Size reduction factor:
as.numeric(object.size(all_data)/object.size(small_data))

## -----------------------------------------------------------------------------
all_data <- grabMSdata(data_files, grab_what = c("MS1", "MS2"), verbosity = "very")

## ---- eval=FALSE--------------------------------------------------------------
#  ## Not run:
#  data_files <- list.files(msdata_dir, pattern = "mzML", full.names = TRUE)
#  
#  library(parallel)
#  output_data <- mclapply(data_files, grabMzmlData, mc.cores = detectCores()-1)
#  
#  library(foreach)
#  library(doParallel)
#  registerDoParallel(numCores)
#  output_data <- foreach (i=data_files) %dopar% {
#    grabMzmlData(i)
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

