R-based access to Mass-Spec data (RaMS)
================

<!-- badges: start -->

[![R build
status](https://github.com/wkumler/RaMS/workflows/R-CMD-check/badge.svg)](https://github.com/wkumler/RaMS/actions/)
<!-- badges: end -->

## Overview

`RaMS` is a lightweight package that provides rapid and tidy access to
mass-spectrometry data. This package is *lightweight* because it’s built
from the ground up rather than relying on an extensive network of
external libraries. No Rcpp, no Bioconductor, no long load times and
strange startup warnings. Just XML parsing provided by `xml2` and data
handling provided by `data.table`. Access is *rapid* because an absolute
minimum of data processing occurs. Unlike other packages, `RaMS` makes
no assumptions about what you’d like to do with the data and is simply
providing access to the encoded information in an intuitive and
R-friendly way. Finally, the access is *tidy* in the philosophy of [tidy
data](https://r4ds.had.co.nz/tidy-data.html). Tidy data neatly resolves
the ragged arrays that mass spectrometers produce and plays nicely with
other [tidy data packages](https://www.tidyverse.org/).

## Installation

Until `RaMS` is on CRAN, the easiest way to install is via `devtools`:

``` r
devtools::install_github("wkumler/RaMS")

library(RaMS)
```

## Usage

(For more usage examples, see the
[vignette](vignettes/RaMS-vignette.pdf).)

There’s only one main function in `RaMS`: the aptly named `grabMSdata`.
This function accepts the names of mass-spectrometry files as well as
the data you’d like to extract (e.g. MS1, MS2, BPC, etc.) and produces a
list of data tables. Each table is intuitively named within the list and
formatted tidily:

``` r
demo_dir <- system.file("extdata", package = "RaMS")
msdata_files <- list.files(demo_dir, pattern = "mzML", full.names = TRUE)
output <- grabMSdata(files = msdata_files, grab_what = c("TIC", "MS1", "MS2"))
```

``` r
knitr::kable(head(output$TIC, 3))
```

|      rt |       int | filename          |
|--------:|----------:|:------------------|
| 240.051 | 291632610 | DDApos\_2.mzML.gz |
| 240.393 | 304016350 | DDApos\_2.mzML.gz |
| 240.966 | 299228260 | DDApos\_2.mzML.gz |

``` r
knitr::kable(head(output$MS1, 3))
```

|      rt |       mz |       int | filename          |
|--------:|---------:|----------:|:------------------|
| 240.051 | 80.05009 | 12057.776 | DDApos\_2.mzML.gz |
| 240.051 | 80.26269 |  8178.767 | DDApos\_2.mzML.gz |
| 240.051 | 80.94841 | 19075.213 | DDApos\_2.mzML.gz |

``` r
knitr::kable(head(output$MS2, 3))
```

|      rt |    premz |   fragmz |       int | voltage | filename          |
|--------:|---------:|---------:|----------:|--------:|:------------------|
| 240.184 | 127.0325 | 56.05023 | 17019.613 |      35 | DDApos\_2.mzML.gz |
| 240.184 | 127.0325 | 59.46344 |  1003.812 |      35 | DDApos\_2.mzML.gz |
| 240.184 | 127.0325 | 70.06580 |  1521.526 |      35 | DDApos\_2.mzML.gz |

This means that basic R functions work exactly as we expect them to, no
new functionality necessary:

``` r
# Outputs are data.tables so we can use their intuitive indexing on column name
first_file_data <- output$TIC[filename==basename(msdata_files[2])]
plot(first_file_data$rt, first_file_data$int, type = "l")
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

------------------------------------------------------------------------

And of course, the tidy data format means that it plays nicely with
every other tidy data package.

### Chromatograms with ggplot2

``` r
output <- grabMSdata(files = msdata_files[-1], grab_what = c("TIC", "MS1"))
```

``` r
library(ggplot2)
ggplot(output$TIC) + geom_line(aes(x = rt, y=int, color=filename)) + theme(legend.position="top")
```

![](man/figures/README-unnamed-chunk-6-1.png)<!-- -->

``` r
ggplot(output$TIC) + geom_line(aes(x = rt, y=int)) +
  facet_wrap(~filename, scales = "free_y", ncol = 1) +
  labs(x="Retention time (min)", y="Intensity")
```

![](man/figures/README-unnamed-chunk-6-2.png)<!-- -->

### Interactive MSMS with dplyr, stringr, and plotly

``` r
output <- grabMSdata(files = msdata_files, grab_what = c("EIC", "EIC_MS2"),
                     mz=118.0865, ppm=5)

library(dplyr)
library(stringr)
library(plotly)

clean_EIC_MS2 <- output$EIC_MS2 %>% 
  group_by(rt) %>%
  arrange(desc(int)) %>%
  summarise(frags=paste(
    paste(round(fragmz, digits = 3), round(int), sep = ": "), collapse = "\n")
  )
output$EIC %>% 
  filter(!str_detect(filename, "DDA")) %>%
  plot_ly() %>%
  add_trace(type="scatter", mode="lines", x=~rt, y=~int, color=~filename,
            hoverinfo="none") %>%
  add_trace(type="scatter", mode="markers", x=~rt, y=0,
            text=~frags, hoverinfo="text", showlegend=FALSE,
            marker=list(color="black"), data = clean_EIC_MS2) %>%
  layout(annotations=list(x=min(clean_EIC_MS2$rt), y=0, 
                          text="Mouse over to see\nMSMS fragments"),
         title="(See vignette for interactive version)")
```

![](man/figures/plotlyplot.png)

For more usage examples, see [the
vignette](vignettes/RaMS-vignette.pdf).

## File types

RaMS is currently limited to the modern **mzML** data format and the
slightly older **mzXML** format. Tools to convert data from other
formats are available through
[Proteowizard](http://proteowizard.sourceforge.net/tools.shtml)’s
`msconvert` tool. Data can, however, be gzip compressed (file ending
.gz) and this compression actually speeds up data retrieval
significantly as well as reducing file sizes.

Currently, `RaMS` also handles only MS<sup>1</sup> and MS<sup>2</sup>
data. This should be easy enough to expand in the future, but right now
I haven’t observed a demonstrated need for higher fragmentation level
data collection.

Additionally, note that files can be streamed from the internet directly
if a URL is provided to `grabMSdata` and the `check_exists` argument is
set to `FALSE`, although this will usually take longer than reading a
file from disk:

``` r
## Not run:
# Find a file with a web browser:
browseURL("https://www.ebi.ac.uk/metabolights/MTBLS703/files")

# Copy link address by right-clicking "download" button:
sample_url <- paste0("https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS703/",
                     "download/acefcd61-a634-4f35-9c3c-c572ade5acf3?file=",
                     "161024_Smp_LB12HL_AB_pos.mzXML")
file_data <- grabMSdata(sample_url, grab_what="everything",
                       check_exists=FALSE, verbosity="very")
file_data$metadata
```
