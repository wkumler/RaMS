Rapid access to Mass-Spec data (RaMS)
================

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

There’s only one main function in `RaMS`: the aptly named `grabMSdata`.
This function accepts the names of mass-spectrometry files as well as
the data you’d like to extract (e.g. MS1, MS2, BPC, etc.) and produces a
list of data tables. Each table is intuitively named within the list and
formatted tidily:

``` r
msdata_files <- list.files("inst/extdata", pattern = "mzML", full.names = TRUE)
output <- grabMSdata(files = msdata_files, grab_what = c("TIC", "MS1"))
```

``` r
knitr::kable(head(output$MS1))
```

|       rt |       mz |       int | filename                |
| -------: | -------: | --------: | :---------------------- |
| 2.004194 |  90.0555 |  25501.37 | FK180310\_Full1.mzML.gz |
| 2.004194 | 104.0710 | 572123.50 | FK180310\_Full1.mzML.gz |
| 2.004194 | 112.0509 |  50395.21 | FK180310\_Full1.mzML.gz |
| 2.004194 | 116.0709 |  60353.67 | FK180310\_Full1.mzML.gz |
| 2.004194 | 118.0865 | 483736.62 | FK180310\_Full1.mzML.gz |
| 2.004194 | 122.0270 |  11133.76 | FK180310\_Full1.mzML.gz |

``` r
knitr::kable(head(output$TIC))
```

|       rt |     int | filename                |
| -------: | ------: | :---------------------- |
| 2.004194 | 2594068 | FK180310\_Full1.mzML.gz |
| 2.021873 | 2396701 | FK180310\_Full1.mzML.gz |
| 2.037056 | 3098733 | FK180310\_Full1.mzML.gz |
| 2.052177 | 3416233 | FK180310\_Full1.mzML.gz |
| 2.067364 | 3536198 | FK180310\_Full1.mzML.gz |
| 2.082489 | 3772849 | FK180310\_Full1.mzML.gz |

This means that the basic R functions work exactly as we expect them to,
no new method definition or memorization necessary:

``` r
str(output)
```

    ## List of 2
    ##  $ MS1:Classes 'data.table' and 'data.frame':    58565 obs. of  4 variables:
    ##   ..$ rt      : num [1:58565] 2 2 2 2 2 ...
    ##   ..$ mz      : num [1:58565] 90.1 104.1 112.1 116.1 118.1 ...
    ##   ..$ int     : num [1:58565] 25501 572124 50395 60354 483737 ...
    ##   ..$ filename: chr [1:58565] "FK180310_Full1.mzML.gz" "FK180310_Full1.mzML.gz" "FK180310_Full1.mzML.gz" "FK180310_Full1.mzML.gz" ...
    ##   ..- attr(*, ".internal.selfref")=<externalptr> 
    ##  $ TIC:Classes 'data.table' and 'data.frame':    1573 obs. of  3 variables:
    ##   ..$ rt      : num [1:1573] 2 2.02 2.04 2.05 2.07 ...
    ##   ..$ int     : num [1:1573] 2594068 2396701 3098733 3416233 3536198 ...
    ##   ..$ filename: chr [1:1573] "FK180310_Full1.mzML.gz" "FK180310_Full1.mzML.gz" "FK180310_Full1.mzML.gz" "FK180310_Full1.mzML.gz" ...
    ##   ..- attr(*, ".internal.selfref")=<externalptr>

``` r
summary(output$MS1)
```

    ##        rt              mz              int              filename        
    ##  Min.   :2.004   Min.   : 60.08   Min.   :     5683   Length:58565      
    ##  1st Qu.:4.045   1st Qu.:127.05   1st Qu.:    21087   Class :character  
    ##  Median :6.054   Median :143.08   Median :    47865   Mode  :character  
    ##  Mean   :5.945   Mean   :152.03   Mean   :   997549                     
    ##  3rd Qu.:7.707   3rd Qu.:166.09   3rd Qu.:   162405                     
    ##  Max.   :9.998   Max.   :425.18   Max.   :294877088

``` r
# Outputs are data.tables so we can use their intuitive indexing on column name
first_file_data <- output$TIC[filename==basename(msdata_files[1])]
plot(first_file_data$rt, first_file_data$int, type = "l")
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

And of course, the tidy data format means that it plays nicely with
every other tidy data package.

``` r
# Plot some chromatograms with ggplot2
library(ggplot2)
ggplot(output$TIC) + geom_line(aes(x = rt, y=int, color=filename)) + theme(legend.position="top")
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggplot(output$TIC) + geom_line(aes(x = rt, y=int)) +
  facet_wrap(~filename, scales = "free_y", ncol = 1) +
  labs(x="Retention time (min)", y="Intensity")
```

![](README_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
# Manually calculate a BPC with dplyr
suppressMessages(library(dplyr))
manual_TIC <- output$MS1 %>%
  group_by(rt, filename) %>%
  summarise(manual_int=sum(int), .groups="drop") %>%
  left_join(output$TIC, by = c("rt", "filename"))

# Return to base to plot TIC comparison
par(pty="s")
plot(manual_TIC$int, manual_TIC$manual_int, 
     ylab = "Manually calculated TIC",
     xlab = "TIC read from file")
grid()
corr <- cor(manual_TIC$int, manual_TIC$manual_int)
legend("topleft", legend=paste("Correlation =", corr))
```

![](README_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

RaMS is currently limited to the modern **mzML** data format and the
slightly older **mzXML** format. Tools to convert data from other
formats are available through
[Proteowizard](http://proteowizard.sourceforge.net/tools.shtml)’s
`msconvert` tool. Data can, however, be gzip compressed (file ending
.gz) and this compression actually speeds up data retrieval
significantly as well as reducing file sizes.

## Motivation
