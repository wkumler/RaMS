R-based access to Mass-Spec data (RaMS)
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R build
status](https://github.com/wkumler/RaMS/workflows/R-CMD-check/badge.svg)](https://github.com/wkumler/RaMS/actions/)
[![Codecov test
coverage](https://codecov.io/gh/wkumler/RaMS/branch/master/graph/badge.svg)](https://codecov.io/gh/wkumler/RaMS)
<!-- badges: end -->

**Table of contents:**
[Overview](https://github.com/wkumler/RaMS#overview) -
[Installation](https://github.com/wkumler/RaMS#installation) -
[Usage](https://github.com/wkumler/RaMS#usage) - [File
types](https://github.com/wkumler/RaMS#file-types) -
[Contact](https://github.com/wkumler/RaMS#contact)

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

To install the stable version on CRAN:

``` r
install.packages('RaMS')
```

To install the current development version:

``` r
devtools::install_github("wkumler/RaMS")
```

Finally, load RaMS like every other package:

``` r
library(RaMS)
```

## Usage

There’s only one main function in `RaMS`: the aptly named `grabMSdata`.
This function accepts the names of mass-spectrometry files as well as
the data you’d like to extract (e.g. MS1, MS2, BPC, etc.) and produces a
list of data tables. Each table is intuitively named within the list and
formatted tidily:

``` r
msdata_dir <- system.file("extdata", package = "RaMS")
msdata_files <- list.files(msdata_dir, pattern = "mzML", full.names=TRUE)

msdata <- grabMSdata(files = msdata_files[2:4], grab_what = c("BPC", "MS1"))
```

#### BPC/TIC data:

Base peak chromatograms (BPCs) and total ion chromatograms (TICs) have
three columns, making them super-simple to plot with either base R or
the popular \[ggplot2\] library:

``` r
knitr::kable(head(msdata$BPC, 3))
```

|       rt |      int | filename           |
|---------:|---------:|:-------------------|
| 4.009000 | 11141859 | LB12HL\_AB.mzML.gz |
| 4.024533 |  9982309 | LB12HL\_AB.mzML.gz |
| 4.040133 | 10653922 | LB12HL\_AB.mzML.gz |

``` r
plot(msdata$BPC$rt, msdata$BPC$int, type = "l")
```

![](man/figures/README-showbaseplot-1.png)<!-- -->

``` r
library(ggplot2)
ggplot(msdata$BPC) + geom_line(aes(x = rt, y=int, color=filename)) +
  facet_wrap(~filename, scales = "free_y", ncol = 1) +
  labs(x="Retention time (min)", y="Intensity", color="File name: ") +
  theme(legend.position="top")
```

![](man/figures/README-showggplot-1.png)<!-- -->

#### MS1 data:

MS<sup>1</sup> data includes an additional dimension, the *m/z* of each
ion measured, and has multiple entries per retention time:

``` r
knitr::kable(head(msdata$MS1, 3))
```

|    rt |       mz |        int | filename           |
|------:|---------:|-----------:|:-------------------|
| 4.009 | 104.0710 | 1297755.00 | LB12HL\_AB.mzML.gz |
| 4.009 | 104.1075 |  140668.12 | LB12HL\_AB.mzML.gz |
| 4.009 | 112.0509 |   67452.86 | LB12HL\_AB.mzML.gz |

This tidy format means that it plays nicely with other tidy data
packages. Here, we use \[data.table\] and a few other tidyverse packages
to compare a molecule’s 13C and 15N peak areas to that of the base peak,
giving us some clue as to its molecular formula.

``` r
library(data.table)
library(tidyverse)

M <- 118.0865
M_13C <- M + 1.003355
M_15N <- M + 0.997035

iso_data <- imap_dfr(lst(M, M_13C, M_15N), function(mass, isotope){
  peak_data <- msdata$MS1[mz%between%pmppm(mass) & rt%between%c(7.6, 8.2)]
  cbind(peak_data, isotope)
})

iso_data %>%
  group_by(filename, isotope) %>%
  summarise(area=sum(int)) %>%
  pivot_wider(names_from = isotope, values_from = area) %>%
  mutate(ratio_13C_12C = M_13C/M) %>%
  mutate(ratio_15N_14N = M_15N/M) %>%
  select(filename, contains("ratio")) %>%
  pivot_longer(cols = contains("ratio"), names_to = "isotope") %>%
  group_by(isotope) %>%
  summarize(avg_ratio = mean(value), sd_ratio = sd(value), .groups="drop") %>%
  mutate(isotope=str_extract(isotope, "(?<=_).*(?=_)")) %>%
  knitr::kable()
```

| isotope | avg\_ratio | sd\_ratio |
|:--------|-----------:|----------:|
| 13C     |  0.0543929 | 0.0006015 |
| 15N     |  0.0033375 | 0.0001846 |

With [natural
abundances](https://en.wikipedia.org/wiki/Natural_abundance) for
<sup>13</sup>C and <sup>15</sup>N of 1.11% and 0.36%, respectively, we
can conclude that this molecule likely has five carbons and a single
nitrogen.

Of course, it’s always a good idea to plot the peaks and perform a
manual check of data quality:

``` r
ggplot(iso_data) +
  geom_line(aes(x=rt, y=int, color=filename)) +
  facet_wrap(~isotope, scales = "free_y", ncol = 1)
```

![](man/figures/README-isoexampleplot-1.png)<!-- -->

#### MS2 data:

DDA (fragmentation) data can also be extracted, allowing rapid and
intuitive searches for fragments or neutral losses:

``` r
msdata <- grabMSdata(files = msdata_files[1], grab_what = "MS2")
```

For example, we may be interested in the major fragments of a specific
molecule:

``` r
msdata$MS2[premz%between%pmppm(118.0865) & int>mean(int)] %>%
  plot(int~fragmz, type="h", data=., ylab="Intensity", xlab="Fragment m/z")
```

![](man/figures/README-plotfragdata-1.png)<!-- -->

Or want to search for a specific neutral loss:

``` r
msdata$MS2[, neutral_loss:=premz-fragmz] %>%
  filter(neutral_loss%between%pmppm(60.02064, 5))
```

    ##           rt    premz    fragmz          int voltage         filename
    ##  1: 4.182333 118.0864  58.06590   390179.500      35 DDApos_2.mzML.gz
    ##  2: 4.276100 116.0709  56.05036     1093.988      35 DDApos_2.mzML.gz
    ##  3: 4.521367 118.0864  58.06589   343084.000      35 DDApos_2.mzML.gz
    ##  4: 4.649867 170.0810 110.06034     4792.479      35 DDApos_2.mzML.gz
    ##  5: 4.857983 118.0865  58.06590   314075.312      35 DDApos_2.mzML.gz
    ##  6: 5.195617 118.0865  58.06590   282611.688      35 DDApos_2.mzML.gz
    ##  7: 5.536383 118.0865  58.06592   300432.906      35 DDApos_2.mzML.gz
    ##  8: 5.642417 116.0709  56.05035     1985.602      35 DDApos_2.mzML.gz
    ##  9: 5.719433 138.0550  78.03460     1749.604      35 DDApos_2.mzML.gz
    ## 10: 5.876633 118.0865  58.06593   336346.625      35 DDApos_2.mzML.gz
    ## 11: 6.063700 138.0549  78.03448  4783814.000      35 DDApos_2.mzML.gz
    ## 12: 6.197033 119.0899  59.06926     1459.779      35 DDApos_2.mzML.gz
    ## 13: 6.214217 118.0865  58.06592    87062.281      35 DDApos_2.mzML.gz
    ## 14: 6.399933 138.0549  78.03451  4311163.000      35 DDApos_2.mzML.gz
    ## 15: 6.552300 118.0865  58.06591   188614.250      35 DDApos_2.mzML.gz
    ## 16: 6.569983 174.1125 114.09184     4416.066      35 DDApos_2.mzML.gz
    ## 17: 6.648017 119.0899  59.06929     7744.978      35 DDApos_2.mzML.gz
    ## 18: 6.739650 138.0550  78.03453    45146.910      35 DDApos_2.mzML.gz
    ## 19: 6.889233 118.0865  58.06590   247183.109      35 DDApos_2.mzML.gz
    ## 20: 7.080050 138.0550  78.03455    15764.113      35 DDApos_2.mzML.gz
    ## 21: 7.232333 118.0865  58.06590   317423.344      35 DDApos_2.mzML.gz
    ## 22: 7.418867 138.0550  78.03444     5292.677      35 DDApos_2.mzML.gz
    ## 23: 7.571183 118.0865  58.06591  1110205.500      35 DDApos_2.mzML.gz
    ## 24: 7.913200 118.0863  58.06591 53963628.000      35 DDApos_2.mzML.gz
    ## 25: 8.094183 148.0425  88.02216     1632.094      35 DDApos_2.mzML.gz
    ## 26: 8.197183 162.1124 102.09158     1861.552      35 DDApos_2.mzML.gz
    ## 27: 8.250383 118.0864  58.06590  3601956.250      35 DDApos_2.mzML.gz
    ## 28: 8.262433 119.0835  59.06298    17410.162      35 DDApos_2.mzML.gz
    ## 29: 8.576400 130.0863  70.06579     7045.373      35 DDApos_2.mzML.gz
    ## 30: 8.585383 118.0864  58.06588  1340395.125      35 DDApos_2.mzML.gz
    ## 31: 8.925200 118.0864  58.06589   920915.938      35 DDApos_2.mzML.gz
    ##           rt    premz    fragmz          int voltage         filename
    ##     neutral_loss
    ##  1:     60.02055
    ##  2:     60.02050
    ##  3:     60.02056
    ##  4:     60.02070
    ##  5:     60.02057
    ##  6:     60.02057
    ##  7:     60.02060
    ##  8:     60.02058
    ##  9:     60.02043
    ## 10:     60.02057
    ## 11:     60.02037
    ## 12:     60.02066
    ## 13:     60.02054
    ## 14:     60.02039
    ## 15:     60.02060
    ## 16:     60.02068
    ## 17:     60.02063
    ## 18:     60.02046
    ## 19:     60.02060
    ## 20:     60.02043
    ## 21:     60.02058
    ## 22:     60.02058
    ## 23:     60.02056
    ## 24:     60.02042
    ## 25:     60.02038
    ## 26:     60.02086
    ## 27:     60.02052
    ## 28:     60.02054
    ## 29:     60.02052
    ## 30:     60.02055
    ## 31:     60.02052
    ##     neutral_loss

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

## Contact

Feel free to submit questions, bugs, or feature requests on the [GitHub
Issues page](https://github.com/wkumler/RaMS/issues).

------------------------------------------------------------------------

README last built on 2021-03-17
