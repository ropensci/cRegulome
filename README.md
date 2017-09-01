[![Build Status](https://travis-ci.org/MahShaaban/cRegulome.svg?branch=master)](https://travis-ci.org/MahShaaban/cRegulome)
[![Coverage Status](https://img.shields.io/codecov/c/github/MahShaaban/cRegulome/master.svg)](https://codecov.io/github/MahShaaban/cRegulome?branch=master)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/MahShaaban/cRegulome?branch=master&svg=true)](https://ci.appveyor.com/project/MahShaaban/cRegulome)

# cRegulome
## Overview  
Transcription factors and microRNAs are importing for regulating the gene
expression in normal physiology and pathological conditions. Many
bioinformatics tools were built to predict and identify transcription
factors and microRNA targets and their role in development of diseases
including cancers. The availability of public access high-throughput data
allowed for data-driven validations and discoveries of these predictions.
Here, we build on that kind of tools and integrative analysis to provide a
tool to access, manage and visualize data from open source databases.
cRegulome provides a programmatic access to the regulome (microRNA and
transcription factor) correlations with target genes in cancer. The package
obtains a local instance of 
[Cistrome Cancer](http://cistrome.org/CistromeCancer/) and 
[miRCancerdb](https://mahshaaban.shinyapps.io/miRCancerdb/) databases and
provides objects and methods to interact with and visualize the correlation
data.  

## Getting started  
To get starting with cRegulome we show a very quick example. We first start
by downloading a small test database file, make a simple query and convert
the output to a cRegulome object to print and visualize.  


```
# install package from github
devtools::install_github('MahShaaban/cRegulome')
```
```r
# load required libraries
library(cRegulome)
library(R.utils)
library(DBI)
library(RSQLite)
library(dplyr)
library(dbplyr)
library(tidyr)
library(ggplot2)
```

```r
#  download and decompress the database file
if(!file.exists('cRegulome.db')) {
    get_db(test = TRUE)
    gunzip('cRegulome.db.gz')
}
# connect to the database file
conn <- dbConnect(SQLite(), 'cRegulome.db')
```

```r
# query the database
dat <- get_mir(conn,
               mir = c('hsa-let-7b', 'hsa-mir-134'),
               study = c('ACC', 'BLCA'),
               min_cor = .3,
               targets_only = TRUE)

# make a cmicroRNA object               
ob <- cmicroRNA(dat)
dbDisconnect(conn)
```

```r
print(ob)
```

```r
plot(ob, study = 'ACC')
```
## More
Using cRegulome
Case study

## Citation  
Please cite: 
