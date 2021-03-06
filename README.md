---
title: "Ravel 2012 R data package"
author: "Binu Mathew and Leo Lahti"
date: "April 24, 2018"
output: html_document
---

# R data package: ravel2012

This package contains one phyloseq object of data from Ravel et al. (2012). This data set contains high-throughput profiling data of the vaginal microbiota collected from reproductive-age women (12-45 years). Exclusion criteria included pregnancy, irregular menstruations, sexual activity, and antibiotic or antimycotic compounds in the past 30 days. See the original publication by Ravel et al. (2012) for further details.


# Installation and use 

Install the package in R:

```{r}
library(devtools)
install_github("binmat/ravel2012")
```

Load the package

```{r}
library(ravel2012)
```


# Description

We have organized the original public data set into a phyloseq object and prepared the R data package to facilitate the reuse of this data resource. In particular, we expect that it will be useful for scientists who develop methods to model the dynamics of microbial community profiling data.

The data processing code that converts the original published data set into the R phyloseq object is available in the inst/exdata/ folder. The phyloseq object contains the following elements:

   1. otu_table: OTU abundance (count) table
   2. tax_table: Taxonomy table created with Greengenes reference database
   3. sam_data: Sample metadata (phenotype data)


For further information regarding the study, see: Ravel J, Gajer P, Abdo Z, Schneider GM, Koenig SS, McCulle SL, et al. Vaginal microbiome of reproductive-age women. Proc Natl Acad Sci U S A (2011) 108(Suppl 1):4680–7 [doi:10.1073/pnas.1002611107](https://doi.org/10.1073/pnas.1002611107).



