---
title: "README"
author: "Binu Mathew"
date: "April 24, 2018"
output: html_document
---
# ravel2012 data package

This package contains one phyloseq object of data from Ravel et al.(2012). It contains data from following sample:

    Consists high-throughput profiling data collected from vaginal microbiota samples collected from women, who were not pregnant of reproductive age, ranging from 12 to 45 years, regularly menstruating, with a history of sexual activity, and had not taken any antibiotic or antimycotic compounds in the past 30 days.
    
    
Data processing code is available in inst/exdata.

Further information regarding the study can be obtained from: Ravel J, Gajer P, Abdo Z, Schneider GM, Koenig SS, McCulle SL, et al. Vaginal microbiome of reproductive-age women. Proc Natl Acad Sci U S A (2011) 108(Suppl 1):4680–7. doi:10.1073/pnas.1002611107

This phyloseq object contains:

   1. otu_table: which contains OTU counts.
   2. tax_table: created by using Greengenes database as reference
   3. sam_data: it contains metadata



# Installing package 


```{r}
install_github("binmat/ravel2012")
```

