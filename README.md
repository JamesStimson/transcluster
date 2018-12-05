# **transcluster**

### Introduction

This is **transcluster**, a package for inferring and viewing transmission clusters from sequence alignments and sample dates.

The package takes two kinds of input data

* SNP distance matrix
* Sample dates 

which can be provided directly or from fasta and csv files respectively.

There is a simple starter vignette in the vignettes folder.

A formal description of the methods is written up in [Beyond the SNP threshold: identifying outbreak clusters using inferred transmissions](https://www.biorxiv.org/content/early/2018/12/03/319707), by 
James Stimson, Jennifer Leigh Gardy, Barun Mathema, Valeriu Crudu, Theodore Cohen, and Caroline Colijn.

### Installation

You can install **transcluster** in **R** using the following command:
```{r}
devtools::install_github("JamesStimson/transcluster", build_vignettes = TRUE)
```

There are some example input files that come with the installation. To find out where they are on your system, use system.file() like this:
```{r}
system.file("extdata", "bc_snp_matrix_R.csv", package = "transcluster", mustWork = TRUE)
```

You will see something like this in response:
```{r}
[1] "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/transcluster/extdata/bc_snp_matrix_R.csv"
```

### Getting help

To view the vignette once installed, run
```{r}
vignette("how-to-guide", package = "transcluster")
```

Alternatively, you can run the R markdown *vignettes/how-to-guide.Rmd* yourself.

If you need further assistance using **transcluster**, you can get in touch by emailing 
```{r}
james.stimson16@imperial.ac.uk
```
