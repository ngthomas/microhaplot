# microhaplot   

`microhaplot` generates visual summaries of microhaplotypes found in short read alignments. All you need are alignment SAM files and variant call VCF file.

This software exists as a an R package `microhaplot` that includes within it the code to set up and 
establish an Rstudio/Shiny server to visualize and manipulate the data.  There are two key steps in 
the `microhaplot` worflow:

1. The first step is to summarize alignment and variation data into a single data frame that is 
easily operated upon.  This is done using the function `microhaplot::runHaplot`.  You must supply a 
VCF file that includes variants that you are interested in extracting, and as many SAM files 
(one for each individual) that you want to extract read information from at each of the variants. 
The function `microhaplot::runHaplot` makes a call
to PERL to parse the CIGAR strings in the SAM files to extract the variant information at each read
and store this information into a data frame which gets saved with the installed Shiny app (see below)
for later use.  Depending on the size of the data set, this can take a few minutes.  

2. The second step is to run the haPLOType Shiny app to visualize the sequence information, call genotypes using
simple read-depth based filtering criteria, and curate the loci. haPLOType is suitable for quick assesement
and quality control of haplotype generated from library runs. Plot summaries include read depth, fraction of callable haplotypes, Hardy-Weinberg
equilibrium plots, and more. 


See the **Example Data** section to learn about how to run each of these steps on the example data that are provided
with the package.  

   
### Installation and Quick Start

You can either clone the repository and build the `microhaplot` package yourself, or, more easily, you can
install it using  [devtools](https://github.com/hadley/devtools). You can get `devtools` by `install.packages("devtools")`.

Once you have `devtools` available in R, you can get `microhaplot` this way:
```r
devtools::install_github("ngthomas/microhaplot", build_vignettes = TRUE)
```

Once you have installed the `microhaplot` R package with devtools there you need to use the `microhaplot::mvHaplotype`
to establish the haPLOType Shiny App in a convenient location on your system. The following line
creates the directory `Shiny` in my home directory and then within that it creates the 
directory `haPLOType` and fills it with the Shiny app as well as the example data that go 
along with that.  

```r
microhaplot::mvHaplotype("~/Shiny") # provide a directory path to host the haPLOType app
```
To start familiarizing yourself with haPLOType using the provided example data.  We recommend
going through our first vignette.  Call it up with:
```r
vignette("haPLOType-walkthrough")
```

Now, having done that, we can launch haPLOType on the example data:
```r
library(microhaplot)
app.path <- "~/Shiny/microhaplot"
runHaplotype(app.path)
```

### Quick Guide to use microhaplot to parse out SAM and VCF files

This microhaplot package comes with a small customized sample data drawn from an actual run 
of short read sequencing run on Rockfish species. The sample data
contains sequences of eight genomic loci for four populations of five individuals each, 
with a total of twenty individuals. 

First you need to create a tab-separate **label** file with 3 info columns: path to SAM file name, individual ID, and group label (in this particular order). If you do not want assign any group label for the individuals, you can just leave it as "NA". It is recommended that you have all of the SAM files under one directory to make this labeling task easier.

The `label` file looks like this:
```txt
s6.sam  s6      gopher
s11.sam s11     copper
s13.sam s13     kelp
s14.sam s14     kelp
s18.sam s18     gopher
``` 

Once you have the label file in place, you can run `runHaplot`, a R function that generates tables of microhaplotype, by providing the following:
 * a label to display in haPLOType
 * path to the directory with all SAM files 
 * path to the `label` file you just created
 * path to the VCF file  
  
```R
library(microhaplot)

run.label <- "sebastes"
sam.path <- system.file("extdata","." , package="microhaplot")
label.path <- system.file("extdata", "label.txt", package = "microhaplot")
vcf.path <- system.file("extdata", "sebastes.vcf", package = "microhaplot")
app.path <- "~/Shiny/microhaplot"
# -or- app.path <- system.file("shiny","microhaplot" , package="microhaplot")

haplo.read.tbl <- runHaplot(run.label = run.label,
           sam.path=sam.path,
           label.path=label.path,
           vcf.path=vcf.path,
           app.path=app.path)
           
runHaplotype(app.path)
```


### Suggestions
- SAM files: For pair-ended experiment, both directional reads should be flashed into one.

- VCF: `SrMicroHap` might have trouble infering individual's true haplotype if no reads are aligned to the variant site.

