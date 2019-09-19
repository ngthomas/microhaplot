# microhaplot   

`microhaplot` generates visual summaries of microhaplotypes found in short read alignments. All you need are alignment SAM 
files and a variant call VCF file. (The latter tells `microhaplot` which SNPs to include into microhaplotypes).  It was 
designed for extracting and visualized haplotypes from high-quality amplicon sequencing data.  We have used it extensively
to process amplicon sequencing data (with 100 to 500 amplicons) from rockfish and Chinook salmon, generated on an Illumina 
MiSeq sequencer.  It should be extensible to sequences from capture arrays, like RAPTURE data.

This software exists as an R package `microhaplot` that includes within it the code to set up and 
establish an Rstudio/Shiny server to visualize and manipulate the data.  There are two key steps in 
the `microhaplot` workflow:

1. The first step is to summarize alignment and variant (SNP) data into a single data frame that is 
easily operated upon.  This is done using the function `microhaplot::prepHaplotFiles`.  You must supply a 
VCF file that includes variants that you are interested in extracting, and as many SAM files 
(one for each individual) that you want to extract read information from at each of the variants. 
The function `microhaplot::prepHaplotFiles` makes a call
to PERL to parse the CIGAR strings in the SAM files to extract the variant information at each read
and store this information into a data frame which gets saved with the installed Shiny app (see below)
for later use.  Depending on the size of the data set, this can take a few minutes.  

2. The second step is to run the microhaplot Shiny app to visualize the sequence information, call genotypes using
simple read-depth based filtering criteria, and curate the loci. microhaplot is suitable for quick assessment
and quality control of haplotypes generated from library runs. Plot summaries include read depth, fraction of callable haplotypes, Hardy-Weinberg
equilibrium plots, and more. 


See the **Example Data** section to learn about how to run each of these steps on the example data that are provided
with the package.  

   
### Installation and Quick Start

#### required Perl dependencies:
You need to have Perl (version >5.014) installed in your OS in order to run Microhaplot.  
For Window users, we recommend install it via http://strawberryperl.com/.  
For Mac and Linux users, Perl can be downloaded from https://www.perl.org/get.html  

You can either clone the repository and build the `microhaplot` package yourself, or, more easily, you can
install it using  [devtools](https://github.com/hadley/devtools). You can get `devtools` by `install.packages("devtools")`.
  
**To mac user: remember to install [XQuartz](https://www.xquartz.org/), when upgrading your macOS to a new major version.**   
 
Once you have `devtools` available in R, you can get `microhaplot` this way:
```r
devtools::install_github("ngthomas/microhaplot", build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```

Once you have installed the `microhaplot` R package with devtools there you need to use the `microhaplot::mvHaplotype`
to establish the microhaplot Shiny App in a convenient location on your system. The following line
creates the directory `Shiny` in my home directory and then within that it creates the 
directory `microhaplot` and fills it with the Shiny app as well as the example data that go 
along with that.  

```r
microhaplot::mvShinyHaplot("~/Shiny") # provide a directory path to host the microhaplot app
```
To start familiarizing yourself with microhaplot using the provided example data.  We recommend
going through our first vignette.  Call it up with:
```r
browseVignettes("microhaplot")
```
and check out `microhaplot-walkthrough`.

Now, having done that, we can launch Shiny microhaplot on the example data:
```r
library(microhaplot)
app.path <- "~/Shiny/microhaplot"
runShinyHaplot(app.path)
```

### Quick Guide to use microhaplot to parse out SAM and VCF files

This microhaplot package comes with a small customized sample data drawn from an actual run 
of short read sequencing run on Rockfish species. The sample data
contains sequences of eight genomic loci for four populations of five individuals each, 
with a total of twenty individuals. 

First you need to create a tab-separate **label** file with 3 info columns: path to SAM file name, individual ID, and group label (in this particular order). If you do not want assign any group label for the individuals, you can just leave it as "NA". It is recommended that you have all of the SAM files under one directory to make this labeling task easier.

The `label` file looks like this:
```txt
s6.sam  s6      copper
s11.sam s11     copper
s13.sam s13     gold
s14.sam s14     kelp
s18.sam s18     gold
``` 

Once you have the label file in place, you can run `prepHaplotFiles`, a R function that generates tables of microhaplotype, by providing the following:
 * a label to display in haPLOType
 * path to the directory with all SAM files 
 * path to the `label` file you just created
 * path to the VCF file  
 * optional number of threads (for non-Windows user); recommend 2 * # of processors 
 
```R
library(microhaplot)

# to access package sample case study dataset of rockfish
run.label <- "sebastes"

sam.path <- tempdir()
untar(system.file("extdata",
                  "sebastes_sam.tar.gz",
                  package="microhaplot"),
      exdir = sam.path)
      
label.path <- file.path(sam.path, "label.txt")
vcf.path <- file.path(sam.path, "sebastes.vcf")
out.path <- tempdir()
app.path <- "~/Shiny/microhaplot"

# for your dataset: customize the following paths
# sam.path <- "~/microhaplot/extdata/"
# label.path <- "~/microhaplot/extdata/label.txt"
# vcf.path <- "~/microhaplot/extdata/sebastes.vcf"
# app.path <- "~/Shiny/microhaplot"

haplo.read.tbl <- prepHaplotFiles(run.label = run.label,
                            sam.path = sam.path,
                            out.path = out.path,
                            label.path = label.path,
                            vcf.path = vcf.path,
                            app.path = app.path,
                            n.jobs = 4) # assume running on dual core
                            

runShinyHaplot(app.path)
```


### Suggestions
- SAM files: For pair-ended experiment, both directional reads should be flashed into one.


