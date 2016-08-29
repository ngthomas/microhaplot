# haplot   

The goal of `haplot` is to generate visual summary of microhaplotype found in short read alignments.

The process of using `haplot` is quick and straight-forward. It takes two function calls to extract, process and display haplotype, of which can be completed within minutes. `haplot` is suitable for carrying quick assesement and quality control of haplotype generated from library runs. Plot summaries include read depth, fraction of calleable haplotype, Hardy-Weinberg equilibrium plot, and more.   

### Installation

You will need [devtools](https://github.com/hadley/devtools) to install `haplot`. You can get `devtools` by `install.packages("devtools")`.

Once you have `devtools` available in R, you can get `haplot` by this command:
```r
# sudo R
devtools::install_github("ngthomas/callBayes/haplot")

haplot::mvHaplotype("~/bin/haPLOType") #provide a directory path to host haPLOType app
```


### Quick Guide to use Haplot

To upload your alignment files to shiny App `haPLOType`, you will need to generate a tab-separate **label** file with 3 info columns: path to SAM file name, individual ID, and group label (in this particular order). 

If you do not want assign any group label for the individuals, you can just leave it as "NA". 

NOTE: It is recommended that you have all of the SAM files under one directory to make this labeling task easier.

An example of the `label` file:
```txt
satro_flashed_s1_aln.sam        s1      black
satro_flashed_s2_aln.sam        s2      black
satro_flashed_s3_aln.sam        s3      black
satro_flashed_s4_aln.sam        s4      black
satro_flashed_s5_aln.sam        s5      copper
``` 
  
  
Now you can proceed with running `runHaplot`. You will need to provide:

 * a label 
 * path to the directory with all SAM files 
 * path to the `label` file you just created
 * path to the VCF file  
  
  
```R
library(haplot)

# ---- edit ---------
run.label <- "example 1"
sam.path <- "data/satro_sample"
label.path <- "data/satro_sample/sample_label.txt"
vcf.path <- "data/satro_sample/sebastes.vcf"
app.path <- "~/bin/haPLOType" 
# -------------------------

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

