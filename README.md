# lncRna <img src="img/lncRna_logo_small.png" align="right" height = 150/>

is the R package for lncRNA identification i a simple way.

This package allows to carry out the lncRNA identification pipeline in a few very simple steps. In general the pipieline consist of a few stages:
<ul>
  <b>I. Annotated features:</b>
  <li>extraction of known lncRNA from reference GTF,</li>
  <li>filtering potential (unknown) lncRNA based on the features such as sequence length, exon number,</li>
</ul>
<ul>  
  <b>II. Coding potential</b>
  <li>lncRNA classification based on the coding potential</li>
  <li>filtering transcripts basing on the simmilarity with known nucleotide and peptide sequences (homologues in databases),</li>
</ul>
<ul>  
  <b>III. Functions</b>
  <li>detection of cis and trans acting genes,</li>
  <li>functional analyses</li>
</ul>
<ul>  
  <b>IV. Structure</b>
  <li>under construction</li>
</ul>


## How to install
First you need to install and load the '[devtools](https://github.com/r-lib/devtools)' package in R. 
```
# Install devtools from CRAN
install.packages("devtools")
# load the library
library(devtools)
```
Next you should download the 'lncRna' library from GitHub and load it in your environment:
```
install_github("prodakt/lncRna")
library("lncRna")
```

## rapid pipeline (first fast run)

1. First you need to install/load the '[rtracklayer](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)' library to import any GTF and GFF files in a acurate format.
```
install.packages("rtracklayer")
library("rtracklayer")
```
2. Run `??lncRna` to list the available functions in the library.
3. Run the functions of 'lncRna' library:

I. the first stage is to read annotation files (reference GTF - prefered downloaded from ENSEMBL database and GTF generated during or after mapping reads - prefered GTF file merged by Stringtie)
```
# reading GTF files
stringtieGTF <- import.gff("stringtie/stringtie_merged.gtf")
refGTF <- import.gff("reference.gtf") # prefered downloaded from ENSEMBL

# extracting from GTF merged with stringtie the basic features for lncRNA identification
tab1 <- strGTF2stat(stringtieGTF = stringtieGTF)
head(tab1)
pot_lncRNA <- tab1[tab1$exons >1 & tab1$trans_length > 200,]$transcript_id 
head(pot_lncRNA)


# extraction of known lncRNA
known_lncRNA <- refBiotypes(refGTF = refGTF)
# what kind and quantity of biotypes do we have in the reference genome
table(known_lncRNA$transcript_biotype) 
# the list of known lncRNA transcripts ID 
known_lncRNA <- known_lncRNA[known_lncRNA$transcript_biotype %in% "lncRNA",]$transcript_id # lista id lncRNA
head(known_lncRNA)
```

IIa. the second stage is reading the output files generated by codding potential tools:
```
# you can read all results separately
CPC2 <- read.CPC2(".../CPC2.out")
FEELnc <- read.FEELnc(".../feelnc_codpot_out/FEELnc_RF.txt")

CPATcutoff <- CPATcutoff(".../CPATout.feature.xls") # 
CPAT <- read.CPAT(CPAT_outfile = ".../CPAT.out", CPAT_cutoff = CPATcutoff)

PLEK <- read.PLEK(".../PLEK.out")

# or we can read and convert to table all results using one function
tbl2 <- CodPot2tbl(CPC2_outfile = ".../CPC2.out",
                   FEELnc_outfile = ".../feelnc_codpot_out/FEELnc_RF.txt",
                   CPAT_outfile = ".../CPAT.out",
                   CPAT_cutoff = 0.78, # the CPAT_cutor is predefined for some organisms (0.78 for human) or you can calculate it using 'CPATcutoff()' function
                   PLEK_outfile = ".../PLEK.out")

head(tbl2)
```
IIb. in this stage you can read the Pfam scanning output file to filter out transcripts containing protein domains
```
# pfam
pfam <- read.pfam(pfam_outfile = ".../PFAM.out", eval_cutoff = 0.001)
head(pfam)
```

III. the third stage is to predict or estimate some functions and functional connections of identified lncRNAs
```
# cis acting genes
lncRNA_transcripts <- tbl2[SELECT_ROWS,]$seqIDs
cis <- filter.CisAct(is.best = T, FEELnc.classes = ".../FEELnc_classes.txt", lncRNAs = lncRNA_transcripts)  

```
