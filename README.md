# lncRna <img src="img/lncRna_logo.png" align="right" height = 150/>

![GitHub](https://img.shields.io/github/license/prodakt/lncRna?style=plastic)
![GitHub top language](https://img.shields.io/github/languages/top/prodakt/lncRna?style=plastic)
![GitHub R package version](https://img.shields.io/github/r-package/v/prodakt/lncRna?color=gree&label=ver&style=plastic)

![GitHub repo file count](https://img.shields.io/github/directory-file-count/prodakt/lncRna?color=purple&style=plastic)
![GitHub repo size](https://img.shields.io/github/repo-size/prodakt/lncRna?color=yellow&label=size)
![GitHub commit activity](https://img.shields.io/github/commit-activity/y/prodakt/lncRna?color=white&label=activ)
![GitHub last commit](https://img.shields.io/github/last-commit/prodakt/lncRna?color=%23FF0000%20&label=last)


## Introduction
The lncRna is an R package designed to simplify the process of lncRNA identification, it allows researchers to set up and perform the analysis in a straightforward manner. This package provides a user-friendly interface for processing transcriptomic data and implementing the lncRNA identification pipeline. Moreover, it offers the capability to estimate the accuracy and precision of the prepared pipelines, facilitating the selection of the most reliable approach. With its simplicity and comprehensive functionality, this package aims to enhance the efficiency and accuracy of lncRNA identification, ultimately contributing to advancements in the field of non-coding RNA research.


## General description
We are pleased to introduce our innovative package that simplifies the process of lncRNA identification in just a few easy steps. Our package is designed to provide researchers with a user-friendly and efficient tool to carry out the lncRNA identification, while also allowing for accuracy and precision validation to select the most appropriate pipeline for their specific needs.

With the package, researchers can streamline their lncRNA identification process, saving valuable time and resources. The package offers a comprehensive set of tools and algorithms that have been specifically developed to accurately identify lncRNAs from large-scale transcriptomic data. By following a series of simple steps, researchers can easily preprocess their data, perform the identification analysis, and obtain significant results.

One of the key features of our package is ability to validate the accuracy and precision of the prepared pipelines. Researchers can evaluate the performance of different pipelines using various metrics and statistical analyses. This allows for an objective assessment of the results and facilitates the selection of the most appropriate pipeline for specific research objectives.

This package stands out for its user-friendly interface, which ensures that researchers, regardless of their computational background, can easily navigate through the process. The intuitive design and clear instructions make it accessible to both beginners and experienced users, enhancing the acknowledgement of the package across the scientific community.

In summary, lncRna package offers a streamlined, efficient, and accurate solution for researchers working in the field. Its simple step-by-step process, combined with the ability to validate and select the most suitable pipeline, ensures that researchers can confidently and successfully identify lncRNAs from their transcriptomic data. We believe that our package will greatly contribute to advancing lncRNA research and foster new discoveries in this exciting field.

In general the pipeline consist of a few stages:
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
  <li>under construction</li>
  <!---
  <li>detection of cis and trans acting genes,</li>
  <li>functional analyses</li>
-->
</ul>
<ul>  
  <b>IV. Structure</b>
  <li>under construction</li>
</ul>


## Core functionality

Our package provides a comprehensive set of functionalities to simplify the lncRNA identification process and enable researchers to select the most accurate pipeline for predicting lncRNAs. The key features of the package include:

- **Easy Setup of Individual Steps**:<br/> 
The package allows researchers to set up individual steps of the lncRNA identification process into pipelines in a simple and intuitive way. Users can easily define and customize each step, such as data preprocessing, feature extraction, and coding potential analysis, based on their specific requirements.
- **Combination of Results**: <br/> 
The package enables merging of the results from coding potential analyses in all possible combinations. By integrating multiple algorithms or strategies for coding potential assessment, researchers can leverage the strength of different approaches and increase the accuracy and reliability of lncRNA predictions.
- **Statistical Analysis**:<br/> 
Our package incorporates advanced statistical analysis capabilities to compute key metrics such as accuracy and precision. These statistics provide quantitative measures to assess the performance of different pipelines and combinations of steps. Researchers can objectively evaluate the effectiveness of each pipeline and select the most appropriate combination for their specific research goals.
- **Most Accurate Pipeline Execution**:</br> 
Based on the statistical analysis results, the package allows researchers to run the most accurate pipeline for predicting lncRNAs in the most reliable manner. By considering the performance metrics obtained during the evaluation phase, researchers can confidently apply the selected pipeline to their transcriptomic data and obtain highly accurate predictions of lncRNA transcripts.

By providing a user-friendly interface and comprehensive functionalities, our package empowers researchers to streamline the lncRNA identification process, optimize pipeline configurations, and achieve highly accurate predictions. This package serves as a valuable tool for advancing non-coding RNA research and facilitating new discoveries in the field.

-------------------------------------------------------------
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
-------------------------------------------------------------
## Rapid pipeline (first fast run)

1. First you need to install/load the '[rtracklayer](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)' library to import any GTF and GFF files in a acurate format.
```
install.packages("rtracklayer")
library("rtracklayer")
```
2. Run `??lncRna` to list the available functions in the library.
3. Run the functions of 'lncRna' library:

### Ia. Annotated features 
The first stage is to read annotation files (reference GTF - preferably downloaded from ENSEMBL database and GTF generated during or after mapping reads - preferably GTF file merged by Stringtie)
- [Used reference GTF annotation file](https://ftp.ensembl.org/pub/release-107/gtf/mus_musculus/Mus_musculus.GRCm39.107.gtf.gz)
```
# reading GTF files
stringtieGTF <- import.gff("stringtie/stringtie_merged.gtf")
refGTF <- import.gff("reference.gtf") # prefered downloaded from ENSEMBL
# refGTF <- import.gff("Mus_musculus.GRCm39.107.gtf")

# extraction of known lncRNA
known_biotypes <- refBiotypes(refGTF = refGTF)
# what kind and quantity of biotypes do we have in the reference genome
table(known_biotypes$transcript_biotype) 
# the list of known lncRNA transcripts ID 
known_lncRNA <- known_biotypes[known_biotypes$transcript_biotype %in% "lncRNA",]$transcript_id # the list of lncRNA ID's
#
# extraction of known protein coding transcripts
known_pcRNA <- known_biotypes[known_biotypes$transcript_biotype %in% "protein_coding",]$transcript_id # the list of protein coding RNA ID's

# extracting from GTF merged with stringtie the basic features for lncRNA identification
tab1 <- strGTF2stat(stringtieGTF = stringtieGTF)
head(tab1)

       transcript_id exons trans_length
1 ENSMUST00000000001     9         3262
2 ENSMUST00000000003     7          902
3 ENSMUST00000000010     2         2574
4 ENSMUST00000000028     9         2158
5 ENSMUST00000000033     4         3708
6 ENSMUST00000000049     8         1190
```

### Ib. Expression level
To filter out low expressed transcripts you can use 'count_matrix' file:

```
transcripts_counts <- read.table("../transcript_count_matrix.csv", header = T, sep = ",", row.names = 1)
expr <- rownames(transcripts_counts[rowSums(transcripts_counts) > 10,])
```

You can also use the 'strtie2expr()' function, but this function only works if you have used the count option from [Stringtie](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

```
strtie2expr(strdir = "stringtie_folder/", expunit = "FPKM" or "TPM")
```

### Ic. Selection of transcripts for coding potential analyses
The final step in tne first stage is to filter out these transcripts which are definitely NOT lncRNA.

```
# remove known protein coding transcripts
pot_lncRNA <- tab1[tab1$transcript_id %!in% known_pcRNA,] 

# remove transcripts shorter than 200 nucleotides and consist of only 1 exon
pot_lncRNA <- pot_lncRNA[pot_lncRNA$exons > "1" & pot_lncRNA$trans_length >= "200",]$transcript_id 

# optionally you can limit tne list of transcripts removing low expressed transcripts
pot_lncRNA <- pot_lncRNA[pot_lncRNA %in% rownames(transcripts_counts[rowSums(transcripts_counts) > 10,])]

length(pot_lncRNA)
[1] 27019
```

The list of ENSEMBL ID from 'pot_lncRNA' is avaliable [here - the "pot_lncRNA.tab.tar.gz" file in the "test/mouse" folder](https://github.com/prodakt/lncRna/blob/main/test/mouse/pot_lncRNA.tab.tar.gz).

In addition, you can prepare a set of sequences for accuracy analysis of methods used for predicting coding potential (see IIb). For this purpose, non-coding, coding reference sequences should be loaded and divided into a set for training, testing and proper prediction of coding potential.

- [Used reference CDS sequences](https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz)
- [Used reference nc sequences](https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz)

```
# loading reference sequences
cds <- read.fasta("Mus_musculus.GRCm39.cds.all.fa", seqtype = "DNA", as.string = T, set.attributes = F)
nc <- read.fasta("Mus_musculus.GRCm39.ncrna.fa", seqtype = "DNA", as.string = T, set.attributes = F)

# Preparation of the nucleotide sequences files for external programs for the prediction of coding potential
# Here you are preparing 5 groups of sequences in 3 files:
# - known non coding sequences to train the models in external tools (the "nc2train.fa" file),
# - known protein coding sequences to train the models in external tools (the "cds2train.fa" file),
# - known non coding sequences to test the results (part of the "seqs2predict.fa" file),
# - known protein coding sequences to test the results (part of the "seqs2predict.fa" file),
# - unknown/potentially non-coding sequences to predict the coding potential (part of the "seqs2predict.fa" file).

# here you can use the list of sequences previously filtered and stored in the "pot_lncRNA" variable
seqs2predict <- pot_lncRNA
# or load assembled merged transcriptome (the whole transcriptome without the first filtering)
seqs2predict <- read.fasta("lncRNA/merged_transtriptome.fa", seqtype = "DNA", as.string = T, set.attributes = F)

# splitting into training and test sets
set.seed(12345) # to make analysis repeatable
cds_tt <- test.train.cds(cds.fa = cds, percent_train = 0.6)
nc_tt <- test.train.nc(nc.fa = nc, percent_train = 0.6)


# preparing sets of sequences for: training, testing and predicting
seqs2predict <- c(seqs2predict, cds[names(cds) %in% cds_tt$cds.test], nc[names(nc) %in% nc_tt$nc.test])
cds2train <- cds[names(cds) %in% cds_tt$cds.train]
nc2train <- nc[names(nc) %in% nc_tt$nc.train]

write.fasta(seqs2predict, names(seqs2predict), "seqs2predict.fa", nbchar = 80, as.string = T)
write.fasta(cds2train, names(cds2train), "cds2train.fa", nbchar = 80, as.string = T)
write.fasta(nc2train, names(nc2train), "nc2train.fa", nbchar = 80, as.string = T)
```

To get the 'cds_tt' and 'nc_tt' variables prepared as part of this tutorial, you can use the 'readRDS()' function and load the shared '[cds_tt.rds](https://github.com/prodakt/lncRna/blob/main/test/mouse/cds_tt.rds)' and '[nc_tt.rds](https://github.com/prodakt/lncRna/blob/main/test/mouse/nc_tt.rds)' files.

```
cds_tt <- readRDS("cds_tt.rds")
nc_tt <- readRDS("nc_tt.rds")
```

### IIa. Coding potential
The second stage is reading the output files generated by coding potential tools:

```
# you can read all results separately
CPC2 <- read.CPC2(".../CPC2.out")
FEELnc <- read.FEELnc(".../feelnc_codpot_out/FEELnc_RF.txt")
PLEK <- read.PLEK(".../PLEK.out")
# and so on ...

# CPAT requires calculation of cutoff value
CPATcutoff <- CPATcutoff(".../CPATout.feature.xls") # 
CPAT <- read.CPAT(CPAT_outfile = ".../CPAT.out", CPAT_cutoff = CPATcutoff)

# or you can read and convert to table all results using one function
tbl2 <- CodPot2tbl(CPC2_outfile = ".../CPC2.out",
                   FEELnc_outfile = ".../feelnc_codpot_out/FEELnc_RF.txt",
                   CPAT_outfile = ".../CPAT.out",
                   CPAT_cutoff = 0.64, # the CPAT_cutor is predefined for some organisms (0.64 for mouse) or you can calculate it using 'CPATcutoff()' function
                   PLEK_outfile = ".../PLEK.out",
                   LncFinder_outfile = ".../LncFinder_results.csv",
                   lncRNA_Mdeep_outfile = ".../Mdeep.out")
                   

head(tbl2)


              seqIDs CPC2 PLEK FEELnc CPAT CNCI LncFinder
1 ENSMUST00000193812    1    1      0    1    1         0
2 ENSMUST00000082908    1    0      0    1    1         0
3 ENSMUST00000192857    1    1      0    1    1         1
4 ENSMUST00000195335    1    1      0    1    1         0
5 ENSMUST00000192336    1    1      0    1    1         1
6 ENSMUST00000194099    1    0      0    1    1         1
```
It is not necessary to read all mentioned files. If you don't have some of the results you can simply omit arguments, i.e.:
```
tbl2 <- CodPot2tbl(CPC2_outfile = ".../CPC2.out", PLEK_outfile = ".../PLEK.out")
```

Compressed folders containing selected results of the coding potential analysis are available [here. (lncRna/test/mouse/lncCodPot_MMus.rar and lncRna/test/mouse/lncCodPot_MMus.tar.bz2)](https://github.com/prodakt/lncRna/tree/main/test/mouse).

You can easily draw the venn diagram for all used methods:
```
venn.CodPot(CodPot = tbl2)
```
![venn_all](https://github.com/prodakt/lncRna/blob/main/img/vennFull.png)

or you can selet which results you wish to include in the venn diagram:
```
venn.CodPot(CodPot = tbl2, selmet = c(1,1,0,1,1,0))
```
![venn_4](https://github.com/prodakt/lncRna/blob/main/img/vennSel.png)


### IIb. Best accuracy analysis
At this stage, you can do an accuracy analysis of all the tools used to analyze the coding potential. In addition, you can combine the results of individual methods (according to the Venn diagram) or choose another way to combine the results and draw consensus information. For all combinations, you can evaluate accuracy, precision and other statistics to decide which path is the most optimal for analyzing your data. You can generate a confusion matrix (numbers of false positives, false negatives, true positives and true negatives) and generate charts to help you make decisions.

#### Accuracy, sensitivity and specificity analysis of coding potential prediction tools
```
# Preparation of statistics
BestPat1 <- SumSingleTools(CodPot.tbl2 = tbl2, nc_test = nc_tt$nc.test, cds_test = cds_tt$cds.test)
head(BestPat1)

                    seqIDs CPC2 PLEK FEELnc CPAT CNCI LncFinder type isNC sums
79712 ENSMUST00000179520.2    1    0      0    1    1         0  cds    0    3
79713 ENSMUST00000179883.2    1    0      0    1    1         0  cds    0    3
79717 ENSMUST00000177646.2    1    0      0    1    1         0  cds    0    3
79718 ENSMUST00000178230.2    1    0      0    1    1         0  cds    0    3
79723 ENSMUST00000180266.2    1    0      0    1    1         0  cds    0    3
79725 ENSMUST00000103270.4    1    1      0    0    0         0  cds    0    2



# Error analysis for selected tools
selectedTools <- c("CPC2", "PLEK", "FEELnc", "CPAT", "CNCI", "LncFinder")
BP.cpt <- BestTool(BestPat = BestPat1, tools = selectedTools)
BP.cpt

                 CPC2    PLEK FEELnc   CPAT   CNCI LncFinder
Accuracy       0.6692  0.4207 0.5247 0.9096 0.7761    0.5285
Kappa          0.3521 -0.1412 0.0000 0.8187 0.5575    0.0170
AccuracyLower  0.6600  0.4110 0.5149 0.9039 0.7678    0.5187
AccuracyUpper  0.6784  0.4303 0.5345 0.9151 0.7842    0.5382
AccuracyNull   0.5247  0.5247 0.5247 0.5247 0.5247    0.5247
AccuracyPValue 0.0000  1.0000 0.5040 0.0000 0.0000    0.2275
McnemarPValue  0.0000  0.0000 0.0000 0.3541 0.0000    0.0000

```

#### Error analysis of combining multiple methods assumes that a given sequence is classified as non-coding if it has been classified as such by at least "n" tools.
```
# Preparation of statistics for sequences meeting the "at least n tools" criterion.
BestPat2 <- SumAtLeast(BestPat = BestPat1, tools = selectedTools)
head(BestPat2)

                    seqIDs CPC2 PLEK FEELnc CPAT CNCI LncFinder type isNC sums atl1 atl2 atl3 atl4 atl5 atl6
79712 ENSMUST00000179520.2    1    0      0    1    1         0  cds    0    3    1    1    1    0    0    0
79713 ENSMUST00000179883.2    1    0      0    1    1         0  cds    0    3    1    1    1    0    0    0
79717 ENSMUST00000177646.2    1    0      0    1    1         0  cds    0    3    1    1    1    0    0    0
79718 ENSMUST00000178230.2    1    0      0    1    1         0  cds    0    3    1    1    1    0    0    0
79723 ENSMUST00000180266.2    1    0      0    1    1         0  cds    0    3    1    1    1    0    0    0
79725 ENSMUST00000103270.4    1    1      0    0    0         0  cds    0    2    1    1    0    0    0    0


# Error analysis for all above criteria
BP.atl <- BestTool.atleast(BestPat = BestPat2)
BP.atl

                 atl1   atl2   atl3   atl4   atl5   atl6
Accuracy       0.5247 0.6727 0.8415 0.7787 0.5358 0.5247
Kappa          0.0000 0.3606 0.6839 0.5472 0.0245 0.0000
AccuracyLower  0.5149 0.6635 0.8342 0.7705 0.5260 0.5149
AccuracyUpper  0.5345 0.6819 0.8486 0.7868 0.5456 0.5345
AccuracyNull   0.5247 0.5247 0.5247 0.5247 0.5247 0.5247
AccuracyPValue 0.5040 0.0000 0.0000 0.0000 0.0131 0.5040
McnemarPValue  0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
```

#### Error analysis of all combinations of selected tools (according to the Venn diagram)
```
# Preparation os statistics for all combinations
BestPat3 <- SumCombTools(BestPat = BestPat2, selectedTools = selectedTools)

# Error analysis for all selected combinations
selectedCombinations <- colnames(BestPat3)[17:ncol(BestPat3)]
BP.cmb <- BestTool.comb(BestPat = BestPat3, selectComb = selectedCombinations)
BP.cmb
BP.cmb[which(BP.cmb[1,] > 0.8)]

               CPC2+CPAT CPC2+CNCI CPAT+CNCI CPC2+CPAT+CNCI
Accuracy          0.8964    0.8040    0.9269         0.9143
Kappa             0.7920    0.6104    0.8530         0.8272
AccuracyLower     0.8903    0.7961    0.9217         0.9086
AccuracyUpper     0.9023    0.8117    0.9319         0.9197
AccuracyNull      0.5247    0.5247    0.5247         0.5247
AccuracyPValue    0.0000    0.0000    0.0000         0.0000
McnemarPValue     0.0000    0.0000    0.0000         0.0000

```
The functions to perform error analysis (BestTool.comb(), BestTool.atleast() and BestTool()) were coded based on the confisionmatrix() function of the ['CARET' library](https://cran.r-project.org/web/packages/caret). Therefore, the output of the created functions is compatible with the format required by 'CARET' library functions. The arguments of this function and the tests used are described in detail in the [CARET library manual](https://cran.r-project.org/web/packages/caret/caret.pdf) and the individual functions:

"The overall accuracy and unweighted Kappa statistic are calculated. The p-value from the McNemar test is also calculated using mcnemar.test (which can generate NA values with sparse tables).
The overall accuracy rate is calculated along with a 95 percent confidence interval for that rate (using binom.test) and a one-sided test to see if the accuracy is better than the "missing information rate," which is taken as the largest percentage of class in the data. "

As one example of this use of the function, you can also select and display or plot more detailed data of the selected result
```
selected <- "CPAT+CNCI"
cm <- caret::confusionMatrix(data = as.factor(BestPat3[,selected]),
                      reference = as.factor(BestPat3$isNC),
                      positive = "1")
cm

Confusion Matrix and Statistics

          Reference
Prediction    0    1
         0 5123  561
         1  177 4240
                                          
               Accuracy : 0.9269          
                 95% CI : (0.9217, 0.9319)
    No Information Rate : 0.5247          
    P-Value [Acc > NIR] : < 2.2e-16       
                                          
                  Kappa : 0.853           
                                          
 Mcnemar's Test P-Value : < 2.2e-16       
                                          
              Precision : 0.9599          
                 Recall : 0.8831          
                     F1 : 0.9199          
             Prevalence : 0.4753          
         Detection Rate : 0.4198          
   Detection Prevalence : 0.4373          
      Balanced Accuracy : 0.9249
      

```

### IIc. Filtering by simmilarity
in this stage you can read the Pfam scanning output file to filter out transcripts containing protein domains
```
# pfam
pfam <- read.pfam(pfam_outfile = ".../PFAM.out", eval_cutoff = 0.001)
head(pfam)
```
The 'read.rfam()' function allows you to read Rfam output to tabular format.
```
rfam <- read.rfam(rfam_outfile = ".../Rfam.out")
```


### Summary of stages I and II
The final list of predicted lncRNA should contain transcripts which:
- are not protein coding
- are longer than 200 nucleotides (but you can change it)
- passed the potential coding tests (it is recommended to run a few independent methods to strengthen this step)

Additionally it is recommended:
- to select transcripts consisting of more than 1 exon
- to filter out low expressed transcripts
- filter out transcripts which have significant similarity with protein domains.


<!--
```
# if the coding potential tests were done on the whole transcriptome you should run:
predicted_lncRNA <- tbl2[tbl2$seqIDs %in% pot_lncRNA,]

# non-coding potential prooved by at least 5 methods
predicted_lncRNA <- predicted_lncRNA[rowSums(predicted_lncRNA[2:ncol(predicted_lncRNA)]) >= 5,] 

# filter out pfam records
predicted_lncRNA <- predicted_lncRNA[predicted_lncRNA$seqIDs %!in% pfam,]

# extract only the ID list
predicted_lncRNA <- predicted_lncRNA$seqIDs
```
-->


### III. Functional annotation
<i>under construction</i>
 <!--
the third stage is to predict or estimate some functions and functional connections of identified lncRNAs
```
# combine tle list of both predicted and known lncRNA's
lncRNA_transcripts <- unique(c(predicted_lncRNA, known_lncRNA))

# cis acting genes
cis <- cisInter(is.best = T, FEELnc.classes = ".../FEELnc_classes.txt", lncRNAs = lncRNA_transcripts, lncRNA.level = "transcript", mRNA.level = "gene", max.dist = 10000)  
head(cis)

     isBest lncRNA_gene lncRNA_transcript    partnerRNA_gene partnerRNA_transcript direction       type distance     subtype location
276       1 MSTRG.13837     MSTRG.13837.9 ENSMGAG00000000862    ENSMGAT00000000882 antisense intergenic    13338   divergent upstream
433       1 MSTRG.37986     MSTRG.37986.9 ENSMGAG00000008533    ENSMGAT00000009584     sense intergenic    51659 same_strand upstream
1112      1 MSTRG.37986    MSTRG.37986.14 ENSMGAG00000008533    ENSMGAT00000009584     sense intergenic    74766 same_strand upstream
1430      1 MSTRG.26440     MSTRG.26440.4 ENSMGAG00000002759    ENSMGAT00000003192     sense      genic        0      nested intronic
1518      1 MSTRG.43740     MSTRG.43740.4 ENSMGAG00000015916    ENSMGAT00000017829     sense      genic        0  containing intronic
1561      1 MSTRG.37986    MSTRG.37986.11 ENSMGAG00000008533    ENSMGAT00000009584     sense intergenic    51659 same_strand upstream

# trans acting genes
trans <- TransAct(expr.matrix = expr, rval = 0.9, pval = 0.05, lncRNA.list = lncRNA_genes, tarRNA.list = pcRNA_genes)
head(trans)

      lncRNA.id targetRNA.id   r.value      p.value
193 MSTRG.14959  MSTRG.10009 0.9129956 3.387791e-05
333 MSTRG.18834  MSTRG.10009 0.9011705 6.277291e-05
397 MSTRG.20790  MSTRG.10009 0.9126770 3.448375e-05
591 MSTRG.26440  MSTRG.10009 0.9556676 1.251677e-06
690 MSTRG.28330  MSTRG.10009 0.9243804 1.713508e-05
819 MSTRG.31990  MSTRG.10009 0.9151334 3.002583e-05


```
-->

### IV. Structural analysis
<i>under construction</i>
