# First you need to install and load the 'devtools' package in R.
# Install devtools from CRAN
install.packages("devtools")
# load the library
library("devtools")

# Next you should download the 'lncRna' library from GitHub and load it in your environment
install_github("prodakt/lncRna@test", force = T) # if got error use : options(timeout=400)
#options(timeout=400) 
# If you want you can load all packages by yourself on the beginning

library("devtools")
library("lncRna")
library("rtracklayer")
library(seqinr)
library(Polychrome)
library(plotly)
library(gprofiler2)
library(tidyr)
library(fsmb)
library(patchwork)
# or omit above and use function to simply load all needed packages

load_packages(additional_packages = c("dplyr", "patchwork"), "plotly" )

# You can run ??lncRna to list the available functions in the library.

# %!in% function to run
"%!in%" <- Negate("%in%") 

##########################################################################

# You can read all used files now, or later on tutorial:

#reading gff files
stringtieGTF <- import.gff("data/Mus_musculus.GRCm39.107.gtf.gz")
head(stringtieGTF)
refGTF <- import.gff("data/Mus_musculus.GRCm39.107.gtf.gz")
head(refGTF)

# tab 1 read
tab1 <- read.csv2("data/tab1.csv", header = T, sep = ";" ) #zaczytywany tylko jako porównanie czy ten fragment działa poprawnie 

# tbl2 read
tbl2 <- read.csv2("data/tbl2.csv", header = T, sep = ";" )

#reading pot_lncRNA filtered transcripts which are definately NOT lncRNA
pot_lncRNA <- read.table("data/pot_lncRNA.tab")

#reading cds and nc _tt for coding potential analysis

cds_tt <- readRDS("data/cds_tt.rds")
nc_tt <- readRDS("data/nc_tt.rds")

#In this tutorial, we will be working with gene counts data where row names might contain extra information after a pipe symbol ("|"). To simplify our analysis and ensure clean gene identifiers for downstream steps, we need to **clean the row names of our `genes_counts` data frame**. We will achieve this by removing everything after the first pipe symbol, if present.
genes_counts <- read.table("data/gene_count_matrix.csv", header = T, sep = ",", row.names = 1)
head(genes_counts)
rownames(genes_counts) <- gsub("\\|.*", "", rownames(genes_counts))
head(genes_counts)

#DEGs read
DEGs <- read.csv2("data/DEGs.csv", header = T, row.names = 1)

#DELs read
DELs <- read.csv2("data/DELs.csv", row.names = 1)

##########################################################################
# Tutorial start

# Ia: Annotated features
# The first stage is to read annotation files (reference GTF - preferably downloaded from ENSEMBL database 
# and GTF generated during or after mapping reads - preferably GTF file merged by Stringtie)

# reading GTF files ----------# skip if loaded before ------------ 
stringtieGTF <- import.gff("data/stringtie_merged.gtf")
#refGTF <- import.gff("reference.gtf") # prefered downloaded from ENSEMBL
refGTF <- import.gff("data/Mus_musculus.GRCm39.107.gtf.gz")

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

###########################################################################
# Ib: Expression level

#To filter out low expressed transcripts you can use 'count_matrix' file:
transcripts_counts <- read.table("data/transcript_count_matrix.csv", header = T, sep = ",", row.names = 1)
expr <- rownames(transcripts_counts[rowSums(transcripts_counts) > 10,])

# The final step in tne first stage is to filter out these transcripts which are definitely NOT lncRNA
# remove known protein coding transcripts
pot_lncRNA <- tab1[tab1$transcript_id %!in% known_pcRNA,] 

# remove transcripts shorter than 200 nucleotides and consist of only 1 exon
pot_lncRNA <- pot_lncRNA[pot_lncRNA$exons > "1" & pot_lncRNA$trans_length >= "200",]$transcript_id 

# optionally you can limit tne list of transcripts removing low expressed transcripts
pot_lncRNA <- pot_lncRNA[pot_lncRNA %in% rownames(transcripts_counts[rowSums(transcripts_counts) > 10,])] #transcript_count_matrix.csv otherwise skip

length(pot_lncRNA)

###########################################################################
# loading reference sequences (can be downloaded from tutorial by link)
cds <- read.fasta("data/Mus_musculus.GRCm39.cds.all.fa.gz", seqtype = "DNA", as.string = T, set.attributes = F)
nc <- read.fasta("data/Mus_musculus.GRCm39.ncrna.fa.gz", seqtype = "DNA", as.string = T, set.attributes = F)

# here you can use the list of sequences previously filtered and stored in the "pot_lncRNA" variable
seqs2predict <- pot_lncRNA
# or load assembled merged transcriptome (the whole transcriptome without the first filtering)
seqs2predict <- read.fasta(".../merged_transtriptome.fa", seqtype = "DNA", as.string = T, set.attributes = F)

# splitting into training and test sets can be done by yourself or loaded from test file
set.seed(12345) # to make analysis repeatable 
cds_tt <- test.train.cds(cds.fa = cds, percent_train = 0.6)
head(cds_tt$cds.test)
set.seed(12345)
nc_tt <- test.train.nc(nc.fa = nc, percent_train = 0.6)
head(nc_tt$nc.test)

# preparing sets of sequences for: training, testing and predicting
seqs2predict <- c(seqs2predict, cds[names(cds) %in% cds_tt$cds.test], nc[names(nc) %in% nc_tt$nc.test])
cds2train <- cds[names(cds) %in% cds_tt$cds.train]
nc2train <- nc[names(nc) %in% nc_tt$nc.train]

#zrobić listę przykładowych programów do analizy potencjału kodowania
# w przyszłości odpalanie programów prosto z R-ki

write.fasta(seqs2predict, names(seqs2predict), "seqs2predict.fa", nbchar = 80, as.string = T)
write.fasta(cds2train, names(cds2train), "cds2train.fa", nbchar = 80, as.string = T)
write.fasta(nc2train, names(nc2train), "nc2train.fa", nbchar = 80, as.string = T)

###########################################################################
# IIa:Coding potential

# you can read all results separately
CPC2 <- read.CPC2(".../CPC2.out")
FEELnc <- read.FEELnc(".../feelnc_codpot_out/FEELnc_RF.txt")
PLEK <- read.PLEK(".../PLEK.out")
# and so on ...

# CPAT requires calculation of cutoff value
CPATcutoff <- CPATcutoff(".../CPATout.feature.xls") # 
CPAT <- read.CPAT(CPAT_outfile = ".../CPAT.out", CPAT_cutoff = CPATcutoff)

# or you can read and convert to table all results using one function
tbl2 <- CodPot2tbl(CPC2_outfile = "data/lncCodPot_MMus/CPC2_Mm_lnc.txt.txt",
                   FEELnc_outfile = "data/lncCodPot_MMus/FEELnc_codpot_RF.txt",
                   CNCI_outfile = "data/lncCodPot_MMus/CNCI.index",
                   CPAT_outfile = "data/lncCodPot_MMus/CPAT_Mm_lnc.txt",
                   CPAT_cutoff = 0.64, # the CPAT_cutor is predefined for some organisms (0.64 for mouse) or you can calculate it using 'CPATcutoff()' function
                   PLEK_outfile = "data/lncCodPot_MMus/PLEK_Mm_lnc.txt",
                   LncFinder_outfile = "data/lncCodPot_MMus/LncFinder_results_more5.csv")
#lncRNA_Mdeep_outfile = "data/lncCodPot_MMus/Mm107_out.feature.xls") Not sure that is a correct file

head(tbl2)

venn.CodPot(CodPot = tbl2)
venn.CodPot(CodPot = tbl2, selmet = c(1,1,0,1,1,0))

###########################################################################
# IIb:Best accuracy analysis

# Accuracy, sensitivity and specificity analysis of coding potential prediction tools
# Preparation of statistics
BestPat1 <- SumSingleTools(CodPot.tbl2 = tbl2, nc_test = nc_tt$nc.test, cds_test = cds_tt$cds.test)
head(BestPat1)

# Error analysis for selected tools
selectedTools <- c("CPC2", "PLEK", "FEELnc", "CPAT", "CNCI", "LncFinder")
BP.cpt <- BestTool(BestPat = BestPat1, tools = selectedTools)
BP.cpt 

#Error analysis of combining multiple methods assumes that a given sequence is classified as non-coding if it has been classified as such by at least "n" tools.
# Preparation of statistics for sequences meeting the "at least n tools" criterion.
BestPat2 <- SumAtLeast(BestPat = BestPat1, tools = selectedTools)
head(BestPat2)

# Error analysis for all above criteria
BP.atl <- BestTool.atleast(BestPat = BestPat2)
BP.atl


#Error analysis of all combinations of selected tools (according to the Venn diagram)

# Preparation os statistics for all combinations
BestPat3 <- SumCombTools(BestPat = BestPat2, selectedTools = selectedTools)
BestPat3

# Error analysis for all selected combinations
selectedCombinations <- colnames(BestPat3)[17:ncol(BestPat3)]
BP.cmb <- BestTool.comb(BestPat = BestPat3, selectComb = selectedCombinations)
BP.cmb
BP.cmb[which(BP.cmb[1,] > 0.8)]

# Creating list of Confusion Matrics from selected methods and metrics 
# For all methods and simple run
all_cms <- calculate_cm(bp_cmb_data = BP.cmb, best_pat3_data = BestPat3)

# with show of metric treshold
all_cms <- calculate_cm(bp_cmb_data = BP.cmb, best_pat3_data = BestPat3, print_metric_threshold_methods = T)

# with selected treshold value default = 0,8
all_cms <- calculate_cm(bp_cmb_data = BP.cmb, best_pat3_data = BestPat3, print_metric_threshold_methods = T, threshold = 0.9 )

# with specific metric selected default = Accuracy
all_cms <- calculate_cm(bp_cmb_data = BP.cmb, best_pat3_data = BestPat3, print_metric_threshold_methods = T, metric_to_extract = "Precision")

# return only treshold confition, for Precision
all_high_precision <- calculate_cm(bp_cmb_data = BP.cmb, best_pat3_data = BestPat3, print_metric_threshold_methods = T, metric_to_extract = "Precision", return_only_high_methods = T)

# Statistics of cm 

# Radar plot
# simple radar plot, without method wil prompt console
radar_plot_cm(cm_list = all_high_precision)

# for specified method
radar_plot_cm(cm_list = all_high_precision, methods = c("CPC2+CPAT", "PLEK+CPAT", "CPAT+CNCI"))

# for specified metrics
radar_plot_cm(cm_list = all_high_precision, metrics = c("Sensitivity", "Recall", "Precision", "Accuracy"))

# for multiple plots and area
radar_plot_cm(cm_list = all_high_precision, layout = "multiple", display_area = T)

# without fill
radar_plot_cm(cm_list = all_high_precision, layout = "multiple", display_area = T, display_fill = F)

# radar plot with saving source data frame
radar_plot_cm(cm_list = all_high_precision, layout = "multiple", display_area = T, display_fill = F, save_data = T, file_name = "Test")

#Clock plot
# simple clock plot
clock_plot_cm(cm_list = all_high_precision)

# draw muliple plots
clock_plot_cm(cm_list = all_high_precision, layout = "multiple")

# draw from selected methods
clock_plot_cm(cm_list = all_high_precision, methods = c("CPC2+CPAT", "PLEK+CPAT", "CPAT+CNCI"))

###########################################################################
# III Functional annotation:

# if the coding potential tests were done on the whole transcriptome you should run:
predicted_lncRNA <- tbl2[tbl2$seqIDs %in% pot_lncRNA,]

# non-coding potential prooved by at least 5 methods
predicted_lncRNA <- predicted_lncRNA[rowSums(predicted_lncRNA[2:ncol(predicted_lncRNA)]) >= 5,] 

# filter out pfam records
#predicted_lncRNA <- predicted_lncRNA[predicted_lncRNA$seqIDs %!in% pfam,] 

# extract only the ID list
predicted_lncRNA <- predicted_lncRNA$seqIDs
head(predicted_lncRNA)

###########################################################################
# In this stage is to predict or estimate some functions and functional connections of identified lncRNAs

# combine the list of both predicted and known lncRNA's
lncRNA_transcripts <- unique(c(predicted_lncRNA, known_lncRNA))
head(lncRNA_transcripts)

#combined list od known lncRNA genes
lncRNA_genes <- unique(stringtieGTF[stringtieGTF$transcript_id %in% lncRNA_transcripts,]$gene_id)
head(lncRNA_genes)

#combined list of known protein coding genes
pcRNA_genes <- unique(refGTF[refGTF$gene_biotype %in% "protein_coding",]$gene_id)
head(pcRNA_genes)

###########################################################################
# Creating list of differantialy expressed protein coding genes that showed cis and trans interactions with lncRNA

#Cis
cis_table <- cisInter(is.best = T, FEELnc.classes = "data/lncCodPot_MMus/FEELnc_Mm_classes.txt", lncRNAs = lncRNA_transcripts, lncRNA.level = "transcript", mRNA.level = "gene", max.dist = 10000, mRNAs = rownames(DEGs))  

#Trans !!!!!!!!!  In this tutorial, we will be working with gene counts data where row names might contain extra information after a pipe symbol ("|"). To simplify our analysis and ensure clean gene identifiers for downstream steps, we need to **clean the row names of our `genes_counts` data frame**. We will achieve this by removing everything after the first pipe symbol, if present.
trans_table <- TransAct(expr.matrix = genes_counts, rval = 0.95, pval = 0.05, lncRNA.list = row.names(DELs), tarRNA.list = rownames(DEGs)) #see all files on top for genes_counts

# Having a list of cis and trans you can now run gProfiler analysis to functionally annotate genes

cis_gostres <- gost(query = cis_table$partnerRNA_gene, 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)

# You can find  an explanation of all the options of the above function through the command "?gost"
#We need to access the result list

cis_gostres <- cis_gostres$result
head(cis_gostres)

# You need to do the same for the second list of trans interacting genes

trans_gostres <- gost(query = trans_table$targetRNA.id, 
                      organism = "mmusculus", ordered_query = FALSE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = TRUE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = NULL, as_short_link = FALSE)

trans_gostres <- trans_gostres$result
head(trans_gostres)

# lncRNA-mRNA
LncTar_interact <- read.table(file = "data/DELs_vs_DEGs_LncTar.txt", header = T)

LncTar_gostres <- gost(query = LncTar_interact$Target, 
                       organism = "mmusculus", ordered_query = FALSE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = TRUE, 
                       user_threshold = 0.05, correction_method = "g_SCS", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE)

LncTar_gostres <- LncTar_gostres$result
head(trans_gostres)

#
LncTar <- interactions_merge(gprof = LncTar_gostres, interaction_table = LncTar_interact, type = "LncTar")

# LION interactions lncRNA-protein
LION <- read.csv2(file = "data/LION_results_part.csv", header = T)

nrow(LION[LION$RPISeq_retrain_pred %in% "Interact", ])
nrow(LION[LION$rpiCOOL_retrain_pred %in% "Interact", ])
nrow(LION[LION$LION_pred %in% "Interact", ])

LION[LION$LION_pred %in% "Interact", ]$Pro_Name

LIONGOin1 <- LION[LION$LION_pred %in% "Interact", ]$Pro_Name
LIONGOin2 <- LION[LION$rpiCOOL_retrain_pred %in% "Interact", ]$Pro_Name
LIONGOin3 <- LION[LION$RPISeq_retrain_pred %in% "Interact", ]$Pro_Name
  
LION_gostres <- gost(query = c(LIONGOin1,LIONGOin2,LIONGOin3), 
                     organism = "mmusculus", ordered_query = FALSE, 
                     multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                     measure_underrepresentation = FALSE, evcodes = TRUE, 
                     user_threshold = 0.05, correction_method = "g_SCS", 
                     domain_scope = "annotated", custom_bg = NULL, 
                     numeric_ns = "", sources = NULL, as_short_link = FALSE)

LION_gostres <- LION_gostres$result
head(LION_gostres)

# Process interactions
Test <- process_interactions(gprof = trans_gostres, interaction_table = trans_table, type = "banana")

# Create a table of cis interactions of protein coding genes with lncRNAs 
cis <- cis_interactions(cis_gprof = cis_gostres, cis_table = cis_table)
head(cis)

# Create a table of trans interactions of protein coding genes with lncRNAs 
trans <- trans_interactions(trans_gprof = trans_gostres,trans_table = trans_table)
head(trans)

# Combine tables
combined_table <- rbind(trans,cis)

# You can use different functions for creating plot
#plot_by_action
plot_by_action(data = combined_table, cis = T, trans = F, label = T)

#plot_by_lnc
plot_by_lnc(data = combined_table, select_lnc = c("ENSMUSG00000106858", "ENSMUSG00000002769"), label = T)

#plot_by_target
plot_by_target(data = combined_table, select_target = c("ENSMUSG00000000731","ENSMUSG00000000732"), label = T)

#plot_by_terms
plot_by_terms(data = combined_table, select_terms = "response to stress", label = T)
plot_by_terms(data = Test, label = T)


plot_by_type(data = Test, label = T)
plot_by_type(data = Test, type = "cis", label = T)
plot_by_type(data = Test, type = c("cis", "LncTar"), label = T)
