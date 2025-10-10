#' A CodPot2tbl function
#'
#' This function reads several outputs from coding potential programs and join all tesults in table, where 0 is coding and 1 - noncoding
#' @param CPC2_outfile the output file from the CPC2 tool
#' @param PLEK_outfile the output file from the PLEK tool
#' @param FEELnc_outfile the output file from the FEELnc tool
#' @param CPAT_outfile the output file from the CPAT tool
#' @param CPAT_cutoff the cutoff value for CPAT output
#' @param CNCI_outfile the output file from the CNCI tool
#' @param LncFinder_outfile the output file from the LncFinder tool
#' @param lncRNA_Mdeep_outfile the output file from the lncRNA_Mdeep tool
#' @keywords coding potential lncRNA
#' @export
#' @examples
#' CodPot2tbl()
#'
CodPot2tbl <- function(CPC2_outfile = NULL, PLEK_outfile = NULL, FEELnc_outfile = NULL, CPAT_outfile = NULL, CPAT_cutoff = 0.364, CNCI_outfile = NULL, LncFinder_outfile = NULL, lncRNA_Mdeep_outfile = NULL){
  if (!is.null(CPC2_outfile))                  CPC2 <- read.CPC2(CPC2_outfile) else CPC2 <- ""
  if (!is.null(PLEK_outfile))                  PLEK <- read.PLEK(PLEK_outfile) else PLEK <- ""
  if (!is.null(FEELnc_outfile))                FEELnc <- read.FEELnc(FEELnc_outfile) else FEELnc <- ""
  if (!is.null(CPAT_outfile))                  CPAT <- read.CPAT(CPAT_outfile, CPAT_cutoff) else CPAT <- ""
  if (!is.null(CNCI_outfile))                  CNCI <- read.CNCI(CNCI_outfile) else CNCI <- ""
  if (!is.null(LncFinder_outfile))             LncFinder <- read.LncFinder(LncFinder_outfile) else LncFinder <- ""
  if (!is.null(lncRNA_Mdeep_outfile))          lncRNA_Mdeep <- read.lncRNA_Mdeep(lncRNA_Mdeep_outfile) else lncRNA_Mdeep <- ""

  seqID <- unique(c(CPC2, PLEK, CPAT, FEELnc, CNCI, LncFinder, lncRNA_Mdeep))
  #seqID <- seqID[!is.na(seqID)]
  #seqID <- seqID[seqID != ""]

  CodPotTbl <- data.frame(seqIDs = seqID,
                          CPC2 = 0,
                          PLEK = 0,
                          FEELnc = 0,
                          CPAT = 0,
                          CNCI = 0,
                          LncFinder = 0,
                          lncRNA_Mdeep = 0)

  CodPotTbl[CodPotTbl$seqIDs %in% CPC2,]$CPC2 <- 1
  CodPotTbl[CodPotTbl$seqIDs %in% PLEK,]$PLEK <- 1
  CodPotTbl[CodPotTbl$seqIDs %in% FEELnc,]$FEELnc <- 1
  CodPotTbl[CodPotTbl$seqIDs %in% CPAT,]$CPAT <- 1
  CodPotTbl[CodPotTbl$seqIDs %in% CNCI,]$CNCI <- 1
  CodPotTbl[CodPotTbl$seqIDs %in% LncFinder,]$LncFinder <- 1
  CodPotTbl[CodPotTbl$seqIDs %in% lncRNA_Mdeep,]$lncRNA_Mdeep <- 1

  CodPotTbl <- CodPotTbl[CodPotTbl$seqIDs != "", ]
  CodPotTbl <- CodPotTbl[!is.na(CodPotTbl$seqIDs), ]
  CodPotTbl <- CodPotTbl[,c(1978,colSums(CodPotTbl[,-1])) > 0]

return(CodPotTbl)
}



#' A read.CPC2 function
#'
#' This function reads CPC2 output and list all noncoding transcripts IDs
#' @param CPC2_outfile is the CPC2 output file localization including filename
#' @keywords CPC2 lncRNA
#' @export
#' @examples
#' read.CPC2()
#'
read.CPC2 <- function(CPC2_outfile){
  CPC2 <- read.table(CPC2_outfile)
  CPC2 <- CPC2[CPC2$V8 %in% "noncoding",]$V1
  CPC2 <- unique(as.character(CPC2))
return(CPC2)
}

#' A read.PLEK function
#'
#' This function reads PLEK output and list all noncoding transcripts IDs
#' @param PLEK_outfile is the PLEK output file localization including filename
#' @keywords PLEK lncRNA
#' @export
#' @examples
#' read.PLEK()
#'
read.PLEK <- function(PLEK_outfile){
  PLEK <- read.table(PLEK_outfile)
  PLEK$V3 <- gsub(">", "", PLEK$V3)
  PLEK <- PLEK[PLEK$V1 %in% "Non-coding",]$V3
  PLEK <- unique(as.character(PLEK))
return(PLEK)
}

#' A read.FEELnc function
#'
#' This function reads FEELnc output and list all noncoding transcripts IDs
#' @param FEELnc_outfile is the FEELnc output file (table containing 'lable' collumn inticating coding = 1 and noncoding = 0 transcripts) localization including filename
#' @keywords FEELnc lncRNA
#' @export
#' @examples
#' read.FEELnc()
#'
read.FEELnc <- function(FEELnc_outfile){
  FEELnc <- read.table(FEELnc_outfile, header = T)
  FEELnc <- FEELnc[FEELnc$label %in% 0,]$name
  FEELnc <- unique(as.character(FEELnc))
return(FEELnc)
}



#' A read.CPAT function
#'
#' This function reads CPAT output and list all noncoding transcripts IDs
#' @param CPAT_outfile is the CPAT output file (table containing 'lable' collumn inticating coding = 1 and noncoding = 0 transcripts) localization including filename
#' @param CPAT_cutoff is the cutoff value predefined (https://cpat.readthedocs.io/en/latest/#how-to-choose-cutoff) or computed i.e. using 'CPATcutoff()' function.
#' The default value is equal to predefined for human 0.364
#' @keywords CPAT lncRNA
#' @export
#' @examples
#' read.CPAT()
#'
read.CPAT <- function(CPAT_outfile, CPAT_cutoff = 0.364){
  CPAT <- read.table(CPAT_outfile, header = T)
  # WAZNE!!! 10Fold_CrossValidation_lncJJ.r tam jest estymacja progu mRNA/ncRNA
  # http://rna-cpat.sourceforge.net/#how-to-choose-cutoff
  #CPAT_cutoff <- 0.785 # i ta wartosc wynika z grafiki Figure3 z 10Fold_CrossValidation_lncJJ.r
  CPAT <- rownames(CPAT[CPAT$coding_prob < CPAT_cutoff,])
return(CPAT)
}



#' A read.CNCI function
#'
#' This function reads CNCI output and list all noncoding transcripts IDs
#' @param CNCI_outfile is the CNCI output file localization including filename
#' @keywords CNCI lncRNA
#' @export
#' @examples
#' read.CNCI()
#'
read.CNCI <- function(CNCI_outfile){
  CNCI <- read.table(CNCI_outfile, header = T, sep = "\t")
  CNCI <- CNCI[CNCI$index %in% "noncoding",1]
 # CNCI <- CNCI[CNCI$index %in% "noncoding",]$Transcript # it generates ERROR
  CNCI <- unique(as.character(CNCI))
  return(CNCI)
}

#' A read.LncFinder function
#'
#' This function reads LncFinder output and list all noncoding transcripts IDs
#' @param LncFinder_outfile is the LncFinder output file localization including filename
#' @keywords LncFinder lncRNA
#' @export
#' @examples
#' read.LncFinder()
#'
read.LncFinder <- function(LncFinder_outfile){
  LncFinder <- read.csv2(LncFinder_outfile, header = T, row.names = 1)
  LncFinder <- rownames(LncFinder[LncFinder$Pred %in% "NonCoding",])
  LncFinder <- unique(as.character(LncFinder))
  return(LncFinder)
}

#' A read.lncRNA_Mdeep function
#'
#' This function reads lncRNA_Mdeep output and list all noncoding transcripts IDs
#' @param lncRNA_Mdeep_outfile is the lncRNA_Mdeep output file localization including filename
#' @keywords lncRNA_Mdeep lncRNA
#' @export
#' @examples
#' read.lncRNA_Mdeep()
#'
read.lncRNA_Mdeep <- function(lncRNA_Mdeep_outfile){
  lncRNA_Mdeep <- read.table(lncRNA_Mdeep_outfile)
  lncRNA_Mdeep <- lncRNA_Mdeep[lncRNA_Mdeep$V2 %in% "noncoding",]$V1
  lncRNA_Mdeep <- unique(as.character(lncRNA_Mdeep))
  return(lncRNA_Mdeep)
}

#' A venn.CodPot function
#'
#' This function allows to draw Venn diagram od coding potential outputs
#' @param CodPot is the table generated usig CodPot2tbl() function
#' @param selmet is the vector argument indicatin which methods you select to venn diagram - i.e.: c(1,1,1,1) or c(1,0,1,0).
#' @param venncolors the vector of colors
#' @keywords venn.CodPot lncRNA
#' @export
#' @examples
#' venn.CodPot(CodPot = CodPot_table, selmet = c(1,1,0,0,1))
#'
venn.CodPot <- function(CodPot, venncolors = c("green", "red", "blue", "yellow", "magenta", "black", "orange", "white", "darkgreen"), selmet=NULL){

  if(is.null(selmet)) selmet <- rep(1, ncol(CodPot)-1)

 lista <- list()
 index <- 1

 for (i in which(selmet == 1)+1) {
   lista[[index]] <- CodPot[CodPot[,i] %in% "1",]$seqIDs
   index <- index+1
 }

 venn::venn(lista,
            zcolor = venncolors[1:sum(selmet)],
            snames = colnames(CodPot[which(selmet == 1)+1]),
            ilabels = T
 )
}

