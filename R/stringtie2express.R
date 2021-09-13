#' strtie2expr_units function
#'
#' This function allows you to list the possible variables from stringtie output file.
#' @param
#' @keywords expression lncRNA
#' @export
#' @examples
#' strtie2expr_units()

strtie2expr_units <- function(strdir = "stringtie/", pattern = "gtab"){
  pliki <- list.files(paste0(strdir), pattern = pattern)
  tab <- read.table(paste0(strdir,pliki[1]), sep = "\t", header = T)
  return(colnames(tab))
}


#' strtie2expr function
#'
#' This function allows you to extract expression values from the stringtie output file.
#' @param
#' @keywords expression lncRNA
#' @export
#' @examples
#' strtie2expr()
#' FPKM_strtie <- strtie2expr()
#' FPKM_strtie <- strtie2expr(strdir = "stringtie/", pattern = "gtab", expunit = "FPKM")
#' TPM_strtie <- strtie2expr(expunit = "TPM")
strtie2expr <- function(strdir = "stringtie/", pattern = "gtab", expunit = "FPKM"){
  pliki <- list.files(paste0(strdir), pattern = pattern)
  tab <- read.table(paste0(strdir,pliki[1]), sep = "\t", header = T)
  tab <- tab[c("Gene.ID", "Gene.Name", "Reference", "Strand", "Start", "End")]
  for (f in pliki) {
    tabtmp <- read.table(paste0(strdir,f), sep = "\t", header = T)
    tabtmp <- tabtmp[,c("Gene.ID", paste0(expunit))]
    colnames(tabtmp)[2] <- gsub(pattern = ".gtab", replacement = "", x = f)
    tab <- merge(tab, tabtmp, by = "Gene.ID", all = T)
  }
  return(tab)
}


#' strtie2glen function
#'
#' This function allows you to compute genes length from stringtie output file.
#' @param
#' @keywords expression lncRNA
#' @export
#' @examples
#' strtie2glen()
strtie2glen <- function(strdir = "stringtie/", pattern = "gtab"){
  pliki <- list.files(paste0(strdir), pattern = pattern)
  message(paste0("Files: "))
  message(paste0(pliki, ", "))
  tab <- read.table(paste0(strdir,pliki[1]), sep = "\t", header = T)
  tab <- data.frame(Gene.ID = tab$Gene.ID)
  for (f in pliki) {
    tabtmp <- read.table(paste0(strdir,f), sep = "\t", header = T)
    tabtmp <- tabtmp[,c("Gene.ID", "Start", "End")]
    tabtmp$Length <- abs(tabtmp$Start - tabtmp$End)+1
    tabtmp <- tabtmp[,c("Gene.ID", "Length")]
    colnames(tabtmp)[2] <- gsub(pattern = ".gtab", replacement = "", x = f)
    tab <- merge(tab, tabtmp, by = "Gene.ID", all = T)
  }
  tab$Length <- tab[,2]
  tab <- tab[,c("Gene.ID", "Length")]
  # rownames(tab) <- tab$Gene.ID
  # tab <- tab[,-1]

  return(tab)
}

