#' A strGTF2stat function
#'
#' This function allows you to extract some statistics from GTF file
#' @param stringtieGTF is the GTF merged by stringtie and imported with 'importGFF()' function from 'rtracklayer' package
#' @keywords GTF stringtie lncRNA
#' @export
#' @examples
#' strGTF2stat('GTF from stringtie')


require(dplyr)

strGTF2TransLen <- function(stringtieGTF){
  strGTF <- as.data.frame(stringtieGTF)
  strGTF_ex <- strGTF[strGTF$type %in% "exon",]
  trans_len <- strGTF_ex %>%
    group_by(transcript_id) %>%
    summarise(trans_length = sum(width))
  trans_len <- as.data.frame(trans_len)
return(trans_len)
}


strGTF2ExonsN <- function(stringtieGTF){
  strGTF <- as.data.frame(stringtieGTF)
  exons_n <- strGTF_ex %>%
    group_by(transcript_id) %>%
    summarise(exons = max(exon_number))
  exons_n <- as.data.frame(exons_n)
  return(exons_n)
}


strGTF2stat <- function(stringtieGTF){
  stats <- merge(strGTF2ExonsN(stringtieGTF), strGTF2TransLen(stringtieGTF), by="transcript_id", all=T)
  return(stats)
}

