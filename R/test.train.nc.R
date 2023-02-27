#' test.train.nc function
#'
#' This function allows you to split a sequences (from FASTA file) into two sets: test and train. As a result, the function outputs a list consisting of two elements "test" and "train" containing sequence names.
#' @param nc.fa the variable containing sequences in fasta format (from read.fasta() function from seqinr library)
#' @param percent_train the ratio value to divide the input set (he default is 0.6, which means that 60% of the sequence is to be allocated to the training set)
#' @keywords training test proportion lncRNA
#' @export
#' @examples
#' nc <- read.fasta("lncRNA/Mus_musculus.GRCm39.ncrna.fa", seqtype = "DNA", as.string = T, set.attributes = F)
#' nc_tt <- test.train.nc(nc.fa = nc, percent_train = 0.6)

test.train.nc <- function(nc.fa, percent_train = 0.6){
  nc <- names(nc)
  nc_train_n <- sample(1:length(nc), round(length(nc)*percent_train), replace = F)
  nc_train <- nc[nc_train_n]
  nc_test <- nc[-(nc_train_n)]
  return(list(nc.train = nc_train, nc.test = nc_test))
}
