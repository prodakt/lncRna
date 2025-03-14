#' test.train.cds function
#'
#' This function allows you to split a fast format sequence file into two sets: test and train.
#' As a result, the function outputs a list consisting of two elements "test" and "train" containing sequence names.
#'
#' @param cds the variable containing sequences in fasta format (from read.fasta() function from seqinr library)
#' @param percent_train the ratio value to divide the input set (he default is 0.6, which means that 60\% of the sequence is to be allocated to the training set)
#'
#' @examples
#' # cds <- read.fasta("lncRNA/Mus_musculus.GRCm39.cds.all.fa", seqtype = "DNA", as.string = T, set.attributes = F)
#' # cds_tt <- test.train.cds(cds = cds, percent_train = 0.6)
#' @export

test.train.cds <- function(cds, percent_train = 0.6){
  cds <- names(cds)
  cds_train_n <- sample(1:length(cds), round(length(cds)*percent_train), replace = F)
  cds_train <- cds[cds_train_n]
  cds_test <- cds[-(cds_train_n)]
  return(list(cds.train = cds_train, cds.test = cds_test))
}
