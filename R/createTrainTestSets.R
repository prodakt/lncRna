#' Split a Set of Sequence Names into Training and Test Sets
#'
#' This function partitions a vector of sequence names into training and testing
#' subsets based on a specified percentage. It uses random sampling, so for
#' reproducible splits, set a seed with `set.seed()` before calling it.
#' This single function replaces the previous, separate `test.train.cds` and
#' `test.train.nc` functions.
#'
#' @param sequences A character vector of sequence names, or an object with a
#'   `names` attribute (like the list returned by `seqinr::read.fasta`).
#' @param percentTrain A numeric value between 0 and 1 indicating the proportion
#'   of sequences to allocate to the training set (default: 0.6).
#' @param prefix A character string used as a prefix for the names of the
#'   elements in the returned list (default: "set"). For example, `prefix = "cds"`
#'   will result in list elements named "cds.train" and "cds.test".
#'
#' @return A list containing two character vectors. The names of the list
#'   elements are constructed using the `prefix` argument (e.g., `cds.train`
#'   and `cds.test`).
#'
#' @export
#'
#' @examples
#' # --- Example 1: Splitting CDS sequences (replaces test.train.cds) ---
#' all_cds_names <- paste0("cds_seq_", 1:200)
#' set.seed(123) # for a reproducible split
#' cds_split <- createTrainTestSets(
#'   sequences = all_cds_names,
#'   percentTrain = 0.7,
#'   prefix = "cds"
#' )
#' names(cds_split)
#' length(cds_split$cds.train)
#' length(cds_split$cds.test)
#'
#' # --- Example 2: Splitting non-coding sequences (replaces test.train.nc) ---
#' # Input can also be a list with names
#' nc_fasta_like <- as.list(paste0("nc_seq_", 1:100))
#' names(nc_fasta_like) <- paste0("nc_seq_", 1:100)
#' set.seed(456)
#' nc_split <- createTrainTestSets(
#'   sequences = nc_fasta_like,
#'   percentTrain = 0.5,
#'   prefix = "nc"
#' )
#' names(nc_split)
#' length(nc_split$nc.train)
#' length(nc_split$nc.test)
#'
createTrainTestSets <- function(sequences, percentTrain = 0.6, prefix = "set") {
    # --- Input Validation and Name Extraction ---
    if (is.character(sequences)) {
        sequence_names <- sequences
    } else if (!is.null(names(sequences))) {
        sequence_names <- names(sequences)
    } else {
        stop(
            "'sequences' must be a character vector or an object with names."
        )
    }
    
    if (length(sequence_names) == 0) {
        warning("No sequence names found in 'sequences' input.")
        train_name <- paste0(prefix, ".train")
        test_name <- paste0(prefix, ".test")
        result_list <- list(character(0), character(0))
        names(result_list) <- c(train_name, test_name)
        return(result_list)
    }
    
    if (!is.numeric(percentTrain) || percentTrain < 0 || percentTrain > 1) {
        stop("'percentTrain' must be a number between 0 and 1.")
    }
    
    # --- Splitting Logic ---
    n_total <- length(sequence_names)
    n_train <- round(n_total * percentTrain)
    
    # Randomly sample indices for the training set
    train_indices <- sample(seq_len(n_total), size = n_train)
    
    train_names <- sequence_names[train_indices]
    test_names <- sequence_names[-train_indices]
    
    # --- Format Output ---
    train_name <- paste0(prefix, ".train")
    test_name <- paste0(prefix, ".test")
    
    result_list <- list(train_names, test_names)
    names(result_list) <- c(train_name, test_name)
    
    return(result_list)
}