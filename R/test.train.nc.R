#' Split Non-Coding Sequences into Test and Train Sets
#'
#' This function splits a set of non-coding sequences into training and testing
#' subsets based on a specified percentage. It uses random sampling. For
#' reproducible splits, set the random seed using `set.seed()` before calling
#' this function.
#'
#' @param nc Either a character vector of non-coding sequence names, or a list/vector
#'           (like the output of `seqinr::read.fasta`) where the names attribute
#'           contains the sequence IDs.
#' @param percent_train The proportion (0 to 1) of sequences to allocate to the
#'                      training set (default: 0.6 or 60%).
#'
#' @return A list containing two character vectors:
#'         \item{nc.train}{Sequence names allocated to the training set.}
#'         \item{nc.test}{Sequence names allocated to the testing set.}
#' @keywords split data partitioning non-coding lncRNA machine learning sampling
#' @export
#' @examples
#' # Example 1: Input is a character vector of names
#' nc_names_vec <- paste0("nc_seq_", 1:100)
#' set.seed(123) # for reproducible example
#' nc_split1 <- test.train.nc(nc = nc_names_vec, percent_train = 0.7)
#' length(nc_split1$nc.train)
#' length(nc_split1$nc.test)
#'
#' # Example 2: Input is a list with names (like seqinr output)
#' # fasta_like_list <- list(nc_seq_1 = "ACGT", nc_seq_2 = "GGTA", nc_seq_3 = "TTAC")
#' # set.seed(123)
#' # nc_split2 <- test.train.nc(nc = fasta_like_list, percent_train = 0.6)
#' # print(nc_split2)
#'
test.train.nc <- function(nc, percent_train = 0.6) {
  # --- Input Validation and Name Extraction ---
  if (is.character(nc)) {
    # Input is already a character vector of names
    nc_names <- nc
  } else if (!is.null(names(nc))) {
    # Input is likely a list or vector with names (e.g., from read.fasta)
    nc_names <- names(nc)
    if (is.null(nc_names)) { # Double check if names() returned NULL
      stop("'nc' input is not a character vector and does not have names.")
    }
  } else {
    stop("'nc' must be a character vector of sequence names or an object with names (e.g., list from read.fasta).")
  }

  # Check if sequence names were successfully extracted
  if (length(nc_names) == 0) {
    warning("No sequence names found in 'nc' input.")
    return(list(nc.train = character(0), nc.test = character(0)))
  }

  # Validate percent_train
  if (!is.numeric(percent_train) || percent_train < 0 || percent_train > 1) {
    stop("'percent_train' must be a number between 0 and 1.")
  }

  # --- Splitting Logic ---
  n_total <- length(nc_names)
  n_train <- round(n_total * percent_train)

  # Handle edge cases where n_train might be 0 or n_total
  if (n_train <= 0) {
    nc_train_names <- character(0)
    nc_test_names <- nc_names
  } else if (n_train >= n_total) {
    nc_train_names <- nc_names
    nc_test_names <- character(0)
  } else {
    # Randomly sample indices for the training set without replacement
    train_indices <- sample(1:n_total, n_train, replace = FALSE)
    # Select names for training and testing sets based on indices
    nc_train_names <- nc_names[train_indices]
    nc_test_names <- nc_names[-train_indices] # Select elements *not* in train_indices
  }

  # Return the lists
  return(list(nc.train = nc_train_names, nc.test = nc_test_names))
}
