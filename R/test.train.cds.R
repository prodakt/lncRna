#' Split Protein-Coding Sequences into Test and Train Sets
#'
#' This function splits a set of protein-coding sequences into training and
#' testing subsets based on a specified percentage. It uses random sampling. For
#' reproducible splits, set the random seed using `set.seed()` before calling
#' this function.
#'
#' @param cds Either a character vector of protein-coding sequence names, or a list/vector
#'            (like the output of `seqinr::read.fasta`) where the names attribute
#'            contains the sequence IDs.
#' @param percent_train The proportion (0 to 1) of sequences to allocate to the
#'                      training set (default: 0.6 or 60%).
#'
#' @return A list containing two character vectors:
#'         \item{cds.train}{Sequence names allocated to the training set.}
#'         \item{cds.test}{Sequence names allocated to the testing set.}
#' @keywords split data partitioning coding cds machine learning sampling
#' @export
#' @examples
#' # Example 1: Input is a character vector of names
#' cds_names_vec <- paste0("cds_seq_", 1:200)
#' set.seed(456) # for reproducible example
#' cds_split1 <- test.train.cds(cds = cds_names_vec, percent_train = 0.5)
#' length(cds_split1$cds.train)
#' length(cds_split1$cds.test)
#'
#' # Example 2: Input is a list with names (like seqinr output)
#' # cds_fasta_like <- list(cds_seq_A = "ATGC", cds_seq_B = "GGCC", cds_seq_C = "TTAA")
#' # set.seed(456)
#' # cds_split2 <- test.train.cds(cds = cds_fasta_like, percent_train = 0.6)
#' # print(cds_split2)
#'
test.train.cds <- function(cds, percent_train = 0.6) {
  # --- Input Validation and Name Extraction ---
  if (is.character(cds)) {
    # Input is already a character vector of names
    cds_names <- cds
  } else if (!is.null(names(cds))) {
    # Input is likely a list or vector with names (e.g., from read.fasta)
    cds_names <- names(cds)
    if (is.null(cds_names)) { # Double check if names() returned NULL
      stop("'cds' input is not a character vector and does not have names.")
    }
  } else {
    stop("'cds' must be a character vector of sequence names or an object with names (e.g., list from read.fasta).")
  }

  # Check if sequence names were successfully extracted
  if (length(cds_names) == 0) {
    warning("No sequence names found in 'cds' input.")
    return(list(cds.train = character(0), cds.test = character(0)))
  }

  # Validate percent_train
  if (!is.numeric(percent_train) || percent_train < 0 || percent_train > 1) {
    stop("'percent_train' must be a number between 0 and 1.")
  }

  # --- Splitting Logic ---
  n_total <- length(cds_names)
  n_train <- round(n_total * percent_train)

  # Handle edge cases where n_train might be 0 or n_total
  if (n_train <= 0) {
    cds_train_names <- character(0)
    cds_test_names <- cds_names
  } else if (n_train >= n_total) {
    cds_train_names <- cds_names
    cds_test_names <- character(0)
  } else {
    # Randomly sample indices for the training set without replacement
    train_indices <- sample(1:n_total, n_train, replace = FALSE)
    # Select names for training and testing sets based on indices
    cds_train_names <- cds_names[train_indices]
    cds_test_names <- cds_names[-train_indices] # Select elements *not* in train_indices
  }


  # Return the lists
  return(list(cds.train = cds_train_names, cds.test = cds_test_names))
}
