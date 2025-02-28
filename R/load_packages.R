#' Load R packages, with options to include default, CRAN, and GitHub packages.
#'
#' This function loads a set of default R packages, and optionally installs and loads additional packages from CRAN and GitHub.
#' It handles package installation and loading, providing informative messages about the process and any errors encountered.
#'
#' @param additional_packages A character vector of package names to install and load from CRAN, in addition to the default packages. Defaults to \code{NULL}.
#' @param additional_github_packages A character vector of GitHub package paths in the format "user/repo", to install and load, in addition to the default GitHub package ('lncRna'). Defaults to \code{NULL}.
#'
#' @examples
#' # Example 1: Load default packages only
#' # load_packages()
#'
#' # Example 2: Load default CRAN packages PLUS 'dplyr' and 'tidyr'
#' # load_packages(additional_packages = c("dplyr", "tidyr"))
#'
#' # Example 3: Load default packages PLUS additional GitHub packages
#' # load_packages(additional_github_packages = c("your_github_user/your_github_repo_package1", "another_github_user/another_package"))
#' # Note: For GitHub packages, you need to provide the full repo path "user/repo"
#'
#' # Example 4: Load default packages PLUS additional CRAN and GitHub packages
#' # load_packages(additional_packages = c("dplyr", "randomForest"), additional_github_packages = c("your_github_user/your_github_repo_package1"))
#'
#' # Example 5: Load default packages PLUS attempt to load a non-existent package from GitHub (error handling demo)
#' # load_packages(additional_github_packages = c("nonexistent_github_user/nonexistent_package"))
#'
#' # Example 6: Load default packages PLUS valid and invalid GitHub packages
#' # load_packages(additional_github_packages = c("prodakt/lncRna", "nonexistent_github_user/nonexistent_package")) # Note: "prodakt/lncRna" is already in default, but it will handle it correctly
#'
#' @export
load_packages <- function(additional_packages = NULL, additional_github_packages = NULL) {
  # Default packages
  default_packages <- c("devtools", "lncRna", "rtracklayer", "seqinr", "Polychrome", 
                        "plotly", "gprofiler2", "tidyr", "fmsb", "caret", "venn", 
                        "scales", "ggplot2")
  github_packages <- c("lncRna")  # Default GitHub package
  
  # Combine and categorize packages
  all_github <- unique(c(github_packages, additional_github_packages))
  all_cran <- unique(c(default_packages[!default_packages %in% all_github], additional_packages))
  
  # Load packages and track results
  loaded <- character(0)
  not_loaded <- character(0)
  
  # Process GitHub packages
  for (pkg in all_github) {
    if (!is.na(pkg) && pkg != "") {
      result <- install_and_load_package(pkg, "GitHub")
      if (result$loaded) loaded <- c(loaded, pkg) else not_loaded <- c(not_loaded, pkg)
    }
  }
  
  # Process CRAN packages
  for (pkg in all_cran) {
    if (!is.na(pkg) && pkg != "") {
      result <- install_and_load_package(pkg, "CRAN")
      if (result$loaded) loaded <- c(loaded, pkg) else not_loaded <- c(not_loaded, pkg)
    }
  }
  
  # Print results
  if (length(not_loaded) > 0) {
    cat("\nNot Loaded:", paste(not_loaded, collapse = ", "), "\nCheck errors above.\n")
  }
  if (length(loaded) > 0) {
    cat("\nLoaded:", paste(loaded, collapse = ", "), "\n")
  } else if (length(not_loaded) == 0) {
    cat("\nNo packages loaded. Check errors.\n")
  }
}

#######################################
# --- Helper Function to Install and Load a Single Package ---
#' Install and load a single R package from CRAN or GitHub.
#'
#' This helper function attempts to install and then load a specified R package.
#' It supports installation from both CRAN and GitHub, and provides detailed output
#' to the console regarding the installation and loading process, including any errors.
#'
#' @param pkg_name A character string naming the package to install and load.
#' @param source A character string specifying the source of the package, either "CRAN" or "GitHub". Defaults to "CRAN".
#' @return A list containing the status of the package loading attempt:
#'   \itemize{
#'     \item{\code{loaded}}{: A logical value indicating \code{TRUE} if the package was successfully loaded, \code{FALSE} otherwise.}
#'     \item{\code{error}}{: If an error occurred during installation or loading, this contains the error message as a character string. Otherwise, it is \code{NULL}.}
#'   }
#'
#' @keywords internal
#' @noRd
install_and_load_package <- function(pkg_name, source = c("CRAN", "GitHub")) {
  source <- match.arg(source) # Validate source argument
  
  package_status <- list(loaded = FALSE, error = NULL) # Initialize status list
  
  is_installed <- requireNamespace(pkg_name, quietly = TRUE)
  
  if (!is_installed) {
    install_message <- paste("Package", pkg_name, "is not installed. Installing from", source, "...")
    if (source == "GitHub") {
      install_message <- paste(install_message, "(prodakt/", pkg_name, ")..")
    }
    cat(install_message, "\n")
    
    install_attempt <- tryCatch({
      if (source == "CRAN") {
        install.packages(pkg_name, quiet = FALSE) # Show install output for CRAN
      } else if (source == "GitHub") {
        devtools::install_github(paste("prodakt/", pkg_name, sep=""), force = TRUE, quiet = FALSE) # Show install output for GitHub
      }
      TRUE # Return TRUE if installation succeeds
    }, error = function(e) {
      package_status$error <- e$message # Store error message
      cat(paste("Error installing package", pkg_name, "from", source, ".\n"))
      FALSE # Return FALSE if installation fails
    })
    
    if (install_attempt) {
      load_attempt <- tryCatch({
        library(pkg_name, character.only = TRUE, quietly = FALSE) # Show load output
        TRUE # Return TRUE if loading succeeds
      }, error = function(e) {
        package_status$error <- e$message # Store error message
        cat(paste("Error loading package", pkg_name, "after installation.\n"))
        FALSE # Return FALSE if loading fails
      })
      if (load_attempt) {
        package_status$loaded <- TRUE
        load_success_message <- paste("Package", pkg_name, "installed and loaded successfully from", source, ".")
        if (source == "GitHub") {
          load_success_message <- paste(load_success_message, "(GitHub).")
        }
        cat(load_success_message, "\n")
      }
    }
  } else {
    load_attempt <- tryCatch({
      library(pkg_name, character.only = TRUE, quietly = FALSE) # Show load output
      TRUE
    }, error = function(e) {
      package_status$error <- e$message
      cat(paste("Error loading package", pkg_name, ".\n"))
      FALSE
    })
    if (load_attempt) {
      package_status$loaded <- TRUE
      cat(paste("Package", pkg_name, "already installed and loaded.\n"))
    }
  }
  return(package_status) # Return status list
}

