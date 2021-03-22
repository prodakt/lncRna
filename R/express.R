#' ExpSums function
#'
#' This function allows you to sum expression values.
#' @param ExprTab is the table with expression valules
#' @keywords expression lncRNA
#' @export
#' @examples
#' ExprTab()

ExpSums <- function(ExprTab){
  sums <- rowSums(ExprTab)
  names(sums) <- rownames(ExprTab)
return(sums)
}


#' ExpMeans function
#'
#' This function allows you to sum expression values.
#' @param ExprTab is the table with expression valules
#' @keywords expression lncRNA
#' @export
#' @examples
#' ExprTab()

ExpMeans <- function(ExprTab){
  means <- rowMeans(ExprTab)
  names(means) <- rownames(ExprTab)
  return(means)
}
