#' ExpSums function
#'
#' This function allows you to sum expression values.
#' @param ExprTab is the table with expression valules
#' @keywords expression lncRNA
#' @export
#' @examples
#' ExpSums()

ExpSums <- function(ExprTab){
  sums <- rowSums(ExprTab)
  names(sums) <- rownames(ExprTab)
return(sums)
}


#' ExpMeans function
#'
#' This function allows you to calculate the mean values of the expression.
#' @param ExprTab is the table with expression valules
#' @keywords expression lncRNA
#' @export
#' @examples
#' ExpMeans()

ExpMeans <- function(ExprTab){
  means <- rowMeans(ExprTab)
  names(means) <- rownames(ExprTab)
  return(means)
}
