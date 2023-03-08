#' \%!in\% function
#'
#' This function is negation of "contains.
#' @param x vector of values to be filtered
#' @param y vector of values to eliminate
#' @keywords negate
#' @export
#' @examples
#'

'%!in%' <- function(x,y)!('%in%'(x,y))

#' h.t function
#'
#' Display head and tail together.
#' @param table table
#' @param n number of rows to display
#' @keywords table head tail
#' @export
#' @examples
#' h.t(table)
#'
h.t <- function(table, n=3) {
  tbl.tmp <- rbind(head(table,n),"..." = rep("...", ncol(table)) ,tail(table,n))
  return(tbl.tmp)
}
