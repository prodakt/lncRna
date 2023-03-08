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

#' NOTin function
#'
#' negation of IN.
#' @param x vector
#' @param y vector
#' @keywords not in
#' @export
#' @examples
#' A NotIn B

NotIn <- function(x,y)!('%in%'(x,y))
