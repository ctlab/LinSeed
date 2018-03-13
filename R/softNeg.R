#' Soft negative score
#'
#' @param Y matrix
#' @param tau coefficient to penalize for volume
#'
#' @return
#'
#' @examples
softNeg <- function(Y, tau) {
    z <- pmax(abs(Y + tau / 2) - tau / 2, 0)
    z <- z / (z + tau / 2) * (Y + tau / 2)
    return(z)
}
