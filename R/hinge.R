#' Hinge function
#'
#' @param Y matrix describing genes (possibly in reduced space)
#'
#' @return matrix -Y with negative values replaced with zeroes
hinge <- function(Y) {
    return(pmax(-Y, 0))
}
