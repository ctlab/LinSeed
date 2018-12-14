
#' Noise estimation
#'
#' Estimates noise using multiple regression approach. Implements method described in
#' J. M. Bioucas-Dias and J. M. P. Nascimento, "Hyperspectral Subspace Identification," in IEEE Transactions on Geoscience and Remote Sensing, vol. 46, no. 8, pp. 2435-2445, Aug. 2008.
#'
#' Based on MATLAB original code from http://www.lx.it.pt/~bioucas/code.htm
#'
#' @param Y normalized gene expression data matrix, columns are genes and rows are samples
#' @param noiseType character, describing noise type. Two possible values are "additive" and "possion"
#' @param verbose logical, default value is FALSE
#'
#' @return list with two elements, w -- estimated noise and Rw estimated noise correlation matrix
#'
#' @examples
estimateNoise <- function(Y, noiseType="additive", verbose=FALSE) {
    L <- nrow(Y)
    N <- ncol(Y)

    if (L < 2) stop("Too few samples in dataset")
    if (!noiseType %in% c("additive", "poisson"))
        stop("Unknown noise model, accepted values are \"additive\" (default) and \"poisson\"")
    if (verbose) message("Noise estimation started:")

    if (noiseType == "poisson") {
        sqY <- sqrt(Y * (Y > 0))
        noiseEst <- estimateAdditiveNoise(sqY, verbose)
        x <- (sqY - noiseEst$w)^2
        w <- sqrt(x) * noiseEst$w * 2
        Rw <- w %*% t(w) / N
        return(list(w=w, Rw=Rw))
    } else {
        return(estimateAdditiveNoise(Y, verbose))
    }


}

#' Estimate additive noise
#'
#' Additive noise estimation subroutine
#'
#' @param Y normalized gene expression matrix
#' @param verbose verbosity
#'
#' @return
#'
#' @examples
estimateAdditiveNoise <- function(Y, verbose) {
    small <- 1e-6;
    L <- nrow(Y)
    N <- ncol(Y)

    w <- matrix(0, nrow=L, ncol=N)
    if (verbose) message("Computing the sample correlation matrix and its inverse")
    RR <- Y %*% t(Y)
    RRi <- solve(RR + small * diag(nrow=L))
    for (i in 1:L) {
        XX <- RRi - RRi[, i, drop=F] %*% RRi[i, , drop=F] / RRi[i, i]
        RRa <- RR[, i, drop=F]
        RRa[i, ] <- 0
        beta <- XX %*% RRa
        beta[i, ] <- 0
        w[i, ] = Y[i, ] - t(beta) %*% Y
    }

    if (verbose) message("Computing correlation matrix")
    Rw = diag(diag(w %*% t(w) / N))
    return(list(w=w, Rw=Rw))
}
