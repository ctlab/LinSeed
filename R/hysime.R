#' HySime
#'
#' Evalutes number of cell types presented in mixture using (HySime) hyperspectral signal identification by minimum error
#'
#' Original paper is J. M. Bioucas-Dias and J. M. P. Nascimento, "Hyperspectral Subspace Identification," in IEEE Transactions on Geoscience and Remote Sensing, vol. 46, no. 8, pp. 2435-2445, Aug. 2008.
#'
#' Based on MATLAB original code from http://www.lx.it.pt/~bioucas/code.htm
#'
#' @param Y gene expression matrix, columns are genes, rows are samples
#' @param W estimated noise matrix
#' @param Rn estimated noise correlation matrix
#' @param verbose verbosity, default valie is FALSE
#'
#' @return list
#' @export
#'
#' @examples
hysime <- function(Y, W, Rn, verbose=FALSE) {
    L <- nrow(Y)
    N <- ncol(Y)

    Lw <- nrow(W)
    Nw <- ncol(W)

    d1 <- nrow(Rn)
    d2 <- ncol(Rn)

    if (Lw != L || Nw != N) {
        stop("Noise matrix size are not in agreement")
    } else if (d1 != d2 || d1 != L) {
        stop("Bad correlation matrix")
    }

    if (verbose) message("Computing the correlation matrices")
    X <- Y - W
    Ry <- Y %*% t(Y) / N
    Rx <- X %*% t(X) / N

    if (verbose) message("Computing the eigen vectors of the signal correlation matrix")

    svdRes <- svd(Rx)
    dx <- svdRes$d
    E <- svdRes$u

    if (verbose) message("Estimating number of endmembers")

    Rn <- Rn + (sum(diag(Rx)) / L / 10^10) * diag(nrow=L)

    Py <- diag(t(E) %*% Ry %*% E)
    Pn <- diag(t(E) %*% Rn %*% E)

    costF <- -Py + 2 * Pn

    kf <- sum(costF < 0)
    ascOrder <- order(costF)
    Ek <- E[, ascOrder[1:kf]]

    if (verbose) message(sprintf("Estimated signal subspace dimension in k = %d", kf))

    PySort <- sum(diag(Ry)) - cumsum(Py[ascOrder])
    PnSort <- 2 * cumsum(Pn[ascOrder])
    costFSort <- PySort + PnSort

    if (verbose) {
        require(ggplot2)
        ind <- 1:(L - 1)
        toPlot <- data.frame(indice = rep(ind, 3),
                             error = c(costFSort[ind], PySort[ind], PnSort[ind]),
                             error_type = c(rep("MSE", L-1), rep("Proj Error", L-1), rep("Noise Power", L-1)))
        plot <- ggplot2::ggplot(data=toPlot, aes(x=indice, y=error, group=error_type, color=error_type)) +
            ggplot2::geom_point() + ggplot2::geom_line() +
            ggplot2::scale_y_log10() + ggplot2::theme_bw()
        return(list(k=kf, E=Ek, plot=plot))
    }
    return(list(k=kf, E=Ek, plot=NULL))
}
