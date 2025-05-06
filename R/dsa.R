#' Fast DSA algorithm implementation
#'
#' Runs DSA implementation using .fccnls for solving least-squares with multiple right-hand-sides
#'
#' @param dataset gene expression matrix
#' @param genes list with putative signatures for DSA algorithm
#' @import NMF
#' @import corpcor
#' @import BH
#' 
#' @return deconvolution results, list with H and W matrices
fastDSA <- function(dataset, genes) {
    eigengenes <- do.call(rbind, lapply(genes, function(geneSet) colMeans(dataset[geneSet,
        ])))
    eigenMultiplier <- fcnnls(t(eigengenes), matrix(1, ncol(eigengenes), 1))
    H <- diag(as.numeric(eigenMultiplier$x),
              length(as.numeric(eigenMultiplier$x))) %*% eigengenes
    res <- .fcnnls(t(H), t(dataset), pseudo = TRUE)
    return(list(H = H, W = t(res$coef)))
}

#' DSA algorithm implementation for pure points
#'
#' Runs DSA implementation using fcnnls_c for solving least-squares with multiple right-hand-sides
#'
#' @param dataset gene expression matrix
#' @param pure matrix contains expression of signature genes
#' @import BH
#'
#' @return deconvolution results, list with H and W matrices
pureDsa <- function(dataset, pure) {
  eigenMultiplier <- fcnnls(t(pure), matrix(1, ncol(pure), 1))
  H <- diag(as.numeric(eigenMultiplier$x),
            length(as.numeric(eigenMultiplier$x))) %*% pure
  res <- .fcnnls(t(H), t(dataset), pseudo = TRUE)
  return(list(H = H, W = t(res$coef)))
}
#' Run DSA by clusters
#'
#' Runs DSA with provided clusters as putative signatures
#'
#' @param dataset gene expression matrix
#' @param clustering numeric vector, clustering of the rows
#' @param clusters numeric vector, which clusters use as putative signatures
#'
#' @return deconvolution results, list with H and W matrices
runDSA <- function(dataset, clustering, clusters) {
    genes <- lapply(clusters, function(i) rownames(dataset[clustering == i, ]))
    fastDSA(dataset, genes)
}
