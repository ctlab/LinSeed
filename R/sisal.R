#' SISAL algorithm
#'
#' Sisal alogorithm for simplex endpoints identification
#'
#' Implementation of method
#'
#'
#' @param Y gene expression matrix
#' @param p number of endpoints
#' @param iters number of iterations to perform
#' @param tau noise points penalty coefficient
#' @param mu regularization
#' @param spherize spherize or not
#' @param tol numeric tolerance
#' @param m0 starting points for SISAL algorithm, default points are getting from VCA
#' @param verbose verbosity
#' @param returnPlot logical, is it needed to return dataframe or not
#'
#' @import corpcor
#' @import dplyr
#'
#' @return
#' @export
#'
#' @examples
sisal <- function(Y, p, iters = 80, tau = 1,
                  mu = p * 1000 / ncol(Y),
                  spherize = F, tol = 1e-2, m0 = NULL, verbose=F,
                  returnPlot = T, nonNeg = F) {
    rnames <- rownames(Y)
    L <- nrow(Y)
    N <- ncol(Y)

    if (L < p) stop("Insufficient number of columns in y")

    # local stuff

    slack <- 1e-3
    energyDecreasing <- 0
    fValBack <- Inf
    lamSphe <- 1e-8
    lamQuad <- 1e-6
    ALiters <- 4
    flaged <- 0

    ## At first we are getting the affine set
    Ymean <- apply(Y, 1, mean)
    ym <- matrix(Ymean, ncol=1)
    Y <- Y - Ymean
    svdObj <- fast.svd(Y)
    Up <- svdObj$u[, 1:(p-1)]
    proj <- Up %*% t(Up)
    D <- svdObj$d[1:(p-1)]
    
    Y <- proj %*% Y
    Y <- Y + Ymean
    YmeanOrtho <- ym - proj %*% ym
    Up <- cbind(Up, YmeanOrtho / (sqrt(sum(YmeanOrtho ^ 2))))
    singValues <- D
    
    Y <- t(Up) %*% Y
    
    ## spherizing
    if (spherize) {
      Y <- Up %*% Y
      Y <- Y - Ymean
      C <- diag(1 / sqrt(D + lamSphe))
      IC <- solve(C)
      Y <- C %*% t(Up[, 1:(p-1)]) %*% Y
      Y <- rbind(Y, 1)
      Y <- Y / sqrt(p)
    }

    ## Init

    if (is.null(m0)) {
        Mvca <- vca(Y, p)
        M <- Mvca
        Ym <- apply(M, 1, mean)
        dQ <- M - Ym
        M <- M + p * dQ
    } else {
        M <- m0
        M <- M - Ymean
        M <- Up[, 1:(p-1)] %*% t(Up[, 1:(p - 1)]) %*% M
        M <- M + Ymean
        M <- t(Up) %*% M
        if (spherize) {
            M <- Up %*% M - Ymean
            M <- C %*% t(Up[, 1:(p - 1)]) %*% M
            M[p, ] <- 1
            M <- M / sqrt(p)
        }
    }

    if (returnPlot) {
        toPlot <- data.frame(
            x = Y[1, ],
            y = Y[2, ],
            z = Y[3, ],
            type = "data point",
            iter = NA,
            tau = NA
        )
        starting <- data.frame(
            x = M[1, ],
            y = M[2, ],
            z = M[3, ],
            type = "sisal",
            iter = 0,
            tau = tau
        )
        toPlot <- rbind(toPlot, starting)
    }


    Q0 <- solve(M)
    Q <- Q0

    AAT <- kronecker(Y %*% t(Y), diag(nrow=p))
    B <- kronecker(diag(nrow=p), matrix(1, nrow=1, ncol=p))
    qm <- rowSums(solve(Y %*% t(Y)) %*% Y)
    qm <- matrix(qm, ncol=1)

    H <- lamQuad * diag(nrow=p^2)
    FF <- H + mu * AAT
    IFF <- solve(FF)


    G <- IFF %*% t(B) %*% solve(B %*% IFF %*% t(B))
    qmAux <- G %*% qm
    G <- IFF - G %*% B %*% IFF

    Z <- Q %*% Y
    Bk <- 0 * Z
    
    fmin = Inf
    Qmin = NULL
    
    for (k in 1:iters) {
        IQ <- solve(Q)
        g <- -t(IQ)
        dim(g) <- c(nrow(g) * ncol(g), 1)

        q0 <- Q
        dim(q0) <- c(nrow(Q) * ncol(Q), 1)
        Q0 <- Q

        baux <- H %*% q0 - g

        if (verbose) {
            if (spherize) {
                M <- IQ * sqrt(p)
                M <- M[1:(p-1), ]
                M <- Up[, 1:(p-1)] %*% IC %*% M
                M <- M + Ymean
                M <- t(Up) %*% M
            } else {
                M <- IQ
            }
            message(sprintf("Iteration %d, simplex volume = %4f", k, abs(det(M)) / factorial(nrow(M))))
        }

        if (k == iters) {
            ALiters <- 100
        }

        while (T) {
            q <- Q
            dim(q) <- c(nrow(Q) * ncol(Q), 1)

            f0val <- -log(abs(det(Q))) + tau * sum(hinge(Q %*% Y))
            f0quad <- t(q - q0) %*% g + 0.5 * t(q - q0) %*% H %*% (q - q0) + tau * sum(hinge(Q %*% Y))
            
            for (i in 2:ALiters) {
                dqAux <- Z + Bk
                dtzB <- dqAux %*% t(Y)
                dim(dtzB) <- c(nrow(dtzB) * ncol(dtzB), 1)
                b <- baux + mu * dtzB
                q <- G %*% b + qmAux
                Q <- matrix(q, nrow=p, ncol=p)

                Z <- softNeg(Q %*% Y - Bk, tau / mu)

                Bk <- Bk - (Q %*% Y - Z)

            }
            
            fquad_tmp <- t(q - q0) %*% g + 0.5 * t(q - q0) %*% H %*% (q - q0) + tau * sum(hinge(Q %*% Y))
            fval_tmp <- -log(abs(det(Q))) + tau * sum(hinge(Q %*% Y))
            
            # if (nonNeg) {
            #   qorig <- Up %*% solve(Q)
            #   
            #   if (any(qorig < 0)) {
            #     diff <- Ymean - qorig
            #     
            #     shift <- sapply(1:p, function(i) {
            #       tmp <- qorig[, i] / diff[, i]
            #       ifelse(any(qorig[, i] < 0), min(tmp[qorig[, i] < 0]), 0)
            #     })
            #     
            #     shift <- -shift + (0.00001)
            #     qorig <- qorig + diff %*% diag(shift)
            #     Q <- solve(t(Up) %*% qorig)
            #     q <- Q
            #     dim(q) <- c(nrow(Q) * ncol(Q), 1)
            #   }
            # }
            
            fquad <- t(q - q0) %*% g + 0.5 * t(q - q0) %*% H %*% (q - q0) + tau * sum(hinge(Q %*% Y))
            fval <- -log(abs(det(Q))) + tau * sum(hinge(Q %*% Y))
            
            # if (fval < fmin) {
            #   fmin <- fval
            #   Qmin <- Q
            # }

            # message(sprintf("f0 quad: %4f , f0 val %4f ", f0quad, f0val))
            # message(sprintf("f temp quad: %4f , f tmp val %4f ", fquad_tmp, fval_tmp))
            # message(sprintf("f quad: %4f , f val %4f ", fquad, fval))
            # stop("kek")

            if (f0quad >= fquad) {
                while (f0val - fval < 0) {
                    Q <- (Q + Q0) / 2
                    fval <- -log(abs(det(Q))) + tau * sum(hinge(Q %*% Y))
                }
              
                break
            }

        }

        if (returnPlot) {
            M <- solve(Q)
            toAdd <- data.frame(
                x = M[1, ],
                y = M[2, ],
                z = M[3, ],
                type = "sisal",
                iter = k,
                tau = tau
            )
            toPlot <- rbind(toPlot, toAdd)
            toPlot <- tbl_df(toPlot)
        }

    }
    
    if (nonNeg) {
      qorig <- Up %*% solve(Q)
      
      if (any(qorig < 0)) {
        diff <- Ymean - qorig
        
        shift <- sapply(1:p, function(i) {
          tmp <- qorig[, i] / diff[, i]
          ifelse(any(qorig[, i] < 0), min(tmp[qorig[, i] < 0]), 0)
        })
        
        shift <- -shift + (0.00001)
        qorig <- qorig + diff %*% diag(shift)
        Q <- solve(t(Up) %*% qorig)
        q <- Q
        dim(q) <- c(nrow(Q) * ncol(Q), 1)
      }
    }

    distanceToEndpoints <- apply(M, 2, function(x) {
        shifted <- Y - x
        dds <- sqrt(colSums(shifted^2))
        return(dds)
    })
    
    endpointsProjection <- M

    if (spherize) {
        M <- solve(Q)
        M <- M * sqrt(p)
        M <- M[1:(p - 1), ]
        M <- Up[, 1:(p - 1)] %*% IC %*% M
        M <- M + Ymean
    } else {
        M <- Up %*% solve(Q)
    }

    colnames(M) <- paste0("Pure gene ", 1:p)
    rownames(M) <- rnames

    retVal <- list(
        endpoints = M,
        endpointsProjection = endpointsProjection,
        projection = Up,
        shift = Ymean,
        singValues = singValues,
        distances = distanceToEndpoints,
        reduced = Y,
        q = Q
    )
    if (returnPlot) retVal$plotObj <- toPlot
    return(retVal)

}
