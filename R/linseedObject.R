#' Linseed Object
#' 
#' The Linseed Object Class. 
#' 
#' Provides an interface to perform 
#' collinear network construction,
#' linear subspace identification, 
#' signature gene detection and gene expression deconvolution
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' 
#' @return Object of \code{\link{R6Class}} -- an interface to work with gene expression data.
#' @format \code{\link{R6Class}} object.
#' @examples
#' LinseedObject$new("GSE19830", samples=10:42, topGenes=12000)
#' 
#' @field exp List of two elements raw and normalized gene expression dataset
#' @field name Character, optional, dataset name
#' @field cellTypeNumber Identified cell type number, required for projection, 
#' corner detection and deconvolution
#' @field projection Projection of genes into space lower-dimensionality (presumably simplex)
#' @field endpoints Simplex corners (in normalized, non-reduced space)
#' @field endpointsProjection Simplex corners (in reduced space)
#' @field distances Stores distances for every gene to each corner in reduced space
#' @field clusters List that stores signatures genes for deconvolution, can be set manually or can be obtained by \code{selectGenes(k)}
#' @field signatures Deconvolution signature matrix
#' @field proportions Deconvolution proportion matrix
#' @field pairwise Calculated pairwise collinearity measure
#' 
#' @section Methods:
#' \describe{
#' 
#' }
#' @name LinseedObject
#' @import R6
#' @import CellMix
#' @import dplyr
#' @import ggplot2
#' @examples
#' 
LinseedObject <- R6Class("LinseedObject",
  public = list(
    exp = list(raw=NULL, norm=NULL),
    name = NULL,
    
    cellTypeNumber = NULL,
    
    projection = NULL,
    endpoints = NULL,
    endpointsProjection = NULL,
    distances = NULL,
    
    
    
    clusters = NULL,
    
    signatures = NULL,
    proportions = NULL,
    removedGenes = NULL,
    pairwise = NULL,
    
    projectiveProjection = NULL,
    Q = NULL,
    
    initialize = function(...) {
      args <- list(...)
      dataset <- args[[1]]
      if (inherits(dataset, "character") && grepl("^GSE", dataset)) {
          self$name <- dataset
          dataset <- preprocessGSE(...)
      } else {
          dataset <- preprocessDataset(...)
      }
      self$exp$raw <- dataset
      self$exp$norm <- dataset / rowSums(dataset)
      
      if (!is.null(args[["rplNorm"]]) && args[["rplNorm"]] == TRUE) {
        message("removing RPL/RPS signal")
        rpl_mean <- colMeans(self$exp$norm[grepl("RPL|RPS", 
                                                 rownames(self$exp$norm)), ])
        total_mean <- colMeans(self$exp$norm)
        rplNorm <- rpl_mean - total_mean
        rplNorm <- rplNorm / sqrt(sum(rplNorm^2))
        
        tmp <- self$exp$norm - self$exp$norm %*% t(t(rplNorm)) %*% t(rplNorm)
        tmp <- tmp[rowMin(tmp) > 0, ]
        self$exp$norm <- tmp
      }
      
    },
    
    svdPlot = function(dataset="norm") {
      requireNamespace("ggplot2")
      selected <- get(dataset, self$exp)
      vars <- svd(selected)$d^2
      vars <- vars / sum(vars)
      df <- data.frame(dimension=1:length(vars), variance=vars)
      ggplot2::ggplot(data=df, ggplot2::aes(x=dimension, y=variance)) +
        ggplot2::geom_segment(ggplot2::aes(xend=dimension, yend=min(vars), y=variance), size=2, color="black") +
        ggplot2::scale_y_log10() + ggplot2::theme_minimal(base_size = 17) + 
        ggplot2::theme(axis.line.x = ggplot2::element_line(colour = 'black', size=0.5, linetype='solid'),
                       axis.line.y = ggplot2::element_line(colour = 'black', size=0.5, linetype='solid'),
                       legend.position = "none")
      
    },
    
    
    tauSearch = function(taus=2^seq(-20, 10, 1),
                         iters=300,
                         genesToSelect=c(10, 50, 100, 200, 300),
                         methods=c("dsa", "ssFrobenius"),
                         dataset="raw") {
      
      requireNamespace("CellMix")
      requireNamespace("dplyr")
      
      k <- length(taus)
      l <- length(genesToSelect)

      endpoints <- list()
      endpointsProj <- list()
      
      cornerAns <- matrix(0, nrow=k, ncol=1)
      ans <- list()
      for (method in methods) {
        ans[[method]] <- matrix(0, nrow=k, ncol=l)
      }
      
      for (i in 1:k) {
        tau <- taus[i]
        
        self$sisalCorners(tau=tau, nonNeg = T, iter=iters)
        
        endpoints[[i]] <- self$endpoints
        endpointsProj[[i]] <- as.data.frame(t(self$endpointsProjection))
        endpointsProj[[i]]$log_tau <- log2(tau)
        
        self$deconvolveByEndpoints(dataset=dataset)
        cornerAns[i, 1] <- self$deconvolutionError(dataset=dataset)
        
        for (j in 1:l) {
          genes <- genesToSelect[j]
          self$selectGenes(genes)
          
          for (method in methods) {
            self$deconvolve(dataset=dataset, method=method, maxIter=500)
            ans[[method]][i, j] <- self$deconvolutionError(dataset=dataset)
          }
        }
        
      }
      
      rownames(cornerAns) <- log2(taus)
      colnames(cornerAns) <- c(0)
      cornerAns <- melt(cornerAns)
      cornerAns <- as.data.frame(cornerAns)
      cornerAns$type <- "corner"
      
      ans <- lapply(ans, function(an){
        rownames(an) <- log2(taus)
        colnames(an) <- genesToSelect
        an
      })
      
      ans <- lapply(ans, function(an) as.data.frame(melt(an)))
      ans <- mapply(function(an, nam){
        an$type <- nam
        an
      }, ans, methods, SIMPLIFY = F)
      
      ans <- do.call(rbind, ans)
      ans <- rbind(ans, cornerAns)

      return(list(errors=ans,
                  endpointsProj=do.call(rbind, endpointsProj)))
    },
    
    hysime = function(set=FALSE, dataset="norm") {
      selected <- get(dataset, self$exp)
      Y <- t(selected)
      noise <- estimateNoise(Y, verbose = T)
      hysimeRes <- hysime(Y, noise$w, noise$Rw, verbose=T)
      if (set) {
        self$cellTypeNumber <- hysimeRes$k
      }
      return(hysimeRes)
    },
    
    setCellTypeNumber = function(k) {
      self$cellTypeNumber <- k
    },
    
    project = function() {
      Y <- t(self$exp$norm)
      if (is.null(self$cellTypeNumber)) stop("Set cell type number first")
      self$projection <- projectiveProjection(Y, self$cellTypeNumber)
    },
    
    projectionPlot = function(dims=1:2) {
      requireNamespace("ggplot2")
      requireNamespace("dplyr")
      
      if (is.null(self$projection)) {
        stop("You have to call 'project' method first")
      }
      
      toPlot <- t(self$projection[dims, ])
      toPlot <- as.data.frame(toPlot)
      colnames(toPlot) <- c("x", "y")
      toPlot$type <- "gene"
      
      if (!is.null(self$endpointsProjection)) {
        toPlotE <- t(self$endpointsProjection[dims, ])
        toPlotE <- as.data.frame(toPlotE)
        colnames(toPlotE) <- c("x", "y")
        toPlotE$type <- "identified corners"
        toPlot <- rbind(toPlot, toPlotE)
      }
      
      
      if (!is.null(self$clusters)) {
        toPlot$cluster <- NA
        for (i in 1:length(self$clusters)) {
          toPlot[self$clusters[[i]], "cluster"] <- i
        }
        toPlot$cluster <- as.factor(toPlot$cluster)
      }
      
      
      axisNames <- paste0("Projection ", dims)
      pl <- ggplot2::ggplot(data=toPlot, ggplot2::aes(x=x, y=y)) +
        ggplot2::theme_bw() +
        ggplot2::labs(x=axisNames[1], y=axisNames[2])
      
      if (!is.null(toPlot$cluster)) {
        pl <- pl + ggplot2::geom_point(ggplot2::aes(shape=type, size=type, color=cluster))
      } else {
        pl <- pl + ggplot2::geom_point(ggplot2::aes(shape=type, size=type, color=type))
      }
      
      if (!is.null(self$endpointsProjection)) {
        pl <- pl + 
          scale_shape_manual(values=c(20, 18), labels=c("gene", "identified corners")) +
          scale_size_manual(values=c(1, 3), labels=c("gene", "identified corners")) +
          geom_polygon(data=dplyr::filter(toPlot, type=="identified corners"),
                       mapping=ggplot2::aes(x, y), fill=NA, color="black", lty=2)
      }
      pl
    },
    
    sisalCorners = function(...){
      Y <- t(self$exp$norm)
      p <- self$cellTypeNumber
      sisalRes <- sisal(Y, p, ...)
      self$endpoints = sisalRes$endpoints
      self$endpointsProjection = sisalRes$endpointsProjection
      self$distances = sisalRes$distances
      self$Q = sisalRes$Q
      self$projectiveProjection = sisalRes$projection
    },
    
    selectGenes = function(n) {
      pureGeneSets <- apply(self$distances, 2, function(xx) {
        pure <- rownames(self$distances)[order(xx)[1:n]]
        return(pure)
      })
      pureGeneSets <-  split(pureGeneSets, rep(1:ncol(pureGeneSets), each = nrow(pureGeneSets)))
      self$clusters <- pureGeneSets
    },
    
    selectGenesUnique = function(n) {
      # TODO
    },
    
    deconvolve = function(dataset="norm", method="dsa", ...) {
      if (!method %in% c("dsa", "ssFrobenius")) {
        stop("Supported methods are `dsa` and `ssFrobenius`")
      }
      selected <- get(dataset, self$exp)
      
      if (method == "dsa") {
        dsaRes <- fastDSA(selected, self$clusters)
        self$signatures = dsaRes$W
        self$proportions = dsaRes$H

      } else if (method == "ssFrobenius") {
        requireNamespace("CellMix")
        res <- CellMix::ged(selected, self$cellTypeNumber, 
                            data=self$clusters, seed=1, method="ssFrobenius", ...)
        dsaRes <- pureDsa(selected, CellMix::coef(res))
        self$signatures = dsaRes$W
        self$proportions = dsaRes$H
      }
      
      ctNames <- paste0("Cell type ", 1:length(self$clusters))
      colnames(self$signatures) <- ctNames
      rownames(self$proportions) <- ctNames
      
    },
    
    
    
    deconvolutionError = function(dataset="norm") {
      selected <- get(dataset, self$exp)
      reconstruct <- self$signatures %*% self$proportions
      diff <- selected - reconstruct
      return(norm(diff, "F"))
    },
    
    deconvolveByEndpoints = function(dataset="norm") {
      selected <- get(dataset, self$exp)
      dsaRes <- pureDsa(selected, t(self$endpoints))
      
      self$signatures = dsaRes$W
      self$proportions = dsaRes$H
      
      ctNames <- paste0("Cell type ", 1:self$cellTypeNumber)
      colnames(self$signatures) <- ctNames
      rownames(self$proportions) <- ctNames
    },
    
    
    getGenesByCutOff = function(r2Threshold, k) {
      subnetwork <- r2Filtering(self$pairwise, 0.99, 1)
      subnetwork$filteredGenes
    },
    
    calculatePairwiseLinearity = function(negToZero=F) {
      self$pairwise <- pairwiseR2(t(self$exp$norm))
      colnames(self$pairwise) <- rownames(self$pairwise) <- rownames(self$exp$norm)
      if (negToZero) {
        self$pairwise[self$pairwise < 0] <- 0
      }
    },
    
    linearityPlot = function(rs = seq(0, 1, 0.1),
                             ks=c(1, 3, 5, 10)) {
      prof <- r2Profiler(self$pairwise, rs, ks)
      visualizeProfilerResults(prof)
    }
  )
)