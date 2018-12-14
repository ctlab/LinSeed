#' Linseed Object
#' 
#' The Linseed Object Class. 
#' 
#' Provides an interface to perform 
#' collinear network construction,
#' linear subspace identification, 
#' simplex corner identification and gene expression deconvolution
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @useDynLib linseed
#' @export
#' 
#' @return Object of \code{\link{R6Class}} -- an interface to work with gene expression data.
#' @format \code{\link{R6Class}} object.
#' @examples
#' LinseedObject$new("GSE19830", samples=10:42, topGenes=10000)
#' 
#' @field exp List of two elements raw and normalized gene expression dataset
#' @field name Character, optional, dataset name
#' @field cellTypeNumber Identified cell type number, required for projection, 
#' corner detection and deconvolution
#' @field projection Projection of genes into space lower-dimensionality (presumably simplex)
#' @field endpoints Simplex corners (in normalized, non-reduced space)
#' @field endpointsProjection Simplex corners (in reduced space)
#' @field distances Stores distances for every gene to each corner in reduced space
#' @field markers List that stores signatures genes for deconvolution, can be set manually or can be obtained by \code{selectGenes(k)}
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
#' @import dplyr
#' @import ggplot2
#' @import Matrix
#' @import progress
#' @import Rtsne
#' @examples
#' 
LinseedObject <- R6Class("LinseedObject",
                          public = list(
                           exp = list(full=list(raw=NULL, norm=NULL),
                                  filtered=list(raw=NULL, norm=NULL)),
                           name = NULL,
                           
                           cellTypeNumber = NULL,
                           
                           projection = NULL,
                           endpoints = NULL,
                           endpointsProjection = NULL,
                           distances = NULL,
                           
                           
                           
                           markers = NULL,
                           signatures = NULL,
                           proportions = NULL,
                           pairwise = NULL,
                           spearman = NULL,
                           
                           genes = list(
                             pvals = NULL,
                             powers = NULL,
                             degrees = NULL
                           ),
                           
                           
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
                             self$exp$full$raw <- dataset
                             self$exp$full$norm <- dataset / rowSums(dataset)
                           },
                           
                           svdPlot = function(dataset="norm", components=50) {
                             dataFull <- get(dataset, self$exp$full)
                             dataFiltered <- get(dataset, self$exp$filtered)
                             components <- min(components, ncol(dataFull))
                             
                             if (is.null(dataFull)) {
                               stop("Full dataset appears to be NULL, something is wrong")
                             }
                             
                             vars <- svd(dataFull)$d^2
                             vars <- cumsum(vars / sum(vars))
                             df <- data.frame(dimension=1:length(vars), variance=vars, filtered=FALSE)
                             colors <- 1
                             
                             if (!is.null(dataFiltered)) {
                               vars <- svd(dataFiltered)$d^2
                               vars <- cumsum(vars / sum(vars))
                               dfFiltered <- data.frame(dimension=1:length(vars), variance=vars, filtered=TRUE)
                               df <- rbind(df, dfFiltered)
                               colors <- 2
                             }
                             
                             
                             ggplot(data=df, aes(x=dimension, y=variance)) +
                               geom_point(aes(y=variance, color=filtered), size=0.5, alpha=1) +
                               geom_line(aes(y=variance, color=filtered, group=filtered), 
                                                     size=0.5, alpha=1) +
                               scale_color_manual(values=c("#999999", "#E41A1C")[1:colors]) +
                               theme_minimal(base_size = 8) + 
                               theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                                              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                                              legend.position = c(1, 0), legend.justification = c(1, 0), 
                                              legend.background = element_rect(colour="black", size=0.2),
                                              legend.key.size = unit(0.1, "in")) +
                               scale_x_continuous(minor_breaks = 1:components,
                                                           limits=c(0, components))
                          
                           },
                           
                           hysime = function(dataset="filtered", error="norm", set=FALSE) {
                             
                             data <- get(dataset, self$exp)
                             selected <- get(error, data)
                             
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
                           
                           project = function(dataset) {
                             data <- get(dataset, self$exp)
                             Y <- t(data$norm)
                             if (is.null(self$cellTypeNumber)) stop("Set cell type number first")
                             self$projection <- projectiveProjection(Y, self$cellTypeNumber)
                             self$projectiveProjection <- getProjectiveProjection(Y, self$cellTypeNumber)
                           },
                           
                           projectionPlot = function(dims=1:2, color=NULL) {
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
                               
                               
                               toPlotE$corner <- 1:self$cellTypeNumber
                               toPlotE$type <- "identified corners"
                               toPlot$corner <- NA
                               
                               toPlot <- rbind(toPlot, toPlotE)
                               toPlot$corner <- as.factor(toPlot$corner)
                             }
                             
                             
                             if (!is.null(self$markers)) {
                               toPlot$cluster <- NA
                               for (i in 1:length(self$markers)) {
                                 toPlot[self$markers[[i]], "cluster"] <- i
                               }
                               toPlot$cluster <- as.factor(toPlot$cluster)
                               toPlot <- toPlot[order(toPlot$cluster, decreasing = T), ]
                               toPlot <- toPlot[nrow(toPlot):1, ]
                             }
                             
                             if (!is.null(self$exp$filtered$norm)) {
                               geneSubset <- rownames(self$exp$filtered$norm)
                               toPlot$filtered <- F
                               toPlot[geneSubset, "filtered"] <- T
                             }
                             
                             
                             axisNames <- paste0("Projection ", dims)
                             pl <- ggplot(data=toPlot, aes(x=x, y=y)) +
                               theme_bw() +
                               labs(x=axisNames[1], y=axisNames[2])
                             
                             
                             if (!is.null(color)) {
                               pl <- pl + geom_point(aes_string(shape="type", size="type", color=color))
                             } else {
                               pl <- pl + geom_point(aes(shape=type, size=type))
                             }
                             
                             if (!is.null(self$endpointsProjection)) {
                               pl <- pl + 
                                 scale_shape_manual(values=c(20, 18), labels=c("gene", "identified corners")) +
                                 scale_size_manual(values=c(1, 3), labels=c("gene", "identified corners")) +
                                 geom_polygon(data=dplyr::filter(toPlot, type=="identified corners"),
                                              mapping=aes(x, y), fill=NA, color="black", lty=2)
                             }
                             pl
                           },
                           
                           sisalCorners = function(...){
                             
                             if (is.null(self$exp$filtered$norm)) {
                               message("Could not find filtered dataset, using whole dataset as filtered")
                               self$filterDataset(rownames(self$exp$full$norm))
                             }
                             
                             Y <- t(self$exp$filtered$norm)
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
                             self$markers <- pureGeneSets
                           },
                           
                           deconvolve = function(dataset="filtered", error="norm", method="dsa", ...) {
                             if (!method %in% c("dsa", "ssFrobenius")) {
                               stop("Supported methods are `dsa` and `ssFrobenius`")
                             }
                             data <- get(dataset, self$exp)
                             selected <- get(error, data)
                             
                             if (method == "dsa") {
                               dsaRes <- fastDSA(selected, self$markers)
                               self$signatures = dsaRes$W
                               self$proportions = dsaRes$H
                               
                             } else {
                               if (requireNamespace("CellMix")) {
                                 res <- CellMix::ged(selected, self$cellTypeNumber, 
                                                     data=self$markers, seed=1, method=method, ...)
                                 dsaRes <- pureDsa(selected, CellMix::coef(res))
                                 self$signatures = dsaRes$W
                                 self$proportions = dsaRes$H
                               } else {
                                 stop("Different method than DSA provided, but no CellMix found")
                               }
                               
                             }
                             
                             ctNames <- paste0("Cell type ", 1:length(self$markers))
                             colnames(self$signatures) <- ctNames
                             rownames(self$proportions) <- ctNames
                             
                           },
                           
                           
                           deconvolutionError = function(dataset="filtered", error="norm") {
                             data <- get(dataset, self$exp)
                             selected <- get(error, data)
                             reconstruct <- self$signatures %*% self$proportions
                             diff <- selected - reconstruct
                             return(norm(diff, "F"))
                           },
                           
                           deconvolveByEndpoints = function(dataset="filtered", error="norm") {
                             data <- get(dataset, self$exp)
                             selected <- get(error, data)
                             dsaRes <- pureDsa(selected, t(self$endpoints))
                             
                             self$signatures = dsaRes$W
                             self$proportions = dsaRes$H
                             
                             ctNames <- paste0("Cell type ", 1:self$cellTypeNumber)
                             
                             rownames(self$signatures) <- rownames(selected)
                             colnames(self$proportions) <- colnames(selected)
                             
                             colnames(self$signatures) <- ctNames
                             rownames(self$proportions) <- ctNames
                           },
                           
                           calculateSpearmanCorrelation = function() {
                             self$spearman <- cor(t(self$exp$full$norm), method="spearman")
                           },
                           
                           calculateSignificanceLevel = function(iters=1000,
                                                                 spearmanThreshold=0,
                                                                 retVal=F) {
                             
                             if (is.null(self$pairwise)) stop("call calculatePairwiseLinearity first")
                             if (is.null(self$spearman)) stop("call calculateSpearmanCorrelation first")
                             
                             degrees <- rowSums(self$pairwise > 0 & 
                                                self$spearman > spearmanThreshold) - 1
                             
                             ij <- which(upper.tri(self$pairwise) & 
                                           self$pairwise > 0 & 
                                           self$spearman > spearmanThreshold, arr.ind = T)
                             values <- self$pairwise[upper.tri(self$pairwise) & 
                                                       self$pairwise > 0 & 
                                                       self$spearman > spearmanThreshold] * 
                               self$spearman[upper.tri(self$pairwise) & 
                                               self$pairwise > 0 & 
                                               self$spearman > spearmanThreshold]
                             
                             genes <- nrow(self$pairwise)
                             edges <- length(values)
                             
                             sparsePP <- sparseMatrix(i = ij[, 1], j = ij[, 2], x = values, dims=c(genes, genes))
                             sparsePPs <- sparsePP + t(sparsePP)
                             vertexDegrees <- rowSums(sparsePPs)
                             orderedDegrees <- order(vertexDegrees)
                             
                             results <- numeric(genes)
                             gc(verbose = F)
                             
                             pb <- progress_bar$new(
                               format = "Sampling weights [:bar] :percent eta: :eta",
                               total = iters, clear = FALSE, width= 60)
                             pb$tick(0)
                             
                             for (i in 1:iters) {
                               valuesRandom <- sample(values, replace = T)
                               sparseRandomPP <- sparseMatrix(i = ij[, 1], j = ij[, 2], x = valuesRandom, dims=c(genes, genes))
                               sparseRandomPPs <- sparseRandomPP + t(sparseRandomPP)
                               stats <- rowSums(sparseRandomPPs)
                               results <- results + as.numeric(stats >= vertexDegrees)
                               pb$tick(1)
                             }
                             
                             pvals <- (results + 1) / (iters + 1)
                             names(pvals) <- rownames(self$exp$norm)
                             
                             self$genes$pvals <- pvals
                             self$genes$degrees <- degrees
                             self$genes$powers <- vertexDegrees
                             
                             if (retVal) return(pvals)
                           },
                           
                           calculatePairwiseLinearity = function(negToZero=T) {
                             self$pairwise <- pairwiseR2(t(self$exp$full$norm))
                             colnames(self$pairwise) <- rownames(self$pairwise) <- rownames(self$exp$norm)
                             if (negToZero) {
                               self$pairwise[self$pairwise < 0] <- 0
                             }
                           },
                           
                           smartSearchCorners = function(taus = 2^seq(0, -20, -1),
                                                         sisalIter=100,
                                                         dataset="filtered",
                                                         error="norm") {
                             
                             if (is.null(self$exp$filtered$norm)) {
                               message("Could not find filtered dataset, using whole dataset as filtered")
                               self$filterDataset(rownames(self$exp$full$norm))
                             }
                             
                             Y <- t(self$exp$filtered$norm)
                             samples <- nrow(Y)
                             taus_n <- length(taus)
                             k <- self$cellTypeNumber
                             
                             keks <- lapply(taus, function(i) sisal(Y, k, sisalIter, i, returnPlot = F, nonNeg = T))
                             endpoints_only <- lapply(keks, function(kek) kek$endpoints)
                             
                             Ymean <- rowMeans(Y)
                             ref <- endpoints_only[[1]] - Ymean
                             
                             endpoints_only <- lapply(endpoints_only, function(endpoints) {
                               shifted <- endpoints - Ymean
                               cors <- cor(ref, shifted)
                               return(endpoints[, apply(cors, 1, which.max)])
                             })
                             
                             endpoints_array <- array(0, c(samples, k, taus_n))
                             for (i in 1:taus_n) {
                               endpoints_array[, , i] <- endpoints_only[[i]]
                             }
                             
                             starting_point <- rep(1, k)
                             to_change <- 1
                             unchanged <- 0
                             
                             while (unchanged <= k) {
                               starting_mat <- do.call(cbind, lapply(1:k, function(i) {
                                 endpoints_array[, i, starting_point[i]]
                               }))
                               
                               toCheck <- lapply(1:taus_n, function(i) {
                                 starting_mat_copy <- starting_mat
                                 starting_mat_copy[, to_change] <- endpoints_array[, to_change, i]
                                 starting_mat_copy
                               })
                               
                               
                               errors <- sapply(toCheck, function(props) {
                                 self$endpoints <- props
                                 self$endpointsProjection <- t(self$projectiveProjection) %*% props
                                 self$distances <- apply(self$endpointsProjection, 2, function(x) {
                                   shifted <- self$projection - x
                                   dds <- sqrt(colSums(shifted^2))
                                   return(dds)
                                 })
                                 
                                 self$deconvolveByEndpoints(dataset=dataset, error=error)
                                 self$deconvolutionError(dataset=dataset, error=error)
                               })
                               
                               new_position <- which.min(errors)
                               if (starting_point[to_change] == new_position) {
                                 unchanged <- unchanged + 1
                               } else {
                                 starting_point[to_change] = new_position
                               }
                               
                               to_change <- to_change %% k + 1
                             }
                             
                             message("Final vector is ")
                             message(cat(starting_point))
                             
                             starting_mat <- do.call(cbind, lapply(1:k, function(i) {
                               endpoints_array[, i, starting_point[i]]
                             }))
                             
                             self$endpoints <- starting_mat
                             rownames(self$endpoints) <- colnames(self$exp$filtered$norm)
                             colnames(self$endpoints) <- paste0("Pure gene ", 1:k)
                             self$endpointsProjection <- t(self$projectiveProjection) %*% self$endpoints
                             
                             self$distances <- apply(self$endpointsProjection, 2, function(x) {
                               shifted <- self$projection - x
                               dds <- sqrt(colSums(shifted^2))
                               return(dds)
                             })
                             
                           }, 
                           
                           
                           filterDatasetByPval = function(pval=0.001) {
                             message(sprintf("Total number of genes is %d", nrow(self$exp$full$norm)))
                             geneSubset <- rownames(self$exp$full$norm[self$genes$pvals < pval, ])
                             self$filterDataset(geneSubset)
                             message(sprintf("The number of genes after filtering is %d", nrow(self$exp$filtered$norm)))
                           },

                           filterDataset = function(geneSubset) {
                             self$exp$filtered$raw <- self$exp$full$raw[geneSubset, ]
                             self$exp$filtered$norm <- self$exp$full$norm[geneSubset, ]
                           },
                           
                           tsnePlot = function(dataset="filtered", error="norm") {
                             data <- get(dataset, self$exp)
                             selected <- get(error, data)
                             
                             tsne <- Rtsne(selected, perplexity = 100, max_iter = 2000)
                             toPlot <- data.frame(
                               tSNE1=tsne$Y[, 1],
                               tSNE2=tsne$Y[, 2],
                               marker=NA,
                               row.names = rownames(selected)
                             )
                             
                             if (!is.null(lo$markers)) {
                               for (i in 1:length(lo$markers)) {
                                   toPlot[lo$markers[[i]], "marker"] <- i
                               }
                             }
                             
                             
                             toPlot$marker <- as.factor(toPlot$marker)
                             ggplot(data=toPlot, aes(x=tSNE1, y=tSNE2, color=marker)) +
                               geom_point(size=0.5) + theme_bw(base_size=8) +
                               theme(legend.key.size = unit(0.1, "in"), aspect.ratio = 1) +
                               guides(color=guide_legend(title="Simplex\ncorner"))
                             
                           },
                           
                           significancePlot = function(threshold=0.001) {
                             toPlot <- as.data.frame(self$genes)
                             toPlot$name <- rownames(self$exp$full$norm)
                             toPlot$significant <- toPlot$pvals < threshold
                             
                             ggplot(data=toPlot, aes(x=degrees, y=powers, color=significant)) +
                               geom_point() + theme_bw(base_size=8) + 
                               scale_color_manual(values=c("grey", "red"))
                           }
                           
                           
                         )
)
