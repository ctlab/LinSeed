#' Draw a plot of estimated proportions
#'
#' Draws a plot of estimated proprotions
#' If ggplot2 and reshape2 are installed will use them and return ggplot object
#' Otherwise will use standart R functions
#'
#'@param ... matricies, data frames, NMF objects of estimated proportions or paths to file
#'@param pnames experiment titles
#'@return ggplot object
#'
#'@examples
#' set.seed(31)
#' data('datasetLiverBrainLung')
#' mixed <- datasetLiverBrainLung[, 10:42]
#' preprocessed <- preprocessDataset(mixed, k=5) # 5 clusters
#' accuracy <- clusdecAccuracy(preprocessed, 3) # assuming 3 cell types
#' results <- chooseBest(preprocessed, accuracy) # choose best combination of clusters as putative sigantures
#'
#' plotProportions(results$H) # plotting proportions
#'
#' data('proportionsLiverBrainLung') # or we can compare with actual mixed proportions
#' mixedProportions <- proportionsLiverBrainLung[, 10:42]
#' plotProportions(results$H, mixedProportions, pnames=c('ClusDec estimated', 'Actual'))
#'
#'@export
plotProportions <- function(..., pnames = NULL, point_size=2, line_size=1) {
    requireNamespace("ggplot2")
    requireNamespace("reshape2")


    proportions <- list(...)
    proportions <- lapply(proportions, toMatrix)

    newRowNames <- do.call(function(...) {
        mapply(function(...) {
            dots <- list(...)
            rn <- unlist(dots)
            paste0(rn, collapse = "\n")
        }, ...)
    }, lapply(proportions, rownames))

    proportions <- lapply(proportions, function(p) {
        rownames(p) <- newRowNames
        p
    })

    names(proportions) <- pnames


    cellTypes <- nrow(proportions[[1]])
    results.m <- reshape2::melt(proportions)
    results.m[, 4] <- as.factor(results.m[, 4])

    gplot <- ggplot2::ggplot(results.m,
                             ggplot2::aes(x = as.numeric(Var2),
                                          y = value,
                                          fill = Var1,
                                          color = L1)) +
        ggplot2::geom_line(size=line_size) +
        ggplot2::geom_point(size=point_size) +
        ggplot2::scale_x_discrete(labels = colnames(proportions[[1]])) +
        ggplot2::facet_grid(Var1 ~ .) +
        ggplot2::ylab("proportions") +
        ggplot2::ylim(0, 1.1) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 45,
                                                           hjust = 1)) +
        ggplot2::guides(fill = FALSE)
    if (length(proportions) > 1) {
        gplot <- gplot + ggplot2::theme(legend.title = ggplot2::element_blank(),
            legend.position = "top")

    } else {
        gplot <- gplot + ggplot2::theme(legend.position = "none")
    }
    gplot
}



toMatrix <- function(x) {
    if (class(x) == "data.frame")
        return(as.matrix(x))
    if (inherits(x, "NMF"))
        return(toMatrix(coef(x)))
    if (class(x) == "character")
        return(toMatrix(read.table.mine(x)))
    if (class(x) == "matrix")
        return(x)
    stop("invalid type for plotting")
}

