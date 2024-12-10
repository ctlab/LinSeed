#' Draw a plot of estimated proportions
#'
#' Draws a plot of estimated proprotions
#' If ggplot2 and reshape2 are installed will use them and return ggplot object
#' Otherwise will use standart R functions
#'
#'@param ... matricies, data frames, NMF objects of estimated proportions or paths to file
#'@param point_size point size for plot
#'@param line_size line size for plot
#'@param pnames experiment titles
#'
#'@return ggplot object
#'
#'
#'@import ggplot2
#'@import reshape2
#'
#'@examples
#'
#'@export
plotProportions <- function(..., pnames = NULL, point_size=2, line_size=1) {
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
    results.m <- melt(proportions)
    results.m[, 4] <- as.factor(results.m[, 4])
    
    results.m <- results.m[sample(nrow(results.m)), ]
    
    gplot <- ggplot(results.m,
                             aes(x = as.numeric(Var2),
                                          y = value,
                                          fill = Var1,
                                          color = L1)) +
        geom_line(size=line_size) +
        geom_point(size=point_size) +
        scale_x_discrete(labels = colnames(proportions[[1]])) +
        facet_grid(Var1 ~ .) +
        ylab("proportions") +
        ylim(0, 1.1) +
        theme_bw() +
        theme(axis.title.x = element_blank(),
                       axis.text.x = element_text(angle = 45,
                                                           hjust = 1)) +
        guides(fill = FALSE)
    if (length(proportions) > 1) {
        gplot <- gplot + theme(legend.title = element_blank(),
            legend.position = "top")

    } else {
        gplot <- gplot + theme(legend.position = "none")
    }
    gplot
}

#' Proportions dot plot
#'
#' @param predicted matrix of predicted proportions
#' @param actual matrix of actual proportions
#' @param main plot title
#' @param guess if True will function will try to guess how to reorder rows of predicted proportions to match rows of actual proportions
#' @param showR2 calculate and show R squared statistics
#'
#' @import ggplot2
#' @import reshape2
#'
#' @return
#' @export
#'
#' @examples
dotPlotPropotions <- function(predicted, actual, guess=FALSE, main=NULL, showR2=FALSE) {
  predicted <- as.matrix(predicted)
  actual <- as.matrix(actual)
  
  if (guess) {
    predicted <- predicted[guessOrder(predicted, actual), ]
  }
  
  colnames(predicted) <- colnames(actual)
  rownames(predicted) <- rownames(actual)
  
  xmelt <- melt(predicted)
  ymelt <- melt(actual)
  
  colnames(ymelt) <- c("Cell Type", "Sample", "Actual")
  colnames(xmelt) <- c("Cell Type", "Sample", "Predicted")
  
  total <- cbind(ymelt, xmelt[, 3, drop=F])
  
  pred <- as.numeric(predicted)
  act <- as.numeric(actual)
  
  r2 <- summary(lm(pred ~ act))$adj.r.squared
  
  pl <- ggplot(data=total, aes(x=Actual, y=Predicted, color=`Cell Type`)) +
    geom_point() + theme_bw(base_size=8) +
    theme(aspect.ratio = 1) + geom_abline(slope=1, intercept = 0, lty=2) +
    xlim(c(0, 1)) + ylim(c(0, 1))
  if (!is.null(main)) {
    pl <- pl + labs(title=main)
  }
  if (showR2) {
    subs <- substitute(italic(R)^2~"="~r2, list(r2=r2))
    pl <- pl + annotate("text", label=as.character(as.expression(subs)), parse=T, x = 0.2, y=0.9)
  }
  pl
}


#' guess the order
#'
#' Function tries to guess ordering for rows of predicted proportions to match rows of actual proportions
#'
#'
#' @importFrom combinat permn
#' @param predicted predicted propotions
#' @param actual actual proportions
#'
#' @return numeric, correct order of predicted proportions
#'
#' @examples
guessOrder <- function(predicted, actual) {
  ctn <- nrow(predicted)
  allPerms <- permn(ctn)
  
  vals <- sapply(allPerms, function(perm) {
    sum(diag(cor(t(predicted[perm, ]), t(actual))))
  })
  perm <- allPerms[[which.max(vals)]]
  return(perm)
}


toMatrix <- function(x) {
    if (is.data.frame(x)) {
        # Convert data frame (or tibble) to a plain matrix
        return(as.matrix(x))
    }
    if (inherits(x, "NMF")) {
        # Handle NMF objects
        return(toMatrix(coef(x)))
    }
    if (is.matrix(x)) {
        # Return if already a matrix
        return(x)
    }
    stop("Invalid type for plotting: ", paste(class(x), collapse = ", "))
}
