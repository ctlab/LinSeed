#' Rotation plotter
#' 
#' Useful rotation plotter for high-dimensional simplexes
#'
#' @param reduced 
#' @param dims 
#' @param step 
#' @param colors 
#' @param fileName 
#'
#' @return
#' @export
#'
#' @examples
plotRotation <- function(reduced, dims, step=5, colors=NULL, fileName="rotation.gif") {
    if (is.null(colors)) {
        for (i in seq(0, 180, step)) {
            png(sprintf("tmpRotation%03d.png", i), width=4, height=4, units="in", res=300)
            s3d <- scatterplot3d(as.matrix(reduced[, (dims-2):dims]),
                                 angle = i, pch = 16, scale.y = 0.8, grid=F, cex.symbols = 0.5)
            dev.off()
        }
    } else {
        for (i in seq(0, 180, step)) {
            png(sprintf("tmpRotation%03d.png", i), width=4, height=4, units="in", res=300)
            s3d <- scatterplot3d(as.matrix(reduced[, (dims-2):dims]),
                                 color=colors,
                                 angle = i, pch = 16, scale.y = 0.8, grid=F, cex.symbols = 0.5)
            dev.off()
        }
    }

    system(paste0("convert -delay 20 tmpRotation*.png ", fileName))
    system("rm tmpRotation*")
}






