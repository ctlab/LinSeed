# # utils
# 
# write.table.mine <- function(data, fn, ...) write.table(data, fn, sep = "\t", quote = FALSE,
#     col.names = NA, ...)
# read.table.mine <- function(fn, ...) read.table(fn, sep = "\t", header = 1, row.names = 1,
#     ...)
# norm.length <- function(r) r/sqrt(sum(r^2))
# norm.length.matrix <- function(m) t(apply(m, 1, norm.length))
# norm.relative <- function(r) (r - min(r))/(max(r) - min(r))
# norm.relative.matrix <- function(m) t(apply(m, 1, norm.relative))
# similarity.cosine <- function(x, y) crossprod(x, y)/sqrt(crossprod(x) * crossprod(y))
# space.metric <- function(x, y) 1 - similarity.cosine(x, y)

#' linearizeDataset
#'
#' @param ge gene expression matrix
#'
#' @return gene expression matrix in linear scale
linearizeDataset <- function(ge) {
    if (is_logscale(ge))
        return(2^ge - 1)
    return(ge)
}

#' logDataset
#'
#' @param ge gene expression matrix
#'
#' @return gene expression matrix in log scale
logDataset <- function(ge) {
    if (is_logscale(ge))
        return(ge)
    return(log2(ge + 1))
}


#' is_logscale
#'
#' @param x gene expression matrix
#'
#' @return logical, whether x is in the log scale
is_logscale <- function(x) {
    qx <- quantile(as.numeric(x), na.rm = T)
    if (qx[5] - qx[1] > 100 || qx[5] > 100) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}


# r2Profiler <- function(r2Table,
#                        rCheck = seq(0, 1, 0.1),
#                        kCheck=c(1, 3, 5, 10)) {
#         # rTemp <- (r2Table + t(r2Table)) / 2
#         results <- do.call(cbind, lapply(rCheck, function(rTreshold) {
#             rMask <- (r2Table > rTreshold)
#             rMaskSums <- rowSums(rMask)
#             sapply(kCheck, function(kn) {
#                 sum(rMaskSums > kn)
#             })
#         }))
#         colnames(results) <- rCheck
#         rownames(results) <- kCheck
#         return(results)
# }
# 
# visualizeProfilerResults <- function(results) {
#     requireNamespace("ggplot2")
#     requireNamespace("reshape2")
#     melted <- reshape2::melt(t(results))
#     colnames(melted) <- c("R2", "K", "candidates")
# 
#     pl <- ggplot2::ggplot(data=melted, ggplot2::aes(x=R2, y=log10(candidates),
#                                   group=as.factor(K), color=as.factor(K))) +
#         ggplot2::geom_point() + ggplot2::geom_line() +
#         ggplot2::geom_hline(yintercept = log10(500), color="green") +
#         ggplot2::geom_hline(yintercept = log10(1000), color="red") +
#         ggplot2::theme_bw()
#     pl
# }
# 
# r2Filtering <- function(r2Table, r2val, k) {
#     if (r2val > 1) stop("Value of R^2 to filter should be less or equal")
#     if (r2val < 0) warning("Negative values of R^2 might be not meaningful")
#     if (k <= 0) stop("Number of neighbours should be positive")
#     rSums <- rowSums(r2Table > r2val)
#     filteredGenes <- rownames(r2Table)[rSums > k]
#     filteredOutGenes<- rownames(r2Table)[rSums <= k]
#     subset <- r2Table[filteredGenes, filteredGenes]
#     return(list(filteredGenes=filteredGenes,
#                 filteredOutGenes=filteredOutGenes,
#                 r2Filtered=subset))
# }