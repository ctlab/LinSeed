#' Preprocess Dataset
#'
#' Preprocesses given dataset. Preprocessing consists of 3 major steps:
#' 1) If needed, probes corresponding to the same genes are collapsed, only most expressed probe is taken for further analysis.
#'    It's common technique in microarray data analysis.
#' 2) If needed, only highly expressed genes are taken for further analysis. (Say hello to noize reduction)
#' 3) All genes are clustered with Kmeans using cosine simillarity as distance.
#'
#' @param dataset matrix, data.frame, path to file or GSE accession with expression data
#' @param annotation dataframe, matrix, named vector with annotation to probes
#' @param geneSymbol column from annotation to collapse the genes, deafult value is 'Gene Symbol'
#' @param samples character vector of samples. If column were not in samples, it would be excluded from analysis.
#' Default value is NULL, which takes every sample from dataset
#' @param topGenes integer How many genes include in analysis. We suppose to include only expressed genes. Default value is 10000
#'
#' @return clustered dataset, matrix, first column identifies cluster of the row
#' @import methods
#' @export
preprocessDataset <- function(dataset, annotation = NULL, geneSymbol = "Gene symbol",
    samples = NULL, topGenes = 10000) {
    if (inherits(dataset, "character")) {
        if (file.exists(dataset)) {
            message("File ", dataset, " exists")
            message("Reading dataset from file ", dataset)
	    message("Make sure file is tab-separated and has row and column names")
            dataset <- read.table(dataset, header=1, row.names=1, sep="\t")
	    message("File successfully read")
        } else {
            stop("File does not exist: ", dataset)
        }
    }
    if (inherits(dataset, "data.frame")) {
        dataset <- as.matrix(dataset)
    }

    if (!inherits(dataset, "matrix")) {
        stop("Unsupported type of dataset: please ensure first argument is matrix, data.frame, path to file or GSE accesssion")
    }

    # sample selection
    if (!is.null(samples)) {
        dataset <- dataset[, samples]
    }

    # annotating if necessary
    if (!is.null(annotation)) {
        fdata <- annotation[, geneSymbol, drop = FALSE]
        dataset <- collapseGenes(dataset, fdata)
    }

    # removing zeroes
    dataset <- dataset[!(rowSums(dataset) == 0), ]
    topGenes <- min(topGenes, nrow(dataset))
    dataset <- logDataset(dataset)
    topRows <- order(rowSums(dataset), decreasing = TRUE)[1:topGenes]
    topDataset <- dataset[topRows, ]
    
    # clustering in linear space
    topDataset <- linearizeDataset(topDataset)
    topDataset <- topDataset[!duplicated(topDataset), ]
    # clustered <- clusterCosine(topDataset, k)
    return(topDataset)

}

#' Preprocess GSE Dataset
#'
#' Downloads GSE dataset by GEO accession and performs preprocessing
#'
#' @param geoAccesion e.g 'GSE19830'
#' @param annotate annotate with feature data from provided geo platform
#' @param normalize quantile normalize GEO dataset
#' @param ... arguments further passed to preprocessDataset
#'
#' @return clustered dataset, matrix, first column identifies cluster of the row
#' @import GEOquery
#' @importFrom Biobase exprs
#' @import preprocessCore
#' @export
preprocessGSE <- function(geoAccesion, annotate = TRUE, normalize=TRUE, ...) {
    gse <- getGEO(geoAccesion, AnnotGPL = T)
    if (length(gse) > 1) {
        stop("This GSE has multiple expression sets. It's probably multiseries. Provide single series experiment")
    }
    
    gse <- gse[[1]]
    expressionData <- Biobase::exprs(gse)

    if (normalize) {
        expressionData <- logDataset(expressionData)
        expressionDataCopy <- normalize.quantiles(expressionData)
        colnames(expressionDataCopy) <- colnames(expressionData)
        rownames(expressionDataCopy) <- rownames(expressionData)
        expressionData <- expressionDataCopy
    }

    if (annotate) {
        preprocessDataset(expressionData, annotation = fData(gse), ...)
    } else {
        preprocessDataset(expressionData, ...)
    }


}

