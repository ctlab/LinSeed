# collapsing genes

#' Collapse Genes
#'
#' Collapses given dataset so every gene is presented by the highest expressed probe
#'
#' @param ge matrix, rows are probes, columns are samples
#' @param probes df, matrix or named vector
#'
#' @return collapsed gene expression matrix with probes replaced with corresponding genes
setGeneric("collapseGenes", function(ge, probes) {
    standardGeneric("collapseGenes")
})

setMethod("collapseGenes", c(ge = "ANY", probes = c("data.frame")), function(ge,
    probes) collapseGenes(ge, as.matrix(probes)))
setMethod("collapseGenes", c("ANY", "matrix"), function(ge, probes) {
    probes <- setNames(as.character(probes[, 1]), rownames(probes))
    collapseGenes(ge, probes)
})

setMethod("collapseGenes", c("ANY", "character"), function(ge, probes) {
    step0 <- length(probes)
    # removing NA, empty string and non-unique mappings
    probes <- probes[!is.na(probes)]
    probes <- probes[probes != ""]
    probes <- probes[!grepl("//", probes)]

    # also removing LOC and orf
    probes <- probes[!grepl("^LOC\\d+", probes, ignore.case = TRUE)]
    probes <- probes[!grepl("^C\\w+orf\\d+", probes, ignore.case = TRUE)]

    step1 <- length(probes)

    message(paste0(step0 - step1, " probes were removed while mapping probes to genes as non-mapped probes or non-uniqely mapped probes"))

    geneToProbes <- lapply(split(probes, probes), names)

    logGE <- logDataset(ge)
    # TODO: needs speed up
    bestProbes <- sapply(geneToProbes, function(probes) {
        subset <- logGE[probes, , drop = FALSE]
        probes[which.max(rowMeans(subset))]
    })

    bestGenes <- ge[bestProbes, ]
    rownames(bestGenes) <- names(geneToProbes)

    step2 <- nrow(bestGenes)
    message(paste0(step1 - step2, " genes were collapsed while mapping genes to most expressed probe"))
    message(paste0(step2, " genes left after collapsing"))
    bestGenes
})



