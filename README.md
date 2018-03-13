Linseed tutorial
================

Linseed (LINear Subspace identification for gene Expresion
Deconvolution) is a package that provides tools and interface to explore
gene expression datasets in linear space.

Load the library
----------------

    library(linseed)

Getting started with linseed
----------------------------

To start working with gene expression data, we need to create a new
LinseedObject, in this tutorial we will use GSE19830 (mixture of Liver,
Brain and Lung), we will take only mixed samples (10-42) and will take
only 10000 most expressed genes.

    lo <- LinseedObject$new("GSE19830", samples=10:42, topGenes=10000)

Coolinearity networks
---------------------

To build a coolinearity network we first have to evaluate all pairwise
collinearity coefficients and then select genes that have at least one
gene that is very collinear ( *R*<sup>2</sup> &gt; =0.99 ) with it and
10 genes that are less collinear ( *R*<sup>2</sup> &gt; =0.95 ). We only
leave genes that meet both requirements. Be very carefull
`calculatePairwiseLinearity` is very memory-demanding method (
*N*<sup>2</sup> memory where *N* is the number of genes).

    lo$calculatePairwiseLinearity(negToZero=T)
    genes1 <- lo$getGenesByCutOff(0.99, 1)
    genes2 <- lo$getGenesByCutOff(0.95, 10)
    goodGenes <- intersect(genes1, genes2)
    subnetwork <- lo$pairwise[goodGenes, goodGenes]

To visualiaze network we can use pheatmap package or use `hclust`
itself.

    library(pheatmap)
    library(RColorBrewer)

    pheatmap(subnetwork, clustering_method = "average",
             show_rownames = F, show_colnames = F, 
             color = colorRampPalette(c("blue", "white", "red"))(100),
             border_color = NA)

![](https://www.dropbox.com/s/bn9dku52angl6al/visi-1.png?raw=1)

Complete deconvolution
----------------------

First step is identification of the number of the cell types presented
in the mixture.

    lo$svdPlot(dataset="norm")

![](https://www.dropbox.com/s/82xs3vaacb0czi9/unnamed-chunk-3-1.png?raw=1)

we can suggest from the figure, that since most of the variance is
explained by first three singular vectors of SVD, that our dataset is as
mixture of 3 cell types.

Projection
----------

Once we know number of cell types, we can project the data into a plane.

    lo$setCellTypeNumber(3)
    lo$project()
    lo$projectionPlot()

![](https://www.dropbox.com/s/gma4krs0y4vfmop/unnamed-chunk-4-1.png?raw=1)

Corner identification
---------------------

Once we projected the data into plane, we can identify simplex corner
using SISAL algorithm and then select genes closest to identified
corners.

    set.seed(1)
    lo$sisalCorners(tau=2^-6, nonNeg=T, iters=300)
    lo$selectGenes(100)
    lo$projectionPlot()

![](/https://www.dropbox.com/s/w8vql525f3c42td/unnamed-chunk-5-1.png?raw=1)

Deconvolution
-------------

Once we selected signature genes, we can use them as an input for DSA
algorithm to perform deconvolution.

    set.seed(1)
    lo$deconvolve(dataset="raw")

    data("proportionsLiverBrainLung")
    actualProportions <- proportionsLiverBrainLung[, 10:42]

    plotProportions(lo$proportions, actualProportions[c(2, 1, 3), ],
                    pnames=c("Linseed", "Actual proportions"))

![](https://www.dropbox.com/s/wz2bcp3ycgu77mr/unnamed-chunk-6-1.png?raw=1)
