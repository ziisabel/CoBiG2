#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

"This script produces multiple plots for DE analysis using DESeq2 (PCA, MA, boxplot, dispersion plot)
Usage in bash: Rscript pca.R samplesName dds pca_colors"

# Install and load required packages
#install.packages("pacman")
#pacman::p_load(DESeq2, ggplot2)

DE_plotter <- function(samplesName, dds, transformation, pca_colors) {

    "Parameters:
        vector with the samples names ordered according to the DESeq object
        DESeq object (DE analysis with DESeq2)
        transformation method (must be rlog or vst)
        vector of colors to personalize the PCA plot

    Output: multiple plots with the DE results (PCA, MA, boxplot, dispersion plot)"

    
    # transform data and remove the dependence of the variance on the mean, unbiasialy

    # transform data with blind regularized logarithm (rlog)
    if transformation == 'rlog':
        transf_dds <- rlog(dds, blind=TRUE)
    # alternatively, vst transformation mode can be used for high volumes of data (faster)
    if transformation == 'vst':
        transf_dds <- vst(dds, blind=TRUE)

    # Plot PCA to visualize the overall effect of experimental covariates and batch effects
    pcaData <- plotPCA(transf_dds, intgroup=c("condition"), returnData=TRUE)
    
    # costumize PCA plot
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC2, color=condition)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed() +
        geom_text(aes(label=samplesName), hjust=0.5, vjust=-0.7) +
        scale_color_manual(values=pca_colors)

    # boxplot of the Cook’s distances to perform manual visualization and filtering (diagnostic test for outliers)
    par(mar=c(8,5,2,2))
    boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

    # plot the dispersion estimates
    plotDispEsts(dds)

    # remove the noise associated with log2FC from low count genes
    resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")

    # MA-plot: log2FC attributable to a given variable over the mean of normalized counts 
    drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
    plotMA(resLFC, ylim=c(-2.5,2.5)); drawLines()

}

DE_plotter(args[1], args[2], args[3])
