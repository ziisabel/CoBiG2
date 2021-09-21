#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

"This script performs Differential Expression (DE) analysis with NOISeq
Usage in bash: Rscript noiseq_DE_analysis.R gene_counts.txt factor treatment control"

# Install and load required packages
#install.packages("pacman")
#pacman::p_load(NOISeq)

noiseq_DE_analysis <- function(filename, factor, treatment, control) {
  
    "Parameters:
        gene counts filename
        vector of variables indicating the experimental group for each sample
        treatment sample name
        control sample name
    
    Output:
        text file (.csv) with significant DE results (FDR < 0.01)
        .pdf file with MA and volcano plots with significant DE results (FDR < 0.01)"
    
    #
    mycounts = read.table(filename, header=T, row.names=1, com='')
    
    "Prepare data"
    # define factors
    # the order of the elements must coincide with the order of the data columns
    # multiple conditions can be added
    myfactors = data.frame(condition=rownames(mycounts))
    
    # converte data into a NOISeq object
    # additional annotations can be added (e.g. chromosome, lenght, biotype, etc.)
    mydata <- readData(data=mycounts, factors=myfactors)

    "Computing probability od differential expression when no replicates are available"
    # perform TMM normalization (norm="tmm") and lenght correction (lc=1)
    # replace zero values with the midpoint between 0 and the next non-zero value in the matrix (k=NULL)
    # determine the number of the simulated technical replicates (nss=5)
    # determine the size (percentage of the sequencing depth) of the simulated samples (pnr=0.2)
    # determine the variability of the simulated sample size (v=0.02)
    myresults <- noiseq(mydata, factor=factor, conditions=c(treatment, control), k = NULL, norm = "tmm", pnr = 0.2,
                        nss = 5, v = 0.02, lc = 1, replicates = "no")
    
    "Computing DE"
    # select the DE features for a given threshold q (0.9 is recommendable when no replicates are available)
    # M can be set to "down" or "up" to select only down- or up-regulated genes, respectively
    myresults.deg = degenes(myresults, q = 0.9, M = NULL)

    "Plot the results"
    # plot the summary of the expression values in both conditions (DE highlighted in red)
    DE.plot(myresults, q = 0.9, graphic = "expr", log.scale = TRUE)
    # plot the log2FC (M) and the absolute value of the DE between conditions (D)
    DE.plot(myresults, q = 0.9, graphic = "MD")

    return(myresults.deg)
  
}

noiseq_DE_analysis(args[1], args[2], args[3], args[4])