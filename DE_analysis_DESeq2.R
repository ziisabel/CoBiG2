#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

"This function performs Differential Expression (DE) analysis with DESeq2
Usage in bash: Rscript deseq2_DE_analysis.R input_directory/ treatment_name treatment_n control_name control_n output_directory/"

# Install and load required packages
#install.packages("pacman")
#pacman::p_load(DESeq2)

deseq2_DE_analysis <- function(indir, treatment, n_T, control, n_C, outdir) {

  "Parameters:
      gene counts directory
      treament sample name (must be of format [0-9][A-E])
      number of treatment samples (int)
      control sample name (must be of format [0-9][A-E])
      number of control samples (int)
      output directory for DE results

  Output: a text file (.csv) with significant DE results (FDR < 0.01)"

  # get treated and control sample filenames and save to a string
  T_grep <- paste(treatmet, "[A-E]", sep="")
  C_grep <- paste(control, "[A-E]", sep="")
  T_C_grep <- paste(T_grep, C_grep, sep="|")

  # specify which files to use based on the sample filenames string
  sampleFiles <- grep(T_C_grep, list.files(directory), value=TRUE)

  # if sample treatment is at the end of the string, reverse vector
  if (grepl(T_grep, sampleFiles[-1])) {
    sampleFiles <- rev(sampleFiles)
  }

  # obtain the condition status (treatment/control)
  sampleCondition <- c(rep(treatment, n_T), rep(control, n_C))
  sampleTable <- data.frame(sampleName = sampleFiles,
                            fileName = sampleFiles,
                            condition = sampleCondition)
  sampleTable$condition <- factor(sampleTable$condition)

  # build the DESeqDataSet (for files produced by HTSeq)
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = directory,
                                        design= ~ condition)

  # differential expression analysis (treatment vs control)
  de <- DESeq(ddsHTSeq)

  # generate results table
  # filter results based on the mean of normalized counts for each gene (FDR < 0.01)
  res01 <- results(de, alpha=0.01)

  # filter DE significant results and order by ascending p value (FDR < 0.01)
  resOrdered <- res01[order(res01$padj),]
  resSig <- subset(resOrdered, padj < 0.01)

  # write results into a .csv file
  write.csv(as.data.frame(resSig), file=paste(outdir, paste(treatment, control, sep="vs"), "_DE_results.csv", sep=""))

  return(de, sampleFiles)

}

deseq2_DE_analysis(args[1], args[2], args[3], args[4], args[5], args[6])