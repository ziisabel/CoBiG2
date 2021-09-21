#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

"This function performs Differential Expression (DE) analysis with edgeR
Usage in bash: Rscript edgeR_DE_analysis.R gene_counts.txt treatment_name treatment_n control_name control_n output_directory/"

# Install and load required packages and scripts
#install.packages("pacman")
#pacman::p_load(edgeR)

edgeR_DE_analysis <- function(filename, treatment, n_T, control, n_C, outdir) {

  "Parameters:
        gene counts filename
        treament sample name
        number of treatment samples
        control sample name
        number of control samples
        output directory for DE results
 
  Output:
        text file (.csv) with significant DE results (FDR < 0.01)
        .pdf file with MA and volcano plots with significant DE results (FDR < 0.01)"

    # read file with gene counts (tab-delimited file)
    data = read.table(filename, header=T, row.names=1, com='')

    # optional filter columns with only the samples in intererest
    #data = data[,col_ordering]

    # convert data to integer counts
    data = round(data)

    # filter absente or lowly expressed genes
    # use counts-per-million (cpm) to account for differences in library sizes between samples
    data = data[rowSums(cpm(data) > 1) >= 2,]

    # join conditions names according to the number of replicates
    conditions_reps <- c(rep(treatment, n_T), rep(control, n_C))
    conditions = factor(conditions_reps)

    # store data in a simple list-based data using a grouping factor
    exp_study = DGEList(counts=data, group=conditions)

    # perform TMM normalization to estimate effective library sizes
    exp_study = calcNormFactors(exp_study)

    # define the design matrix based on the experimental design
    condition <- factor(conditions_reps)
    data.frame(Sample=colnames(exp_study),condition)
    design <- model.matrix(~condition)
    rownames(design) <- colnames(exp_study)

    # check if the samples include replicates
    # estimate the quantile-adjusted conditional maximum likelihood (qCML) common dispersion accordingly
    if (n_T > 1 & n_C > 1) {
    dispersion <- estimateGLMCommonDisp(exp_study, design, verbose=TRUE);
    }
    else {
    design0 <- matrix(1,2,1)
    dispersion <- estimateGLMCommonDisp(exp_study,design0,method="deviance",robust=TRUE);
    }

    # Test for DE genes based on the qCML method (applicable to experiments with a single factor)
    conditions <- unique(conditions_reps)
    et <- exactTest(exp_study, pair=conditions, dispersion=dispersion)

    # save all DE results into a table
    tTags <- topTags(et,n=NULL)
    result_table = tTags$table
    result_table = data.frame(sampleA=conditions[1], sampleB=conditions[2], result_table)

    # if the samples include replicares, filter significant results (FDR < 0.01)
    if (replicates == TRUE) {
    result_table_filtered = result_table[result_table$FDR < 0.01, ];
    }
    # else filter by log2FC
    else {
    result_table_filtered = result_table[result_table$logFC > 2 | result_table$logFC < -2, ]
    }

    # name the output files (DE results) based on the treated and control sample names
    output_filename <- paste(treatment, control, sep='_vs_')

    # write signifcant results in a text file (.tab)
    write.table(result_table_filtered, file= paste0(paste(outdir, output_filename, sep = '/'), '.tab'), sep='\t', quote=F, row.names=T)

    # output an MA and a volcano plot with significant results (.pdf)
    pdf(paste0(paste(outdir, output_filename, sep = '/'), '.pdf'))
    plot_MA_and_Volcano(rownames(result_table_filtered), result_table_filtered$logCPM, result_table_filtered$logFC, result_table_filtered$FDR)
    dev.off()

}

edgeR_DE_analysis(args[1], args[2], args[3], args[4], args[5], args[6])