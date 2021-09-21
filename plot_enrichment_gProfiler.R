#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

" This script plots enrichment results (ORA from gProfiler, filtered by revigo and processed with clean_gProfiler.py and clean_revigo.py scripts)
Usage in bash: Rscript plot_enrichment.R filename x_labels"

# Install and load required packages
#install.packages("pacman")
#pacman::p_load(devtools, ggplot2)


plotEnrichment <- function(filename, x_labels) {

  "Parameter:
      clean_revigo.py output filename
      vector with the labels conversion to show on the x axis [e.g. c(a=treatment1, b=treatment2)]

  Output: enrichment dotplot
  "

  # read enrichment results as dataframe
  df <- transform(data.frame(as.matrix(read.table(filename, sep="\t",
                                                  header=TRUE, quote="", fill=FALSE))),
                  padj = as.numeric(padj),
                  Counts = as.numeric(Counts),
                  Treatment = as.character(Treatment))

  # optional filters (others may apply)
  df <- filter(df, padj < 0.01)
  df <- filter(df, Counts > 20)

  # order unique terms according to the input file
  df$Terms <- factor(df$Terms,levels=rev(unique(df$Terms)))
  
  # plot results grouped by GO category (BP, MF, CC) [horizontaly] and treatment [verticaly]
  sp2 <- ggplot(df, split=Category, aes(x = Genotype, y = Terms, colour = Regulation, 
                        size = Counts)) +
      geom_point() +
      facet_grid(Category~Treatment, scale="free", space='free')
  
  # determining regulation color coding (down- and up-regulated DEGs)
  sp2+scale_color_manual(values=c('blue', 'red'))+
      # optional renames treatment labels
      scale_x_discrete(labels=x_labels)+ 
      # determine multiple legend order
      guides(colour = guide_legend(order = 3), size = guide_legend(order = 2)) +
      # personalize text elements design
      theme_bw() +
      theme(text = element_text(size=12),
            axis.text.x = element_text(size=12, face='bold'),
            axis.text.y = element_text(size=12))
}


plotEnrichment(args[1, args[2])