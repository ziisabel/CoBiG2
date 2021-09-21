#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

"This function plots a personalized heatmaps with DE analysis results
Usage in bash: Rscript heatmap.R DE_filename new_colnames rowColorLabs colColorLabs"

# Install and load required packages
#install.packages("pacman")
#pacman::p_load(gplots, RColorBrewer, dendextend)


heatmap_plotter <- function(filename, new_colnames, rowColorLabs, colColorLabs) {

       "Parameters:
              DE filename (.tab)
              vector with heatmap column labels
              vector with heatmap row color labels
              vector with heatmap column color labels
       
       Output: Heatmap plot"

       # read DE file (.tab) as a matrix
       data <- read.table(filename, header = TRUE, stringsAsFactors=FALSE, sep = "\t", dec = ".", row.names = 1, na.strings="unknown")

       # optional edit matrix column names
       colnames(data)<-c(new_colnames)

       # heatmap color code
       coul <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(25))

       # plot heatmap
       hmap2 <- heatmap.2(as.matrix(data),
                            col = coul,
                            symkey=FALSE,
                            density.info="none",
                            trace="none",
                            scale = "row", # scale by row for easier visualization
                            srtCol = 0,
                            ColSideColors=colColorLabs,
                            RowSideColors=rowColorLabs,
                            cexCol=1,
                            offsetCol=2,
                            key = TRUE,
                            keysize = 1.5,
                            key.xlab = NA,
                            key.par=list(mar=c(10.5,6,2,1)),
                            adjCol = c(0.5,0.5)
       )

       # add row color labels legend
       legend(x="topleft",
              legend = unique(data$Terms),
              col = unique(rowColorLabs),
              lty= 1,             
              lwd = 5,           
              cex=.6,
              bty = "o",
              inset = c(0.01,0.01)
       )

       # add column color labels legend
       legend(x="topleft",
              legend = new_colnames,
              col = unique(colColorLabs),
              lty= 1,             
              lwd = 5,           
              cex=.6,
              bty = "o",
              inset = c(0.01,0.14)
       )

}

heatmap_plotter(args[1])