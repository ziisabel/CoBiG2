format_df <- function(df){
  id <- rownames(df)
  output <- cbind(id, df)
  row.names(output) <- 1:length(output[,1])
  return(output)
}

#################
### NaCl = 200mM
#################

# load edgeR DE_results file and convert them to a data frame
edgeR_df_N_2 <- as.data.frame(read.table("Trinity_trans.gene.counts.matrix.N-0_vs_N-200.edgeR.DE_results"))

# process data frame
edgeR_df_N_2 <- edgeR_df_N_2[c(3:6)]
edgeR_df_N_2 <- format_df(edgeR_df_N_2)

# filter and order by PValue most significantly differentially expressed genes
edgeR_N2_P005 <- edgeR_df_N_2[edgeR_df_N_2$PValue < 0.05,]
head(edgeR_N2_P005[order(edgeR_N2_P005$PValue),])

# filter most strongly down-regulated of the significant genes
edgeR_N2_down <- edgeR_N2_P005[edgeR_N2_P005$logFC > 2,]
head(edgeR_N2_down[order(-edgeR_N2_down$logFC),])

# filter most strongly up-regulated of the significant genes
edgeR_N2_up <- edgeR_N2_P005[edgeR_N2_P005$logFC < -2,]
head(edgeR_N2_up[order(edgeR_N2_up$logFC),])

# Make a basic volcano plot
with(edgeR_df_N_2, plot(logFC, -log10(PValue), pch=20, main="edgeR: NOD 0mM vs. 200mM Volcano plot ", xlim=c(-10,10)))
# Add colored points: blue if PValue<0.05, green if log2FC>2 and PValue<0.01)
with(subset(edgeR_df_N_2, PValue<.05 ), points(logFC, -log10(PValue), pch=20, col="royalblue3"))
with(subset(edgeR_df_N_2, PValue<.05 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="red3"))


#################
### NaCl = 400mM
#################

# load edgeR DE_results file and convert them to a data frame
edgeR_df_N_4 <- as.data.frame(read.table("Trinity_trans.gene.counts.matrix.N-0_vs_N-400.edgeR.DE_results"))

# process data frame
edgeR_df_N_4 <- edgeR_df_N_4[c(3:6)]
edgeR_df_N_4 <- format_df(edgeR_df_N_4)

# filter and order by PValue most significantly differentially expressed genes
edgeR_N4_P005 <- edgeR_df_N_4[edgeR_df_N_4$PValue < 0.05,]
head(edgeR_N4_P005[order(edgeR_N4_P005$PValue),])

# filter most strongly down-regulated of the significant genes
edgeR_N4_down <- edgeR_N4_P005[edgeR_N4_P005$logFC > 2,]
head(edgeR_N4_down[order(edgeR_N4_down$logFC),])

# filter most strongly up-regulated of the significant genes
edgeR_N4_up <- edgeR_N4_P005[edgeR_N4_P005$logFC < -2,]
head(edgeR_N4_up[order(-edgeR_N4_up$logFC),])

# Make a basic volcano plot
with(edgeR_df_N_4, plot(logFC, -log10(PValue), pch=20, main="edgeR: NOD 0mM vs. 400mM Volcano plot ", xlim=c(-10,10)))
# Add colored points: blue if PValue<0.05, green if log2FC>2 and PValue<0.01)
with(subset(edgeR_df_N_4, PValue<.05 ), points(logFC, -log10(PValue), pch=20, col="royalblue3"))
with(subset(edgeR_df_N_4, PValue<.05 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="red3"))


#################
### NaCl = 600mM
#################

# load edgeR DE_results file and convert them to a data frame
edgeR_df_N_6 <- as.data.frame(read.table("Trinity_trans.gene.counts.matrix.N-0_vs_N-600.edgeR.DE_results"))

# process data frame
edgeR_df_N_6 <- edgeR_df_N_6[c(3:6)]
edgeR_df_N_6 <- format_df(edgeR_df_N_6)

# filter and order by PValue most significantly differentially expressed genes
edgeR_N6_P005 <- edgeR_df_N_6[edgeR_df_N_6$PValue < 0.05,]
head(edgeR_N6_P005[order(edgeR_N6_P005$PValue),])

# filter most strongly down-regulated of the significant genes
edgeR_N6_down <- edgeR_N6_P005[edgeR_N6_P005$logFC > 2,]
head(edgeR_N6_down[order(edgeR_N6_down$logFC),])

# filter most strongly up-regulated of the significant genes
edgeR_N6_up <- edgeR_N6_P005[edgeR_N6_P005$logFC < -2,]
head(edgeR_N6_up[order(-edgeR_N6_up$logFC),])

# Make a basic volcano plot
with(edgeR_df_N_6, plot(-logFC, -log10(PValue), pch=20, main="edgeR: NOD 0mM vs. 600mM Volcano plot ", xlim=c(-10,10)))
# Add colored points: blue if PValue<0.05, green if log2FC>2 and PValue<0.01)
with(subset(edgeR_df_N_6, PValue<.05 ), points(-logFC, -log10(PValue), pch=20, col="royalblue3"))
with(subset(edgeR_df_N_6, PValue<.05 & abs(logFC)>2), points(-logFC, -log10(PValue), pch=20, col="red3"))

library(dplyr)
# differentially expressed genes for both 200mM and 400mM
edgeR_join_N2_N4 <- merge(edgeR_N2_P005, edgeR_N4_P005, by = "id")
edgeR_up_N2_N4 <- merge(edgeR_N2_up, edgeR_N4_up, by = "id")
edgeR_down_N2_N4 <- merge(edgeR_N2_down, edgeR_N4_down, by = "id")

# differentially expressed genes for both 200mM and 600mM
edgeR_join_N2_N6 <- merge(edgeR_N2_P005, edgeR_N6_P005, by = "id")
edgeR_up_N2_N6 <- merge(edgeR_N2_up, edgeR_N6_up, by = "id")
edgeR_down_N2_N6 <- merge(edgeR_N2_down, edgeR_N6_down, by = "id")

# differentially expressed genes for both 400mM and 600mM
edgeR_join_N4_N6 <- merge(edgeR_N4_P005, edgeR_N6_P005, by = "id")
edgeR_up_N4_N6 <- merge(edgeR_N4_up, edgeR_N6_up, by = "id")
edgeR_down_N4_N6 <- merge(edgeR_N4_down, edgeR_N6_down, by = "id")

# differentially expressed genes for all salinity conditions
edgeR_join_N2_N4_N6 <- merge(edgeR_join_N2_N4, edgeR_N6_P005, by = "id")
edgeR_up_N2_N4_N6 <- merge(edgeR_up_N2_N4, edgeR_N6_up, by = "id")
edgeR_down_N2_N4_N6 <- merge(edgeR_down_N2_N4, edgeR_N6_down, by = "id")

# differentially expressed genes for 200mM and/or 400mM
edgeR_out_join_N2_N4 <- merge(edgeR_N2_P005, edgeR_N4_P005, by="id", all=TRUE)
edgeR_out_join_N2_N4_N6 <- merge(edgeR_out_join_N2_N4, edgeR_N6_P005, by="id", all=TRUE)

edgeR_out_up_N2_N4 <- merge(edgeR_N2_up, edgeR_N4_up, by="id", all=TRUE)
edgeR_out_up_N2_N4_N6 <- merge(edgeR_out_up_N2_N4, edgeR_N6_up, by="id", all=TRUE)

edgeR_out_down_N2_N4 <- merge(edgeR_N2_down, edgeR_N4_down, by="id", all=TRUE)
edgeR_out_down_N2_N4_N6 <- merge(edgeR_out_down_N2_N4, edgeR_N6_down, by="id", all=TRUE)


# genes differentially expressed only with 600mM
edgeR_unique_N6 <- anti_join(edgeR_N6_P005['id'], edgeR_out_join_N2_N4['id'])
edgeR_up_unique_N6 <- anti_join(edgeR_N6_up['id'], edgeR_out_up_N2_N4['id'])
edgeR_down_unique_N6 <- anti_join(edgeR_N6_down['id'], edgeR_out_down_N2_N4['id'])

BiocManager::install("VennDiagram")
BiocManager::install("RColorBrewer")

library(VennDiagram)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = c(as.list(edgeR_N2_P005['id']), as.list(edgeR_N4_P005['id']), as.list(edgeR_N6_P005['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'testing005_edgeR_NOD.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1920 , 
  width = 1920 , 
  resolution = 1080,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


venn.diagram(
  x = c(as.list(edgeR_N2_up['id']), as.list(edgeR_N4_up['id']), as.list(edgeR_N6_up['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'testing005_edgeR_NOD_up.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1920 , 
  width = 1920 , 
  resolution = 1080,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = c(as.list(edgeR_N2_down['id']), as.list(edgeR_N4_down['id']), as.list(edgeR_N6_down['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'testing005_edgeR_NOD_down.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1920 , 
  width = 1920 , 
  resolution = 1080,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
