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
edgeR_df_K_2 <- as.data.frame(read.table("Trinity_trans.gene.counts.matrix.K-0_vs_K-200.edgeR.DE_results"))

# process data frame
edgeR_df_K_2 <- edgeR_df_K_2[c(3:6)]
edgeR_df_K_2 <- format_df(edgeR_df_K_2)

# filter and order by PValue most significantly differentially expressed genes
edgeR_K2_P005 <- edgeR_df_K_2[edgeR_df_K_2$PValue < 0.05,]
head(edgeR_K2_P005[order(edgeR_K2_P005$PValue),])

# filter most strongly down-regulated of the significant genes
edgeR_K2_down <- edgeR_K2_P005[edgeR_K2_P005$logFC < -2,]
head(edgeR_K2_down[order(-edgeR_K2_down$logFC),])

# filter most strongly up-regulated of the significant genes
edgeR_K2_up <- edgeR_K2_P005[edgeR_K2_P005$logFC > 2,]
head(edgeR_K2_up[order(edgeR_K2_up$logFC),])

# Make a basic volcano plot
with(edgeR_df_K_2, plot(-logFC, -log10(PValue), pch=20, main="edgeR: KNO3 0mM vs. 200mM Volcano plot ", xlim=c(-10,10)))
# Add colored points: blue if PValue<0.05, red if log2FC>2 and PValue<0.05)
with(subset(edgeR_df_K_2, PValue<.05 ), points(-logFC, -log10(PValue), pch=20, col="royalblue3"))
with(subset(edgeR_df_K_2, PValue<.05 & abs(logFC)>2), points(-logFC, -log10(PValue), pch=20, col="red3"))


#################
### NaCl = 400mM
#################

# load edgeR DE_results file and convert them to a data frame
edgeR_df_K_4 <- as.data.frame(read.table("Trinity_trans.gene.counts.matrix.K-0_vs_K-400.edgeR.DE_results"))

# process data frame
edgeR_df_K_4 <- edgeR_df_K_4[c(3:6)]
edgeR_df_K_4 <- format_df(edgeR_df_K_4)

# filter and order by PValue most significantly differentially expressed genes
edgeR_K4_P005 <- edgeR_df_K_4[edgeR_df_K_4$PValue < 0.05,]
head(edgeR_K4_P005[order(edgeR_K4_P005$PValue),])

# filter most strongly down-regulated of the significant genes
edgeR_K4_down <- edgeR_K4_P005[edgeR_K4_P005$logFC < -2,]
head(edgeR_K4_down[order(-edgeR_K4_down$logFC),])

# filter most strongly up-regulated of the significant genes
edgeR_K4_up <- edgeR_K4_P005[edgeR_K4_P005$logFC > 2,]
head(edgeR_K4_up[order(edgeR_K4_up$logFC),])

# Make a basic volcano plot
with(edgeR_df_K_4, plot(-logFC, -log10(PValue), pch=20, main="edgeR: KNO3 0mM vs. 400mM Volcano plot ", xlim=c(-10,10)))
# Add colored points: blue if PValue<0.05, green if log2FC>2 and PValue<0.01)
with(subset(edgeR_df_K_4, PValue<.05 ), points(-logFC, -log10(PValue), pch=20, col="royalblue3"))
with(subset(edgeR_df_K_4, PValue<.05 & abs(logFC)>2), points(-logFC, -log10(PValue), pch=20, col="red3"))


#################
### NaCl = 600mM
#################

# load edgeR DE_results file and convert them to a data frame
edgeR_df_K_6 <- as.data.frame(read.table("Trinity_trans.gene.counts.matrix.K-0_vs_K-600.edgeR.DE_results"))

# process data frame
edgeR_df_K_6 <- edgeR_df_K_6[c(3:6)]
edgeR_df_K_6 <- format_df(edgeR_df_K_6)

# filter and order by PValue most significantly differentially expressed genes
edgeR_K6_P005 <- edgeR_df_K_6[edgeR_df_K_6$PValue < 0.05,]
head(edgeR_K6_P005[order(edgeR_K6_P005$PValue),])

# filter most strongly down-regulated of the significant genes (in first condition)
edgeR_K6_down <- edgeR_K6_P005[edgeR_K6_P005$logFC > 2,]
head(edgeR_K6_down[order(-edgeR_K6_down$logFC),])

# filter most strongly up-regulated of the significant genes (in first condition)
edgeR_K6_up <- edgeR_K6_P005[edgeR_K6_P005$logFC < -2,]
head(edgeR_K6_up[order(edgeR_K6_up$logFC),])

# Make a basic volcano plot
with(edgeR_df_K_6, plot(-logFC, -log10(PValue), pch=20, main="edgeR: KNO3 0mM vs. 600mM Volcano plot ", xlim=c(-10,10)))
# Add colored points: blue if PValue<0.05, green if log2FC>2 and PValue<0.01)
with(subset(edgeR_df_K_6, PValue<.05 ), points(-logFC, -log10(PValue), pch=20, col="royalblue3"))
with(subset(edgeR_df_K_6, PValue<.05 & abs(logFC)>2), points(-logFC, -log10(PValue), pch=20, col="red3"))

BiocManager::install("dplyr")
library(dplyr)

# differentially expressed genes for both 200mM and 400mM
edgeR_join_K2_K4 <- merge(edgeR_K2_P005, edgeR_K4_P005, by = "id")
edgeR_up_K2_K4 <- merge(edgeR_K2_up, edgeR_K4_up, by = "id")
edgeR_down_K2_K4 <- merge(edgeR_K2_down, edgeR_K4_down, by = "id")

# differentially expressed genes for both 200mM and 600mM
edgeR_join_K2_K6 <- merge(edgeR_K2_P005, edgeR_K6_P005, by = "id")
edgeR_up_K2_K6 <- merge(edgeR_K2_up, edgeR_K6_up, by = "id")
edgeR_down_K2_K6 <- merge(edgeR_K2_down, edgeR_K6_down, by = "id")

# differentially expressed genes for both 400mM and 600mM
edgeR_join_K4_K6 <- merge(edgeR_K4_P005, edgeR_K6_P005, by = "id")
edgeR_up_K4_K6 <- merge(edgeR_K4_up, edgeR_K6_up, by = "id")
edgeR_down_K4_K6 <- merge(edgeR_K4_down, edgeR_K6_down, by = "id")

# differentially expressed genes for all salinity conditions
edgeR_join_K2_K4_K6 <- merge(edgeR_join_K2_K4, edgeR_join_K4_K6, by = "id")

all_NK_edgeR <- merge(edgeR_join_N2_N4_N6, edgeR_join_K2_K4_K6, by = "id")

all_K4_K6_edgeR <- merge(edgeR_K4_P005, edgeR_K6_P005, by="id", all=TRUE)
all_N4_N6_edgeR <- merge(edgeR_N4_P005, edgeR_N6_P005, by="id", all=TRUE)
all_N46_K46_edgeR <- merge(all_K4_K6_edgeR, all_N4_N6_edgeR, by="id")

edgeR_up_K2_K4_K6 <- merge(edgeR_up_K2_K4, edgeR_K6_up, by = "id")
edgeR_down_K2_K4_K6 <- merge(edgeR_down_K2_K4, edgeR_K6_down, by = "id")

count(edgeR_K2_P005)

# differentially expressed genes for 200mM and/or 400mM
edgeR_out_join_K2_K4 <- merge(edgeR_K2_P005, edgeR_K4_P005, by="id", all=TRUE)
edgeR_out_join_K2_K4_K6 <- merge(edgeR_out_join_K2_K4, edgeR_K6_P005, by="id", all=TRUE)
all_join_edgeR <- merge(edgeR_out_join_K2_K4_K6, edgeR_out_join_N2_N4_N6, by="id")

edgeR_out_up_K2_K4 <- merge(edgeR_K2_up, edgeR_K4_up, by="id", all=TRUE)
edgeR_out_up_K2_K4_K6 <- merge(edgeR_out_up_K2_K4, edgeR_K6_up, by="id", all=TRUE)
all_join_up_edgeR <- merge(edgeR_out_up_K2_K4_K6, edgeR_out_up_N2_N4_N6, by="id", all=TRUE)

edgeR_out_down_K2_K4 <- merge(edgeR_K2_down, edgeR_K4_down, by="id", all=TRUE)
edgeR_out_down_K2_K4_K6 <- merge(edgeR_out_down_K2_K4, edgeR_K6_down, by="id", all=TRUE)
all_join_down_edgeR <- merge(edgeR_out_down_K2_K4_K6, edgeR_out_down_N2_N4_N6, by="id", all=TRUE)


# genes differentially expressed only with 600mM
edgeR_unique_K6 <- anti_join(edgeR_K6_P005['id'], edgeR_out_join_K2_K4['id'])
edgeR_up_unique_K6 <- anti_join(edgeR_K6_up['id'], edgeR_out_up_K2_K4['id'])
edgeR_down_unique_K6 <- anti_join(edgeR_K6_down['id'], edgeR_out_down_K2_K4['id'])

library(VennDiagram)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = c(as.list(edgeR_K2_P005['id']), as.list(edgeR_K4_P005['id']), as.list(edgeR_K6_P005['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'testing005_edgeR_KNO3.png',
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
  x = c(as.list(edgeR_K2_up['id']), as.list(edgeR_K4_up['id']), as.list(edgeR_K6_up['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'testing005_edgeR_KNO3_up.png',
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
  x = c(as.list(edgeR_K2_down['id']), as.list(edgeR_K4_down['id']), as.list(edgeR_K6_down['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'testing005_edgeR_KNO3_down.png',
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

#####################
# 600mM KNO3 vs NOD
#####################

# DE upregulated genes with 600mM for both KNO3 and NOD genotypes
edgeR_join_K6_N6.up <- merge(edgeR_K6_up, edgeR_N6_up, by = "id")

# Chart
venn.diagram(
  x = c(as.list(edgeR_K6_up['id']), as.list(edgeR_N6_up['id'])),
  category.names = c("KNO3" , "NOD"),
  filename = 'testing005_edgeR_up_600mM_KNO3_vs_NOD.png',
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
  fill = c("royalblue3", "seagreen"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans"
)


# DE upregulated genes with 600mM for both KNO3 and NOD genotypes
edgeR_join_K6_N6.down <- merge(edgeR_K6_down, edgeR_N6_down, by = "id")

# Chart
venn.diagram(
  x = c(as.list(edgeR_K6_down['id']), as.list(edgeR_N6_down['id'])),
  category.names = c("KNO3" , "NOD"),
  filename = 'testing005_edgeR_down_600mM_KNO3_vs_NOD.png',
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
  fill = c("royalblue3", "seagreen"),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans"
)


