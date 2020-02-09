###################
### DESeq Analysis
###################

# load count file and convert to integer matrix
genes_matrix_N <- as.matrix(read.table("NOD_genes.counts.matrix"))
forceMatrixToInteger <- function(m){
  apply (m, c (1, 2), function (x) {
    (as.integer(x))
  })
}
matrix_int_N <- forceMatrixToInteger(genes_matrix_N)
# eliminate all entries undiferentially expressed
matrix_int_N <- matrix_int_N[rowSums(matrix_int_N[, -1])>0, ]

# change column names
colnames(matrix_int_N) <- c("N_0", "N_200", "N_400", "N_600")

# define features and coldata
condition_N <- factor(c("0mM", "200mM", "400mM", "600mM"))

# DESeq analysis
library(DESeq)
cds_N <- newCountDataSet(matrix_int_N,condition_N)
cds_N <- estimateSizeFactors(cds_N)
sizeFactors(cds_N)
head(counts(cds_N,normalized=TRUE))
cds_N <- estimateDispersions(cds_N, method="blind", sharingMode="fit-only")
plotDispEsts(cds_N, main="DESeq: NOD per-gene dispersion estimates")
str(fitInfo(cds_N))
head(fData(cds_N))

dispTable(cds_N)

library(gplots)
library(RColorBrewer)
# Heatmap
vsd_N <- varianceStabilizingTransformation(cds_N)
select_N <- order(rowMeans(counts(cds_N)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsd_N)[select_N,], col = hmcol, trace="none", margin=c(10, 10))

# PCA plot
plotPCA(vsd_N)

#################
### NaCl = 200mM
#################

resN2 <- nbinomTest(cds_N, "0mM", "200mM")
# filter lines with only NA entries
resN2 <- na.omit(resN2)
plotMA(resN2)
hist(resN2$pval, breaks=100, col="skyblue", border="slateblue", main="DESeq: NOD 200mM p-value frequency across genes")

#most significantly differentially expressed genes
resSigN2 <- resN2[resN2$pval < 0.1,]
head(resSigN2[order(resSigN2$pval),])

#most strongly down-regulated of the significant genes
downN2 <- resSigN2[resSigN2$log2FoldChange < -2,]
head(downN2[order(downN2$log2FoldChange),])

#most strongly up-regulated ones
upN2 <- resSigN2[resSigN2$log2FoldChange > 2,]
head(upN2[order(-upN2$log2FoldChange),])

# Make a basic volcano plot
with(resN2, plot(log2FoldChange, -log10(pval), pch=20, main="NOD 0mM vs. 200mM Volcano plot", xlim=c(-10,10)))
# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.05)
with(subset(resN2, pval<.05 ), points(log2FoldChange, -log10(pval), pch=20, col="royalblue3"))
with(subset(resN2, pval<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pval), pch=20, col="red3"))


#################
### NaCl = 400mM
#################

resN4 <- nbinomTest(cds_N, "0mM", "400mM")
# filter lines with only NA entries
resN4 <- na.omit(resN4)
plotMA(resN4)
hist(resN4$pval, breaks=100, col="skyblue", border="slateblue", main="DESeq: NOD 400mM p-value frequency across genes")

#most significantly differentially expressed genes
resSigN4 <- resN4[resN4$pval < 0.1,]
head(resSigN4[order(resSigN4$pval),])

#most strongly down-regulated of the significant genes
downN4 <- resSigN4[resSigN4$log2FoldChange < -2,]
head(downN4[order(downN4$log2FoldChange),])

#most strongly up-regulated ones
upN4 <- resSigN4[resSigN4$log2FoldChange > 2,]
head(upN4[order(-upN4$log2FoldChange),])

# Make a basic volcano plot
with(resN4, plot(log2FoldChange, -log10(pval), pch=20, main="NOD 0mM vs. 400mM Volcano plot", xlim=c(-10,10)))
# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.05)
with(subset(resN4, pval<.05 ), points(log2FoldChange, -log10(pval), pch=20, col="royalblue3"))
with(subset(resN4, pval<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pval), pch=20, col="red3"))


#################
### NaCl = 600mM
#################

resN6 <- nbinomTest(cds_N, "0mM", "600mM")
# filter lines with only NA entries
resN6 <- na.omit(resN6)
plotMA(resN6)
hist(resN6$pval, breaks=100, col="skyblue", border="slateblue", main="DESeq: NOD 600mM p-value frequency across genes")

#most significantly differentially expressed genes
resSigN6 <- resN6[resN6$pval < 0.1,]
head(resSigN6[order(resSigN6$pval),])

#most strongly down-regulated of the significant genes
downN6 <- resSigN6[resSigN6$log2FoldChange < -2,]
head(downN6[order(downN6$log2FoldChange),])

#most strongly up-regulated of the significant genes
upN6 <- resSigN6[resSigN6$log2FoldChange > 2,]
head(upN6[order(-upN6$log2FoldChange),])

# Make a basic volcano plot
with(resN6, plot(log2FoldChange, -log10(pval), pch=20, main="NOD 0mM vs. 600mM Volcano plot", xlim=c(-10,10)))
# Add colored points: blue if pval<0.01, red if log2FC>2 and pval<0.05)
with(subset(resN6, pval<.05 ), points(log2FoldChange, -log10(pval), pch=20, col="royalblue3"))
with(subset(resN6, pval<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pval), pch=20, col="red3"))

# differentially expressed genes for both 200mM and 400mM
join_N2_N4 <- merge(resSigN2, resSigN4, by = "id")
up_N2_N4 <- merge(upN2, upN4, by = "id")
down_N2_N4 <- merge(downN2, downN4, by = "id")

# differentially expressed genes for both 200mM and 600mM
join_N2_N6 <- merge(resSigN2, resSigN6, by = "id")
up_N2_N6 <- merge(upN2, upN6, by = "id")
down_N2_N6 <- merge(downN2, downN6, by = "id")

# differentially expressed genes for both 400mM and 600mM
join_N4_N6 <- merge(resSigN4, resSigN6, by = "id")
up_N4_N6 <- merge(upN4, upN6, by = "id")
down_N4_N6 <- merge(downN4, downN6, by = "id")

# differentially expressed genes for all salinity conditions
join_N2_N4_N6 <- merge(join_N2_N4, resSigN6, by = "id")

all_NK <- merge(join_N2_N4_N6, join_K2_K4_K6, by = "id")

all_N4_N6 <- merge(resSigN4, resSigN6, by="id", all=TRUE)
all_K4_K6 <- merge(resSigK4, resSigK6, by="id", all=TRUE)
all_N46_K46 <- merge(all_N4_N6, all_K4_K6, by="id")

up_N2_N4_N6 <- merge(up_N2_N4, upN6, by = "id")
down_N2_N4_N6 <- merge(down_N2_N4, downN6, by = "id")

# differentially expressed genes for 200mM and/or 400mM
out_join_N2_N4 <-merge(resSigN2, resSigN4, by="id", all=TRUE)
out_join_N2_N4_N6 <-merge(out_join_N2_N4, resSigN6, by="id", all=TRUE)
join_all <- merge(out_join_N2_N4_N6, out_join_K2_K4_K6, by="id")
out_join_all <- merge(out_join_N2_N4_N6, out_join_K2_K4_K6, by="id", all=TRUE)

up_out_N2_N4 <-merge(upN2, upN4, by="id", all=TRUE)
up_out_N2_N4_N6 <-merge(up_out_N2_N4, upN6, by="id", all=TRUE)
join_all_up <- merge(up_out_N2_N4_N6, up_all_K2_K4_K6, by="id", all=TRUE)

down_out_N2_N4 <-merge(downN2, downN4, by="id", all=TRUE)
down_out_N2_N4_N6 <-merge(down_out_N2_N4, downN6, by="id", all=TRUE)
join_all_down <- merge(down_out_N2_N4_N6, down_all_K2_K4_K6, by="id", all=TRUE)

library(dplyr)
# genes differentially expressed only with 600mM
unique_N6 <- anti_join(resSigN6['id'], out_join_N2_N4['id'])
unique_N6 <- anti_join(resSigN6['id'], out_join_N2_N4['id'])
up_unique_N6 <- anti_join(upN6['id'], up_out_N2_N4['id'])
down_unique_N6 <- anti_join(downN6['id'], down_out_N2_N4['id'])

library(VennDiagram)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")


# Chart
venn.diagram(
  x = c(as.list(resSigN2['id']), as.list(resSigN4['id']), as.list(resSigN6['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'deseq_NOD.png',
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


# Chart
venn.diagram(
  x = c(as.list(upN2['id']), as.list(upN4['id']), as.list(upN6['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'deseq_NOD_up.png',
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

# Chart
venn.diagram(
  x = c(as.list(downN2['id']), as.list(downN4['id']), as.list(downN6['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'deseq_NOD_down.png',
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
