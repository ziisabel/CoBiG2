###################
### DESeq Analysis
###################

# load count file and convert to integer matrix
genes_matrix_K <- as.matrix(read.table("KNO3_genes.counts.matrix"))
forceMatrixToInteger <- function(m){
  apply (m, c (1, 2), function (x) {
    (as.integer(x))
  })
}
matrix_int_K <- forceMatrixToInteger(genes_matrix_K)
# eliminate all entries undiferentially expressed
matrix_int_K <- matrix_int_K[rowSums(matrix_int_K[, -1])>0, ]

# change column names
colnames(matrix_int_K) <- c("K_0", "K_200", "K_400", "K_600")

# define features and coldata
condition_K <- factor(c("0mM", "200mM", "400mM", "600mM"))

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("DESeq")
library(DESeq)

# DESeq analysis
cds_K <- newCountDataSet(matrix_int_K,condition_K)
cds_K <- estimateSizeFactors(cds_K)
sizeFactors(cds_K)
head(counts(cds_K,normalized=TRUE))
cds_K <- estimateDispersions(cds_K, method="blind", sharingMode="fit-only")
plotDispEsts(cds_K, main="DESeq: KNO3 Per-gene dispersion estimates")
str(fitInfo(cds_K))
head(fData(cds_K))
dispTable(cds_K)

library(gplots)
library(RColorBrewer)
# Heatmap
vsd_K <- varianceStabilizingTransformation(cds_K)
select_K <- order(rowMeans(counts(cds_K)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
heatmap.2(exprs(vsd_K)[select_K,], col = hmcol, trace="none", margin=c(6, 10), main="DESeq: KNO3 Heatmap")
heatmap.2(exprs(vsd_KN)[select_KN,], col = hmcol, trace="none", margin=c(6, 10), Rowv=NA, key.title = "Color key and histogram")

a <- merge("TRINITY_DN2177_c0_g1", resK2["id"])

exprs(vsd_KN)[select_KN,]

a <- merge(upK6, upN6, by="id")
a <- a[order(a$log2FoldChange.x),]
a <- a[1:10,]

a0 <- merge(downK6, downN6, by="id")
a0 <- a0[order(a0$log2FoldChange.x),]
a0 <- a0[1:10,]
  
# PCA plot
plotPCA(vsd_K)

#################
### NaCl = 200mM
#################

resK2 <- nbinomTest(cds_K, "0mM", "200mM")
# filter lines with only NA entries
resK2 <- na.omit(resK2)
plotMA(resK2)
hist(resK2$pval, breaks=100, col="skyblue", border="slateblue", main="DESeq: KNO3 200mM p-value frequency across genes")

# most significantly differentially expressed genes
resSigK2 <- resK2[resK2$pval < 0.05,]
head(resSigK2[order(resSigK2$pval),])

# most strongly down-regulated of the significant genes
downK2 <- resSigK2[resSigK2$log2FoldChange < -2,]
head(downK2[order(downK2$log2FoldChange),])

# most strongly up-regulated ones
upK2 <- resSigK2[resSigK2$log2FoldChange > 2,]
head(upK2[order(-upK2$log2FoldChange),])

# Make a basic volcano plot
with(resK2, plot(log2FoldChange, -log10(pval), pch=20, main="KNO3 0mM vs. 200mM Volcano plot ", xlim=c(-10,10)))
# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.05)
with(subset(resK2, pval<.05 ), points(log2FoldChange, -log10(pval), pch=20, col="royalblue3"))
with(subset(resK2, pval<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pval), pch=20, col="red3"))

#################
### NaCl = 400mM
#################

resK4 <- nbinomTest(cds_K, "0mM", "400mM")
# filter lines with only NA entries
resK4 <- na.omit(resK4)
plotMA(resK4)
hist(resK4$pval, breaks=100, col="skyblue", border="slateblue", main="DESeq: KNO3 400mM p-value frequency across genes")

#most significantly differentially expressed genes
resSigK4 <- resK4[resK4$pval < 0.05,]
head(resSigK4[order(resSigK4$pval),])

#most strongly down-regulated of the significant genes
downK4 <- resSigK4[resSigK4$log2FoldChange < -2,]
head(downK4[order(downK4$log2FoldChange),])

#most strongly up-regulated ones
upK4 <- resSigK4[resSigK4$log2FoldChange > 2,]
head(upK4[order(-upK4$log2FoldChange),])

# Make a basic volcano plot
with(resK4, plot(log2FoldChange, -log10(pval), pch=20, main="DESeq: KNO3 400mM Volcano plot", xlim=c(-10,10)))
# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.05)
with(subset(resK4, pval<.05 ), points(log2FoldChange, -log10(pval), pch=20, col="royalblue3"))
with(subset(resK4, pval<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pval), pch=20, col="red3"))


#################
### NaCl = 600mM
#################

resK6 <- nbinomTest(cds_K, "0mM", "600mM")
# filter lines with only NA entries
resK6 <- na.omit(resK6)
plotMA(resK6)
hist(resK6$pval, breaks=100, col="skyblue", border="slateblue", main="DESeq: KNO3 600mM p-value frequency across genes")

#most significantly differentially expressed genes
resSigK6 <- resK6[resK6$pval < 0.05,]
head(resSigK6[order(resSigK6$pval),])

#most strongly down-regulated of the significant genes (in first condition)
downK6 <- resSigK6[resSigK6$log2FoldChange < -2,]
head(downK6[order(downK6$log2FoldChange),])

#most strongly up-regulated of the significant genes (in first condition)
upK6 <- resSigK6[resSigK6$log2FoldChange > 2,]
head(upK6[order(-upK6$log2FoldChange),])

# Make a basic volcano plot
with(resK6, plot(log2FoldChange, -log10(pval), pch=20, main="DESeq: KNO3 600mM Volcano plot", xlim=c(-10,10)))
# Add colored points: blue if pval<0.01, red if log2FC>2 and pval<0.05)
with(subset(resK6, pval<.05 ), points(log2FoldChange, -log10(pval), pch=20, col="royalblue3"))
with(subset(resK6, pval<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pval), pch=20, col="red3"))

# differentially expressed genes (up and down) for both 200mM and 400mM
join_K2_K4 <- merge(resSigK2, resSigK4, by = "id")
up_K2_K4 <- merge(upK2, upK4, by = "id")
down_K2_K4 <- merge(downK2, downK4, by = "id")

# differentially expressed genes (up and down)  for both 200mM and 600mM
join_K2_K6 <- merge(resSigK2, resSigK6, by = "id")
up_K2_K6 <- merge(upK2, upK6, by = "id")
down_K2_K6 <- merge(downK2, downK6, by = "id")

# differentially expressed genes (up and down) for both 400mM and 600mM
join_K4_K6 <- merge(resSigK4, resSigK6, by = "id")
up_K4_K6 <- merge(upK4, upK6, by = "id")
down_K4_K6 <- merge(downK4, downK6, by = "id")

# differentially expressed genes (up and down) for all salinity conditions
join_K2_K4_K6 <- merge(join_K2_K4, join_K4_K6, by = "id")
up_K2_K4_K6 <- merge(up_K2_K4, upK6, by = "id")
down_K2_K4_K6 <- merge(down_K2_K4, downK6, by = "id")

# differentially expressed genes (up and down) for 200mM and/or 400mM
out_join_K2_K4 <- merge(resSigK2, resSigK4, by="id", all=TRUE)
out_join_K2_K4_K6 <- merge(out_join_K2_K4, resSigK6, by="id", all=TRUE)

up_all_K2_K4 <- merge(upK2, upK4, by="id", all=TRUE)
up_all_K2_K4_K6 <- merge(up_all_K2_K4, upK6, by="id", all=TRUE)

down_all_K2_K4 <- merge(downK2, downK4, by="id", all=TRUE)
down_all_K2_K4_K6 <- merge(down_all_K2_K4, downK6, by="id", all=TRUE)

# genes differentially expressed (up and down) only with 600mM
unique_K6 <- anti_join(resSigK6['id'], out_join_K2_K4['id'])
up_unique_K6 <- anti_join(upK6['id'], up_all_K2_K4['id'])
down_unique_K6 <- anti_join(downK6['id'], down_all_K2_K4['id'])

library(VennDiagram)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# venn diagram between the 3 salinity conditions
venn.diagram(
  x = c(as.list(resSigK2['id']), as.list(resSigK4['id']), as.list(resSigK6['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'DESEq_KNO3.png',
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



# venn diagram between the 3 salinity conditions and up regulated genes
venn.diagram(
  x = c(as.list(upK2['id']), as.list(upK4['id']), as.list(upK6['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'DESEq_KNO3_up.png',
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


# venn diagram between the 3 salinity conditions and down-regulated genes
venn.diagram(
  x = c(as.list(downK2['id']), as.list(downK4['id']), as.list(downK6['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'DESEq_KNO3_down.png',
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
join_K6_N6.up <- merge(upK6, upN6, by = "id")

# Chart
venn.diagram(
  x = c(as.list(upK6['id']), as.list(upN6['id'])),
  category.names = c("KNO3" , "NOD"),
  filename = 'DESEq_up_600mM_KNO3_vs_NOD.png',
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

# DE downregulated genes with 600mM for both KNO3 and NOD genotypes
join_K6_N6.down <- merge(downK6, downN6, by = "id")

# Chart
venn.diagram(
  x = c(as.list(downK6['id']), as.list(downN6['id'])),
  category.names = c("KNO3" , "NOD"),
  filename = 'DESEq_down_600mM_KNO3_vs_NOD.png',
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



####################
# PCA KNO3 + NOD
###################

count_set <- forceMatrixToInteger(cbind(genes_matrix_K, genes_matrix_N))
colnames(count_set) <- c("K_0", "K_200", "K_400", "K_600", "N_0", "N_200", "N_400", "N_600")
condition_KN <- factor(c("0mM", "200mM", "400mM", "600mM", "0mM", "200mM", "400mM", "600mM"))
genotype_KN <- factor(c("KNO3", "KNO3", "KNO3", "KNO3", "NOD", "NOD", "NOD", "NOD"))
coldata_KN <- data.frame(row.names=colnames(count_set), condition_KN, genotype_KN)

cds_KN <- newCountDataSet(count_set, coldata_KN)
cds_KN <- estimateSizeFactors(cds_KN)
sizeFactors(cds_KN)
head(counts(cds_KN,normalized=TRUE))
cds_KN <- estimateDispersions(cds_KN, method="blind", sharingMode="fit-only")
plotDispEsts(cds_KN, main="DESeq: Per-gene dispersion estimates")
str(fitInfo(cds_KN))
head(fData(cds_KN))
dispTable(cds_KN)

vsd_KN <- varianceStabilizingTransformation(cds_KN)

# PCA plot
plotPCA(vsd_KN, intgroup=c("condition_KN", "genotype_KN"))
