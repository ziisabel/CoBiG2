#BiocManager::install("NOISeq")
library(NOISeq)

factors_K <- data.frame(row.names=colnames(matrix_int_K), condition_K)

data_K <- readData(data = matrix_int_K, factors = factors_K)
#str(data_K)
#head(assayData(data_K)$exprs)
#head(pData(data_K))


#################
### NaCl = 200mM
#################

# computing probability od differential expression
K2_noiseq <- noiseq(data_K, factor="condition_K",conditions=c("200mM", "0mM"), k=NULL, norm="n", pnr=0.2, nss=5, v=0.02, lc=1, replicates="no")

# computing DE 
K2_noiseq.dg = degenes(K2_noiseq, q = 0.9, M = NULL)
# filtering down-regulated genes
K2_noiseq.dg.down = degenes(K2_noiseq, q = 0.9, M = "down")
# filtering DE up-regulated genes
K2_noiseq.dg.up = degenes(K2_noiseq, q = 0.9, M = "up")

# ploting expression
DE.plot(K2_noiseq, q = 0.9, graphic = "expr", log.scale = TRUE)
# ploting MD
DE.plot(K2_noiseq, q = 0.9, graphic = "MD")


#################
### NaCl = 400mM
#################

# computing probability od differential expression
K4_noiseq <- noiseq(data_K, factor="condition_K",conditions=c("400mM", "0mM"), k=NULL, norm="n", pnr=0.2, nss=5, v=0.02, lc=1, replicates="no")

# computing DE 
K4_noiseq.dg = degenes(K4_noiseq, q = 0.9, M = NULL)
# filtering down-regulated genes
K4_noiseq.dg.down = degenes(K4_noiseq, q = 0.9, M = "down")
# filtering up-regulated genes
K4_noiseq.dg.up = degenes(K4_noiseq, q = 0.9, M = "up")

# ploting expression
DE.plot(K4_noiseq, q = 0.9, graphic = "expr", log.scale = TRUE)
# ploting MD
DE.plot(K4_noiseq, q = 0.9, graphic = "MD")


#################
### NaCl = 600mM
#################

# computing probability od differential expression
K6_noiseq <- noiseq(data_K, factor="condition_K",conditions=c("600mM", "0mM"), k=NULL, norm="n", pnr=0.2, nss=5, v=0.02, lc=1, replicates="no")

# computing DE 
K6_noiseq.dg = degenes(K6_noiseq, q = 0.9, M = NULL)
# filtering down-regulated genes
K6_noiseq.dg.down = degenes(K6_noiseq, q = 0.9, M = "down")
# filtering up-regulated genes
K6_noiseq.dg.up = degenes(K6_noiseq, q = 0.9, M = "up")

# ploting expression
DE.plot(K6_noiseq, q = 0.9, graphic = "expr", log.scale = TRUE)
# ploting MD
DE.plot(K6_noiseq, q = 0.9, graphic = "MD")

format_df <- function(df){
  id <- rownames(df)
  output <- cbind(id, df)
  row.names(output) <- 1:length(output[,1])
  return(output)
}

# convert df to matrices with gene names as the first column (instead of row names)
K2_noiseq_1 <- format_df(K2_noiseq.dg)
K4_noiseq_1 <- format_df(K4_noiseq.dg)
K6_noiseq_1 <- format_df(K6_noiseq.dg)
K2_noiseq.dg.down_1 <- format_df(K2_noiseq.dg.down)
K2_noiseq.dg.up_1 <- format_df(K2_noiseq.dg.up)
K4_noiseq.dg.down_1 <- format_df(K4_noiseq.dg.down)
K4_noiseq.dg.up_1 <- format_df(K4_noiseq.dg.up)
K6_noiseq.dg.down_1 <- format_df(K6_noiseq.dg.down)
K6_noiseq.dg.up_1 <- format_df(K6_noiseq.dg.up)

library(dplyr)
# differentially expressed genes for both 200mM and 400mM
noiseq_join_K2_K4 <- merge(K2_noiseq_1, K4_noiseq_1, by = "id")
up_noiseq_K2_K4 <- merge(K2_noiseq.dg.up_1, K4_noiseq.dg.up_1, by = "id")
down_noiseq_K2_K4 <- merge(K2_noiseq.dg.down_1, K4_noiseq.dg.up_1, by = "id")

# differentially expressed genes for both 200mM and 600mM
noiseq_join_K2_K6 <- merge(K2_noiseq_1, K6_noiseq_1, by = "id")
up_noiseq_K2_K6 <- merge(K2_noiseq.dg.up_1, K6_noiseq.dg.up_1, by = "id")
down_noiseq_K2_K6 <- merge(K2_noiseq.dg.down_1, K6_noiseq.dg.up_1, by = "id")

# differentially expressed genes for both 400mM and 600mM
noiseq_join_K4_K6 <- merge(K4_noiseq_1, K6_noiseq_1, by = "id")
up_noiseq_K4_K6 <- merge(K4_noiseq.dg.up_1, K6_noiseq.dg.up_1, by = "id")
down_noiseq_K4_K6 <- merge(K4_noiseq.dg.down_1, K6_noiseq.dg.up_1, by = "id")

# differentially expressed genes for all salinity conditions
noiseq_join_K2_K4_K6 <- merge(noiseq_join_K2_K4, noiseq_join_K4_K6, by = "id")
up_noiseq_K2_K4_K6 <- merge(up_noiseq_K2_K4, K6_noiseq.dg.up_1, by = "id")
down_noiseq_K2_K4_K6 <- merge(down_noiseq_K2_K4, K6_noiseq.dg.down_1, by = "id")

# differentially expressed genes for 200mM and/or 400mM
noiseq_out_join_K2_K4 <- merge(K2_noiseq_1, K4_noiseq_1, by="id", all=TRUE)
noiseq_out_join_K2_K4_K6 <- merge(noiseq_out_join_K2_K4, K6_noiseq_1, by="id", all=TRUE)

up_noiseq_out_K2_K4 <- merge(K2_noiseq.dg.up_1, K4_noiseq.dg.up_1, by="id", all=TRUE)
up_noiseq_out_K2_K4_K6 <- merge(up_noiseq_out_K2_K4, K6_noiseq.dg.up_1, by="id", all=TRUE)

down_noiseq_out_K2_K4 <- merge(K2_noiseq.dg.down_1, K4_noiseq.dg.down_1, by="id", all=TRUE)
down_noiseq_out_K2_K4_K6 <- merge(down_noiseq_out_K2_K4, K6_noiseq.dg.down_1, by="id", all=TRUE)

# genes differentially expressed only with 600mM
noiseq_unique_K6 <- anti_join(K6_noiseq_1['id'], noiseq_out_join_K2_K4['id'])
up_noiseq_unique_K6 <- anti_join(K6_noiseq.dg.up_1, up_noiseq_out_K2_K4, by = "id")
down_noiseq_unique_K6 <- anti_join(K6_noiseq.dg.down_1, down_noiseq_out_K2_K4, by = "id")

library(VennDiagram)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = c(as.list(K2_noiseq_1['id']), as.list(K4_noiseq_1['id']), as.list(K6_noiseq_1['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'NOISeq_KNO3.png',
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
  x = c(as.list(K2_noiseq.dg.up_1['id']), as.list(K4_noiseq.dg.up_1['id']), as.list(K6_noiseq.dg.up_1['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'NOISeq_KNO3_up.png',
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
  x = c(as.list(K2_noiseq.dg.down_1['id']), as.list(K4_noiseq.dg.down_1['id']), as.list(K6_noiseq.dg.down_1['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'NOISeq_KNO3_down.png',
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
