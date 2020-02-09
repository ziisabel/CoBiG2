#BiocManager::install("NOISeq")
library(NOISeq)

factors_N <- data.frame(row.names=colnames(matrix_int_N), condition_N)

data_N <- readData(data = matrix_int_N, factors = factors_N)
#str(data_N)
#head(assayData(data_N)$exprs)
#head(pData(data_N))


#################
### NaCl = 200mM
#################

# computing probability od differential expression
N2_noiseq <- noiseq(data_N, factor="condition_N",conditions=c("200mM", "0mM"), k=NULL, norm="n", pnr=0.2, nss=5, v=0.02, lc=1, replicates="no")

# computing DE 
N2_noiseq.dg = degenes(N2_noiseq, q = 0.9, M = NULL)
# filtering down-regulated genes
N2_noiseq.dg.down = degenes(N2_noiseq, q = 0.9, M = "down")
# filtering DE up-regulated genes
N2_noiseq.dg.up = degenes(N2_noiseq, q = 0.9, M = "up")

# ploting expression
DE.plot(N2_noiseq, q = 0.9, graphic = "expr", log.scale = TRUE)
# ploting MD
DE.plot(N2_noiseq, q = 0.9, graphic = "MD")


#################
### NaCl = 400mM
#################

# computing probability od differential expression
N4_noiseq <- noiseq(data_N, factor="condition_N",conditions=c("400mM", "0mM"), k=NULL, norm="n", pnr=0.2, nss=5, v=0.02, lc=1, replicates="no")

# computing DE 
N4_noiseq.dg = degenes(N4_noiseq, q = 0.9, M = NULL)
# filtering down-regulated genes
N4_noiseq.dg.down = degenes(N4_noiseq, q = 0.9, M = "down")
# filtering up-regulated genes
N4_noiseq.dg.up = degenes(N4_noiseq, q = 0.9, M = "up")

# ploting expression
DE.plot(N4_noiseq, q = 0.9, graphic = "expr", log.scale = TRUE)
# ploting MD
DE.plot(N4_noiseq, q = 0.9, graphic = "MD")


#################
### NaCl = 600mM
#################

# computing probability od differential expression
N6_noiseq <- noiseq(data_N, factor="condition_N",conditions=c("600mM", "0mM"), k=NULL, norm="n", pnr=0.2, nss=5, v=0.02, lc=1, replicates="no")

# computing DE 
N6_noiseq.dg = degenes(N6_noiseq, q = 0.9, M = NULL)
# filtering down-regulated genes
N6_noiseq.dg.down = degenes(N6_noiseq, q = 0.9, M = "down")
# filtering up-regulated genes
N6_noiseq.dg.up = degenes(N6_noiseq, q = 0.9, M = "up")

# ploting expression
DE.plot(N6_noiseq, q = 0.9, graphic = "expr", log.scale = TRUE)
# ploting MD
DE.plot(N6_noiseq, q = 0.9, graphic = "MD")

format_df <- function(df){
  id <- rownames(df)
  output <- cbind(id, df)
  row.names(output) <- 1:length(output[,1])
  return(output)
}

# convert df to matrices with gene names as the first column (instead of row names)
N2_noiseq_1 <- format_df(N2_noiseq.dg)
N4_noiseq_1 <- format_df(N4_noiseq.dg)
N6_noiseq_1 <- format_df(N6_noiseq.dg)
N2_noiseq.dg.down_1 <- format_df(N2_noiseq.dg.down)
N2_noiseq.dg.up_1 <- format_df(N2_noiseq.dg.up)
N4_noiseq.dg.down_1 <- format_df(N4_noiseq.dg.down)
N4_noiseq.dg.up_1 <- format_df(N4_noiseq.dg.up)
N6_noiseq.dg.down_1 <- format_df(N6_noiseq.dg.down)
N6_noiseq.dg.up_1 <- format_df(N6_noiseq.dg.up)

library(dplyr)
# differentially expressed genes for both 200mM and 400mM
noiseq_N2_N4 <- merge(N2_noiseq_1, N4_noiseq_1, by = "id")
up_noiseq_N2_N4 <- merge(N2_noiseq.dg.up_1, N4_noiseq.dg.up_1, by = "id")
down_noiseq_N2_N4 <- merge(N2_noiseq.dg.down_1, N4_noiseq.dg.up_1, by = "id")

# differentially expressed genes for both 200mM and 600mM
noiseq_N2_N6 <- merge(N2_noiseq_1, N6_noiseq_1, by = "id")
up_noiseq_N2_N6 <- merge(N2_noiseq.dg.up_1, N6_noiseq.dg.up_1, by = "id")
down_noiseq_N2_N6 <- merge(N2_noiseq.dg.down_1, N6_noiseq.dg.up_1, by = "id")

# differentially expressed genes for both 400mM and 600mM
noiseq_N4_N6 <- merge(N4_noiseq_1, N6_noiseq_1, by = "id")

noiseq_all_N4_N6 <- merge(N4_noiseq_1, N6_noiseq_1, by="id", all=TRUE)
noiseq_all_K4_K6 <- merge(K4_noiseq_1, K6_noiseq_1, by="id", all=TRUE)
noiseq_all_N46_K46 <- merge(noiseq_all_K4_K6, noiseq_all_N4_N6, by="id")

up_noiseq_N4_N6 <- merge(N4_noiseq.dg.up_1, N6_noiseq.dg.up_1, by = "id")
down_noiseq_N4_N6 <- merge(N4_noiseq.dg.down_1, N6_noiseq.dg.up_1, by = "id")

# differentially expressed genes for all salinity conditions
noiseq_N2_N4_N6 <- merge(noiseq_N2_N4, N6_noiseq_1, by = "id")
noiseq_join_all <- merge(noiseq_out_N2_N4_N6, noiseq_out_join_K2_K4_K6, by = "id")

all_NK_noiseq <- merge(noiseq_N2_N4_N6, noiseq_join_K2_K4_K6, by = "id")
up_noiseq_N2_N4_N6 <- merge(up_noiseq_N2_N4, N6_noiseq.dg.up_1, by = "id")
down_noiseq_N2_N4_N6 <- merge(down_noiseq_N2_N4, N6_noiseq.dg.down_1, by = "id")

# differentially expressed genes for 200mM and/or 400mM
noiseq_out_N2_N4 <- merge(N2_noiseq_1, N4_noiseq_1, by="id", all=TRUE)
noiseq_out_N2_N4_N6 <- merge(noiseq_out_N2_N4, N6_noiseq_1, by="id", all=TRUE)
all_join_noiseq <- merge(noiseq_out_N2_N4_N6, noiseq_out_join_K2_K4_K6, by="id", all=TRUE)
diff_all <- anti_join(all_join_noiseq, all_NK_noiseq, by = "id")

up_noiseq_out_N2_N4 <- merge(N2_noiseq.dg.up_1, N4_noiseq.dg.up_1, by="id", all=TRUE)
up_noiseq_out_N2_N4_N6 <- merge(up_noiseq_out_N2_N4, N6_noiseq.dg.up_1, by="id", all=TRUE)
all_join_up_noiseq <- merge(up_noiseq_out_N2_N4_N6, up_noiseq_out_K2_K4_K6, by="id", all=TRUE)

down_noiseq_out_N2_N4 <- merge(N2_noiseq.dg.down_1, N4_noiseq.dg.down_1, by="id", all=TRUE)
down_noiseq_out_N2_N4_N6 <- merge(down_noiseq_out_N2_N4, N6_noiseq.dg.down_1, by="id", all=TRUE)
all_join_down_noiseq <- merge(down_noiseq_out_N2_N4_N6, down_noiseq_out_K2_K4_K6, by="id", all=TRUE)

# genes differentially expressed only with 600mM
noiseq_unique_N6 <- anti_join(N6_noiseq_1, noiseq_out_N2_N4, by = "id")
up_noiseq_unique_N6 <- anti_join(N6_noiseq.dg.up_1, up_noiseq_out_N2_N4, by = "id")
down_noiseq_unique_N6 <- anti_join(N6_noiseq.dg.down_1, down_noiseq_out_N2_N4, by = "id")

library(VennDiagram)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = c(as.list(N2_noiseq_1['id']), as.list(N4_noiseq_1['id']), as.list(N6_noiseq_1['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'NOISeq_NOD.png',
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
  x = c(as.list(N2_noiseq.dg.up_1['id']), as.list(N4_noiseq.dg.up_1['id']), as.list(N6_noiseq.dg.up_1['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'NOISeq_NOD_up.png',
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
  x = c(as.list(N2_noiseq.dg.down_1['id']), as.list(N4_noiseq.dg.down_1['id']), as.list(N6_noiseq.dg.down_1['id'])),
  category.names = c("200mM" , "400mM" , "600mM"),
  filename = 'NOISeq_NOD_down.png',
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
