library(dplyr)

######################
### genotype = KNO3
#####################

##### 200mM 

# differentially expressed genes
edgeR_DESeq_K2 <- merge(edgeR_K2_P01, resSigK2, by = "id")
edgeR_DESeq_NOISeq_K2 <- merge(edgeR_DESeq_K2, K2_noiseq_1, by = 'id')

# up-regulated genes
up_NOISeqK2_DESeqK2 <- merge(upK2, K2_noiseq.dg.up_1, by = "id")
up_K_alltools_200 <- merge(up_NOISeqK2_DESeqK2, edgeR_K2_up, by = "id")

# down-regulated genes
down_NOISeqK2_DESeqK2 <- merge(downK2, K2_noiseq.dg.down_1, by = "id")
down_K_alltools_200 <- merge(down_NOISeqK2_DESeqK2, edgeR_K2_down, by = "id")

#‘XML’, ‘annotate’, ‘genefilter’, ‘geneplotter’


##### 400mM 

# differentially expressed genes
edgeR_DESeq_K4 <- merge(edgeR_K4_P01, resSigK4, by = "id")
edgeR_DESeq_NOISeq_K4 <- merge(edgeR_DESeq_K4, K4_noiseq_1, by = 'id')

# up-regulated genes
up_NOISeqK4_DESeqK4 <- merge(upK4, K4_noiseq.dg.up_1, by = "id")
up_K_alltools_400 <- merge(up_NOISeqK4_DESeqK4, edgeR_K4_up, by = "id")

# down-regulated genes
down_NOISeqK4_DESeqK4 <- merge(downK4, K4_noiseq.dg.down_1, by = "id")
down_K_alltools_400 <- merge(down_NOISeqK4_DESeqK4, edgeR_K4_down, by = "id")


##### 600mM 

# differentially expressed genes in K_600 for all tools
edgeR_DESeq_K6 <- merge(edgeR_K6_P01, resSigK6, by = "id")
edgeR_DESeq_NOISeq_K6 <- merge(edgeR_DESeq_K6, K6_noiseq_1, by = 'id')

# up-regulated genes
up_NOISeqK6_DESeqK6 <- merge(upK6, K6_noiseq.dg.up_1, by = "id")
up_K_alltools_600 <- merge(up_NOISeqK6_DESeqK6, edgeR_K6_up, by = "id")
up_K_alltools_600 <- up_K_alltools_600[order(up_K_alltools_600$log2FoldChange),]
topup_alltools_600K <- up_K_alltools_600[1:10,1]
topup_alltools_600K <- gsub("_", " ", topup_alltools_600K)

write.csv(topup_alltools_600K, file = "topup_alltools_600K.csv", row.names=FALSE)

# down-regulated genes
down_NOISeqK6_DESeqK6 <- merge(downK6, K6_noiseq.dg.down_1, by = "id")
down_K_alltools_600 <- merge(down_NOISeqK6_DESeqK6, edgeR_K6_down, by = "id")
down_K_alltools_600 <- down_K_alltools_600[order(down_K_alltools_600$log2FoldChange),]
topdown_alltools_600K <- down_K_alltools_600[1:10,1]
topdown_alltools_600K <- gsub("_", " ", topdown_alltools_600K)
write.csv(topdown_alltools_600K, file = "topdown_alltools_600K.csv", row.names=FALSE)
  
# differentially expressed genes for both 200mM and 400mM
edgeR_DESeq_NOISeq_join_K2_K4 <- merge(edgeR_DESeq_NOISeq_K2, edgeR_DESeq_NOISeq_K4, by = "id")
# differentially expressed genes for both 200mM and 600mM
edgeR_DESeq_NOISeq_join_K2_K6 <- merge(edgeR_DESeq_NOISeq_K2, edgeR_DESeq_NOISeq_K6, by = "id")
# differentially expressed genes for both 400mM and 600mM
edgeR_DESeq_NOISeq_join_K4_K6 <- merge(edgeR_DESeq_NOISeq_K4, edgeR_DESeq_NOISeq_K6, by = "id")
# differentially expressed genes for all salinity conditions
edgeR_DESeq_NOISeq_join_K2_K4_K6 <- merge(edgeR_join_K2_K4, edgeR_join_K4_K6, by = "id")
# differentially expressed genes for 200mM and/or 400mM
edgeR_DESeq_NOISeq_out_join_K2_K4 <- merge(edgeR_DESeq_NOISeq_K2, edgeR_DESeq_NOISeq_K4, by="id", all=TRUE)
# genes differentially expressed only with 600mM
edgeR_DESeq_NOISeq_unique_K6 <- anti_join(edgeR_DESeq_NOISeq_K6['id'], edgeR_DESeq_NOISeq_out_join_K2_K4['id'])



venn.diagram(
  x = c(as.list(edgeR_DESeq_NOISeq_K2['id']), as.list(edgeR_DESeq_NOISeq_K4['id']), as.list(edgeR_DESeq_NOISeq_K6['id'])),
  category.names = c("200mM" , "400mM", "600mM"),
  filename = 'alltools_KNO3.png',
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
  cat.pos = c(-27, 27, 130),
  cat.dist = c(0.055, 0.055, 0.070),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = c(as.list(up_K_alltools_200['id']), as.list(up_K_alltools_400['id']), as.list(up_K_alltools_600['id'])),
  category.names = c("200mM" , "400mM", "600mM"),
  filename = 'up_alltools_KNO3.png',
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
  cat.pos = c(-27, 27, 130),
  cat.dist = c(0.055, 0.055, 0.070),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = c(as.list(down_K_alltools_200['id']), as.list(down_K_alltools_400['id']), as.list(down_K_alltools_600['id'])),
  category.names = c("200mM" , "400mM", "600mM"),
  filename = 'down_alltools_KNO3.png',
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
  cat.pos = c(-27, 27, 130),
  cat.dist = c(0.055, 0.055, 0.070),
  cat.fontfamily = "sans",
  rotation = 1
)


######################
### genotype = NOD
#####################

##### 200mM 

# differentially expressed genes
edgeR_DESeq_N2 <- merge(edgeR_N2_P01, resSigN2, by = "id")
edgeR_DESeq_NOISeq_N2 <- merge(edgeR_DESeq_N2, N2_noiseq_1, by = 'id')

# up-regulated genes
up_NOISeqN2_DESeqN2 <- merge(upN2, N2_noiseq.dg.up_1, by = "id")
up_N_alltools_200 <- merge(up_NOISeqN2_DESeqN2, edgeR_N2_up, by = "id")

# down-regulated genes
down_NOISeqN2_DESeqN2 <- merge(downN2, N2_noiseq.dg.down_1, by = "id")
down_N_alltools_200 <- merge(down_NOISeqN2_DESeqN2, edgeR_N2_down, by = "id")


##### 400mM 

# differentially expressed genes
edgeR_DESeq_N4 <- merge(edgeR_N4_P01, resSigN4, by = "id")
edgeR_DESeq_NOISeq_N4 <- merge(edgeR_DESeq_N4, N4_noiseq_1, by = 'id')

# up-regulated genes
up_NOISeqN4_DESeqN4 <- merge(upN4, N4_noiseq.dg.up_1, by = "id")
up_N_alltools_400 <- merge(up_NOISeqN4_DESeqN4, edgeR_N4_up, by = "id")

# down-regulated genes
down_NOISeqN4_DESeqN4 <- merge(downN4, N4_noiseq.dg.down_1, by = "id")
down_N_alltools_400 <- merge(down_NOISeqN4_DESeqN4, edgeR_N4_down, by = "id")


##### 600mM 

# differentially expressed genes
edgeR_DESeq_N6 <- merge(edgeR_N6_P01, resSigN6, by = "id")
edgeR_DESeq_NOISeq_N6 <- merge(edgeR_DESeq_N6, N6_noiseq_1, by = 'id')
edgeR_DESeq_NOISeq_K6 <- merge(edgeR_DESeq_K6, K6_noiseq_1, by = 'id')
all_deg_600 <- merge(edgeR_DESeq_NOISeq_N6, edgeR_DESeq_NOISeq_K6, by="id", all=TRUE)
all_deg_600 <- all_deg_600[order(all_deg_600$id),]
all_deg_600 <- gsub("_", " ", all_deg_600[,1])
write.csv(all_deg_600[,1], file = "all_deg_600.csv", row.names=FALSE)

# up-regulated genes
up_NOISeqN6_DESeqN6 <- merge(upN6, N6_noiseq.dg.up_1, by = "id")
up_N_alltools_600 <- merge(up_NOISeqN6_DESeqN6, edgeR_N6_up, by = "id")
up_N_alltools_600 <- up_N_alltools_600[order(up_N_alltools_600$log2FoldChange),]
topup_alltools_600N <- up_N_alltools_600[1:10,1]
topup_alltools_600N <- gsub("_", " ", topup_alltools_600N)
write.csv(topup_alltools_600N, file = "topup_alltools_600N.csv", row.names=FALSE)

# down-regulated genes
down_NOISeqN6_DESeqN6 <- merge(downN6, N6_noiseq.dg.down_1, by = "id")
down_N_alltools_600 <- merge(down_NOISeqN6_DESeqN6, edgeR_N6_down, by = "id")
down_N_alltools_600 <- down_N_alltools_600[order(down_N_alltools_600$log2FoldChange),]
topdown_alltools_600N <- down_N_alltools_600[1:10,1]
topdown_alltools_600N <- gsub("_", " ", topdown_alltools_600N)

write.csv(topdown_alltools_600N, file = "topdown_alltools_600N.csv", row.names=FALSE)

# differentially expressed genes for all salinity conditions
edgeR_DESeq_NOISeq_join_N2_N4_N6 <- merge(edgeR_join_N2_N4, edgeR_join_N4_N6, by = "id")
# differentially expressed genes for 200mM and/or 400mM
edgeR_DESeq_NOISeq_out_join_N2_N4 <- merge(edgeR_DESeq_NOISeq_N2, edgeR_DESeq_NOISeq_N4, by="id", all=TRUE)
# genes differentially expressed only with 600mM
edgeR_DESeq_NOISeq_unique_N6 <- anti_join(edgeR_DESeq_NOISeq_N6['id'], edgeR_DESeq_NOISeq_out_join_N2_N4['id'])

# differentially expressed genes for both 200mM and 400mM
edgeR_DESeq_NOISeq_join_N2_N4 <- merge(edgeR_DESeq_NOISeq_N2, edgeR_DESeq_NOISeq_N4, by = "id")
# differentially expressed genes for both 200mM and 600mM
edgeR_DESeq_NOISeq_join_N2_N6 <- merge(edgeR_DESeq_NOISeq_N2, edgeR_DESeq_NOISeq_N6, by = "id")
# differentially expressed genes for both 400mM and 600mM
edgeR_DESeq_NOISeq_join_N4_N6 <- merge(edgeR_DESeq_NOISeq_N4, edgeR_DESeq_NOISeq_N6, by = "id")
# differentially expressed genes for all salinity conditions
edgeR_DESeq_NOISeq_join_N2_N4_N6 <- merge(edgeR_join_N2_N4, edgeR_join_N4_N6, by = "id")
# genes differentially expressed only with 600mM
edgeR_DESeq_unique_N6 <- merge(unique_N6, edgeR_unique_N6, by = "id")
edgeR_DESeq_NOISeq_unique_N6 <- merge(edgeR_DESeq_unique_N6['id'], noiseq_unique_N6['id'])
# up regulated - unique 600
up_edgeR_DESeq_unique_N6 <- merge(up_unique_N6, edgeR_up_unique_N6, by = "id")
up_alltools_unique_N6 <- merge(up_edgeR_DESeq_unique_N6['id'], up_noiseq_unique_N6['id'])
# downregulated - unique 600
down_edgeR_DESeq_unique_N6 <- merge(down_unique_N6, edgeR_down_unique_N6, by = "id")
down_alltools_unique_N6 <- merge(down_edgeR_DESeq_unique_N6['id'], down_noiseq_unique_N6['id'])

edgeR_up_down <- merge(edgeR_N6_up['id'], edgeR_N6_down['id'], all=TRUE)

###### BEST BEST!!!!!!
# intersecao entre os up e down regulated de noiseq com todos os deg de deseq
N2_noiseq_up_down <- merge(N2_noiseq.dg.up_1['id'], N2_noiseq.dg.down_1['id'], all=TRUE)
N4_noiseq_up_down <- merge(N4_noiseq.dg.up_1['id'], N4_noiseq.dg.down_1['id'], all=TRUE)
N6_noiseq_up_down <- merge(N6_noiseq.dg.up_1['id'], N6_noiseq.dg.down_1['id'], all=TRUE)
N2_alldeseq_vs_updown_noiseq <- merge(resSigN2['id'], N2_noiseq_up_down['id'])
# dos 120 deg no deseq 120 sao os mais expressos (up e down) do noiseq
N4_alldeseq_vs_updown_noiseq <- merge(resSigN4['id'], N4_noiseq_up_down['id'])
# dos 698 deg no deseq 695 sao os mais expressos (up e down) do noiseq
N6_alldeseq_vs_updown_noiseq <- merge(resSigN6['id'], N6_noiseq_up_down['id'])
# dos 939 deg no deseq 927 sao os mais expressos (up e down) do noiseq
# deseq e menos sensivel mas podera ter menos falsos positivos pq deteta apenas
# as maiores variacoes de expressao, o que pode ser util neste caso em que nao 
# ha replicados

# intersecao entre os up e down regulated de noiseq com todos os deg de deseq
K2_noiseq_up_down <- merge(K2_noiseq.dg.up_1['id'], K2_noiseq.dg.down_1['id'], all=TRUE)
K4_noiseq_up_down <- merge(K4_noiseq.dg.up_1['id'], K4_noiseq.dg.down_1['id'], all=TRUE)
K6_noiseq_up_down <- merge(K6_noiseq.dg.up_1['id'], K6_noiseq.dg.down_1['id'], all=TRUE)
K2_alldeseq_vs_updown_noiseq <- merge(resSigK2['id'], K2_noiseq_up_down['id'])
# dos 165 deg no deseq 160 sao os mais expressos (up e down) do noiseq
K4_alldeseq_vs_updown_noiseq <- merge(resSigK4['id'], K4_noiseq_up_down['id'])
# dos 880 deg no deseq 879 sao os mais expressos (up e down) do noiseq
K6_alldeseq_vs_updown_noiseq <- merge(resSigK6['id'], K6_noiseq_up_down['id'])
# dos 1121 deg no deseq 1121 sao os mais expressos (up e down) do noiseq
# deseq e menos sensivel mas podera ter menos falsos positivos pq deteta apenas
# as maiores variacoes de expressao, o que pode ser util neste caso em que nao 
# ha replicados
K2_deseq_updw_noiseq <- merge(resK2['id'], K2_noiseq_up_down['id'])
# 958 / 958
K4_deseq_updw_noiseq <- merge(resK4['id'], K4_noiseq_up_down['id'])
# 2237 / 2237
K6_deseq_updw_noiseq <- merge(resK6['id'], K6_noiseq_up_down['id'])
# 2601 / 2601

######## BEST!!!!!
edgeR_up_down_K2 <- merge(edgeR_K2_up, edgeR_K2_down, by="id", all=TRUE)
resSigK2_edgeR_updown <- merge(resSigK2, edgeR_up_down_K2, by="id")
# 93 / 165 = 56%
edgeR_up_down_N2 <- merge(edgeR_N2_up, edgeR_N2_down, by="id", all=TRUE)
resSigN2_edgeR_updown <- merge(resSigN2, edgeR_up_down_N2, by="id")
# 71 / 120 = 59%

edgeR_up_down_K4 <- merge(edgeR_K4_up, edgeR_K4_down, by="id", all=TRUE)
resSigK4_edgeR_updown <- merge(resSigK4['id'], edgeR_up_down_K4['id'])
# 609/880 = 69%
edgeR_up_down_N4 <- merge(edgeR_N4_up, edgeR_N4_down, by="id", all=TRUE)
resSigN4_edgeR_updown <- merge(resSigN4['id'], edgeR_up_down_N4['id'])
# 409/698 = 59%

edgeR_up_down_K6 <- merge(edgeR_K6_up, edgeR_K6_down, by="id", all=TRUE)
resSigK6_edgeR_updown <- merge(resSigK6['id'], edgeR_up_down_K6['id'])
# 779/1121 = 70%
edgeR_up_down_N6 <- merge(edgeR_N6_up, edgeR_N6_down, by="id", all=TRUE)
resSigN6_edgeR_updown <- merge(resSigN6['id'], edgeR_up_down_N6['id'])
# 613/939 = 65%

resN2_deseq_df_edgeR <- merge(resSigN2['id'], edgeR_df_N_2['id'])
# 120 / 120
resN4_deseq_df_edgeR <- merge(resSigN4['id'], edgeR_df_N_4['id'])
# 698 / 697
resN6_deseq_df_edgeR <- merge(resSigN6['id'], edgeR_df_N_6['id'])
# 939 / 935

resK2_updw_edgeR <- merge(resK2['id'], K2_edgeR_up_down['id'])
# 931 / 799
resK4_updw_edgeR <- merge(resK4['id'], K4_edgeR_up_down['id'])
# 1696 / 1565
resK6_updw_edgeR <- merge(resK6['id'], K6_edgeR_up_down['id'])
# 1989 / 1858
resN2_updw_edgeR <- merge(resN2['id'], N2_edgeR_up_down['id'])
# 664 / 769
resN4_updw_edgeR <- merge(resN4['id'], N4_edgeR_up_down['id'])
# 1213 / 1318
resN6_updw_edgeR <- merge(resN6['id'], N6_edgeR_up_down['id'])
# 1345 / 1407

edgeR_up_down <- merge(edgeR_N6_up['id'], edgeR_N6_down['id'], all=TRUE)

# intersecao entre os up e down regulated de noiseq com todos os deg de deseq
N2_edgeR_up_down <- merge(edgeR_N2_up['id'], edgeR_N2_down['id'], all=TRUE)
N4_edgeR_up_down <- merge(edgeR_N4_up['id'], edgeR_N4_down['id'], all=TRUE)
N6_edgeR_up_down <- merge(edgeR_N6_up['id'], edgeR_N6_down['id'], all=TRUE)
N2_deseq_vs_updown_edgeR <- merge(resSigN2['id'], N2_edgeR_up_down['id'])
# dos 120 deg no deseq 71 sao os mais expressos (up e down) do edgeR
N4_deseq_vs_updown_edgeR <- merge(resSigN4['id'], N4_edgeR_up_down['id'])
# dos 698 deg no deseq 409 sao os mais expressos (up e down) do edgeR
N6_deseq_vs_updown_edgeR <- merge(resSigN6['id'], N6_edgeR_up_down['id'])
# dos 939 deg no deseq 613 sao os mais expressos (up e down) do edgeR
# deseq e menos sensivel mas podera ter menos falsos positivos pq deteta apenas
# as maiores variacoes de expressao, o que pode ser util neste caso em que nao 
# ha replicados

# intersecao entre os up e down regulated de noiseq com todos os deg de deseq
K2_edgeR_up_down <- merge(edgeR_K2_up['id'], edgeR_K2_down['id'], all=TRUE)
K4_edgeR_up_down <- merge(edgeR_K4_up['id'], edgeR_K4_down['id'], all=TRUE)
K6_edgeR_up_down <- merge(edgeR_K6_up['id'], edgeR_K6_down['id'], all=TRUE)
K2_deseq_vs_updown_edgeR <- merge(resSigK2['id'], K2_edgeR_up_down['id'])
# dos 165 deg no deseq 93 sao os mais expressos (up e down) do edgeR
K4_deseq_vs_updown_edgeR <- merge(resSigK4['id'], K4_edgeR_up_down['id'])
# dos 880 deg no deseq 609 sao os mais expressos (up e down) do edgeR
K6_deseq_vs_updown_edgeR <- merge(resSigK6['id'], K6_edgeR_up_down['id'])
# dos 1121 deg no deseq 779 sao os mais expressos (up e down) do edgeR
# deseq e menos sensivel mas podera ter menos falsos positivos pq deteta apenas
# as maiores variacoes de expressao, o que pode ser util neste caso em que nao 
# ha replicados


resK2_deseq_updw_edgeR <- merge(resK2['id'], edgeR_df_K_2['id'])
# 931 / 799 (resK2)
# 93 dos 165
# 16923 edgeR df == 16090 resK2
resK4_deseq_updw_edgeR <- merge(resK4['id'], edgeR_df_K_4['id'])
# 1565 / 1696
# 609 dos 880
# 16411 == 15578 resK4
resK6_deseq_updw_edgeR <- merge(resK6['id'], edgeR_df_K_6['id'])
# 1989 / 1858
# 779 dos 1121
# 16870 == 16041


resN2_deseq_updw_edgeR <- merge(resN2['id'], edgeR_df_N_2['id'])
# 664 / 769
# 71 dos 120
# 16159 == 15332
resN4_deseq_updw_edgeR <- merge(resN4['id'], edgeR_df_N_4['id'])
# 1213 / 1318
# 409 dos 698
# 17001 == 16177
resN6_deseq_updw_edgeR <- merge(resN6['id'], edgeR_df_N_6['id'])
# 1345 / 1407
# 613 dos 939
# 16360 == 15538

count(resN6)
count(edgeR_N6_P01)
count(resSigN6)
count(N2_edgeR_up_down)
count(N4_edgeR_up_down)
count(N6_edgeR_up_down)
count(resSigN6)
count(N2_noiseq_up_down)
count(N4_noiseq_up_down)
count(N6_noiseq_up_down)

venn.diagram(
  x = c(as.list(edgeR_DESeq_NOISeq_N2['id']), as.list(edgeR_DESeq_NOISeq_N4['id']), as.list(edgeR_DESeq_NOISeq_N6['id'])),
  category.names = c("200mM" , "400mM", "600mM"),
  filename = 'alltools_NOD.png',
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
  cat.pos = c(-27, 27, 130),
  cat.dist = c(0.055, 0.055, 0.070),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = c(as.list(down_N_alltools_200['id']), as.list(down_N_alltools_400['id']), as.list(down_N_alltools_600['id'])),
  category.names = c("200mM" , "400mM", "600mM"),
  filename = 'down_alltools_NOD.png',
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
  cat.pos = c(-27, 27, 130),
  cat.dist = c(0.055, 0.055, 0.070),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = c(as.list(up_N_alltools_200['id']), as.list(up_N_alltools_400['id']), as.list(up_N_alltools_600['id'])),
  category.names = c("200mM" , "400mM", "600mM"),
  filename = 'up_alltools_NOD.png',
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
  cat.pos = c(-27, 27, 130),
  cat.dist = c(0.055, 0.055, 0.070),
  cat.fontfamily = "sans",
  rotation = 1
)

library(VennDiagram)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")



#################
### NaCl = 200mM
#################

venn.diagram(
  x = c(as.list(edgeR_K2_P01['id']), as.list(resSigK2['id']), as.list(K2_noiseq_1['id'])),
  category.names = c("edgeR" , "DESeq", "NOISeq"),
  filename = 'down_edgeR_DESeq_NOISeq_K2.png',
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
  cat.pos = c(-27, 27, 105),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


#################
### NaCl = 400mM
#################

venn.diagram(
  x = c(as.list(edgeR_K4_P01['id']), as.list(resSigK4['id']), as.list(K4_noiseq_1['id'])),
  category.names = c("edgeR" , "DESeq", "NOISeq"),
  filename = 'edgeR_DESeq_NOISeq_K4.png',
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
  cat.pos = c(-27, 27, 130),
  cat.dist = c(0.055, 0.055, 0.070),
  cat.fontfamily = "sans",
  rotation = 1
)


#################
### NaCl = 600mM
#################

venn.diagram(
  x = c(as.list(edgeR_K6_P01['id']), as.list(resSigK6['id']), as.list(K6_noiseq_1['id'])),
  category.names = c("edgeR" , "DESeq", "NOISeq"),
  filename = 'edgeR_DESeq_NOISeq_K6.png',
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


##############################
# NOD vs. KNO3 in all tools
#############################

K6_noiseq.dg.up_1 <- format_df(K6_noiseq.dg.up)
K6_noiseq.dg.down_1 <- format_df(K6_noiseq.dg.down)
N6_noiseq.dg.up_1 <- format_df(N6_noiseq.dg.up)
N6_noiseq.dg.down_1 <- format_df(N6_noiseq.dg.down)

# DE upregulated genes with 600mM for all tools and KNO3 genotype
DESeq_NOISeq_join_K6.up <- merge(upK6, K6_noiseq.dg.up_1, by = "id") # no matches
edgeR_DESeq_NOISeq_join_K6.up <- merge(edgeR_K6_up, DESeq_NOISeq_join_K6.up, by = "id")
#a_edgeR_DESeq_NOISeq_join_K6.up <- merge(edgeR_K6_up, DESeq_NOISeq_join_K6.up, by = "id")

# DE downregulated genes with 600mM for all tools and KNO3 genotype
DESeq_NOISeq_join_K6.down <- merge(downK6, K6_noiseq.dg.down_1, by = "id") 
edgeR_DESeq_NOISeq_join_K6.down <- merge(edgeR_K6_down, DESeq_NOISeq_join_K6.down, by = "id")

# DE upregulated genes with 600mM for all tools and NOD genotype
DESeq_NOISeq_join_N6.up <- merge(upN6, N6_noiseq.dg.up_1, by = "id") 
edgeR_DESeq_NOISeq_join_N6.up <- merge(edgeR_N6_up, DESeq_NOISeq_join_N6.up, by = "id")

# DE downregulated genes with 600mM for all tools and both KNO3 and NOD genotypes
DESeq_NOISeq_join_N6.down <- merge(downN6, N6_noiseq.dg.down_1, by = "id")
edgeR_DESeq_NOISeq_join_N6.down <- merge(edgeR_N6_down, DESeq_NOISeq_join_N6.down, by = "id")


venn.diagram(
  x = c(as.list(edgeR_DESeq_NOISeq_join_K6.down['id']), as.list(edgeR_DESeq_NOISeq_join_N6.down['id'])),
  category.names = c("KNO3" , "NOD"),
  filename = '1alltools_600down_KNO3_vs_NOD.png',
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
  fill = c('royalblue3', 'seagreen'),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans"
)

venn.diagram(
  x = c(as.list(edgeR_DESeq_NOISeq_join_K6.up['id']), as.list(edgeR_DESeq_NOISeq_join_N6.up['id'])),
  category.names = c("KNO3" , "NOD"),
  filename = '1alltools_600up_KNO3_vs_NOD.png',
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
  fill = c('royalblue3', 'seagreen'),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans"
)
