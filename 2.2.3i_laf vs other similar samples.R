#### Part 3:3i/LAF vs other iPS

# c2:samples indicating cluster2 of corelation plot
### Expression levels of gene sets in each group of condition -------

## get samples expression matrix

c2.tpm <- all.tpm[,c(1,2,8:16,21,22,37:47,54:59,62:69)]
c2.data <- all.data[,c(1,2,8:16,21,22,37:47,54:59,62:69)]
c2.info <- geneinfo.new[match(colnames(c2.tpm),geneinfo.new$Sample), c(1,7:11)]
rownames(c2.info) <- c2.info$Sample

## pluripotency genes set
# tpm
c2.plu.tpm <- c2.tpm[key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)],]
colSums(c2.plu.tpm)

ha_top.col <- list(Cell_type = new.col[c(7:14)], 
                   Author = c(pal_simpsons("springfield")(6)[3:6],pal_simpsons("springfield")(2)), 
                   System = new.col[20:24])
names(ha_top.col$Cell_type) <- factor(c('EDSC','EPSC','ESC','iPS',"pgEpiSC", "NGR.iPS", "pEF.iPS", "pEF.pMX.iPS"))
names(ha_top.col$Author) <- factor(c('Kinoshita','Gao','Choi','Yoshimatsu','Zhi','Zhu'))
names(ha_top.col$System) <- factor(c('AFX','ALCIW','AFCI','AFTI','3i_LAF'))

pheatmap::pheatmap(log2(c2.plu.tpm+1),
                   border_color = NA,
                   clustering_method = 'complete',
                   color = viridis_pal(alpha = 1, begin = .2)(256),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'plu.genes.heatmap',
                   breaks = unique(c(seq(0,12, length=256))),
                   legend_breaks = c(0,4,8,12),
                   fontsize_row = 8,
                   fontsize_col = 5)

# rank
c2.tpm.rank <- data.frame(apply(c2.tpm, 2, function(x){rank(nrow(c2.tpm)-rank(x,ties.method = 'first')+1)}))
c2.tpm.rank.bin <- data.frame(bin = 1:nrow(c2.tpm.rank))
for (i in colnames(c2.tpm.rank)) {
  bin <- data.frame(i = ifelse(c2.tpm.rank[,i]>0 & c2.tpm.rank[,i] <= 10, '1', 
                               ifelse(c2.tpm.rank[,i]>10 & c2.tpm.rank[,i] <= 100,'2',
                                      ifelse(c2.tpm.rank[,i]>50 & c2.tpm.rank[,i] <= 100,'3',
                                             ifelse(c2.tpm.rank[,i]>100 & c2.tpm.rank[,i] <= 500, '4',
                                                    ifelse(c2.tpm.rank[,i]>500 & c2.tpm.rank[,i] <= 1000, '5',
                                                           ifelse(c2.tpm.rank[,i]>1000 & c2.tpm.rank[,i] <= 2000, '6',
                                                                  ifelse(c2.tpm.rank[,i]>2000 & c2.tpm.rank[,i] <= 4000,'7',
                                                                         ifelse(c2.tpm.rank[,i]>4000 & c2.tpm.rank[,i] <= 6000,'8',
                                                                                ifelse(c2.tpm.rank[,i]>6000 & c2.tpm.rank[,i] <= 8000,'9',
                                                                                       ifelse(c2.tpm.rank[,i]>8000 & c2.tpm.rank[,i] <= 10000,'10',
                                                                                              ifelse(c2.tpm.rank[,i]>10000 & c2.tpm.rank[,i] <= 12000,'11',
                                                                                                     ifelse(c2.tpm.rank[,i]>12000 & c2.tpm.rank[,i] <= 14000,'12',
                                                                                                            ifelse(c2.tpm.rank[,i]>14000  & c2.tpm.rank[,i] <= 16000,'13',
                                                                                                                   ifelse(c2.tpm.rank[,i]>16000 & c2.tpm.rank[,i] <= 18000,'14',
                                                                                                                          ifelse(c2.tpm.rank[,i]>18000 & c2.tpm.rank[,i] <= 20000,'15',
                                                                                                                                 ifelse(c2.tpm.rank[,i]>20000 & c2.tpm.rank[,i] <= 22000,'16','17')))))))))))))))) )
  colnames(bin) <- i
  c2.tpm.rank.bin <- cbind(c2.tpm.rank.bin,bin)
}
c2.tpm.rank.bin <- c2.tpm.rank.bin[,-1] 
rownames(c2.tpm.rank.bin) <- rownames(c2.tpm.rank)
c2.tpm.rank.bin <- data.frame(lapply(c2.tpm.rank.bin,as.numeric))
c2.plu.tpm.rank.bin <- c2.tpm.rank.bin[rownames(c2.plu.tpm),]

pheatmap::pheatmap(c2.plu.tpm.rank.bin,
                   border_color = 'NA',
                   clustering_method = 'complete',
                   color = rev(viridis_pal(alpha = 1, begin = .1)(17)),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'plu.genes.rank.bin.heatmap',
                   breaks = unique(c(seq(1,16, length=17))),
                   legend_breaks = c(1,8,13,16),
                   legend_labels = c('bin1','bin8','bin13','bin16'),
                   fontsize_row = 8,
                   display_numbers = T,
                   number_format = "%.0f",
                   fontsize_number = 5,
                   number_color = 'gray90',
                   fontsize_col = 5)

## Ectoderm genes set
c2.ecto.tpm <- c2.tpm[key.genes$Ectoderm$Ectoderm,] %>% na.omit()

pheatmap::pheatmap(log2(c2.ecto.tpm+1),
                   border_color = NA,
                   clustering_method = 'complete',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'plu.genes.heatmap',
                   breaks = unique(c(seq(0,8, length=256))),
                   legend_breaks = c(0,4,8),
                   fontsize_row = 8,
                   fontsize_col = 5)

c2.ecto.tpm.rank.bin <- c2.tpm.rank.bin[rownames(c2.ecto.tpm),]

pheatmap::pheatmap(c2.ecto.tpm.rank.bin,
                   border_color = 'NA',
                   clustering_method = 'complete',
                   color = rev(viridis_pal(alpha = 1, begin = .1)(17)),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'plu.genes.rank.bin.heatmap',
                   breaks = unique(c(seq(1,16, length=17))),
                   legend_breaks = c(1,4,8,13,16),
                   legend_labels = c('bin1','bin4','bin8','bin13','bin16'),
                   fontsize_row = 8,
                   display_numbers = T,
                   number_format = "%.0f",
                   fontsize_number = 6,
                   number_color = 'gray90',
                   fontsize_col = 5)

## Mesoderm genes set
c2.meso.tpm <- c2.tpm[key.genes$Mesoderm$Mesoderm,] %>% na.omit()

pheatmap::pheatmap(log2(c2.meso.tpm+1),
                   border_color = NA,
                   clustering_method = 'complete',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'plu.genes.heatmap',
                   breaks = unique(c(seq(0,7, length=256))),
                   legend_breaks = c(0,3,7),
                   fontsize_row = 8,
                   fontsize_col = 5)

c2.meso.tpm.rank.bin <- c2.tpm.rank.bin[rownames(c2.meso.tpm),]

pheatmap::pheatmap(c2.meso.tpm.rank.bin,
                   border_color = 'NA',
                   clustering_method = 'complete',
                   color = rev(viridis_pal(alpha = 1, begin = .1)(17)),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'plu.genes.rank.bin.heatmap',
                   # breaks = unique(c(seq(1,16, length=17))),
                   # legend_breaks = c(1,7,16),
                   # legend_labels = c('bin1','bin7','bin16'),
                   fontsize_row = 8,
                   display_numbers = T,
                   number_format = "%.0f",
                   fontsize_number = 6,
                   number_color = 'gray90',
                   fontsize_col = 5)

## Endoderm genes set
c2.endo.tpm <- c2.tpm[key.genes$Endoderm$Endoderm,] %>% na.omit()

pheatmap::pheatmap(log2(c2.endo.tpm+1),
                   border_color = NA,
                   clustering_method = 'complete',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'plu.genes.heatmap',
                   breaks = unique(c(seq(0,10, length=256))),
                   legend_breaks = c(0,5,10),
                   fontsize_row = 8,
                   fontsize_col = 5)

c2.endo.tpm.rank.bin <- c2.tpm.rank.bin[rownames(c2.endo.tpm),]

pheatmap::pheatmap(c2.endo.tpm.rank.bin,
                   border_color = 'NA',
                   clustering_method = 'complete',
                   color = rev(viridis_pal(alpha = 1, begin = .1)(17)),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'endo.genes.rank.bin.heatmap',
                   breaks = unique(c(seq(1,16, length=17))),
                   legend_breaks = c(1,4,16),
                   legend_labels = c('bin1','bin4','bin16'),
                   fontsize_row = 8,
                   display_numbers = T,
                   number_format = "%.0f",
                   fontsize_number = 6,
                   number_color = 'gray90',
                   fontsize_col = 5)

## PE genes set
c2.PE.tpm <- c2.tpm[key.genes$PE$PE,] %>% na.omit()

pheatmap::pheatmap(log2(c2.PE.tpm+1),
                   border_color = NA,
                   clustering_method = 'complete',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'PE.genes.heatmap',
                   breaks = unique(c(seq(0,12, length=256))),
                   legend_breaks = c(0,6,12),
                   fontsize_row = 8,
                   fontsize_col = 5)

c2.PE.tpm.rank.bin <- c2.tpm.rank.bin[rownames(c2.PE.tpm),]

pheatmap::pheatmap(c2.PE.tpm.rank.bin,
                   border_color = 'NA',
                   clustering_method = 'complete',
                   color = rev(viridis_pal(alpha = 1, begin = .1)(17)),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'PE.genes.rank.bin.heatmap',
                   breaks = unique(c(seq(1,16, length=17))),
                   legend_breaks = c(1,8,16),
                   legend_labels = c('bin1','bin8','bin16'),
                   fontsize_row = 8,
                   display_numbers = T,
                   number_format = "%.0f",
                   fontsize_number = 6,
                   number_color = 'gray90',
                   fontsize_col = 5)
## TE
c2.TE.tpm <- c2.tpm[unique(key.genes$TE$TE),] %>% na.omit()

pheatmap::pheatmap(log2(c2.TE.tpm+1),
                   border_color = NA,
                   clustering_method = 'ward.D',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'TE.genes.heatmap',
                   breaks = unique(c(seq(0,12, length=256))),
                   legend_breaks = c(0,5,12),
                   fontsize_row = 8,
                   fontsize_col = 5)

c2.TE.tpm.rank.bin <- c2.tpm.rank.bin[rownames(c2.TE.tpm),]

pheatmap::pheatmap(c2.TE.tpm.rank.bin,
                   border_color = 'NA',
                   clustering_method = 'ward.D',
                   color = rev(viridis_pal(alpha = 1, begin = .1)(17)),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'TE.genes.rank.bin.heatmap',
                   breaks = unique(c(seq(1,16, length=17))),
                   legend_breaks = c(1,8,16),
                   legend_labels = c('bin1','bin8','bin16'),
                   fontsize_row = 8,
                   display_numbers = T,
                   number_format = "%.0f",
                   fontsize_number = 6,
                   number_color = 'gray90',
                   fontsize_col = 5)

# wnt genes set
library(KEGGREST)
wnt.genes <- keggGet("ssc04310")
wnt.genes <- sapply(seq(2,334,2), function(x){
  gns = unlist(strsplit(wnt.genes[[1]][["GENE"]][x],";"))[1]
})

c2.wnt.tpm <- all.tpm.clean.pick[unique(wnt.genes),colnames(c2.tpm)] %>% na.omit()
c2.wnt.tpm <- c2.wnt.tpm[-107,]

pheatmap::pheatmap(log2(c2.wnt.tpm+1),
                   border_color = NA,
                   clustering_method = 'ward.D',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'TE.genes.heatmap',
                   breaks = unique(c(seq(0,11, length=256))),
                   legend_breaks = c(0,5,11),
                   fontsize_row = 5,
                   fontsize_col = 5)


c2.wnt.tpm.rank.bin <- c2.tpm.rank.bin[rownames(c2.wnt.tpm),]

pheatmap::pheatmap(c2.wnt.tpm.rank.bin,
                   border_color = 'NA',
                   clustering_method = 'ward.D',
                   color = rev(viridis_pal(alpha = 1, begin = .1)(17)),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'wnt.genes.rank.bin.heatmap',
                   breaks = unique(c(seq(1,16, length=17))),
                   legend_breaks = c(1,8,16),
                   legend_labels = c('bin1','bin8','bin16'),
                   fontsize_row = 5,
                   display_numbers = T,
                   number_format = "%.0f",
                   fontsize_number = 5,
                   number_color = 'gray90',
                   fontsize_col = 5)

### WNT GSEA----------
library(clusterProfiler)
library(enrichplot)
library(GSEABase)
library(org.Ss.eg.db)
library(org.Hs.eg.db)

## only pig genes annotation & biomart ---- wnt pathway just located in kegg database

library(biomaRt)
ensmebl <- useDataset("sscrofa_gene_ensembl", useMart("ensembl"))
pig.mart <- useEnsembl('ensembl',dataset = "sscrofa_gene_ensembl")

gene_fc <- DEG12[,c(1,3)] %>% 
  left_join(y = pig.gene.info,by=c('Symbol'='genename'))
genes_mart = getBM(attributes = c('ensembl_gene_id','entrezgene_id'), values = gene_fc$geneid, mart = pig.mart, filters = 'ensembl_gene_id') %>% na.omit()
genes_mart <- genes_mart[!duplicated(genes_mart$entrezgene_id),]

gene_fc <- gene_fc %>% 
  left_join(y = genes_mart,by = c('geneid'='ensembl_gene_id')) %>% 
  na.omit() %>% 
  dplyr::arrange(desc(log2FoldChange))

gene_for_gsea <- gene_fc$log2FoldChange 
names(gene_for_gsea) <- gene_fc$entrezgene_id

gsea_k <- gseKEGG(
  gene_for_gsea,
  organism = "ssc",
  keyType = "kegg",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea"
)
gsea_k_order <- gsea_k[order(gsea_k$enrichmentScore, decreasing = T),]
gseaplot2(gsea_k, "ssc04310", color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table = T)


## orth 2 human genes annotation & biomart ---- wnt pathway just located in kegg database
human.mart <- useEnsembl('ensembl',dataset = "hsapiens_gene_ensembl")

gene_fc <- DEG11[,c(1,3)] %>% 
  left_join(y = orth[,2:3],by=c('Symbol'='genename'))
genes_mart = getBM(attributes = c('ensembl_gene_id','entrezgene_id'), values = gene_fc$human.geneid, mart = human.mart, filters = 'ensembl_gene_id') %>% na.omit()
genes_mart <- genes_mart[!duplicated(genes_mart$entrezgene_id),]

gene_fc <- gene_fc %>% 
  left_join(y = genes_mart,by = c('human.geneid'='ensembl_gene_id')) %>% 
  na.omit() %>% 
  dplyr::arrange(desc(log2FoldChange))

gene_for_gsea <- gene_fc$log2FoldChange 
names(gene_for_gsea) <- gene_fc$entrezgene_id

gsea_k <- gseKEGG(
  gene_for_gsea,
  organism = "hsa",
  keyType = "kegg",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea"
)
gsea_k_order <- gsea_k[order(gsea_k$enrichmentScore, decreasing = T),]
gseaplot2(gsea_k, "hsa04310", color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table = T)

# JAK-STAT GSEA-----
gene_fc <- DEG1[,c(1,3)] %>% 
  left_join(y = pig.gene.info,by=c('Symbol'='genename'))
genes_mart = getBM(attributes = c('ensembl_gene_id','entrezgene_id'), values = gene_fc$geneid, mart = pig.mart, filters = 'ensembl_gene_id') %>% na.omit()
genes_mart <- genes_mart[!duplicated(genes_mart$entrezgene_id),]

gene_fc <- gene_fc %>% 
  left_join(y = genes_mart,by = c('geneid'='ensembl_gene_id')) %>% 
  na.omit() %>% 
  dplyr::arrange(desc(log2FoldChange))

gene_for_gsea <- gene_fc$log2FoldChange 
names(gene_for_gsea) <- gene_fc$entrezgene_id

gsea_k <- gseKEGG(
  gene_for_gsea,
  organism = "ssc",
  keyType = "kegg",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea"
)
gsea_k_order <- gsea_k[order(gsea_k$enrichmentScore, decreasing = T),]
gseaplot2(gsea_k, "ssc04630", color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table = T)

# EMT GSEA-----
gene_fc <- DEG12[,c(1,3)] %>% 
  left_join(y = pig.gene.info,by=c('Symbol'='genename'))
genes_mart = getBM(attributes = c('ensembl_gene_id','entrezgene_id'), values = gene_fc$geneid, mart = pig.mart, filters = 'ensembl_gene_id') %>% na.omit()
genes_mart <- genes_mart[!duplicated(genes_mart$entrezgene_id),]

gene_fc <- gene_fc %>% 
  left_join(y = genes_mart,by = c('geneid'='ensembl_gene_id')) %>% 
  na.omit() %>% 
  dplyr::arrange(desc(log2FoldChange))

gene_for_gsea <- gene_fc$log2FoldChange 
names(gene_for_gsea) <- gene_fc$entrezgene_id

gsea_k <- gseKEGG(
  gene_for_gsea,
  organism = "ssc",
  keyType = "kegg",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea"
)
gsea_k_order <- gsea_k[order(gsea_k$enrichmentScore, decreasing = T),]
gseaplot2(gsea_k, "ssc01521", color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table = T)

# MAPK GSEA-----
gene_fc <- DEG1[,c(1,3)] %>% 
  left_join(y = pig.gene.info,by=c('Symbol'='genename'))
genes_mart = getBM(attributes = c('ensembl_gene_id','entrezgene_id'), values = gene_fc$geneid, mart = pig.mart, filters = 'ensembl_gene_id') %>% na.omit()
genes_mart <- genes_mart[!duplicated(genes_mart$entrezgene_id),]

gene_fc <- gene_fc %>% 
  left_join(y = genes_mart,by = c('geneid'='ensembl_gene_id')) %>% 
  na.omit() %>% 
  dplyr::arrange(desc(log2FoldChange))

gene_for_gsea <- gene_fc$log2FoldChange 
names(gene_for_gsea) <- gene_fc$entrezgene_id

gsea_k <- gseKEGG(
  gene_for_gsea,
  organism = "ssc",
  keyType = "kegg",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea"
)
gsea_k_order <- gsea_k[order(gsea_k$enrichmentScore, decreasing = T),]
gseaplot2(gsea_k, "ssc04010", color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table = T)

### combine scRNA&RNA seq data--------

## embryo dataset
library(Seurat)
library(DeconRNASeq)
library(edgeR)
pig.em.seu <- readRDS('../pig.zml.em.seu.rds')
pig.epi.seu <- subset(pig.em.seu,cells=c(rownames(pig.em.seu@meta.data[pig.em.seu$Cell_type=='ICM' | pig.em.seu$Cell_type=='Epiblast' | pig.em.seu$Cell_type == 'Ectoderm',])))

# hvg 1000 or 2000? 

pig.epi.seu <- pig.epi.seu %>%
  NormalizeData( normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures( selection.method = "vst", nfeatures = 1500) %>% 
  ScaleData( features = rownames(pig.epi.seu), vars.to.regress = 'mitoRatio') %>% 
  RunUMAP( dims = 1:10)

Idents(pig.epi.seu) <- pig.epi.seu$Embryonic_Day
ed.markers <- FindAllMarkers(pig.epi.seu,only.pos = T,logfc.threshold = 0.25) %>% 
  filter(p_val_adj < 0.05) %>% 
  group_by(cluster) %>% 
  arrange(avg_log2FC)

ed.markers.top100 <- ed.markers %>% 
  group_by(cluster) %>% 
  top_n(n=100,wt=avg_log2FC)

DimPlot(pig.epi.seu,cols = new.col,group.by = 'Embryonic_Day',label = T)
DimPlot(pig.epi.seu,reduction = 'pca',cols = new.col,group.by = 'Embryonic_Day',label = T)

## SF score
library(KEGGREST)
sf_gene <- keggGet("ssc03040")
sf_gene <- sapply(seq(2,250,2), function(x){
  gns = unlist(strsplit(sf_gene[[1]][["GENE"]][x],";"))[1]
})

DefaultAssay(pig.epi.seu) <- 'RNA'
pig.epi.seu <- pig.epi.seu %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)

sf.gene.ave <- data.frame(AverageExpression(pig.epi.seu,features = sf_gene, group.by = 'Embryonic_Day',assays = 'RNA'))

pheatmap::pheatmap(log2(sf.gene.ave+1),
                   border_color = 'NA',
                   clustering_method = 'ward.D',
                   color = (viridis_pal(alpha = 1, begin = .1)(256)),
                   angle_col = '315',
                   # annotation_col = c2.info[,c(2,3,5)],
                   # show_rownames = T,
                   # annotation_colors = ha_top.col,
                   # main = 'wnt.genes.rank.bin.heatmap',
                   breaks = unique(c(seq(1,3, length=256))),
                   # legend_breaks = c(1,8,16),
                   # legend_labels = c('bin1','bin8','bin16'),
                   fontsize_row = 5,
                   # display_numbers = T,
                   number_format = "%.0f",
                   cluster_cols = F,
                   # fontsize_number = 5,
                   # number_color = 'gray90',
                   fontsize_col = 15)


pig.epi.seu$sf_score <- PercentageFeatureSet(pig.epi.seu, features = rownames(sf.gene.ave))
view(pig.epi.seu)

FeaturePlot(pig.epi.seu, features = 'sf_score',cols = c('grey90','red3'))
VlnPlot(pig.epi.seu, features = 'sf_score')+scale_color_simpsons()

VlnPlot(pig.epi.seu, features = c('SRSF1','SRSF2','SRSF3','SRSF5','SRSF10','TPA2B','U2AF1','QKI','HNRNPA1','HNRNPA2','HNRNPH','HNRNPK','HNRNPL','HNRNPM','PTBP1','RBM5','RBM10','RBM39','ESPR1','ESPR2','RBFOX2'), stack = T,flip = T) 


# HVGS
hvgs <- VariableFeatures(pig.epi.seu)

pig.epi.counts <- data.frame(GetAssayData(pig.epi.seu,slot = 'count',assay = 'RNA'))
pig.epi.counts.cpm <- cpm(pig.epi.counts)

pig.epi.cpm.hvgs <- data.frame(pig.epi.counts.cpm[hvgs,]) 
pig.epi.cpm.hvgs <- round(pig.epi.cpm.hvgs,2)

## part3 data cpm
c2.cpm <- data.frame(cpm(c2.data))
c2.cpm.hvgs <- data.frame(c2.cpm[hvgs,])
c2.cpm.hvgs <- round(c2.cpm.hvgs,2)

## method1:pca & predict cos2
d_p <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 40),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               mean.point=F,
               axes = c(1, 2),
               # repel = T,
               label = "none",
               geom.ind = c("point",'text'),
               fill.ind = tpm.t$group_list,
               palette = new.col,
               legend.title = "Embryonic day",
               pointsize = 2,
               pointshape = 21,
               col.ind = "black",
               title = 'Combine sc&bulk data',
               # addEllipses = T
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}

# pca plot
group_list <- factor(pig.epi.seu$Embryonic_Day)
d_p(pig.epi.cpm.hvgs,group_list)

# get PCA object
cpm.t = as.data.frame(t(pig.epi.cpm.hvgs))
cpm.t = cbind(cpm.t,group_list)
cpm.pca <- PCA(cpm.t[,-ncol(cpm.t)],graph = FALSE)
fviz_screeplot(cpm.pca, addlabels = TRUE, ylim = c(0, 10),title='Dim choose')

# predict

# all other cells 
ind.sup.coord <- predict(cpm.pca, newdata = t(c2.cpm.hvgs))
ind.sup.coord$coord
ind.sup.coord$cos2

# only our cells
ind.sup.coord <- predict(cpm.pca, newdata = t(c2.cpm.hvgs[,25:ncol(c2.cpm.hvgs)]))
ind.sup.coord$coord
ind.sup.coord$cos2

# add predict result to PCA object
p <- fviz_pca_ind(cpm.pca, 
                  mean.point=F,
                  axes = c(1, 2),
                  label = "none",
                  geom.ind = c("point",'text'),
                  fill.ind = cpm.t$group_list,
                  palette = RColorBrewer::brewer.pal(9,"Paired"),
                  legend.title = "Embryonic day",
                  pointsize = 2,
                  pointshape = 21,
                  col.ind = "gray90",
                  title = 'Combine sc&bulk data',
                  alpha = 0.5
)
fviz_add(p, ind.sup.coord$coord, 
         color = c(rep("#2A82C6",2), rep("#7C4374",9),rep('#B09C85FF',2), rep('#339E55',2), rep('#000376',9), rep('#8C6D36',2), rep('#CB70C0',4), rep('#EBB854',4), rep('#FC8D37',4)),
         repel = T,
         addlabel = F,
         # labelsize =4,
         # geom = c("point", "text"),
         pointsize = 3,
         alpha=0.5)
fviz_add(p, ind.sup.coord$coord, 
         color = c( rep('#8C6D36',2), rep('#CB70C0',4), rep('#EBB854',4), rep('#FC8D37',4)),
         repel = T,
         addlabel = F,
         # labelsize =4,
         # geom = c("point", "text"),
         pointsize = 3,
         alpha=0.5)

show_col(c(rep("#2A82C6",2), rep("#7C4374",9),rep('#B09C85FF',2), rep('#339E55',2), rep('#000376',9), rep('#8C6D36',2), rep('#CB70C0',4), rep('#EBB854',4), rep('#FC8D37',4)),ncol = 1,labels = F)

## method2:DeconRNASeq
pig.epi.counts.mean <- data.frame(AverageExpression(pig.epi.seu,slot = 'count',assays = 'RNA',group.by = 'Embryonic_Day'))
colnames(pig.epi.counts.mean) <- unique(pig.epi.seu$Embryonic_Day)
pig.epi.counts.mean.cpm <- data.frame(cpm(pig.epi.counts.mean))
pig.epi.counts.mean.hvg.cpm <- pig.epi.counts.mean.cpm[hvgs,]
pig.epi.ed.marker.cpm <- pig.epi.counts.mean[ed.markers$gene,]
ed.markers.top100.cpm <- pig.epi.counts.mean[ed.markers.top100$gene,]
# signatures=ed.all.cpm
# DeconRNASeq(c2.cpm, pig.epi.counts.mean.cpm, cell.ratio,checksig=FALSE,
#             known.prop = F, use.scale = TRUE, fig = TRUE)
# 
# # signatures=ed.hvg.cpm
# DeconRNASeq(c2.cpm.hvgs, pig.epi.counts.mean.hvg.cpm[,-c(1:4)],checksig=FALSE,
#             known.prop = F, use.scale = TRUE, fig = TRUE)
# 
# # signatures=ed.markers.cpm
# DeconRNASeq(c2.cpm, pig.epi.ed.marker.cpm[,-4],checksig=FALSE,
#             known.prop = F, use.scale = TRUE, fig = TRUE)
# 
# # signatures=ed.top100.markers.cpm
# DeconRNASeq(c2.cpm, ed.markers.top100.cpm[,-4],checksig=FALSE,
#             known.prop = F, use.scale = TRUE, fig = TRUE)

# signatures=ed.markers.cpm---final choose
cell.fraction <- DeconRNASeq(c2.cpm, pig.epi.ed.marker.cpm[,-c(1:4)],checksig=FALSE,
            known.prop = F, use.scale = TRUE, fig = TRUE)
export(cell.fraction,file = '../3.part3/cell.fraction.xlsx')
# signatures=ed.top100.markers.cpm
# DeconRNASeq(c2.cpm, ed.markers.top100.cpm[,-c(1:4)],checksig=FALSE,
#             known.prop = F, use.scale = TRUE, fig = TRUE)

## ratio
# # signatures=ed.all.cpm
# DeconRNASeq(c2.cpm, pig.epi.counts.mean.cpm, cell.ratio,checksig=FALSE,
#             known.prop = T, use.scale = TRUE, fig = TRUE)
# 
# # signatures=ed.hvg.cpm
# DeconRNASeq(c2.cpm, pig.epi.counts.mean.hvg.cpm,cell.ratio,checksig=FALSE,
#             known.prop = T, use.scale = TRUE, fig = TRUE)
# 
# # signatures=ed.markers.cpm
# DeconRNASeq(c2.cpm, pig.epi.ed.marker.cpm,cell.ratio,checksig=FALSE,
#             known.prop = T, use.scale = TRUE, fig = TRUE)
# 
# # signatures=ed.markers.cpm
# DeconRNASeq(c2.cpm, ed.markers.top100.cpm,checksig=FALSE,
#             known.prop = F, use.scale = TRUE, fig = TRUE)
# plot
cell.fraction.plot <- data.frame(cell.fraction$out.all)
cell.fraction.plot$Samples <- colnames(c2.cpm.hvgs[,25:ncol(c2.cpm.hvgs)])
cell.fraction.plot <- melt(cell.fraction.plot,variable.name = 'Embryonic_Day',id.vars = 'Samples',value.name = 'Fraction')
cell.fraction.plot$Embryonic_Day <- factor(cell.fraction.plot$Embryonic_Day, levels = c('E14','E13','E12','E11','E10'))

ggplot(cell.fraction.plot, aes(fill=Embryonic_Day, y=Fraction, x=Samples)) + 
  geom_bar( stat="identity") + 
  coord_flip() + 
  guides(fill = guide_legend(reverse=TRUE)) +
  theme_bw() +
  scale_fill_manual( values = c("#CAB2D6","#FF7F00", "#FDBF6F",  "#E31A1C", "#FB9A99"))

# WGCNA ---------

### STEP1.MATRIX
# c2.meta <- data.frame(Cell_type = c2.info$System[c(1:24,31:34)])
c2.meta <- data.frame(Cell_type = c2.info$System)
# rownames(c2.meta) <- c2.info$Sample[c(1:24,31:34)]
rownames(c2.meta) <- c2.info$Sample

datExpr <- all.tpm.clean.pick %>% dplyr::select(rownames(c2.meta)) %>% na.omit()
datExpr <- t(datExpr[order(apply(datExpr,1,mad), decreasing = T)[1:10000],])

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

### STEP2:soft-thresholding power--
powers = c(c(1:10), seq(from = 12, to=50, by=2))
sft = pickSoftThreshold(datExpr, RsquaredCut = 0.9,powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col="red");
abline(h=0.9,col="green")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1,col="red")

best_beta=sft$powerEstimate

### STEP3:build co-expresssion network-
net = blockwiseModules(datExpr, power = best_beta,
                       maxBlockSize = 6000,TOMType = "unsigned", 
                       minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "All.tpm-TOM",
                       verbose = 3)

### STEP4:Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
moduleColors=mergedColors
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

datExpr_tree <- hclust(dist(datExpr), method = "ward.D2")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)

### STEP5:module & Traits-

design=model.matrix(~0+ c2.meta$Cell_type)
colnames(design)=c("3i_LAF", 'AFCI', 'AFTI', 'AFX', 'ALCIW')
moduleColors <- labels2colors(net$colors)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors =  blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
par(mar = c(6, 8.5, 3, 3))
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",plotDendrograms = F,
                      marDendro = c(4,4,2,4))

### STEP6:3i_LAF:cayn--
# part1 3ilaf
i3_LAF = as.data.frame(design[,1])
names(i3_LAF) = "3i_LAF"
module = "greenyellow"

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneModuleMembership[1:4,1:4]

geneTraitSignificance = as.data.frame(cor(datExpr, i3_LAF, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(i3_LAF), sep="")
names(GSPvalue) = paste("p.GS.", names(i3_LAF), sep="")

column = match(module, modNames)
moduleGenes = moduleColors==module

par(mfrow = c(2,2))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for i3_LAF",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)

# part2 AFCI
AFCI = as.data.frame(design[,2])
names(AFCI) = "AFCI"
module = "blue"

geneTraitSignificance = as.data.frame(cor(datExpr, AFCI, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(AFCI), sep="")
names(GSPvalue) = paste("p.GS.", names(AFCI), sep="")

column = match(module, modNames)
moduleGenes = moduleColors==module

par(mfrow = c(2,2))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for AFCI",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)

# part3 AFTI
AFTI = as.data.frame(design[,3])
names(AFTI) = "AFTI"
module = "pink"

geneTraitSignificance = as.data.frame(cor(datExpr, AFTI, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(AFTI), sep="")
names(GSPvalue) = paste("p.GS.", names(AFTI), sep="")

column = match(module, modNames)
moduleGenes = moduleColors==module

par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for AFTI",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)

# part4 AFX
AFX = as.data.frame(design[,4])
names(AFX) = "AFX"
module = "brown"

geneTraitSignificance = as.data.frame(cor(datExpr, AFX, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(AFX), sep="")
names(GSPvalue) = paste("p.GS.", names(AFX), sep="")

column = match(module, modNames)
moduleGenes = moduleColors==module

par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for AFX",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)

# part5 ALCIW
ALCIW = as.data.frame(design[,5])
names(ALCIW) = "ALCIW"
module = "red"

geneTraitSignificance = as.data.frame(cor(datExpr, ALCIW, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(ALCIW), sep="")
names(GSPvalue) = paste("p.GS.", names(ALCIW), sep="")

column = match(module, modNames)
moduleGenes = moduleColors==module

par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for ALCIW",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)

### STEP7:network visualization---
# part1:genes network 
geneTree = net$dendrograms[[1]] 
dissTOM = 1-TOMsimilarityFromExpr(as.data.frame(datExpr), power = best_beta)
plotTOM = dissTOM^7
diag(plotTOM) = NA 
# TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

nSelect = 500
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "ward.D2")
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
rm(plotTOM,dissTOM,selectTOM,selectTree,selectColors,select,plotDiss)

# part2:samples_module network 

MET = orderMEs(cbind(MEs, i3_LAF))
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


### STEP8:提取指定模块的基因名--
table(mergedColors)
# PART1:3I_LAF
module = "greenyellow"
probes = colnames(datExpr)
inModule = (moduleColors==module)
greenyellow_3i_LAF = probes[inModule]
greenyellow_3i_LAF_orth <- data.frame(greenyellow_3i_LAF = greenyellow_3i_LAF) %>% left_join(y = orth[,2:4], by = c('greenyellow_3i_LAF' = 'genename'))

# PART2:AFCI
module = "blue"
probes = colnames(datExpr)
inModule = (moduleColors==module)
blue_AFCI = probes[inModule]
blue_AFCI_orth <- data.frame(blue_AFCI = blue_AFCI) %>% left_join(y = orth[,2:4], by = c('blue_AFCI' = 'genename'))

# PART3:AFTI
module = "pink"
probes = colnames(datExpr)
inModule = (moduleColors==module)
pink_AFTI = probes[inModule]
pink_AFTI_orth <- data.frame(pink_AFTI = pink_AFTI) %>% left_join(y = orth[,2:4], by = c('pink_AFTI' = 'genename'))

# PART4:AFX
module = "brown"
probes = colnames(datExpr)
inModule = (moduleColors==module)
brown_AFX = probes[inModule]
brown_AFX_orth <- data.frame(brown_AFX = brown_AFX) %>% left_join(y = orth[,2:4], by = c('brown_AFX' = 'genename'))

# PART5:ALCIW
module = "red"
probes = colnames(datExpr)
inModule = (moduleColors==module)
red_ALCIW = probes[inModule]
red_ALCIW_orth <- data.frame(red_ALCIW = red_ALCIW) %>% left_join(y = orth[,2:4], by = c('red_ALCIW' = 'genename'))

slience.spe.module <- list(greenyellow_3i_LAF_orth = greenyellow_3i_LAF_orth, 
                           blue_AFCI_orth = blue_AFCI_orth,
                           pink_AFTI_orth = pink_AFTI_orth,
                           brown_AFX_orth = brown_AFX_orth,
                           red_ALCIW_orth = red_ALCIW_orth)

export(slience.spe.module, file = '../5.part5/slience.spe.module.genes.xlsx')
# STEP9:模块导出-

TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate)
# Select module
module = "greenyellow"
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
)

vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)

# STEP10:hub genes

# (1) Intramodular connectivity
# moduleColors <- labels2colors(net$colors)
connet=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(connet, moduleColors)
head(Alldegrees1)

# (2) Relationship between gene significance and intramodular connectivity
which.module="greenyellow"

GS1 = as.numeric(cor(i3_LAF,datExpr, use="p"))
GeneSignificance=abs(GS1)
colorlevels=unique(moduleColors)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))){
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

# (3) Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
head(datKME)
##Finding genes with high gene significance and high intramodular connectivity in
# interesting modules
FilterGenes = abs(GS1) > .8 & abs(datKME$MM.greenyellow)>.8
table(FilterGenes)
test <- colnames(datExpr)[FilterGenes] 
test <- intersect(test,greenyellow_3i_LAF)
write.table(test,file = '../5.part5/3i_LAF.hub.genes.GS0.8_MM0.8.txt',quote = F,col.names = F,row.names = F)
test1 <- setdiff(greenyellow_3i_LAF,test)
test <- data.frame(genes = c(test,test1), type = c(rep('hub',26), rep('nohub',76)))
export(test,file = '../5.part5/3i_LAF.hub.genes.GS0.8_MM0.8.xlsx')

chooseTopHubInEachModule(datExpr,moduleColors)

# # (4) other methods
# geneTraitSignificance = as.data.frame(cor(datExpr, i3_LAF, use = "p"))
# geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
# 
# module = "greenyellow"
# column = match(module, modNames)
# moduleGenes = moduleColors==module
# greenyellow_module<-as.data.frame(dimnames(data.frame(datExpr))[[2]][moduleGenes])
# names(greenyellow_module)="genename"
# MM <- abs(geneModuleMembership[moduleGenes,column])
# GS <- abs(geneTraitSignificance[moduleGenes, 1])
# greenyellow_MMGS<-as.data.frame(cbind(MM,GS))
# rownames(greenyellow_MMGS)=greenyellow_module$genename
# hub_b<-abs(greenyellow_MMGS$MM)>0.8&abs(greenyellow_MMGS$GS)>0.8
# table(hub_b)
# rownames(subset(greenyellow_MMGS, abs(greenyellow_MMGS$MM)>0.8&abs(greenyellow_MMGS$GS)>0.8))

# STEP11: 其他分析
plotMEpairs(MEs,y=c2.meta$cellType)

# Diagnostics: heatmap plots of module expression
#par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
# for greenyellow module
sizeGrWindow(8,9)
which.module="greenyellow";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )

sizeGrWindow(8,7);
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="MPP")

# 优化
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))

pheatmap::pheatmap(t(scale(datExpr[,moduleColors==which.module ]) ),
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .8),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'GREENYELLOW.tpm.scale',
                   breaks = unique(c(seq(-3,3, length=256))),
                   legend_breaks = c(-3,0,3),
                   fontsize_row = 5,
                   fontsize_col = 5)

par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="MPP")

save.image()
