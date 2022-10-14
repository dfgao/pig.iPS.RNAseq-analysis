### Part1-2:iPS vs pgEpiSC
# PCA-------
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
               palette = new.col[c(11:14)],
               legend.title = "Groups",
               pointsize = 4,
               pointshape = 21,
               col.ind = "black",
               title = 'g&piPS PCA',
               # addEllipses = T
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}

group_list <- factor(geneinfo.all$ct3[c(54:59,62:69)],levels = c('Zhu_NGR.iPS_3i_LAF','Zhu_pEF.iPS_3i_LAF','Zhu_pEF.pMX.iPS_3i_LAF','Zhi_pgEpiSC_3i_LAF'))
d_p(my.tpm[,-c(7,8)],group_list)

# correlation------
my.cor.p <- cor(my.tpm[,-c(7,8)]) %>% round(digits = 2)
my.cor.s <- cor(my.tpm[,-c(7,8)],method = 'spearman') %>% round(digits = 2)

ha_left.col <- list(Cell_type = new.col[11:14], Author = pal_simpsons("springfield")(2))
names(ha_left.col$Cell_type) <- factor(unique(my.ann[1:14,]$Cell_type),levels = c("pgEpiSC", "NGR.iPS", "pEF.iPS", "pEF.pMX.iPS"))
names(ha_left.col$Author) <- factor(unique(my.ann[1:14,]$Author),levels = c('Zhi','Zhu'))

pheatmap::pheatmap(my.cor.s,
                   border_color = 'gray80',
                   clustering_method = 'ward.D',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                   angle_col = '315',
                   annotation_row = my.ann,
                   annotation_colors = ha_left.col,
                   main = 'ward.D2.tpm.pearson',
                   display_numbers = T,
                   number_color = 'white',
                   # number_format = "%.2f",
                   fontsize_number = 5,
                   # breaks = unique(c(seq(0.9,1, length=256))),
                   # legend_breaks = c(0.9,0.95,1),
                   fontsize = 10)

corrplot::corrplot(my.cor.p,
                   is.corr = F,
                   col = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                   order = "hclust", addrect = 3,type = "lower")

# DEG--------

## DEG5:pgEpiSC--NGR
group = my.data[,c(1:6)]
group_list = c(rep('pgEpiSC',2),rep('NGR.iPS',4)) 
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','NGR.iPS','pgEpiSC'))
DEG5 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

## DEG6:pgEpiSC--pEF-EPI
group = my.data[,c(1:2,9:12)]
group_list = c(rep('pgEpiSC',2),rep('pEF_EPI',4)) 
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pEF_EPI','pgEpiSC'))
DEG6 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

## DEG7:pgEpiSC--pEF_pMX
group = my.data[,c(1:2,13:16)]
group_list = c(rep('pgEpiSC',2),rep('pEF_pMX',4)) 
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pEF_pMX','pgEpiSC'))
DEG7 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

# venn-----

# up
my.up <- list(NGR.iPS = DEG5[DEG5$log2FoldChange > 2 & DEG5$padj < 0.05,]$Symbol,
              pEF.iPS = DEG6[DEG6$log2FoldChange > 2 & DEG6$padj < 0.05,]$Symbol,
              pEF.pMX = DEG7[DEG7$log2FoldChange > 2 & DEG7$padj < 0.05,]$Symbol)
my.up.venn <- Venn(my.up)
Vennerable::plot(my.up.venn,doWeights = F)

venn(my.up, 
     zcolor=new.col[12:14],
     opacity = .7,
     box=F,
     sncs=1.2,
     ilcs=1.2)

# down
my.down <- list(NGR.iPS = DEG5[DEG5$log2FoldChange < -2 & DEG5$padj < 0.05,]$Symbol,
                pEF.iPS = DEG6[DEG6$log2FoldChange < -2 & DEG6$padj < 0.05,]$Symbol,
                pEF.pMX = DEG7[DEG7$log2FoldChange < -2 & DEG7$padj < 0.05,]$Symbol)

my.down.venn <- Venn(my.down)
Vennerable::plot(my.down.venn,doWeights = F)

venn(my.down, 
     zcolor=new.col[12:14],
     opacity = .7,
     box=F,
     sncs=1.2,
     ilcs=1.2)

# deg heatmap
pheatmap::pheatmap(t(scale(t(my.tpm[unique(c(my.down$pgEpiSC,my.down$NGR.iPS,my.down$pEF.iPS,my.down$pEF.pMX, my.up$pgEpiSC,my.up$NGR.iPS,my.up$pEF.iPS,my.up$pEF.pMX)),-c(7:8)]))),
                   border_color = NA,
                   clustering_method = 'average',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .8),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = F,
                   annotation_colors = ha_left.col,
                   main = 'All.DEG.TPM',
                   breaks = unique(c(seq(-2,3, length=256))),
                   clustering_distance_rows = 'minkowski',
                   fontsize = 10)

pheatmap::pheatmap(t(scale(t(my.tpm[unique(c(my.down.venn@IntersectionSets$`111`, my.up$pgEpiSC,my.up$NGR.iPS,my.up$pEF.iPS,my.up$pEF.pMX)),-c(7,8)]))),
                   border_color = NA,
                   clustering_method = 'average',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .8),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = F,
                   annotation_colors = ha_left.col,
                   main = 'co-down.all-up.DEG.TPM',
                   breaks = unique(c(seq(-2,3, length=256))),
                   clustering_distance_rows = 'binary',
                   fontsize = 10)

pheatmap::pheatmap(t(scale(t(my.tpm[my.up.venn@IntersectionSets$`111`,-c(7,8)]))),
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = T,
                   annotation_colors = ha_left.col,
                   main = 'CO.DEG.UP.TPM',
                   breaks = unique(c(seq(-2,3, length=256))),
                   fontsize = 10)

pheatmap::pheatmap(t(scale(t(my.tpm[my.down.venn@IntersectionSets$`111`,-c(7,8)]))),
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = T,
                   annotation_colors = ha_left.col,
                   main = 'CO.DEG.DOWN.TPM',
                   breaks = unique(c(seq(-2,3, length=256))),
                   fontsize = 10)

# orth 
my.up.DEG5.orth <- data.frame(NGR.iPS = my.up$NGR.iPS) %>% left_join(y = orth[,2:4], by = c('NGR.iPS' = 'genename'))
my.up.DEG6.orth <- data.frame(pEF.iPS = my.up$pEF.iPS) %>% left_join(y = orth[,2:4], by = c('pEF.iPS' = 'genename'))
my.up.DEG7.orth <- data.frame(pEF.pMX = my.up$pEF.pMX) %>% left_join(y = orth[,2:4], by = c('pEF.pMX' = 'genename'))
my.up.co <- data.frame(co.up = my.up.venn@IntersectionSets$`111`) %>% left_join(y = orth[,2:4], by = c('co.up' = 'genename'))

my.down.DEG5.orth <- data.frame(NGR.iPS = my.down$NGR.iPS) %>% left_join(y = orth[,2:4], by = c('NGR.iPS' = 'genename'))
my.down.DEG6.orth <- data.frame(pEF.iPS = my.down$pEF.iPS) %>% left_join(y = orth[,2:4], by = c('pEF.iPS' = 'genename'))
my.down.DEG7.orth <- data.frame(pEF.pMX = my.down$pEF.pMX) %>% left_join(y = orth[,2:4], by = c('pEF.pMX' = 'genename'))
my.down.co <- data.frame(co.down = my.down.venn@IntersectionSets$`111`) %>% left_join(y = orth[,2:4], by = c('co.down' = 'genename'))

for.enrich.genes.pg_ips <- list(my.up.DEG5.orth,my.up.DEG6.orth,my.up.DEG7.orth,my.up.co,my.down.DEG5.orth,my.down.DEG6.orth,my.down.DEG7.orth,my.down.co)
export(for.enrich.genes.pg_ips,file = '../1.part1/for.enrich.genes.pg_ips.xlsx')

