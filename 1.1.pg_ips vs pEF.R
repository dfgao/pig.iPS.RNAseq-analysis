### Part1-1:pg_ips vs PEF

key.genes <- import_list('../1.part1/pickgenes.xlsx')

my.data <- all.data[,c(54:69)]
my.tpm <- all.tpm[,c(54:69)]

mydata.genes <- unique(c(rownames(my.tpm[(matrixStats::rowMedians(as.matrix(my.tpm), cols = 1:2) > 0.5) ,]),
                                          rownames(my.tpm[(matrixStats::rowMedians(as.matrix(my.tpm), cols = 3:4) > 0.5) ,]),
                                          rownames(my.tpm[(matrixStats::rowMedians(as.matrix(my.tpm), cols = 5:6) > 0.5) ,]),
                                          rownames(my.tpm[(matrixStats::rowMedians(as.matrix(my.tpm), cols = 7:8) > 0.5) ,]),
                                          rownames(my.tpm[(matrixStats::rowMedians(as.matrix(my.tpm), cols = 9:10) > 0.5) ,]),
                                          rownames(my.tpm[(matrixStats::rowMedians(as.matrix(my.tpm), cols = 11:12) > 0.5) ,]),
                                          rownames(my.tpm[(matrixStats::rowMedians(as.matrix(my.tpm), cols = 13:14) > 0.5) ,]),
                                          rownames(my.tpm[(matrixStats::rowMedians(as.matrix(my.tpm), cols = 15:16) > 0.5) ,])
))
my.data <- my.data[mydata.genes,] # 15808 genes
my.tpm <- my.tpm[mydata.genes,] 

# correlation----
library(ComplexHeatmap)
library(circlize)
library(corrplot)
my.ann <- rbind(geneinfo.new[48:61,c(1,7:8)],geneinfo.all[60:61,1:3])
rownames(my.ann) <- my.ann$Sample
my.ann <- my.ann[,-1]

ha_left.col <- list(Cell_type = new.col[11:15], Author = pal_simpsons("springfield")(3))
names(ha_left.col$Cell_type) <- factor(unique(my.ann$Cell_type),levels = c("pgEpiSC", "NGR.iPS", "pEF.iPS", "pEF.pMX.iPS",'pEF'))
names(ha_left.col$Author) <- factor(unique(my.ann$Author),levels = c('Zhi','Zhu','Xu'))

my.cor.p <- cor(my.tpm) %>% round(digits = 2)
pheatmap::pheatmap(my.cor.p,
                         border_color = 'gray80',
                         clustering_method = 'ward.D2',
                         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                         angle_col = '315',
                         annotation_row = my.ann,
                         annotation_colors = ha_left.col,
                         main = 'ward.D2.tpm.pearson',
                         display_numbers = T,
                         number_color = 'white',
                         # number_format = "%.2f",
                         fontsize_number = 5,
                         breaks = unique(c(seq(0.4,1, length=256))),
                         fontsize = 10)

pheatmap::pheatmap(my.cor.p,
                         border_color = 'gray80',
                         clustering_method = 'ward.D2',
                         color = scales::alpha(colorRampPalette(colors = c('#00509d','gray80','#f35b04'),alpha=T,bias=1)(256),alpha = 1),
                         angle_col = '315',
                         annotation_row = my.ann,
                         annotation_colors = ha_left.col,
                         main = 'ward.D2.tpm.pearson',
                         display_numbers = T,
                         number_color = 'white',
                         # number_format = "%.2f",
                         fontsize_number = 5,
                         breaks = unique(c(seq(0.4,1, length=256))),
                         fontsize = 10)

# PCA-------
library(ggplot2)
library(factoextra)
library(FactoMineR)

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
               palette = new.col[c(11:15)],
               legend.title = "Groups",
               pointsize = 4,
               pointshape = 21,
               col.ind = "black",
               title = 'pEF&pg&piPS PCA',
               # addEllipses = T
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}

group_list <- factor(geneinfo.all$ct3[c(54:69)],levels = c('Zhi_pgEpiSC_3i_LAF','Zhu_NGR.iPS_3i_LAF','Zhu_pEF.iPS_3i_LAF','Zhu_pEF.pMX.iPS_3i_LAF','Xu_pEF_NA'))
d_p(my.tpm,group_list)

# use rgl plot 3D PCA 

tpm.t = as.data.frame(t(my.tpm))
tpm.t = cbind(tpm.t,group_list)
tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
a = tpm.pca[["ind"]][["coord"]] 

open3d() 
par3d(family="serif",
      cex=1.2,
      font=1) 

plot3d(a[,1:3],
       type="p",
       col=new.col[c(11:15)][as.numeric(group_list)],
       lwd = 1,
       size=12,
       ticktype = "detailed",
       #theta = 20,
       #phi = 10, 
       #fov = 60, 
       #zoom = 1 
)


view3d(theta = 20,
       phi = 10, 
       fov = 90, 
       zoom = .8 )
legend3d("right",
         col = new.col,
         legend = levels(group_list),
         pt.bg =  new.col, 
         pch = 16,
         # inset = -0.5, 
         xpd = T, 
         horiz = F,
         bty = 'n',
         cex = 1.5,
         bg = "transparent")

rgl.postscript("../1.part1/rgl.PCA-dim3-3D.pdf", fmt = "pdf", drawText = T)

close3d() 

# loading genes-----

## factoextra loading score
factor.score <- tpm.pca$var$contrib[,1]
factor.score <- sort(factor.score, decreasing = F)
factor.score <- as.data.frame(factor.score)
colnames(factor.score) <- 'PC1'
factor.score$Rank <- seq(1,nrow(factor.score))
top10.genes <- rownames(factor.score)[1:10]

ggplot(factor.score, aes(x=Rank, y=PC1)) +
  geom_point() +
  scale_x_log10()

# prcomp score
a <- prcomp(tpm.t[,-ncol(tpm.t)],scale=TRUE)
pca.data <- data.frame(Sample=rownames(a$x),
                       X=a$x[,1],
                       Y=a$x[,2])
pca.var <- a$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
  theme_bw()+
  ggtitle("My PCA Graph")

# only PC1
prcomp.score <- a$rotation[,1]
prcomp.score <- sort(prcomp.score, decreasing = F)
prcomp.score <- as.data.frame(prcomp.score)
colnames(prcomp.score) <- 'PC1'
prcomp.score$Rank <- seq(1,nrow(prcomp.score))

# prcomp.score <- prcomp.score %>%
#   rownames_to_column(var="Symbol") %>% 
#   mutate(highlight = case_when(
#     str_detect(string = Symbol, pattern = key.genes$Pluripotency$Pluripotency) ~ Symbol,
#     TRUE ~ ""
#   ))

prcomp.score <- prcomp.score %>%
  rownames_to_column(var="Symbol") %>% 
  mutate(highlight = "")

prcomp.score[c(na.omit(match(key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)], prcomp.score$Symbol))), 4] <- key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)]
p1 <- ggplot(prcomp.score, aes(x=Rank, y=PC1)) +
  geom_point() +
  geom_text_repel(
    aes(label = highlight),
    size = 3,
    min.segment.length = 0,
    seed = 100,
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = 800,
    nudge_y = .001,
    color = "grey30"
  ) + 
  geom_rug(col="steelblue",alpha=0.1, size=.1,sides = 'l') +
  theme_bw()
p1
ggsave(filename = '../1.part1/PC1.key.genes.loading.score.pdf',height = 5, width = 6)

# PC1&2
prcomp.score <- data.frame(a$rotation[,1:2])
prcomp.score <- prcomp.score[order(prcomp.score$PC1, decreasing = F),]
prcomp.score$Rank <- seq(1,nrow(prcomp.score))

prcomp.score <- prcomp.score %>%
  rownames_to_column(var="Symbol") %>% 
  mutate(highlight = "",
         showcol = case_when(
           PC1 > 0.012 | PC1 < -0.012 ~ "yes",
           TRUE ~ "no"
         ))
prcomp.score[c(na.omit(match(key.genes$Pluripotency$Pluripotency, prcomp.score$Symbol))), 5] <- key.genes$Pluripotency$Pluripotency[-c(14,16,17)]

prcomp.score$Z.PC1 <- scale(prcomp.score$PC1)
prcomp.score$Z.PC2 <- scale(prcomp.score$PC2)

prcomp.score <- prcomp.score %>%
  mutate(showcol.radius = case_when(
           sqrt(Z.PC1^2 + Z.PC2^2) > 2 ~ "yes",
           TRUE ~ "no"
         ))
table(prcomp.score$showcol.radius)

ggplot(prcomp.score, aes(x=Z.PC1, y=Z.PC2, color = showcol.radius)) +
  geom_point(size = .1) +
  geom_text_repel(
    aes(label = highlight),
    family = "Poppins",
    size = 3,
    min.segment.length = 0,
    seed = 100,
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    # nudge_x = 800,
    # nudge_y = .001,
    color = "grey30"
  ) +
  xlim(-5,5) +
  ylim(-5,5) +
  scale_color_npg()+
  theme_ipsum()
  # scale_fill_manual(values=c("#999999", "#E69F00"),
  #   name = "1",
  #   labels = c('Radius < 2*sd', 'Radius > 2*sd'))
  

# DEG-----
colnames(my.data)

# part1: pg vs PEF
group = my.data[,c(1,2,7,8)]
group_list = c(rep('pgEpiSC',2),rep('pEF',2)) 
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pgEpiSC','pEF'))
DEG1 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

# part2: NGR vs PEF
group = my.data[,c(3:6,7,8)]
group_list = c(rep('NGR.iPS',4),rep('pEF',2))
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','NGR.iPS','pEF'))
DEG2 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

# part3: PEF-EPI vs PEF
group = my.data[,c(9:12,7,8)]
group_list = c(rep('pEF_EPI',4),rep('pEF',2))
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pEF_EPI','pEF'))
DEG3 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

# part4: PEF-pMX vs PEF
group = my.data[,c(13:16,7,8)]
group_list = c(rep('pEF_pMX',4),rep('pEF',2))
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pEF_pMX','pEF'))
DEG4 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

# venn-----
library(VennDiagram) 
library(venn)
# up
my.up <- list(pgEpiSC = DEG1[DEG1$log2FoldChange > 2 & DEG1$padj < 0.05,]$Symbol, 
              NGR.iPS = DEG2[DEG2$log2FoldChange > 2 & DEG2$padj < 0.05,]$Symbol,
              pEF.iPS = DEG3[DEG3$log2FoldChange > 2 & DEG3$padj < 0.05,]$Symbol,
              pEF.pMX = DEG4[DEG4$log2FoldChange > 2 & DEG4$padj < 0.05,]$Symbol)
my.up.venn <- Venn(my.up)
Vennerable::plot(my.up.venn,doWeights = F)

venn(my.up, 
     zcolor=new.col[11:14],
     opacity = .7,
     box=F,
     sncs=1.2,
     ilcs=1.2)

# down
my.down <- list(pgEpiSC = DEG1[DEG1$log2FoldChange < -2 & DEG1$padj < 0.05,]$Symbol, 
              NGR.iPS = DEG2[DEG2$log2FoldChange < -2 & DEG2$padj < 0.05,]$Symbol,
              pEF.iPS = DEG3[DEG3$log2FoldChange < -2 & DEG3$padj < 0.05,]$Symbol,
              pEF.pMX = DEG4[DEG4$log2FoldChange < -2 & DEG4$padj < 0.05,]$Symbol)

my.down.venn <- Venn(my.down)
Vennerable::plot(my.down.venn,doWeights = F)

venn(my.down, 
     zcolor=new.col[11:14],
     opacity = .7,
     box=F,
     sncs=1.2,
     ilcs=1.2)

# deg heatmao
pheatmap::pheatmap(t(scale(t(my.tpm[unique(c(my.down$pgEpiSC,my.down$NGR.iPS,my.down$pEF.iPS,my.down$pEF.pMX, my.up$pgEpiSC,my.up$NGR.iPS,my.up$pEF.iPS,my.up$pEF.pMX)),]))),
                         border_color = NA,
                         clustering_method = 'ward.D2',
                         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .8),
                         angle_col = '315',
                         annotation_col = my.ann,
                         show_rownames = F,
                         annotation_colors = ha_left.col,
                         main = 'All.DEG.TPM',
                   breaks = unique(c(seq(-2,3, length=256))),
                         fontsize = 10)

pheatmap::pheatmap(t(scale(t(my.tpm[unique(c(my.down.venn@IntersectionSets$`1111`, my.up$pgEpiSC,my.up$NGR.iPS,my.up$pEF.iPS,my.up$pEF.pMX)),]))),
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .8),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = F,
                   annotation_colors = ha_left.col,
                   main = 'co-down.all-up.DEG.TPM',
                   breaks = unique(c(seq(-2,3, length=256))),
                   fontsize = 10)

pheatmap::pheatmap(t(scale(t(my.tpm[my.up.venn@IntersectionSets$`1111`,]))),
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = F,
                   annotation_colors = ha_left.col,
                   main = 'CO.DEG.UP.TPM',
                   breaks = unique(c(seq(-2,3, length=256))),
                   fontsize = 10)

pheatmap::pheatmap(t(scale(t(my.tpm[my.down.venn@IntersectionSets$`1111`,]))),
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = F,
                   annotation_colors = ha_left.col,
                   main = 'CO.DEG.DOWN.TPM',
                   breaks = unique(c(seq(-2,3, length=256))),
                   fontsize = 10)
# volcano-----
# need deteach ggtern
p1 <- EnhancedVolcano::EnhancedVolcano(na.omit(DEG1),
                lab = na.omit(DEG1)$Symbol,
                selectLab = key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)],
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 2,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.5,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'pgEpiSC vs pEF',
                subtitle = ''
) + coord_flip()
ggsave(p1,filename = '../1.part1/pgEpiSC vs pEF.volcano.pdf',h=6,w=8)

p1 <- EnhancedVolcano(DEG2,
                lab = DEG2$Symbol,
                selectLab = key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)],
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 2,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.5,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'NGR.iPS vs pEF',
                subtitle = '',
) + coord_flip()
ggsave(p1,filename = '../1.part1/NGR.iPS vs pEF.volcano.pdf',h=6,w=8)

p1 <- EnhancedVolcano(DEG3,
                lab = DEG3$Symbol,
                selectLab = key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)],
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 2,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.5,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'pEF.iPS vs pEF',
                subtitle = '',
) + coord_flip()
ggsave(p1,filename = '../1.part1/pEF.iPS vs pEF.volcano.pdf',h=6,w=8)

p1 <- EnhancedVolcano(DEG4,
                lab = DEG4$Symbol,
                selectLab = key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)],
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 2,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.5,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'pEF.pMX vs pEF',
                subtitle = '',
) + coord_flip()
ggsave(p1,filename = '../1.part1/pEF.pMX vs pEF.volcano.pdf',h=6,w=8)

# kmeans NOT USE IN HERE-----

my.tpm.median <- data.frame(pgEpiSC = matrixStats::rowMedians(as.matrix(my.tpm),cols = 1:2),
                            NGR.iPS = matrixStats::rowMedians(as.matrix(my.tpm),cols = 3:6),
                            pEF.iPS = matrixStats::rowMedians(as.matrix(my.tpm),cols = 9:12),
                            pEF.pMX = matrixStats::rowMedians(as.matrix(my.tpm),cols = 13:16),
                            pEF = matrixStats::rowMedians(as.matrix(my.tpm),cols = 7:8))
rownames(my.tpm.median) <- rownames(my.tpm)
my.tpm.median.scale <- data.frame(t(scale(t(my.tpm.median))))

nk=2:15
set.seed(222)
Wss<-sapply(nk,function(k){
  kmeans(na.omit(my.tpm.median.scale),centers = k,iter.max = 100)$tot.withinss})
plot(nk,Wss,type = "l",xlab="Number of k",ylab="Within sum of squares")
plot(nk,Wss,type = "o",xlab="Number of k",ylab="Within sum of squares",col='red3')
abline(v=7,col='black')
legend("topleft", legend = 'k = 7', 
       col= "red3",
       pch = 15, bty = "n", pt.cex = 2, cex = 1.2,  horiz = F, inset =  0.1)

fit <- kmeans(x = na.omit(my.tpm.median.scale),7,iter.max = 100)

table(fit$cluster)
my.tpm.median.scale <- na.omit(my.tpm.median.scale)
my.tpm.median.scale$kmeans.cluster <- fit$cluster
my.tpm.median.scale$gene <- names(fit[["cluster"]])

kmectocluster <- data.frame(cluster = fit$cluster)

my.tpm.median.log2 <- log2(my.tpm.median+1)
my.tpm.median.log2$gene <- rownames(my.tpm.median.log2)
my.tpm.median.log2 <- my.tpm.median.log2[my.tpm.median.scale$gene,]

kmectocluster <- 
  kmectocluster %>% 
  rownames_to_column(var = 'gene') %>%
  left_join(y = my.tpm.median.log2,
            by = 'gene')
kmectocluster$cluster <- paste0('C_',kmectocluster$cluster)
kmectocluster$cluster <- factor(kmectocluster$cluster,levels = c(paste0('C_',seq(1,7))))
kmectocluster <- kmectocluster %>% group_by(cluster) %>%  mutate(cluster=paste0(cluster,"(",n()," genes)"))
table(kmectocluster$cluster)
kmectocluster$cluster <- factor(kmectocluster$cluster,levels = c('C_1(1090 genes)', 'C_2(1776 genes)', 'C_3(1675 genes)', 'C_4(4103 genes)',  'C_5(3770 genes)', 'C_6(928 genes)',  'C_7(2239 genes)'))

kmectocluster %>% 
  select("pgEpiSC", "NGR.iPS" ,"pEF.iPS" ,"pEF.pMX", "pEF" ,'cluster') %>%
  gather(key = 'stage',value = 'exp', -cluster) %>%
  mutate(stage=fct_relevel(stage, "pgEpiSC", "NGR.iPS" ,"pEF.iPS" ,"pEF.pMX", "pEF" )) %>%
  ggplot( aes(x=stage, y=round(exp,4),group=1,fill=cluster)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               #width=.2, 
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line") +
  labs(x="Stage(cluster base on scale data)", y="Expression level(log2(TPM+1))") +
  theme_bw() +
  theme(legend.position="none") +
  scale_fill_manual(values = rep('red',50),
                    guide=guide_legend(direction="vertical",
                                       label.position="right",
                                       title=NULL,
                                       ncol=6,
                                       label.hjust=0.8))+
  scale_color_manual(values =  rep('red',50),guide = 'none')+
  # geom_smooth()+
  # ylim(0,1)+
  facet_wrap(~cluster,scales = 'free_y',ncol = 4,drop = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("pgEpiSC", "NGR.iPS" ,"pEF.iPS" ,"pEF.pMX", "pEF" ))

# WGCNA NOT USE IN HERE--------
library(WGCNA)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

my.meta <- data.frame(stage = factor(c(rep('pgEpiSC',2), rep('NGR.iPS',4), rep('pEF',2), rep('pEF.iPS',4),rep('pEF.pMX',4)), levels = c('pgEpiSC','NGR.iPS','pEF','pEF.iPS','pEF.pMX')))
rownames(my.meta) <- colnames(my.tpm)

datExpr <- t(my.tpm) 

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# β value
powers = c(c(1:10), seq(from = 12, to=500, by=4))
#设置beta值的取值范围
sft = pickSoftThreshold(datExpr, RsquaredCut = 0.9,powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="green")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1,col="red")
#选择beta值
best_beta=sft$powerEstimate

net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       maxBlockSize = 6000,TOMType = "unsigned", 
                       minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "my.data.TPM-TOM",
                       verbose = 3)


# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
moduleColors=mergedColors
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

datExpr_tree<-hclust(dist(datExpr), method = "ward.D2")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
#针对前面构造的样品矩阵添加对应颜色
sample_colors <- numbers2colors(as.numeric(factor(my.meta$stage)), 
                                colors = c("grey","blue","red","green"),signed = FALSE)
## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目。
#  sample_colors <- numbers2colors( datTraits ,signed = FALSE)
## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，
##当然，这样给的颜色不明显，意义不大。
#构造10个样品的系统聚类树及性状热图
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = rownames(my.meta),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

design=model.matrix(~0+ my.meta$stage)
colnames(design)=unique(my.meta$stage)
# names(design) <- colnames(my.tpm)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


connet=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(connet, moduleColors)
head(Alldegrees1)

# (2) Relationship between gene significance and intramodular connectivity
which.module="black"
EB= as.data.frame(design[,2]); # change specific 
names(EB) = "EB"
GS1 = as.numeric(cor(EB,datExpr, use="p"))
GeneSignificance=abs(GS1)
colorlevels=unique(moduleColors)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}



module = "turquoise";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因
inModule = (moduleColors==module)
modProbes = probes[inModule]
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste('C0_p11', collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste('C0_p11', collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule]
)

# get orth-----
orth <- import('~/project/01.zhugaoxiang/04.pig.RNA.MSC/05.analysis/pigref105_human.unique.xlsx',which = 'Sheet1')
my.up.DEG1.orth <- data.frame(pgEpiSC = my.up$pgEpiSC) %>% left_join(y = orth[,2:4], by = c('pgEpiSC' = 'genename'))
my.up.DEG2.orth <- data.frame(NGR.iPS = my.up$NGR.iPS) %>% left_join(y = orth[,2:4], by = c('NGR.iPS' = 'genename'))
my.up.DEG3.orth <- data.frame(pEF.iPS = my.up$pEF.iPS) %>% left_join(y = orth[,2:4], by = c('pEF.iPS' = 'genename'))
my.up.DEG4.orth <- data.frame(pEF.pMX = my.up$pEF.pMX) %>% left_join(y = orth[,2:4], by = c('pEF.pMX' = 'genename'))
my.up.co <- data.frame(co.up = my.up.venn@IntersectionSets$`1111`) %>% left_join(y = orth[,2:4], by = c('co.up' = 'genename'))

my.down.DEG1.orth <- data.frame(pgEpiSC = my.down$pgEpiSC) %>% left_join(y = orth[,2:4], by = c('pgEpiSC' = 'genename'))
my.down.DEG2.orth <- data.frame(NGR.iPS = my.down$NGR.iPS) %>% left_join(y = orth[,2:4], by = c('NGR.iPS' = 'genename'))
my.down.DEG3.orth <- data.frame(pEF.iPS = my.down$pEF.iPS) %>% left_join(y = orth[,2:4], by = c('pEF.iPS' = 'genename'))
my.down.DEG4.orth <- data.frame(pEF.pMX = my.down$pEF.pMX) %>% left_join(y = orth[,2:4], by = c('pEF.pMX' = 'genename'))
my.down.co <- data.frame(co.down = my.down.venn@IntersectionSets$`1111`) %>% left_join(y = orth[,2:4], by = c('co.down' = 'genename'))

for.enrich.genes <- list(my.up.DEG1.orth,my.up.DEG2.orth,my.up.DEG3.orth,my.up.DEG4.orth,my.up.co,my.down.DEG1.orth,my.down.DEG2.orth,my.down.DEG3.orth,my.down.DEG4.orth,my.down.co)
export(for.enrich.genes,file = '../1.part1/for.enrich.genes.xlsx')

# pluripotency genes--------
plu.my.tpm <- my.tpm[key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)],] %>% na.omit()
pheatmap::pheatmap(t(scale(t(plu.my.tpm))),
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = T,
                   # annotation_colors = ha_left.col,
                   main = 'Pluripotency genes',
                   fontsize = 10)

pheatmap::pheatmap(log10(plu.my.tpm+1),
                   border_color = NA,
                   clustering_method = 'median',
                   color = scales::alpha(colorRampPalette(RColorBrewer::brewer.pal(8,"YlOrRd"),alpha=T,bias=1)(256),alpha = 1),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = T,
                   annotation_colors = ha_left.col,
                   main = 'Pluripotency genes',
                   # breaks = unique(c(seq(0,4, length=256))),
                   fontsize = 10)

pheatmap::pheatmap(log10(plu.my.tpm+1),
                   border_color = NA,
                   clustering_method = 'median',
                   color = rev(heat.colors(256)),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = T,
                   annotation_colors = ha_left.col,
                   main = 'Pluripotency genes',
                   # breaks = unique(c(seq(0,4, length=256))),
                   fontsize = 10)

pheatmap::pheatmap(log10(plu.my.tpm+1),
                   border_color = NA,
                   clustering_method = 'centroid',
                   color = viridis_pal()(256),
                   angle_col = '315',
                   annotation_col = my.ann,
                   show_rownames = T,
                   annotation_colors = ha_left.col,
                   main = 'Pluripotency genes',
                   breaks = unique(c(seq(0,3.5, length=256))),
                   legend_breaks = c(0:3),
                   fontsize = 10)
# ternary--------
library(ggtern)

for.tern.tpm <- data.frame(pgEpiSC = matrixStats::rowMedians(as.matrix(my.tpm),cols = 1:2),
                           iPS = matrixStats::rowMedians(as.matrix(my.tpm),cols = c(3:6,9:16)),
                           pEF = matrixStats::rowMedians(as.matrix(my.tpm),cols = 7:8))
for.tern.tpm$gene <- rownames(my.tpm)
for.tern.tpm$average <- rowMeans(for.tern.tpm[,1:3])
for.tern.tpm$gene.type <- 'Other'
for.tern.tpm[c(na.omit(match(key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)], for.tern.tpm$gene))), 6] <- 'Pluripotency'

for.tern.tpm$gene.type <- factor(for.tern.tpm$gene.type,levels = c('Pluripotency','Other'))
for.tern.tpm <- for.tern.tpm[order(for.tern.tpm$gene.type),]

ggtern(data=for.tern.tpm,aes(x=pgEpiSC,y=iPS,z=pEF))+    
  geom_point(data = for.tern.tpm[for.tern.tpm$gene.type == 'Other',], aes(size = average), color = "grey80",
             alpha = 0.8, show.legend = FALSE) +
  # scale_colour_manual(values = c('red3',"grey90"))+ 
  theme_rgbw(base_size = 12 )+   
  labs(title = "3 cell type genes")+ 
  theme(plot.title = element_text(size=15,hjust = 0.5)) +
  theme_showarrows() + 
  # geom_density_tern(h=2,
  #                   bdl = 0.08,
  #                   expand=0.75,
  #                   alpha=0.5,
  #                   bins=5)+ 
  stat_density_tern(aes(fill=..level.., alpha=..level..),geom='polygon',bdl = 0.08, h = 2,alpha=0.8,expand=0.75,bins = 6) +
  scale_fill_gradient2(high = "#71D0F5FF") +
  geom_point(data = for.tern.tpm[for.tern.tpm$gene.type == 'Pluripotency',], aes(size = average, color = '#B2182B'), 
             show.legend = F)+    
  geom_mask() 
  # geom_text(data = for.tern.tpm[for.tern.tpm$gene.type == 'Pluripotency',],
  #           aes(label = gene),
  #           check_overlap = T,
  #           color = 'red3')
  

# ggtern(data=for.tern.tpm,aes(x=pgEpiSC,y=iPS,z=pEF))+    
#   geom_point( aes(size = average,color = gene.type),
#              alpha = 0.8, show.legend = FALSE) +   
#   scale_colour_manual(values = c('red3',"grey90"))+
#   theme_rgbw(base_size = 12 )+   
#   labs(title = "3 cell type genes")+ 
#   theme(plot.title = element_text(size=15,hjust = 0.5)) +
#   # theme_showarrows() + 
#   geom_density_tern(h=2, 
#                     bdl = 0.08,
#                     expand=0.75,alpha=0.5,bins=5) +
#  geom_text(data = for.tern.tpm[for.tern.tpm$gene.type == 'Pluripotency',],
#                                         aes(label = gene))
