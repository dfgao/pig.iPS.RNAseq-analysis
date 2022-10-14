#### kmeans for two culture system -- e xogenous Gene silence/not silence

### get high expression genes matrix ---------
sys2.data <- all.data[,c(1:7, 17:22, 25:34, 56:59, 62:69)]
sys2.tpm <- all.tpm[,c(1:7, 17:22, 25:34, 56:59, 62:69)]
sys2.info <- geneinfo.all[c(1:7, 17:22, 25:34, 56:59, 62:69), -(2:6)]
sys2.info$Exogenous = c(rep('Silencing',2), rep('Activation',9), rep('Silencing',2), rep('Activation',10), rep('Silencing',12))

sys2.high.genes <- unique(c(rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 1:2) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 3:5) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 6:7) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 8:9) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 10:11) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 12:13) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 14:16) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 17:19) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 20:21) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 22:23) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 24:27) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 28:31) > 5) ,]),
                         rownames(sys2.tpm[(matrixStats::rowMedians(as.matrix(sys2.tpm), cols = 32:35) > 5) ,])
))

sys2.high.tpm <- sys2.tpm[sys2.high.genes,]

### corralation --------
sys.anno <- sys2.info[,c(2,5,7)]
rownames(sys.anno) <- sys2.info$Sample

ha_top.col <- list(Exogenous = new.col[c(16:17)], 
                   Author = c(pal_simpsons("springfield")(16)[7:14]), 
                   System = pal_igv("default")(11))
names(ha_top.col$Exogenous) <- factor(c("Silencing","Activation"))
names(ha_top.col$Author) <- factor(unique(sys2.info$Author))
names(ha_top.col$System) <- factor(unique(sys2.info$System))

sys.cor.s <- cor(sys2.high.tpm,method = 'spearman') %>% round(digits = 2)
pheatmap::pheatmap(sys.cor.s,
                   border_color = 'NA',
                   clustering_method = 'ward.D2',
                   color = scales::alpha(colorRampPalette(colors = c('#00509d','gray80','#f35b04'),alpha=T,bias=1)(256),alpha = 1),
                   angle_col = '315',
                   annotation_col = sys.anno,
                   annotation_colors = ha_top.col,
                   main = 'ward.D2.tpm.pearson',
                   display_numbers = T,
                   number_color = 'white',
                   # number_format = "%.2f",
                   fontsize_number = 4,
                   breaks = unique(c(seq(0.5,1, length=256))),
                   fontsize = 9)
### kmeans for median-------
## for scale data

sys2.high.tpm.median <- data.frame(	Gao_iPS = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 1:2),
                                    Xu_iPS = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 3:5),
                                    Shi_iPS = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 6:7),
                                    Yuan_iPS_FL6i = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 8:9),
                                    Yuan_iPS_FLB2i = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 10:11),
                                    Yoshimatsu_iPS = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 12:13),
                                    Secher_iPS_F = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 14:16),
                                    Secher_iPS_L = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 17:19),
                                    Mao_iPS_L_Rex1_pos = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 20:21),
                                    Mao_iPS_L_Rex1_neg = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 22:23),
                                    Zhu_NGR = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 24:27),
                                    Zhu_EPI = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 28:31),
                                    Zhu_pMX = matrixStats::rowMedians(as.matrix(sys2.high.tpm),cols = 32:35))
rownames(sys2.high.tpm.median) <- rownames(sys2.high.tpm)
sys2.high.tpm.median.scale <- data.frame(t(scale(t(sys2.high.tpm.median))))

nk=2:100
Wss<-sapply(nk,function(k){
  kmeans(sys2.high.tpm.median.scale,centers = k,iter.max = 100)$tot.withinss})
plot(nk,Wss,type = "o",xlab="Number of k",ylab="Within sum of squares",col='red3',)
abline(v=15,col='black')
legend("topright", legend = 'k = 15', 
       col= "red3",
       pch = 15, bty = "n", pt.cex = 2, cex = 1.2,  horiz = F, inset =  0.1)


fit <- kmeans(x = sys2.high.tpm.median.scale,15,iter.max = 100)
# kmeansBIC(fit)
table(fit$cluster)
sys2.high.tpm.median.scale$kmeans.cluster <- fit$cluster
sys2.high.tpm.median.scale$gene <- rownames(sys2.high.tpm.median)

kmectocluster <- data.frame(cluster = fit$cluster)

sys2.high.tpm.median.log2 <- log2(sys2.high.tpm.median+1)
sys2.high.tpm.median.log2$gene <- rownames(sys2.high.tpm.median.log2)

kmectocluster <- 
  kmectocluster %>% 
  rownames_to_column(var = 'gene') %>%
  left_join(y = sys2.high.tpm.median.log2,
            by = 'gene')
kmectocluster$cluster <- paste0('C_',kmectocluster$cluster)
kmectocluster$cluster <- factor(kmectocluster$cluster,levels = c(paste0('C_',seq(1,30))))
kmectocluster <- kmectocluster %>% group_by(cluster) %>%  mutate(cluster=paste0(cluster,"(",n()," genes)"))
table(kmectocluster$cluster)
kmectocluster$cluster <- factor(kmectocluster$cluster,levels = c('C_1(828 genes)', 'C_2(1435 genes)', 'C_3(474 genes)', 'C_4(1083 genes)',  'C_5(768 genes)', 'C_6(859 genes)',  'C_7(644 genes)','C_8(774 genes)',  'C_9(707 genes)','C_10(595 genes)', 'C_11(616 genes)', 'C_12(996 genes)', 'C_13(767 genes)', 'C_14(856 genes)', 'C_15(1578 genes)'))

kmectocluster %>% 
  dplyr::select("Gao_iPS","Xu_iPS","Shi_iPS","Yuan_iPS_FL6i", "Yuan_iPS_FLB2i","Yoshimatsu_iPS","Secher_iPS_F","Secher_iPS_L","Mao_iPS_L_Rex1_pos", "Mao_iPS_L_Rex1_neg", "Zhu_NGR","Zhu_EPI",  "Zhu_pMX","cluster" ) %>%
  gather(key = 'System',value = 'exp', -cluster) %>%
  mutate(System=fct_relevel(System, "Xu_iPS","Shi_iPS","Yuan_iPS_FL6i", "Yuan_iPS_FLB2i", "Secher_iPS_F","Secher_iPS_L","Mao_iPS_L_Rex1_pos", "Mao_iPS_L_Rex1_neg", "Gao_iPS", "Yoshimatsu_iPS", "Zhu_NGR","Zhu_EPI", "Zhu_pMX")) %>%
  ggplot( aes(x=System, y=round(exp,4),group=1,fill=cluster)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               #width=.2, 
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line") +
  labs(x="System(cluster base on scale data)", y="Expression level(log2(TPM+1))") +
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
  facet_wrap(~cluster,scales = 'free_y',ncol = 5,drop = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = c(rep('#00A087FF',8), rep('#E64B35FF',5))))

kmectocluster.orth <- kmectocluster %>% left_join(y = orth[,2:4], by = c('gene' = 'genename'))
export(kmectocluster, file = '../4.part4/high.exp.median.scale.kmeans.cluster15.xlsx')
export(kmectocluster.orth, file = '../4.part4/high.exp.median.scale.kmeans.cluster15.orth.xlsx')

## heatmap for C8/14
plot.data <- kmectocluster %>% select( c("Xu_iPS","Shi_iPS","Yuan_iPS_FL6i", "Yuan_iPS_FLB2i", "Secher_iPS_F","Secher_iPS_L","Mao_iPS_L_Rex1_pos", "Mao_iPS_L_Rex1_neg", "Gao_iPS", "Yoshimatsu_iPS", "Zhu_NGR","Zhu_EPI", "Zhu_pMX"))

cluster_anno <- data.frame(Exogenous = c(rep('Activation', 8), rep('Silencing',5)))
rownames(cluster_anno) <- c("Xu_iPS","Shi_iPS","Yuan_iPS_FL6i", "Yuan_iPS_FLB2i", "Secher_iPS_F","Secher_iPS_L","Mao_iPS_L_Rex1_pos", "Mao_iPS_L_Rex1_neg", "Gao_iPS", "Yoshimatsu_iPS", "Zhu_NGR","Zhu_EPI", "Zhu_pMX")

ha_top.col <- list(Exogenous = new.col[c(16:17)])
names(ha_top.col$Exogenous) <- factor(c("Silencing","Activation"))

pheatmap::pheatmap(log2(plot.data[plot.data$cluster == 'C_8(774 genes)', 2:ncol(plot.data)] +1),
                   border_color = 'NA',
                   clustering_method = 'ward.D2',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = cluster_anno,
                   annotation_colors = ha_top.col,
                   main = 'ward.D2.C8',
                   display_numbers = F,
                   fontsize_number = 4,
                   cluster_cols = F,
                   show_rownames = F,
                   breaks = unique(c(seq(0,4, length=256))),
                   fontsize = 9)

pheatmap::pheatmap(t(scale(t(plot.data[plot.data$cluster == 'C_8(774 genes)', 2:ncol(plot.data)]))),
                   border_color = 'NA',
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .8),
                   angle_col = '315',
                   annotation_col = cluster_anno,
                   annotation_colors = ha_top.col,
                   main = 'ward.D2.C8',
                   display_numbers = F,
                   fontsize_number = 4,
                   cluster_cols = F,
                   show_rownames = F,
                   breaks = unique(c(seq(-3,3, length=256))),
                   fontsize = 9)

pheatmap::pheatmap(log2(plot.data[plot.data$cluster == 'C_14(856 genes)', 2:ncol(plot.data)] +1),
                   border_color = 'NA',
                   clustering_method = 'ward.D2',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = cluster_anno,
                   annotation_colors = ha_top.col,
                   main = 'ward.D2.C14',
                   display_numbers = F,
                   fontsize_number = 4,
                   cluster_cols = F,
                   show_rownames = F,
                   breaks = unique(c(seq(0,4, length=256))),
                   fontsize = 9)

pheatmap::pheatmap(t(scale(t(plot.data[plot.data$cluster == 'C_14(856 genes)', 2:ncol(plot.data)]))),
                   border_color = 'NA',
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .8),
                   angle_col = '315',
                   annotation_col = cluster_anno,
                   annotation_colors = ha_top.col,
                   main = 'ward.D2.C14',
                   display_numbers = F,
                   fontsize_number = 4,
                   cluster_cols = F,
                   show_rownames = F,
                   breaks = unique(c(seq(-2,3, length=256))),
                   fontsize = 9)

### kmeans for rank--fail---------

sys2.high.median.rank <- data.frame(apply(sys2.high.tpm.median, 2, function(x){rank(nrow(sys2.high.tpm.median)-rank(x,ties.method = 'first')+1)}))

nk=2:100
Wss<-sapply(nk,function(k){
  kmeans(sys2.high.median.rank,centers = k,iter.max = 100)$tot.withinss})
plot(nk,Wss,type = "o",xlab="Number of k",ylab="Within sum of squares",col='red3',)
abline(v=10,col='black')
legend("topright", legend = 'k = 10', 
       col= "red3",
       pch = 15, bty = "n", pt.cex = 2, cex = 1.2,  horiz = F, inset =  0.1)


fit.rank <- kmeans(x = sys2.high.median.rank,10,iter.max = 100)
# kmeansBIC(fit)
table(fit.rank$cluster)
sys2.high.median.rank$kmeans.cluster <- fit.rank$cluster
sys2.high.median.rank$gene <- rownames(sys2.high.tpm.median)

kmectocluster.rank <- data.frame(cluster = fit.rank$cluster)

# sys2.high.tpm.median.log2 <- log2(sys2.high.tpm.median+1)
# sys2.high.tpm.median.log2$gene <- rownames(sys2.high.tpm.median.log2)

kmectocluster.rank <- 
  kmectocluster.rank %>% 
  rownames_to_column(var = 'gene') %>%
  left_join(y = sys2.high.median.rank,
            by = 'gene')
kmectocluster.rank$cluster <- paste0('C_',kmectocluster.rank$cluster)
kmectocluster.rank$cluster <- factor(kmectocluster.rank$cluster,levels = c(paste0('C_',seq(1,10))))
kmectocluster.rank <- kmectocluster.rank %>% group_by(cluster) %>%  mutate(cluster=paste0(cluster,"(",n()," genes)"))
table(kmectocluster.rank$cluster)
# kmectocluster.rank$cluster <- factor(kmectocluster.rank$cluster,levels = c('C_1(1232 genes)', 'C_2(1523 genes)', 'C_3(1253 genes)', 'C_4(1333 genes)',  'C_5(1349 genes)', 'C_6(644 genes)',  'C_7(1343 genes)','C_8(1830 genes)',  'C_9(1373 genes)','C_10(2666 genes)'))

kmectocluster.rank %>% 
  dplyr::select("Gao_iPS","Xu_iPS","Shi_iPS","Yuan_iPS_FL6i", "Yuan_iPS_FLB2i","Yoshimatsu_iPS","Secher_iPS_F","Secher_iPS_L","Mao_iPS_L_Rex1_pos", "Mao_iPS_L_Rex1_neg", "Zhu_NGR","Zhu_EPI",  "Zhu_pMX","cluster" ) %>%
  gather(key = 'System',value = 'exp', -cluster) %>%
  mutate(System=fct_relevel(System, "Xu_iPS","Shi_iPS","Yuan_iPS_FL6i", "Yuan_iPS_FLB2i", "Secher_iPS_F","Secher_iPS_L","Mao_iPS_L_Rex1_pos", "Mao_iPS_L_Rex1_neg", "Gao_iPS", "Yoshimatsu_iPS", "Zhu_NGR","Zhu_EPI", "Zhu_pMX")) %>%
  ggplot( aes(x=System, y=round(exp,4),group=1,fill=cluster)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               #width=.2, 
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line") +
  labs(x="System(cluster base on scale data)", y="Expression level(log2(TPM+1))") +
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
  facet_wrap(~cluster,scales = 'free_y',ncol = 5,drop = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 5))

kmectocluster.rank.orth <- kmectocluster.rank %>% left_join(y = orth[,2:4], by = c('gene' = 'genename'))
export(kmectocluster.rank.orth, file = '../4.part4/high.exp.rank.kmeans.cluster10.orth.v3.xlsx')

# export(kmectocluster.rank, file = '../4.part4/high.exp.rank.kmeans.cluster10.xlsx')

kmectocluster.rank <- import('../4.part4/high.exp.rank.kmeans.cluster10.xlsx')

# MAPK----
mapk.genes <- read.table('../4.part4/mapk.txt',sep = ',',header = T)

c2.mapk.tpm <- c2.tpm[colnames(mapk.genes),c(1,2,12,13,27:38)] %>% na.omit()

pheatmap::pheatmap(log2(c2.mapk.tpm+1),
                   border_color = NA,
                   clustering_method = 'complete',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'plu.genes.heatmap',
                   # breaks = unique(c(seq(0,8, length=256))),
                   legend_breaks = c(0,4,8),
                   fontsize_row = 8,
                   fontsize_col = 8)

pheatmap::pheatmap(t(scale(t(c2.mapk.tpm))),
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'mapk.genes.heatmap',
                   # breaks = unique(c(seq(0,8, length=256))),
                   legend_breaks = c(0,4,8),
                   fontsize_row = 8,
                   fontsize_col = 8)
# WNT-------
wnt.genes <- read.table('../4.part4/wnt.txt',sep = ',',header = T)
c2.wnt.tpm <- c2.tpm[colnames(wnt.genes),c(1,2,27:38)] %>% na.omit()

pheatmap::pheatmap(log2(c2.wnt.tpm+1),
                   border_color = NA,
                   clustering_method = 'complete',
                   color = viridis_pal(alpha = 1, begin = .1)(256),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'plu.genes.heatmap',
                   # breaks = unique(c(seq(0,8, length=256))),
                   legend_breaks = c(0,4,8),
                   fontsize_row = 8,
                   fontsize_col = 8)

pheatmap::pheatmap(t(scale(t(c2.wnt.tpm))),
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
                   angle_col = '315',
                   annotation_col = c2.info[,c(2,3,5)],
                   show_rownames = T,
                   annotation_colors = ha_top.col,
                   main = 'wnt.genes.heatmap',
                   breaks = unique(c(seq(-1,3, length=256))),
                   # legend_breaks = c(0,4,8),
                   fontsize_row = 8,
                   fontsize_col = 8)

library(reshape2)
c2.wnt.tpm.rank <- c2.tpm.rank[colnames(wnt.genes),c(1,2,27:38)]
apply(c2.wnt.tpm.rank, 2, median)
c2.wnt.tpm.rank.long <- melt(c2.wnt.tpm.rank)
colnames(c2.wnt.tpm.rank.long) <- c('wnt','rank')

library(ggstatsplot)

ggbetweenstats(
  data = c2.wnt.tpm.rank.long,
  x = wnt,
  y = rank,
  type = 'nonparametric',
  title = 'wnt',
  palette          = "default_igv",
  package          = "ggsci",
  p.adjust.method = 'none',
  pairwise.display = 'significant'
)





