# add 
library(ggplot2)
library(ggsci)
library(scales)

# pluripotency in all samples
clean.plu.tpm <- all.tpm.clean.pick[key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)],]
colSums(clean.plu.tpm)

ann <- data.frame(geneinfo.new[,c(7,11)])
rownames(ann) <- geneinfo.new$Sample

ha_left.col <- list(ct3 = new.col[1:19], Author = pal_simpsons("springfield")(11))
names(ha_left.col$ct3) <- unique(geneinfo.new$ct3)
names(ha_left.col$Author) <- unique(geneinfo.new$Author)

pheatmap::pheatmap(log2(clean.plu.tpm+1),
                   border_color = NA,
                   clustering_method = 'complete',
                   color = viridis_pal(alpha = 1, begin = .2)(256),
                   angle_col = '315',
                   annotation_col = ann,
                   show_rownames = T,
                   annotation_colors = ha_left.col,
                   main = 'plu.genes.heatmap',
                   breaks = unique(c(seq(0,12, length=256))),
                   legend_breaks = c(0,4,8,12),
                   fontsize_row = 8,
                   fontsize_col = 4)

pheatmap::pheatmap(log2(clean.plu.tpm+1),
                   border_color = NA,
                   clustering_method = 'complete',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .9),
                   angle_col = '315',
                   annotation_col = ann,
                   show_rownames = T,
                   annotation_colors = ha_left.col,
                   main = 'plu.genes.heatmap',
                   breaks = unique(c(seq(0,12, length=256))),
                   legend_breaks = c(0,4,8,12),
                   fontsize_row = 8,
                   fontsize_col = 4)

# pluripotency in two culture system samples

sys.plu.tpm <- sys2.tpm[key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)],]

plu.genes.cluster.anno <- data.frame(cluster = kmectocluster[match(key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)], kmectocluster$gene),2])
rownames(plu.genes.cluster.anno) <- key.genes$Pluripotency$Pluripotency[-c(5,24,14,16,17)]
colnames(plu.genes.cluster.anno) <- 'Cluster'

tmp.gray <- colorRampPalette(colors = c('#495057','#343a40'),alpha=T,bias=1)(13)
ha_top.col <- list(Exogenous = new.col[c(16:17)], 
                   Author = c(pal_simpsons("springfield")(16)[7:14]), 
                   System = pal_igv("default")(11),
                   Cluster = c(tmp.gray[1:2],'#eb5e28',tmp.gray[4:6],'#e5383b',tmp.gray[13] ))
names(ha_top.col$Exogenous) <- factor(c("Silencing","Activation"))
names(ha_top.col$Author) <- factor(unique(sys2.info$Author))
names(ha_top.col$System) <- factor(unique(sys2.info$System))
names(ha_top.col$Cluster) <- factor(c( 'C_2(1435 genes)', 'C_3(474 genes)', 'C_4(1083 genes)',  'C_5(768 genes)', 'C_11(616 genes)', 'C_12(996 genes)', 'C_14(856 genes)', 'C_15(1578 genes)'))


ComplexHeatmap::pheatmap(log2(sys.plu.tpm+1),
                   border_color = 'NA',
                   clustering_method = 'ward.D2',
                   color = scales::alpha(colorRampPalette(colors = c('#00509d','gray80','#f35b04'),alpha=T,bias=1)(256),alpha = 1),
                   angle_col = '315',
                   annotation_col = sys.anno,
                   annotation_row = plu.genes.cluster.anno,
                   annotation_colors = ha_top.col,
                   main = 'ward.D2.tpm.plu',
                   breaks = unique(c(seq(0,12, length=256))),
                   legend_breaks = c(0,4,8,12),
                   fontsize_row = 8,
                   fontsize_col = 5)

ComplexHeatmap::pheatmap(log2(sys.plu.tpm+1),
                   border_color = 'NA',
                   clustering_method = 'ward.D2',
                   color = viridis_pal(alpha = 1, begin = .2)(256),
                   angle_col = '315',
                   annotation_col = sys.anno,
                   annotation_colors = ha_top.col,
                   main = 'ward.D2.tpm.pearson',
                   breaks = unique(c(seq(0,12, length=256))),
                   legend_breaks = c(0,4,8,12),
                   fontsize_row = 8,
                   fontsize_col = 5)


