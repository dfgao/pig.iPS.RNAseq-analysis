# part1. global

#----gene info-----
geneinfo <- data.frame(TPM_100 = apply(all.tpm.clean,2,function(x) {table(x>=100)["TRUE"]}),
                       TPM_100_50 = apply(all.tpm.clean,2,function(x) {table(x>=50 & x< 100)["TRUE"]}),
                       TPM_50_5 = apply(all.tpm.clean,2,function(x) {table(x>=5 & x< 50)["TRUE"]}),
                       TPM_5_0.5 = apply(all.tpm.clean,2,function(x) {table(x>=0.5 & x< 5)["TRUE"]}),
                       TPM_0.5 = apply(all.tpm.clean,2,function(x) {table(x< 0.5)["TRUE"]})
)

geneinfo.all <- geneinfo %>% rownames_to_column( var = 'Sample') %>% left_join(y = all.data.info,by = 'Sample')
geneinfo.all$ct2 <- paste0(geneinfo.all$Author,'_',geneinfo.all$Cell_type)
export(geneinfo.all,file = './geneinfo.all.xlsx')
geneinfo.all <- import('./geneinfo+culture medium info.xlsx')
geneinfo.all$ct3 <- paste(geneinfo.all$ct2, geneinfo.all$System,sep = '_')

new.col <- c("#E7959B","#DB5E67","#CFD99F",'#C5EC41','#E5CD94',"#FC0C00","#7C4374","#339E55","#000376","#2A82C6","#8C6D36","#CB70C0",
             "#EBB854",'#FC8D37',"#63753A","#6D2811","#DD9AD2","#68AADE","#3B397C","#9D9AE5","#B8CF6E","#949494","#BF4F8F","#844346")

#------ PCA------

# FactoMineR
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
               legend.title = "Groups",
               pointsize = 4,
               pointshape = 21,
               col.ind = "black",
               title = 'Samples PCA',
               # addEllipses = T
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}

# remove yuanye-fgf4/fgf2 xu-pef liu K3-----BEST 
group_list <- factor(geneinfo.all$ct3[-c(23,24,48,49,60,61,35,36)],levels = unique(geneinfo.all$ct3[-c(23,24,48,49,60,61,35,36)]))
d_p(all.tpm.clean[,-c(23,24,48,49,60,61,35,36)],group_list)
d_p(all.data.clean[,-c(23,24,48,49,60,61,35,36)],group_list)

# plot 3d PCA-scatterplot3d 
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
group_list <- factor(geneinfo.all$ct2[-c(23,24,48,49,60,61,35,36)],levels = unique(geneinfo.all$ct2[-c(23,24,48,49,60,61,35,36)]))
tpm.t = as.data.frame(t(all.tpm.clean[,-c(23,24,48,49,60,61,35,36)]))
tpm.t = cbind(tpm.t,group_list)
tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
a = tpm.pca[["ind"]][["coord"]] 

par(mar=c(1,1,1,3),oma = c(1,1,1,8))
p3 = scatterplot3d(a[,1:3],
                   color = 'black',
                   main="PCA 3D",
                   pch = 21,
                   bg = new.col[as.numeric(group_list)],
                   cex.symbols = 1.2,
                   type="p",
                   # grid = F,
                   box = T)
addgrids3d(a[, 1:3], grid = c("xy", "xz", "yz"))
legend("right",
       col = "black", 
       legend = levels(group_list),
       pt.bg =  new.col, 
       pch = 21,
       inset = -0.5, 
       xpd = T, 
       horiz = F,
       bty = 'n',
       cex = .8,
       bg = "transparent")

# use rgl plot 3D PCA 
open3d() 
par3d(family="serif",
      cex=1.2,
      font=1) 

plot3d(a[,1:3],
       type="p",
       col=new.col[as.numeric(group_list)],
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

rgl.postscript("../rgl.PCA-dim3-3D.pdf", fmt = "pdf", drawText = T)

close3d() 

# save(new.col, geneinfo.all,all.tpm.clean,all.tpm.clean.pick,file = 'for.3D.PCA.Rdata')


# get new data table----------

all.tpm.clean.pick <- all.tpm.clean[-c(23,24,48,49,60,61,35,36)]
dim(all.tpm.clean.pick)
all.tpm.clean.pick.pick.genes <- unique(c(rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 1:2) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 3:5) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 6:7) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 8:16) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 17:18) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 19:20) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 21:22) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 23:25) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 26:28) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 29:32) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 33:34) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 35:37) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 38:40) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 41:43) > 0.5) ,]),
                                          # rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 44:45) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 44:45) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 46:47) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 48:49) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 50:51) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 52:53) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 54:55) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 56:57) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 58:59) > 0.5) ,]),
                                          rownames(all.tpm.clean.pick[(matrixStats::rowMedians(as.matrix(all.tpm.clean.pick), cols = 60:61) > 0.5) ,])
))
all.tpm.clean.pick <- all.tpm.clean.pick[all.tpm.clean.pick.pick.genes,] %>% na.omit()

all.data.clean.pick <- all.data[all.tpm.clean.pick.pick.genes,colnames(all.tpm.clean.pick)] %>% na.omit()
dim(all.data.clean.pick)

# gene info new--------
geneinfo.new <- data.frame(TPM_100 = apply(all.tpm.clean.pick,2,function(x) {table(x>=100)["TRUE"]}),
                       TPM_100_50 = apply(all.tpm.clean.pick,2,function(x) {table(x>=50 & x< 100)["TRUE"]}),
                       TPM_50_5 = apply(all.tpm.clean.pick,2,function(x) {table(x>=5 & x< 50)["TRUE"]}),
                       TPM_5_0.5 = apply(all.tpm.clean.pick,2,function(x) {table(x>=0.5 & x< 5)["TRUE"]}),
                       TPM_0.5 = apply(all.tpm.clean.pick,2,function(x) {table(x< 0.5)["TRUE"]})
)

geneinfo.all <- import('./geneinfo+culture medium info.xlsx')
geneinfo.all <- geneinfo.all[,-c(2:6)]

geneinfo.new <- geneinfo.new %>% rownames_to_column(var = 'Sample') %>% left_join(y = geneinfo.all,by = 'Sample')
geneinfo.new$ct3 <- paste0(geneinfo.new$ct2,"_",geneinfo.new$System)
geneinfo.new$ct3 <- factor(geneinfo.new$ct3,levels = unique(geneinfo.new$ct3))

test <- melt(geneinfo.new[,c(1:7)],id.vars = c('Sample','Author'),variable.name = 'Range',value.name = 'genenum')
test$genenum <- as.integer(test$genenum)
test$Sample <- factor(test$Sample, levels = geneinfo.new$Sample)

ggplot(test) +
  theme_bw() +
  # theme_ipsum()+
  geom_bar(stat = "identity",aes(x=Sample,y=genenum,fill=Range),alpha=.8,position = 'stack') +
  # geom_vline(xintercept = c(3.5,6.5,8.5,11.5,14.5)) +
  scale_fill_brewer(palette = "RdBu",direction = 1,) +
  labs(y='Gene number',x=NULL) +
  theme(axis.text.x = element_text(angle = 45,margin = 
                                     margin(1,0,0,0,'cm'),size = 6),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15)
  )

group_list <- factor(geneinfo.new$ct3,levels = unique(geneinfo.new$ct3))
d_p(all.tpm.clean.pick,group_list)

#-----gene complex-----
library(GGally)
library(grid)
all.tpm.sort <- data.frame(apply(all.tpm.clean.pick, 2, function(x){sort(x,decreasing = T)}))
gene.comp <- all.tpm.sort

for (i in colnames(all.tpm.sort)) {
  gene.comp[,i] <- all.tpm.sort[,i]/colSums(all.tpm.sort)[i]
}

gene.comp <- data.frame(apply(gene.comp, 2, cumsum))

group_list <- factor(geneinfo.all.new$Author,levels = unique(geneinfo.all.new$Author))

plot.data <- data.frame(t(gene.comp),Group=group_list) %>% rownames_to_column(var = 'Sample')
colnames(plot.data) <- c('Sample',(seq(1:16923)),'Group')
plot.data.log <- melt(plot.data,id.vars = c('Group','Sample'),variable.name = 'Genenum',value.name = 'Percent')
plot.data.log$Genenum <- as.numeric(plot.data.log$Genenum)
plot.data.log$Sample <- factor(plot.data.log$Sample, levels = geneinfo.all.new$Sample)
plot.data.log$Percent <- plot.data.log$Percent*100

cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=1)(6))
cols <- as.matrix(apply(data.frame(cols),1,function(x){rep(x,12)}))
cols <- as.vector(cols)
ggplot(data = plot.data.log) +
  theme_bw() +
  geom_line(aes(x=Genenum, y=Percent,color=Sample),size=.7 ) + 
  scale_x_log10() +
  scale_color_manual(values = new.col[as.numeric(group_list)]) +
  labs(x = 'Genes number',y = 'Gene accumulation ratio (%)') +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.key.height=unit(.85,"line"))

# correlation-------
library(ComplexHeatmap)
library(circlize)
ann <- data.frame(geneinfo.new[,c(7,11)])
rownames(ann) <- geneinfo.new$Sample

ha_left.col <- list(ct3 = new.col[1:19], Author = pal_simpsons("springfield")(11))
names(ha_left.col$ct3) <- unique(geneinfo.new$ct3)
names(ha_left.col$Author) <- unique(geneinfo.new$Author)

all.cor.p <- cor(all.tpm.clean.pick)
ComplexHeatmap::pheatmap(all.cor.p,
                   border_color = NA,
                   clustering_method = 'ward.D2',
                   color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .9),
                   angle_col = '315', 
                   annotation_row = ann,
                   annotation_colors = ha_left.col,
                   main = 'average.tpm.pearson',
                   fontsize = 8)


all.cor.s <- cor(all.tpm.clean.pick,method = 'spearman')
ComplexHeatmap::pheatmap(all.cor.s,
                         border_color = NA,
                         clustering_method = 'ward.D',
                         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=.8)(256)),alpha = .9),
                         angle_col = '315', 
                         annotation_row = ann,
                         annotation_colors = ha_left.col,
                         main = 'ward.D2.tpm.spearman',
                         fontsize = 5)



#----Tsne fail--------
library(Rtsne)
set.seed(42)

test <- t(unique(all.tpm.clean.pick))

tsne.all.batch <- Rtsne(test,pca=FALSE,dims=2,
                        perplexity=7,theta=.5) # perplexity 5 best

# head(tsne.all.batch)

tsne_res <- as.data.frame(tsne.all.batch$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)


# 使用ggplot2可视化tSNE降维的结果 __ 最好是去掉PEF再进行批次效应矫正
ggplot(tsne_res,aes(tSNE1,tSNE2,color=geneinfo.all.new$ct2)) + 
  geom_point(size=4) + theme_bw() + 
  scale_color_simpsons() +
  geom_hline(yintercept = 0,lty=2,col="red") + 
  geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "tSNE plot",color="Species")

