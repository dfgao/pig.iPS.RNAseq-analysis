#----load packages-------

library(DESeq2)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(reshape2)
library(ggplot2)
library(Vennerable)
library(rio)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(pheatmap)
library(EnhancedVolcano)
library(patchwork)
library(pca3d)
library(rgl)
library(scatterplot3d)
library(WGCNA)
library(stringr)
library(scales)

#------load gene info-----
pig.gene.info <- import('~/project/99.gaodengfeng/01.multi.embory/06.analysis/00.datainfo/01.pig/01.mydata/sus.105.xlsx',which = 'pig_105_clean')
pig.gene.info <- pig.gene.info[pig.gene.info$biotype == 'protein_coding' | pig.gene.info$biotype == 'lncRNA' | pig.gene.info$biotype == 'pseudogene' | pig.gene.info$biotype == 'processed_pseudogene',]
which(pig.gene.info$genename =='DEFB1')
pig.gene.info[10314,2] <- 'AMELY.X'
pig.gene.info[14647,2] <- 'AMELY.Y'
pig.gene.info[12197,2] <- 'GZMA.1'
pig.gene.info[12198,2] <- 'GZMA.2'
pig.gene.info = pig.gene.info[-c(14831),]


#-----data input------
## public data
public.data <-read.table('../../../02.public.ips/03.quant/01.count/IPS-ESC-ref105.remove.MT.count.txt',header = T,sep = '\t',comment.char = '#',row.names = 1)
meta.data <- public.data[,c(1:5)]
public.data <- public.data[,-c(1:5)]
public.data.info <- data.frame(Author = c('Gao','Gao',rep('Xu',3),'Shi','Shi',rep('Kinoshita',9), rep('Yuan',4), 'Yoshimatsu','Yoshimatsu',rep('Yuan',2), rep('Secher',6), rep('Mao',4), rep('Gao',4), rep('Choi',9), rep('Yuan',6)),
                               Cell_type = c(rep('iPS',7), rep('EDSC',9), rep('iPS',18), rep('EPSC',4), rep('ESC',9), rep('ESCLC',6))
)

rownames(public.data.info) <- colnames(public.data)
public.data.info <- import('../../zqq.ips/public.data.info.xlsx')
colnames(public.data) <- public.data.info$Sample

## ngr data
ngr.data <- read.table('../../../01.NGR_piglet_ear_YH_ips/04.quant/01.count/NGR.IPS.ref105.remove.MT.count.txt',header = T,sep = '\t',comment.char = '#',row.names = 1)
ngr.data <- ngr.data[,-c(1:5)]
ngr.data.info <- data.frame(Sample = c('pg_1','pg_2','NGR.iPS_1_rep1','NGR.iPS_1_rep2','NGR.iPS_3_rep1','NGR.iPS_3_rep2','pEF_1','pEF_2'),
                            Author = c(rep('Zhi',2), rep('Zhu',4), rep('Xu',2)),
                            Cell_type = c(rep('pgEpiSC',2), rep('NGR.iPS',4), rep('pEF',2)))
colnames(ngr.data) <- ngr.data.info$Sample

## pef.ips
pef.ips.data <- read.table('../../../03.PEF.ips/04.quant/PEF-IPS-ref105.remove.MT.count.txt',header = T,sep = '\t',comment.char = '#',row.names = 1)
pef.ips.data <- pef.ips.data[,-c(1:5)]
pef.ips.info <- data.frame(Sample = c('PEF_EPI_5_iPS_rep1','PEF_EPI_5_iPS_rep2','PEF_EPI_8_iPS_rep1','PEF_EPI_8_iPS_rep2','PEF_pMX_10_iPS_rep1','PEF_pMX_10_iPS_rep2','PEF_pMX_6_iPS_rep1','PEF_pMX_6_iPS_rep2'),
                           Author = c(rep('Zhu',8)),
                           Cell_type = c(rep('pEF.iPS',4), rep('pEF.NGR.iPS',4)))
colnames(pef.ips.data) <- pef.ips.info$Sample

## all data
all.data <- cbind(public.data, ngr.data, pef.ips.data)
all.data.info <- rbind(public.data.info, ngr.data.info, pef.ips.info)
all.data <- all.data %>% rownames_to_column(var = 'geneid') %>% left_join( y = pig.gene.info[pig.gene.info$biotype=='protein_coding',1:2], by = 'geneid') %>% na.omit()
rownames(all.data) <- all.data$genename
all.data <- all.data[,-grep('gene', colnames(all.data))]

## all ips
all.ips.info <- all.data.info[all.data.info$Cell_type=='iPS' | all.data.info$Cell_type == 'NGR.iPS' | all.data.info$Cell_type == 'pEF.iPS' | all.data.info$Cell_type == 'pEF.NGR.iPS',]
all.ips.info$Exogenous = c(rep('Silencing',2), rep('Activation',5), rep('Silencing',2), rep('Activation',2), rep('Silencing',4), rep('Activation',10), rep('Silencing',12))
all.ips.info$Medium = c(rep('ALCIW',2), rep('LCDMV',3), rep('2i_LIF',2), rep('FL6i',2), rep('FLB2i',2), rep('AFI',2), rep('FGF',2), rep('FGF',3), rep('LIF',3), rep('VPA_LIF',4), rep('3i_LAF',12))

all.ips.data <- all.data[,all.ips.info$Sample]

#----get TPM------
meta.data.pcg <- meta.data[,c(1,5)] %>% rownames_to_column(var = 'geneid') %>% left_join(y = pig.gene.info[pig.gene.info$biotype=='protein_coding',1:2], by = 'geneid') %>% na.omit()

kb = meta.data.pcg$Length / 1000 
rpk <- all.data / kb
all.tpm <- data.frame(t(t(rpk)/colSums(rpk) * 1000000))
rm(kb,rpk)
all.tpm <- round(all.tpm, 2)

all.ips.tpm <- all.tpm[,all.ips.info$Sample]

#----pick.genes-----
all.pick.genes <- unique(c(rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 1:2) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 3:5) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 6:7) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 8:16) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 17:18) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 19:20) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 21:22) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 23:24) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 25:27) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 28:30) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 31:34) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 35:36) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 37:38) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 39:41) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 42:44) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 45:47) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 48:49) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 50:51) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 52:53) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 54:55) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 56:57) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 58:59) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 60:61) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 62:63) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 64:65) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 66:67) > 0.5) ,]),
                           rownames(all.tpm[(matrixStats::rowMedians(as.matrix(all.tpm), cols = 68:69) > 0.5) ,])
))

all.data.clean <- all.data[all.pick.genes,]
all.tpm.clean <- all.tpm[all.pick.genes,]

all.ips.data.clean <- all.ips.data[all.pick.genes,]
all.ips.tpm.clean <- all.ips.tpm[all.pick.genes,]


