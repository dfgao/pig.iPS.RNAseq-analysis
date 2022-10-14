# WGCNA for samples

library(WGCNA)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

# make data info -----
all.meta <- data.frame(Cell_type = geneinfo.all.new$ct3)
rownames(all.meta) <- geneinfo.all.new$Sample

datExpr <- t(all.tpm.clean.pick) 

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# 确定最佳soft-thresholding power--
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr, RsquaredCut = 0.9,powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col="red");
abline(h=0.90,col="green")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1,col="red")

best_beta=sft$powerEstimate

# build co-expresssion network-
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       maxBlockSize = 6000,TOMType = "unsigned", 
                       minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "All.tpm-TOM",
                       verbose = 3)

# Convert labels to colors for plotting
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

# module & Traits------

design=model.matrix(~0+ all.meta$Cell_type)
colnames(design)=unique(all.meta$Cell_type)
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
               colors =  blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",plotDendrograms = F,
                      marDendro = c(4,4,2,4))

# interesting module-
# names (colors) of the modules
modNames = 'MEdarkmagenta'
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

## 只有连续型性状才能只有计算
## 这里把是否属于 pg 表型这个变量用0,1进行数值化。
pg = as.data.frame(design[,16])
names(pg) = "pgEpiSC"
geneTraitSignificance = as.data.frame(cor(datExpr, pg, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(pg), sep="");
names(GSPvalue) = paste("p.GS.", names(pg), sep="")

module = "MEdarkmagenta"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for pg",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



# network vis-
#首先针对所有基因画热图
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
geneTree = net$dendrograms[[1]]; 
dissTOM = TOMsimilarityFromExpr(as.data.frame(t(na.omit(all.tpm.clean.pick[,c(1:24,31:34)])) ), power = 9); # 设置power=sft$powerEstimate，最佳beta值，此处是3
plotTOM = dissTOM^7; 
diag(plotTOM) = NA; 
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

#然后随机选取部分基因作图
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "ward.D2")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
