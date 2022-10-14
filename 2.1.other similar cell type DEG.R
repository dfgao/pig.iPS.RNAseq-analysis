#

# DEG8: gao iPS vs PEF
group = all.data.clean[,c(1,2,60,61)]
group_list = c(rep('gao.iPS',2),rep('pEF',2)) 
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','gao.iPS','pEF'))
DEG8 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

# DEG9: gao EPSC vs PEF
group = all.data.clean[,c(37,38,60,61)]
group_list = c(rep('gao.EPSC',2),rep('pEF',2)) 
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','gao.EPSC','pEF'))
DEG9 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

# DEG10: EDSC vs PEF
group = all.data.clean[,c(8:10,60,61)]
group_list = c(rep('EDSC',3),rep('pEF',2)) 
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','EDSC','pEF'))
DEG10 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

# DEG11: ESC vs PEF
group = all.data.clean[,c(39:41,60,61)]
group_list = c(rep('ESC',3),rep('pEF',2)) 
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','ESC','pEF'))
DEG11 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

# DEG12: Yoshi-ips vs PEF
group = all.data.clean[,c(21:22,60,61)]
group_list = c(rep('AFI',2),rep('pEF',2))
colData = data.frame(row.names = colnames(group),group_list = group_list)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','AFI','pEF'))
DEG12 = as.data.frame(RES) %>% rownames_to_column(var = 'Symbol') %>% na.omit() %>% arrange(padj)

